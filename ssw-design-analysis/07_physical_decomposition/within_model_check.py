#!/usr/bin/env python3
"""
within_model_check.py
=====================
IS THE HEADLINE R^2 WITHIN-MODEL SKILL, OR BETWEEN-MODEL SIGNAL?

THE PROBLEM, FOUND BY STRESS-TESTING RATHER THAN BY A REFEREE
  `headline_stress_test.py` computed CV R^2 separately inside each ensemble
  member and found it scattered around ZERO (-0.18 to +0.17 across 20 members),
  while the pooled cross-member estimate is +0.115.

  That gap has an innocent reading and a fatal one.

  INNOCENT: within-member n is only ~94 against 38 features, so per-member CV is
  swamped by estimation noise and its near-zero average means nothing.

  FATAL: GroupKFold by member stops the model memorising a member, but it does
  NOT stop it exploiting BETWEEN-MODEL structure. If models differ systematically
  in both their typical vortex anomaly and their typical surface response, a
  regression can learn that cross-model relationship and score well on held-out
  members without any ability to rank events INSIDE one climate. The field's
  question -- given an SSW in the real atmosphere, will this one couple
  downward? -- is a within-climate question. Between-model skill does not answer
  it, and reporting 0.115 as if it did would be wrong.

THE TEST
  Standardise Y and every feature WITHIN each member (subtract that member's mean,
  divide by its sd) before pooling. That removes all between-model level
  differences by construction while leaving every within-model event-to-event
  relationship intact. Then rerun the identical pooled GroupKFold pipeline.

    R^2 holds near 0.115  ->  the skill is WITHIN-model. Headline stands.
    R^2 collapses to ~0   ->  the skill was BETWEEN-model. The headline must be
                              restated as a bound and the paper's central number
                              changes.

  Pseudo-events are carried through the same transformation, so the artefact
  share measured in Attack 1 can be re-read on the same footing.

Output: within_model_check.json
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "6_predictability"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]
SIM = ROOT / "ssw-design-analysis" / "06_simulation_validation"
sys.path.insert(0, str(SIM))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
import ensemble_precursor as EP                     # noqa: E402
import stratifier_law_cmip6 as C6                   # noqa: E402
import predictability_ceiling as P                  # noqa: E402
import headline_stress_test as H                    # noqa: E402

RAW = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "raw" / "cmip6"
OUT_WIN = (8, 52)
SEASON = (11, 12, 1, 2, 3, 4)
SEED = 20260805


def demean_within(X, y, g):
    """Remove each member's own mean and scale, leaving only within-member variation."""
    Xo = X.copy().astype(float)
    yo = y.copy().astype(float)
    for m in np.unique(g):
        i = g == m
        for c in range(Xo.shape[1]):
            v = Xo[i, c]
            s = np.nanstd(v)
            Xo[i, c] = (v - np.nanmean(v)) / (s if s > 1e-12 else 1.0)
        sy = np.nanstd(yo[i])
        yo[i] = (yo[i] - np.nanmean(yo[i])) / (sy if sy > 1e-12 else 1.0)
    return Xo, yo


def main():
    rng = np.random.default_rng(SEED)
    real_rows, real_Y, real_g = [], [], []
    ps_rows, ps_Y, ps_g = [], [], []

    print("rebuilding CMIP6 real and pseudo tables ...")
    for f in sorted(RAW.glob("*_zm.nc")):
        try:
            m = EP.load_member(f)
        except Exception:
            continue
        m2 = m[np.isin(m.index.month, SEASON)].dropna()
        if len(m2) < 2000:
            continue
        on = EP.detect_ssw(m2["u10"].values, m2.index)
        if len(on) < 15:
            continue
        am = m2["am"]
        msk = C6.influence_mask(am.index, on)
        cl = am[~msk].groupby(am[~msk].index.dayofyear).mean()
        cln = am[~msk].index
        dy = np.array([t.dayofyear for t in pd.DatetimeIndex(on)])
        ds = xr.open_dataset(f)
        lat, plev, u = ds["lat"].values, ds["plev"].values, ds["u_zm"].values
        ds.close()
        idx = m.index
        if len(idx) != u.shape[0]:
            continue
        Yr = C6.anom(am, on, cl, OUT_WIN)
        real_rows.append(H.features_for(u, plev, lat, idx, on))
        real_Y.append(Yr)
        real_g.append(np.full(len(Yr), f.stem))
        for _ in range(6):
            p_on = C6.draw_pseudo(cln, dy, rng)
            if len(p_on) < 10:
                continue
            Yp = C6.anom(am, p_on, cl, OUT_WIN)
            ps_rows.append(H.features_for(u, plev, lat, idx, p_on))
            ps_Y.append(Yp)
            ps_g.append(np.full(len(Yp), f.stem))

    Xr = pd.concat(real_rows, ignore_index=True)
    yr = np.concatenate(real_Y)
    gr = np.concatenate(real_g)
    Xp = pd.concat(ps_rows, ignore_index=True)
    yp = np.concatenate(ps_Y)
    gp = np.concatenate(ps_g)
    print(f"real {np.isfinite(yr).sum():,} | pseudo {np.isfinite(yp).sum():,} | "
          f"{len(np.unique(gr))} members")

    res = {}
    print("\n" + "=" * 82)
    print("=== POOLED (as published) vs WITHIN-MODEL STANDARDISED ===")
    print(f"{'tier':<32s} {'pooled':>9s} {'within':>9s} {'retained':>9s} "
          f"{'pseudo-w':>9s}")
    print("-" * 76)
    for tier, lab in ((1, "P1 pre-onset"), (2, "P2 at-onset"),
                      (3, "P3 + post-onset stratosphere")):
        cols = P.tier_cols(Xr.columns, tier) + ["doy_sin", "doy_cos"]
        Xv, Xpv = Xr[cols].values, Xp[cols].values
        pooled = P.cv_r2(Xv, yr, gr, "ridge")
        Xw, yw = demean_within(Xv, yr, gr)
        within = P.cv_r2(Xw, yw, gr, "ridge")
        Xpw, ypw = demean_within(Xpv, yp, gp)
        pseudo_w = P.cv_r2(Xpw, ypw, gp, "ridge")
        ret = within / pooled if pooled else np.nan
        res[lab] = {"pooled_cv_r2": round(pooled, 4),
                    "within_model_cv_r2": round(within, 4),
                    "retained_fraction": None if not np.isfinite(ret) else round(float(ret), 3),
                    "pseudo_within_cv_r2": round(pseudo_w, 4),
                    "event_specific": round(float(within - pseudo_w), 4)}
        print(f"{lab:<32s} {pooled:+9.4f} {within:+9.4f} {ret:9.2f} {pseudo_w:+9.4f}")

    p1 = res["P1 pre-onset"]
    p3 = res["P3 + post-onset stratosphere"]
    print("\n=== READING ===")
    if p1["within_model_cv_r2"] < 0.3 * p1["pooled_cv_r2"]:
        print("  P1 skill COLLAPSES once between-model level differences are removed.")
        print("  The pooled 0.115 was largely BETWEEN-MODEL signal, and the honest")
        print("  within-climate predictable share is far smaller. The headline number")
        print("  must be restated.")
    else:
        print("  P1 skill SURVIVES within-model standardisation, so it is genuine")
        print("  event-to-event predictability inside a single climate, which is the")
        print("  quantity the field's question is about. Headline stands.")
    print(f"\n  event-specific within-model share (within minus pseudo-within):")
    for k, v in res.items():
        print(f"    {k:<32s} {v['event_specific']:+.4f}")

    (RESULTS / "within_model_check.json").write_text(json.dumps(res, indent=2),
                                                  encoding="utf8")
    print("\nSaved -> within_model_check.json")


if __name__ == "__main__":
    main()
