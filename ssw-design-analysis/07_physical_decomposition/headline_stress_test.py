#!/usr/bin/env python3
"""
headline_stress_test.py
=======================
THE THREE ATTACKS THAT DECIDE WHETHER THE HEADLINE SURVIVES REVIEW.

`predictability_ceiling.py` reports pre-onset out-of-sample R^2 = 0.115 against
0.422 with post-onset stratospheric information, and concludes that 73% of
apparent diagnostic skill comes from the classifier window overlapping the
response window. Three objections can each destroy that, and none had been tested.

ATTACK 1 -- "THAT IS MEDIATION, NOT ARTEFACT."
  P3 predicts the SURFACE from post-onset STRATOSPHERIC state. A referee will say
  the gain is real physics: the SSW disturbs the lower stratosphere, which
  genuinely transmits downward, so of course knowing it helps. On that reading
  0.422 is a true causal chain and only its *timing* disqualifies it as a
  forecast -- which is a far weaker claim than the one being made.

  THE TEST: run the identical P3 pipeline on PSEUDO-onsets drawn from
  event-cleaned days, where no SSW exists and no mediation is possible. If P3
  still scores ~0.42 there, the skill is structural stratosphere-surface
  covariance and the artefact reading is right. If it collapses to ~0, the
  mediation reading is right and the headline must be reworded.

ATTACK 2 -- "CMIP6 UNDER-COUPLES, SO 0.115 IS TOO LOW."
  This project's own `FINDING_ensemble_precursor.md` measures an observed/model
  response ratio of 1.74, and CanESM5 and MIROC6 disagree in SIGN on the
  precursor. If weak coupling suppresses predictability, the real atmosphere
  could be far more predictable than 0.115 and the whole argument inverts.

  THE TEST: compute CV R^2 SEPARATELY per model, and regress it on that model's
  own coupling strength (its post-onset composite magnitude). A positive slope
  means predictability scales with coupling and the ensemble estimate is a floor
  that must be extrapolated. A flat slope means R^2 is a structural property that
  does not inherit the coupling bias.

ATTACK 3 -- "YOU COMPARED OUT-OF-SAMPLE AGAINST IN-SAMPLE."
  The published splits' implied R^2 of 0.63 comes from a contrast fitted in
  sample, on an outcome-based split. Comparing it against a cross-validated 0.115
  is not like for like, and the gap could be nothing but the ordinary optimism of
  an unvalidated fit.

  THE TEST: report this model's OWN in-sample R^2 alongside its out-of-sample
  R^2. If in-sample lands near 0.6 while out-of-sample is 0.115, then 0.63 is
  precisely what an unvalidated fit produces on this data, and the comparison
  becomes stronger rather than weaker -- but stated in the correct terms.

Output: headline_stress_test.json
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

from sklearn.linear_model import RidgeCV            # noqa: E402
from sklearn.pipeline import make_pipeline          # noqa: E402
from sklearn.preprocessing import StandardScaler    # noqa: E402

RAW = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "raw" / "cmip6"
OUT_WIN = (8, 52)
SEASON = (11, 12, 1, 2, 3, 4)
SEED = 20260805


def doy_std(s):
    d = s.index.dayofyear
    return (s - s.groupby(d).transform("mean")) / s.groupby(d).transform("std")


def features_for(u, plev, lat, idx, onsets):
    """Identical feature construction to predictability_ceiling, any date list."""
    i60 = int(np.argmin(np.abs(lat - 60.0)))
    cap = lat >= 65.0
    w = np.cos(np.deg2rad(lat[cap]))
    feat = {}
    for lev in P.LEVELS_PA:
        j = int(np.argmin(np.abs(plev - lev)))
        if not np.isclose(plev[j], lev, rtol=0.1):
            continue
        s60 = pd.Series(u[:, j, i60], index=idx)
        scap = pd.Series(np.average(u[:, j, cap], weights=w, axis=1), index=idx)
        for src, tag in ((s60, "60N"), (scap, "cap")):
            z = doy_std(src)
            for wname, win, _ in P.WINDOWS:
                feat[f"{P.LEVEL_LAB[lev]}_{tag}_{wname}"] = P._wm(z, onsets, win)
    o = pd.DatetimeIndex(onsets)
    feat["doy_sin"] = np.sin(2 * np.pi * o.dayofyear / 365.25)
    feat["doy_cos"] = np.cos(2 * np.pi * o.dayofyear / 365.25)
    return pd.DataFrame(feat)


def in_sample_r2(X, y):
    ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
    X, y = X[ok], y[ok]
    m = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 4, 25)))
    m.fit(X, y)
    p = m.predict(X)
    return float(1 - np.sum((y - p) ** 2) / np.sum((y - y.mean()) ** 2))


def main():
    rng = np.random.default_rng(SEED)
    res = {}

    real_rows, real_Y, real_g = [], [], []
    ps_rows, ps_Y, ps_g = [], [], []
    per_model = []

    print("building real and PSEUDO event tables from CMIP6 ...")
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
        real_rows.append(features_for(u, plev, lat, idx, on))
        real_Y.append(Yr)
        real_g.append(np.full(len(Yr), f.stem))
        # coupling strength for ATTACK 2: this model's own post-onset composite
        per_model.append({"member": f.stem, "coupling": float(np.nanmean(Yr)),
                          "n": int(np.isfinite(Yr).sum())})

        # ---- ATTACK 1: pseudo-onsets, no event, matched day-of-year ----
        for _ in range(6):
            p_on = C6.draw_pseudo(cln, dy, rng)
            if len(p_on) < 10:
                continue
            Yp = C6.anom(am, p_on, cl, OUT_WIN)
            ps_rows.append(features_for(u, plev, lat, idx, p_on))
            ps_Y.append(Yp)
            ps_g.append(np.full(len(Yp), f.stem))

    Xr = pd.concat(real_rows, ignore_index=True)
    yr = np.concatenate(real_Y)
    gr = np.concatenate(real_g)
    Xp = pd.concat(ps_rows, ignore_index=True)
    yp = np.concatenate(ps_Y)
    gp = np.concatenate(ps_g)
    print(f"real events {np.isfinite(yr).sum():,} | pseudo-events {np.isfinite(yp).sum():,}")

    # ================= ATTACK 1 =================
    print("\n" + "=" * 76)
    print("=== ATTACK 1: is the P3 gain mediation, or structural covariance? ===")
    print("  If P3 scores ~0.42 on PSEUDO-events (no SSW, no mediation possible),")
    print("  the skill is structural. If it collapses, it is real mediation.\n")
    a1 = {}
    print(f"  {'tier':<32s} {'real':>9s} {'pseudo':>9s}  reading")
    print("  " + "-" * 68)
    for tier, lab in ((1, "P1 pre-onset"), (2, "P2 at-onset"),
                      (3, "P3 + post-onset stratosphere")):
        cols = P.tier_cols(Xr.columns, tier) + ["doy_sin", "doy_cos"]
        rr = P.cv_r2(Xr[cols].values, yr, gr, "ridge")
        pp = P.cv_r2(Xp[cols].values, yp, gp, "ridge")
        a1[lab] = {"real_cv_r2": round(rr, 4), "pseudo_cv_r2": round(pp, 4),
                   "retained_fraction": round(float(pp / rr), 3) if rr else None}
        note = ("structural" if pp > 0.5 * rr else
                "mostly event-related" if pp > 0.15 * rr else "event-related")
        print(f"  {lab:<32s} {rr:+9.4f} {pp:+9.4f}  {note}")
    res["attack1_pseudo_event_skill"] = a1

    # ================= ATTACK 2 =================
    print("\n" + "=" * 76)
    print("=== ATTACK 2: does R^2 scale with a model's coupling strength? ===")
    print("  If yes, 0.115 is a FLOOR and must be extrapolated to the real")
    print("  atmosphere, which couples ~1.74x more strongly than CanESM5.\n")
    pm = pd.DataFrame(per_model)
    cols2 = P.tier_cols(Xr.columns, 2) + ["doy_sin", "doy_cos"]
    rows = []
    for mem in pm.member:
        sel = gr == mem
        if np.isfinite(yr[sel]).sum() < 60:
            continue
        # within-member CV needs groups; split on winter-decade blocks instead
        yy = yr[sel]
        XX = Xr[cols2].values[sel]
        blk = (np.arange(sel.sum()) // 12).astype(str)
        r = P.cv_r2(XX, yy, blk, "ridge", n_splits=5)
        c = float(pm.loc[pm.member == mem, "coupling"].iloc[0])
        rows.append({"member": mem, "coupling": c, "cv_r2": r})
    d2 = pd.DataFrame(rows).dropna()
    print(f"  {'member':<28s} {'coupling':>9s} {'CV R^2':>9s}")
    print("  " + "-" * 50)
    for r_ in d2.itertuples():
        print(f"  {r_.member:<28s} {r_.coupling:+9.3f} {r_.cv_r2:+9.4f}")
    if len(d2) >= 5:
        # coupling is NEGATIVE (a negative AO anomaly); stronger = more negative
        sl, ic = np.polyfit(np.abs(d2.coupling), d2.cv_r2, 1)
        rho = float(np.corrcoef(np.abs(d2.coupling), d2.cv_r2)[0, 1])
        print(f"\n  regression of CV R^2 on |coupling|: slope {sl:+.4f}, r = {rho:+.3f}")
        obs_coupling = 1.026          # observed |post-onset AO| at +15..+29
        pred = ic + sl * obs_coupling
        print(f"  extrapolated to the observed coupling of {obs_coupling:.3f}: "
              f"R^2 ~ {pred:+.4f}")
        res["attack2"] = {"slope": round(float(sl), 5), "r": round(rho, 4),
                          "intercept": round(float(ic), 5),
                          "extrapolated_to_observed": round(float(pred), 4),
                          "per_member": d2.round(4).to_dict("records")}

    # ================= ATTACK 3 =================
    print("\n" + "=" * 76)
    print("=== ATTACK 3: in-sample vs out-of-sample, the like-for-like check ===")
    a3 = {}
    for tier, lab in ((1, "P1 pre-onset"), (2, "P2 at-onset"),
                      (3, "P3 + post-onset stratosphere")):
        cols = P.tier_cols(Xr.columns, tier) + ["doy_sin", "doy_cos"]
        ins = in_sample_r2(Xr[cols].values, yr)
        oos = P.cv_r2(Xr[cols].values, yr, gr, "ridge")
        a3[lab] = {"in_sample_r2": round(ins, 4), "out_of_sample_r2": round(oos, 4),
                   "optimism": round(float(ins - oos), 4)}
        print(f"  {lab:<32s} in-sample {ins:+.4f}   out-of-sample {oos:+.4f}   "
              f"optimism {ins - oos:+.4f}")
    res["attack3_in_vs_out"] = a3
    print("\n  The published splits' implied R^2 of 0.63 is an UNVALIDATED,")
    print("  outcome-based, in-sample quantity. Compare it against the in-sample")
    print("  column, not the out-of-sample one.")

    (RESULTS / "headline_stress_test.json").write_text(json.dumps(res, indent=2),
                                                    encoding="utf8")
    print("\nSaved -> headline_stress_test.json")


if __name__ == "__main__":
    main()
