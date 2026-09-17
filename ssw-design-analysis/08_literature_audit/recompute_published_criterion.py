#!/usr/bin/env python3
"""
recompute_published_criterion.py
================================
THE BODY COUNT. Applies the published classification exactly as specified and
decomposes the number it produces.

THE CRITERION, VERIFIED FROM SOURCE 2026-08-01
  Karpechko et al. (2017), QJRMS 143, 1459 (doi 10.1002/qj.3017), 229 citations,
  still ~28/yr in 2024. An SSW is "downward propagating" (DW) if, over days
  +8..+52 after onset:
    1. the mean NAM index at 1000 hPa is negative;
    2. the fraction of those 45 days with negative 1000 hPa NAM is > 0.5;
    3. the fraction with negative 150 hPa NAM is > 0.7.
  Conditions 1 and 2 are the surface response itself.

  In active use. ACP 26, 3723 (2026) states "The definition of DWs used in this
  paper follows the method by Karpechko et al. (2017)" (moving 1000 -> 850 hPa)
  and reports, for 60 days post-onset in ERA5, mean NAO of -0.762 (BOTH subtype),
  -0.567 (EA), -0.435 (NA) against +0.088 for NDW -- a DW-minus-NDW contrast of
  -0.850. J. Climate 32, 85 (2019) uses the same definition (850/100 hPa) and
  reports DW tropospheric anomalies "of around twice that of the total" composite.

WHAT IS COMPUTED HERE
  The same criterion on the frozen 43-event catalogue and this project's indices,
  then the decomposition established in FINDING_stratifier_law.md:

      DeltaY_observed  =  beta x DeltaS  +  DeltaY_residual

  beta is the OLS slope of the outcome window on the SELECTION window, fitted
  ONLY on event-free pseudo-onsets. DeltaS is the separation in the selection
  variable that the criterion itself creates. The residual is what survives.

  The null (pseudo-onsets on event-cleaned days, true differential effect exactly
  zero) gives the p-value and the bound.

WHAT THIS IS NOT
  This does NOT reproduce ACP 2026's -0.850. They use ERA5 850 hPa NAM and a
  1000 hPa-height NAO; this uses the CPC AO as the 1000 hPa NAM index and the CPC
  NAO, on a different catalogue. The published value is quoted for scale only.
  The claim made here is about the DECOMPOSITION of a contrast produced by this
  criterion, not about reproducing any particular paper's digits.

Output: recompute_published_criterion.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "9_literature"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]
SIM = ROOT / "ssw-design-analysis" / "06_simulation_validation"
sys.path.insert(0, str(SIM))
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
from build_catalogue import load_catalogue          # noqa: E402
import gate3_clean_null as G                        # noqa: E402
import multi_index_event_study as M                 # noqa: E402
import selection_on_outcome as SO                   # noqa: E402

SEL_WIN = (8, 52)          # the classification window, as published
OUT_WIN = (1, 60)          # ACP 2026's reporting window
N_BETA = 4000
N_NULL = 3000
SEED = 20260801

PUBLISHED = {"ACP_26_3723_2026_NAO_DW_BOTH": -0.762,
             "ACP_26_3723_2026_NAO_NDW": 0.088,
             "ACP_26_3723_2026_contrast": -0.850}


def strat_nam_150():
    s = pd.read_parquet(ROOT / "data/processed/atmospheric/ncep_stratosphere.parquet")
    if s.index.tz is not None:
        s.index = s.index.tz_convert("UTC").tz_localize(None)
    z = s["hgt_m_100hPa"]          # nearest available level to the published 150 hPa
    d = z.index.dayofyear
    return -((z - z.groupby(d).transform("mean")) / z.groupby(d).transform("std"))


def win_mean(series, onsets, win, need=20):
    out = []
    for o in onsets:
        lo, hi = o + pd.Timedelta(days=win[0]), o + pd.Timedelta(days=win[1])
        v = series[(series.index >= lo) & (series.index <= hi)]
        out.append(float(v.mean()) if len(v) >= need else np.nan)
    return np.array(out)


def classify(onsets, ao, snam):
    """Karpechko et al. (2017), conditions 1-3, exactly as published."""
    lab = []
    for o in onsets:
        lo, hi = o + pd.Timedelta(days=SEL_WIN[0]), o + pd.Timedelta(days=SEL_WIN[1])
        a = ao[(ao.index >= lo) & (ao.index <= hi)]
        if len(a) < 20:
            lab.append(False)
            continue
        c1 = a.mean() < 0
        c2 = (a < 0).mean() > 0.5
        n = snam[(snam.index >= lo) & (snam.index <= hi)]
        c3 = ((n < 0).mean() > 0.7) if len(n) >= 20 else None
        lab.append(bool(c1 and c2 and (c3 is not False)))
    return np.array(lab)


def main():
    rng = np.random.default_rng(SEED)
    ao = M.load("ao")["y"]
    snam = strat_nam_150()
    outs = {"nao": M.load("nao")["y"], "ao": ao, "pna": M.load("pna")["y"]}

    real = load_catalogue("primary")
    real = real[(real >= ao.index.min()) & (real <= ao.index.max())]
    mask = G.real_influence_mask(ao.index, real)
    clean_idx = ao[~mask].index
    doys = np.array([t.dayofyear for t in pd.DatetimeIndex(real)])
    clims = {k: v[~G.real_influence_mask(v.index, real)].groupby(
        v[~G.real_influence_mask(v.index, real)].index.dayofyear).mean()
        for k, v in outs.items()}

    lab = classify(real, ao, snam)
    print(f"criterion applied to {len(real)} events: "
          f"DW {int(lab.sum())}, NDW {int((~lab).sum())} "
          f"({100 * lab.mean():.0f}% downward)")

    # selection variable = the 1000 hPa NAM over the classification window
    Sr = win_mean(ao, real, SEL_WIN)

    # ---- beta on event-free windows --------------------------------------
    print(f"fitting beta on ~{N_BETA} event-free pseudo-onsets ...")
    Yp = {k: [] for k in outs}
    Sp = []
    drawn = 0
    while drawn < N_BETA:
        f = G.draw_clean(clean_idx, doys, rng)
        if len(f) < 10:
            continue
        Sp.append(win_mean(ao, f, SEL_WIN))
        for k, v in outs.items():
            Yp[k].append(SO.per_event(v, f, clims[k]))
        drawn += len(f)
    Sp = np.concatenate(Sp)
    Yp = {k: np.concatenate(v) for k, v in Yp.items()}
    beta = {}
    for k in outs:
        ok = np.isfinite(Yp[k]) & np.isfinite(Sp)
        beta[k] = float(np.cov(Yp[k][ok], Sp[ok])[0, 1] / np.var(Sp[ok], ddof=1))

    # ---- null: same criterion on pseudo-onsets ---------------------------
    print(f"null: applying the criterion to {N_NULL} pseudo-onset sets ...")
    null = {k: [] for k in outs}
    for i in range(N_NULL):
        f = G.draw_clean(clean_idx, doys, rng)
        if len(f) < 10:
            continue
        fl = classify(f, ao, snam)
        if fl.sum() < 3 or (~fl).sum() < 3:
            continue
        Sf = win_mean(ao, f, SEL_WIN)
        for k, v in outs.items():
            Yf = SO.per_event(v, f, clims[k])
            ok = np.isfinite(Yf) & np.isfinite(Sf)
            if (fl & ok).sum() < 3 or ((~fl) & ok).sum() < 3:
                continue
            dY = float(Yf[fl & ok].mean() - Yf[(~fl) & ok].mean())
            dS = float(Sf[fl & ok].mean() - Sf[(~fl) & ok].mean())
            null[k].append((dY, dS, dY - beta[k] * dS))
        if (i + 1) % 1000 == 0:
            print(f"  {i + 1}/{N_NULL}", flush=True)

    res = {"criterion": "Karpechko et al. 2017 conditions 1-3, verified from source",
           "selection_window": list(SEL_WIN), "outcome_window": list(OUT_WIN),
           "n_events": int(len(real)), "n_DW": int(lab.sum()),
           "rate_DW": round(float(lab.mean()), 3),
           "published_for_scale": PUBLISHED, "outcomes": {}}

    print(f"\n{'outcome':6s} {'DW':>8s} {'NDW':>8s} {'contrast':>9s} {'beta':>7s} "
          f"{'beta*dS':>8s} {'residual':>9s} {'p':>7s} {'%selection':>11s}")
    print("-" * 82)
    for k, v in outs.items():
        Y = SO.per_event(v, real, clims[k])
        ok = np.isfinite(Y) & np.isfinite(Sr)
        dw, ndw = float(Y[lab & ok].mean()), float(Y[(~lab) & ok].mean())
        dY = dw - ndw
        dS = float(Sr[lab & ok].mean() - Sr[(~lab) & ok].mean())
        pb = beta[k] * dS
        resid = dY - pb
        a = np.array(null[k])
        p = float((np.abs(a[:, 2]) >= abs(resid)).mean()) if len(a) else np.nan
        pct = 100 * pb / dY if abs(dY) > 1e-9 else np.nan
        res["outcomes"][k] = {
            "DW_composite": round(dw, 4), "NDW_composite": round(ndw, 4),
            "contrast": round(dY, 4), "dS": round(dS, 4),
            "beta": round(beta[k], 4), "selection_term": round(float(pb), 4),
            "residual": round(float(resid), 4), "p_residual": p,
            "pct_contrast_from_selection": round(float(pct), 1),
            "null_mean_contrast": round(float(a[:, 0].mean()), 4) if len(a) else None,
            "null_residual_sd": round(float(a[:, 2].std()), 4) if len(a) else None}
        print(f"{k:6s} {dw:+8.3f} {ndw:+8.3f} {dY:+9.3f} {beta[k]:+7.3f} "
              f"{pb:+8.3f} {resid:+9.3f} {p:7.4f} {pct:10.0f}%")

    print("\n=== THE NULL: same criterion, pseudo-events, TRUE contrast = 0 ===")
    for k in outs:
        a = np.array(null[k])
        r = res["outcomes"][k]
        print(f"  {k:6s} null DW-minus-NDW contrast = {a[:, 0].mean():+.3f}  "
              f"(observed {r['contrast']:+.3f}; "
              f"{100 * a[:, 0].mean() / r['contrast']:.0f}% reproduced with no effect)")

    nao = res["outcomes"]["nao"]
    res["verdict"] = (
        f"The published criterion produces a DW-minus-NDW NAO contrast of "
        f"{nao['contrast']:+.3f}; {nao['pct_contrast_from_selection']:.0f}% of it is "
        f"the selection term beta x dS. Residual {nao['residual']:+.3f} "
        f"(p={nao['p_residual']:.3f}). For scale, ACP 26, 3723 (2026) publishes "
        f"{PUBLISHED['ACP_26_3723_2026_contrast']:+.3f} from the same criterion.")
    print(f"\n=== VERDICT ===\n  {res['verdict']}")

    (RESULTS / "recompute_published_criterion.json").write_text(
        json.dumps(res, indent=2), encoding="utf8")
    print("\nSaved -> recompute_published_criterion.json")


if __name__ == "__main__":
    main()
