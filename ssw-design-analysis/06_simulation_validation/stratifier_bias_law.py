#!/usr/bin/env python3
"""
stratifier_bias_law.py
======================
THE GENERALISATION, AND THE TEST THAT CAN KILL IT.

WHERE THIS COMES FROM
  selection_on_outcome.py showed that classifying SSWs with a criterion that IS
  the surface response manufactures 30-44% of the reported anomaly under a true
  null. selection_projection_law.py derived the size exactly:

      bias(Y) = beta(Y_window on A_window) x bias(A)                        (1)

  with beta the ORDINARY LEAST-SQUARES SLOPE (R^2 = 0.999 across 5 outcomes,
  0.985 on an unrelated selection variable).

  Both scripts treat A as the OUTCOME ITSELF. That is the maximal case, beta = 1.
  But nothing in (1) requires it. beta is defined for ANY stratifier, and the
  literature audit's "S2 = not outcome-classified" stratum was built on the
  assumption that a stratospheric stratifier carries NO bias.

  That assumption is wrong, and it is wrong for a physical reason. A stratospheric
  stratifier is not independent of the surface -- correlation between the vortex
  and the annular mode is the entire subject of the field. So splitting events on
  v'T', on vortex strength, or on the depth of the reversal ALSO manufactures a
  between-group surface contrast, just with a smaller beta.

  What (1) then buys is not a warning but a CORRECTION:

      DeltaY_observed  =  beta x DeltaS  +  DeltaY_causal                   (2)

  measure the first two, subtract, and what remains is the differential surface
  effect actually attributable to the event. That number has never been reported
  for any SSW stratification because the beta x DeltaS term has never been removed.

WHAT IS TESTED HERE, ON OBSERVATIONS, NOT IN SIMULATION
  For each stratifier S, at the frozen 43-event catalogue:
    DeltaY_obs   observed surface contrast, top-half minus bottom-half by S
    DeltaS       the stratifier separation the split creates
    beta         OLS slope of Y_window on S_window, estimated ONLY on event-free
                 pseudo-onsets, so it carries no event effect by construction
    predicted    beta x DeltaS
    residual     DeltaY_obs - beta x DeltaS, the unbiased differential effect

  THE FALSIFICATION. On pseudo-onsets drawn from event-cleaned data the true
  differential effect is exactly zero, so (2) demands

      DeltaY_null  ==  beta x DeltaS_null      for EVERY stratifier.

  If the null residual is not zero, equation (1) does not hold outside the
  simulation that produced it and this whole line closes. The test is run for
  every stratifier and reported whether it passes or fails.

  SUPPORT CHECK. beta is fitted on event-free windows and applied to SSW-disturbed
  ones. For the post-onset stratospheric stratifiers DeltaS can fall outside the
  fitted range, making (2) an extrapolation. The overlap is measured and reported
  per stratifier; a stratifier with poor support gets its residual flagged, not
  quietly used.

Output: stratifier_bias_law.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "4_selection_bias"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
from build_catalogue import load_catalogue          # noqa: E402
import gate3_clean_null as G                        # noqa: E402
import multi_index_event_study as M                 # noqa: E402
import selection_on_outcome as S                    # noqa: E402

import os                                            # noqa: E402
# SSW_STRAT_1958=1 swaps in the pre-satellite-extended record (n=42 instead of
# n=29). Only legitimate once era_dependence_gate.py reports pooling_allowed.
STRAT = ROOT / "data" / "processed" / "atmospheric" / (
    "ncep_stratosphere_1958.parquet" if os.environ.get("SSW_STRAT_1958")
    else "ncep_stratosphere.parquet")
HF = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "heatflux_daily.parquet"

OUT_WIN = (8, 52)          # the published outcome window, kept identical
N_BETA = 3000              # pseudo-onsets used to fit beta on event-free data
N_NULL = 1500              # pseudo-onset draws for the falsification test
N_BOOT = 1200
SEED = 20260801


# --------------------------------------------------------------------------
# stratifiers.  "pre" = strictly before onset, usable in a forecast.
# "strat" = post-onset but STRATOSPHERIC ONLY, i.e. exactly the class the audit
# protocol called S2 and assumed unbiased.  "outcome" = the published criterion.
# --------------------------------------------------------------------------
STRATIFIERS = [
    ("vT_pre30",   "pre",     "vt",   (-30, -1),  "wave driving, 30 d before onset"),
    ("vT_pre15",   "pre",     "vt",   (-15, -1),  "wave driving, 15 d before onset"),
    ("u10_pre",    "pre",     "u10",  (-45, -16), "vortex preconditioning"),
    ("u10_post30", "strat",   "u10",  (0, 30),    "reversal depth, 10 hPa"),
    ("u10_post45", "strat",   "u10",  (0, 45),    "reversal persistence, 10 hPa"),
    ("z100_post",  "strat",   "z100", (0, 30),    "lower-stratospheric descent (PJO-like)"),
    ("AO_post",    "outcome", "ao",   (8, 52),    "the published surface criterion"),
]


def doy_standardise(x):
    d = x.index.dayofyear
    return (x - x.groupby(d).transform("mean")) / x.groupby(d).transform("std")


def load_predictors():
    s = pd.read_parquet(STRAT)
    if s.index.tz is not None:
        s.index = s.index.tz_convert("UTC").tz_localize(None)
    hf = pd.read_parquet(HF).set_index("date")["vT_100hPa_45_75N"]
    return {"u10": doy_standardise(s["uwnd_ms_10hPa"]),
            "z100": -doy_standardise(s["hgt_m_100hPa"]),   # sign: + = strong vortex
            "vt": doy_standardise(hf)}


def window_mean(series, onsets, win):
    """Plain window mean of an ALREADY day-of-year standardised series."""
    out = []
    for o in onsets:
        lo, hi = o + pd.Timedelta(days=win[0]), o + pd.Timedelta(days=win[1])
        v = series[(series.index >= lo) & (series.index <= hi)]
        out.append(float(v.mean()) if len(v) >= max(6, (win[1] - win[0]) // 3)
                   else np.nan)
    return np.array(out)


def contrast(y, s):
    """Top-half minus bottom-half of y, split at the median of s."""
    ok = np.isfinite(y) & np.isfinite(s)
    if ok.sum() < 10:
        return np.nan, np.nan, 0
    yy, ss = y[ok], s[ok]
    hi = ss >= np.median(ss)
    if hi.sum() < 3 or (~hi).sum() < 3:
        return np.nan, np.nan, int(ok.sum())
    return (float(yy[hi].mean() - yy[~hi].mean()),
            float(ss[hi].mean() - ss[~hi].mean()), int(ok.sum()))


def winter_of(idx):
    t = pd.DatetimeIndex(idx)
    return np.where(t.month >= 11, t.year + 1, t.year)


def main():
    rng = np.random.default_rng(SEED)
    ao = M.load("ao")["y"]
    pred = load_predictors()

    allev = load_catalogue("primary")   # EXCLUSION set: every event, even ones not scored here
    real = allev[(allev >= ao.index.min()) & (allev <= ao.index.max())]
    # every stratifier must be computable, so restrict to the common coverage
    lo_cov = max(v.index.min() for v in pred.values())
    hi_cov = min(v.index.max() for v in pred.values())
    real = real[(real >= lo_cov + pd.Timedelta(days=60))
                & (real <= hi_cov - pd.Timedelta(days=60))]
    print(f"events with full predictor coverage: {len(real)} "
          f"({real.min().date()} .. {real.max().date()})")

    # CLEAN day-of-year climatology for the OUTCOME -- the contamination trap
    # that has already cost this project three times.
    # Mask and pseudo pool exclude EVERY catalogued event, not only the scored
    # ones: passing the coverage-trimmed `real` left the 1958-78 events inside
    # the pool and the climatology (review 2026-09-25, 21% of NDJFM candidates).
    mask = G.real_influence_mask(ao.index, allev)
    clim = ao[~mask].groupby(ao[~mask].index.dayofyear).mean()
    clean_idx = G.zone_free_index(ao.index, allev)
    doys = np.array([t.dayofyear for t in pd.DatetimeIndex(real)])

    Y_real = S.per_event(ao, real, clim)
    print(f"all-event surface composite, AO {OUT_WIN}: {np.nanmean(Y_real):+.3f}")

    # one standardised source series per stratifier, built once
    pred["ao"] = doy_standardise(ao)
    NAMES = [s[0] for s in STRATIFIERS]
    SRC = {name: pred[var] for name, _, var, _, _ in STRATIFIERS}
    WINS = {name: win for name, _, _, win, _ in STRATIFIERS}

    # ---------------- beta, fitted ONLY on event-free windows ---------------
    print(f"\nfitting beta on {N_BETA} event-free pseudo-onsets ...")
    Yp, Sp = [], {k: [] for k in NAMES}
    drawn = 0
    while drawn < N_BETA:
        f = G.draw_clean(clean_idx, doys, rng)
        if len(f) < 10:
            continue
        Yp.append(S.per_event(ao, f, clim))
        for name in NAMES:
            Sp[name].append(window_mean(SRC[name], f, WINS[name]))
        drawn += len(f)
    Yp = np.concatenate(Yp)
    Sp = {k: np.concatenate(v) for k, v in Sp.items()}

    beta, support = {}, {}
    for name in NAMES:
        ok = np.isfinite(Yp) & np.isfinite(Sp[name])
        beta[name] = float(np.cov(Yp[ok], Sp[name][ok])[0, 1]
                           / np.var(Sp[name][ok], ddof=1))
        support[name] = (float(np.percentile(Sp[name][ok], 1)),
                         float(np.percentile(Sp[name][ok], 99)))

    # ---------------- observations -----------------------------------------
    print(f"\n{'stratifier':12s} {'class':8s} {'dY_obs':>8s} {'dS':>7s} "
          f"{'beta':>7s} {'beta*dS':>8s} {'residual':>9s} {'95% CI':>18s} {'supp':>6s}")
    print("-" * 96)
    res = {"n_events": int(len(real)), "outcome_window": list(OUT_WIN),
           "all_event_composite": round(float(np.nanmean(Y_real)), 4),
           "n_beta_pseudo": int(len(Yp)), "stratifiers": {}}

    wid = winter_of(real)
    for name, kind, var, win, desc in STRATIFIERS:
        Sr = window_mean(SRC[name], real, win)
        dY, dS, n = contrast(Y_real, Sr)
        pred_bias = beta[name] * dS
        resid = dY - pred_bias

        # winter-block bootstrap of the residual
        ok = np.isfinite(Y_real) & np.isfinite(Sr)
        uw = np.unique(wid[ok])
        bs = []
        for _ in range(N_BOOT):
            pick = rng.choice(uw, len(uw), replace=True)
            sel = np.concatenate([np.flatnonzero(wid == q) for q in pick])
            sel = sel[np.isfinite(Y_real[sel]) & np.isfinite(Sr[sel])]
            if len(sel) < 10:
                continue
            d1, d2, _ = contrast(Y_real[sel], Sr[sel])
            if np.isfinite(d1):
                bs.append(d1 - beta[name] * d2)
        ci = [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))] if bs else None

        lo, hi = support[name]
        inrange = float(np.mean((Sr[ok] >= lo) & (Sr[ok] <= hi)))
        flag = "" if inrange > 0.8 else "  <- weak support, extrapolated"

        res["stratifiers"][name] = {
            "class": kind, "variable": var, "window": list(win),
            "description": desc,
            "dY_observed": round(dY, 4), "dS": round(dS, 4),
            "beta": round(beta[name], 4),
            "predicted_bias": round(float(pred_bias), 4),
            "residual_causal": round(float(resid), 4),
            "residual_CI95": [round(ci[0], 3), round(ci[1], 3)] if ci else None,
            "frac_pct_of_dY_explained_by_stratifier":
                round(float(100 * pred_bias / dY), 1) if abs(dY) > 1e-9 else None,
            "support_fraction": round(inrange, 3),
            "n": n}
        cis = f"[{ci[0]:+6.2f},{ci[1]:+6.2f}]" if ci else " " * 18
        print(f"{name:12s} {kind:8s} {dY:+8.3f} {dS:+7.3f} {beta[name]:+7.3f} "
              f"{pred_bias:+8.3f} {resid:+9.3f} {cis} {inrange:5.2f}{flag}")

    # ---------------- THE FALSIFICATION TEST -------------------------------
    print(f"\n=== FALSIFICATION: on pseudo-onsets the residual MUST be zero ===")
    print(f"{'stratifier':12s} {'dY_null':>9s} {'beta*dS':>9s} {'residual':>9s} "
          f"{'95% CI':>18s}  verdict")
    print("-" * 78)
    null = {k: [] for k in NAMES}
    for i in range(N_NULL):
        f = G.draw_clean(clean_idx, doys, rng)
        if len(f) < 10:
            continue
        Yf = S.per_event(ao, f, clim)
        for name in NAMES:
            Sf = window_mean(SRC[name], f, WINS[name])
            d1, d2, _ = contrast(Yf, Sf)
            if np.isfinite(d1) and np.isfinite(d2):
                null[name].append((d1, d2, d1 - beta[name] * d2))
        if (i + 1) % 500 == 0:
            print(f"  {i + 1}/{N_NULL} null draws", flush=True)

    n_pass = 0
    for name in NAMES:
        a = np.array(null[name])
        if not len(a):
            continue
        dY_n, pb_n, r_n = a[:, 0].mean(), (beta[name] * a[:, 1]).mean(), a[:, 2]
        ci = [float(np.percentile(r_n, 2.5)), float(np.percentile(r_n, 97.5))]
        ok = ci[0] <= 0.0 <= ci[1]
        n_pass += ok
        res["stratifiers"][name]["null_dY"] = round(float(dY_n), 4)
        res["stratifiers"][name]["null_predicted"] = round(float(pb_n), 4)
        res["stratifiers"][name]["null_residual"] = round(float(r_n.mean()), 4)
        res["stratifiers"][name]["null_residual_CI95"] = [round(ci[0], 3), round(ci[1], 3)]
        res["stratifiers"][name]["law_holds_under_null"] = bool(ok)
        print(f"{name:12s} {dY_n:+9.3f} {pb_n:+9.3f} {r_n.mean():+9.3f} "
              f"[{ci[0]:+6.3f},{ci[1]:+6.3f}]  {'PASS' if ok else 'FAIL'}")

    res["n_stratifiers_passing_null"] = int(n_pass)
    res["n_stratifiers"] = len(STRATIFIERS)
    res["verdict"] = (
        "LAW HOLDS on observations: every stratifier's null residual contains zero, "
        "so beta x dS accounts for the whole non-causal contrast and the residual is "
        "an unbiased differential effect."
        if n_pass == len(STRATIFIERS) else
        f"LAW FAILS for {len(STRATIFIERS) - n_pass}/{len(STRATIFIERS)} stratifiers -- "
        "equation (1) does not generalise beyond the simulation that produced it.")
    print(f"\n=== VERDICT ===\n  {res['verdict']}")

    tag = "_1958" if os.environ.get("SSW_STRAT_1958") else ""
    res["stratosphere_record"] = STRAT.name
    (RESULTS / f"stratifier_bias_law{tag}.json").write_text(json.dumps(res, indent=2),
                                                        encoding="utf8")
    print(f"\nSaved -> stratifier_bias_law{tag}.json")


if __name__ == "__main__":
    main()
