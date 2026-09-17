#!/usr/bin/env python3
"""
design_sensitivity.py
=====================
GATE 8: is the difference between estimators predictable from dynamical
conditions, or is it noise?

This is the piece that turns a statistical correction into a physical result.
If the gap between the conventional and corrected estimates is systematic in
identifiable dynamics, the paper can say *when* a composite will mislead and
*why* -- not merely that it can.

PER-EVENT OUTCOME
  For each of the 39 frozen-catalogue events, the same window (+15..+29 d, the
  bin where the AO response peaks) is estimated two ways:

    conventional  window mean minus the day-of-year-matched mean taken from
                  winters containing NO event  (the field's standard composite)
    corrected     window mean minus that SAME winter's own baseline, restricted
                  to days >75 d from onset and day-of-year adjusted

    design_sensitivity = conventional - corrected

  A negative value means the conventional composite reports a more negative AO
  response than the within-winter comparison does.

PREDICTORS (all fixed before fitting)
  precursor_vT     100 hPa eddy heat flux 45-75N, mean over days -45..-1.
                   The wave driving that produced the warming.
  vortex_u10       10 hPa 60N zonal wind, mean over days -5..+5. Event strength.
  polarcap_z100    100 hPa height, mean over days 0..+30. Downward coupling.
  onset_doy        timing within the season.
  winter_AO_anom   that winter's mean AO minus the all-winter climatology.
                   THE selection variable: if design sensitivity tracks this,
                   the estimator gap is winter selection, measured directly
                   rather than inferred.

The last predictor is the point. Earlier rounds inferred winter selection from
a placebo failing; here it is entered as a regressor and its coefficient is
estimated.

Inference: ordinary least squares across events with heteroskedasticity-robust
errors, plus a permutation test (event labels shuffled) because 39 events and
5 predictors is a small regression and asymptotic p-values are optimistic.

Output: design_sensitivity.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "5_mechanism"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "06_simulation_validation"))
from build_catalogue import load_catalogue          # noqa: E402
import calibrate_seasonality as C                   # noqa: E402

OUT = RESULTS / "design_sensitivity.json"
IDX = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "daily_indices.parquet"
HF = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "heatflux_daily.parquet"
STRAT = ROOT / "data" / "processed" / "atmospheric" / "ncep_stratosphere.parquet"

SEASON = (11, 12, 1, 2, 3, 4)
WINDOW = (15, 29)
BASELINE_GAP = 75
N_PERM = 5000


def winter_of(idx):
    return np.where(idx.month >= 11, idx.year + 1, idx.year)


def load_ao():
    t = pd.read_parquet(IDX)
    s = t[t["index"] == "ao"].set_index("date")["value"].sort_index()
    d = pd.DataFrame({"ao": s.values}, index=s.index)
    d = d[np.isin(d.index.month, SEASON)]
    d["winter"] = winter_of(d.index)
    d["doy"] = d.index.dayofyear
    return d


def load_predictor_series():
    hf = pd.read_parquet(HF).set_index("date")["vT_100hPa_45_75N"].sort_index()
    st = pd.read_parquet(STRAT)
    if st.index.tz is not None:
        st.index = st.index.tz_convert("UTC").tz_localize(None)
    return hf, st


def main():
    ao = load_ao()
    hf, st = load_predictor_series()
    onsets = load_catalogue("primary")
    onsets = onsets[(onsets >= ao.index.min()) & (onsets <= ao.index.max())]

    ssw_winters = set(winter_of(pd.DatetimeIndex(onsets)))
    nonssw = ao[~ao["winter"].isin(ssw_winters)]
    ao_clim_by_winter = ao.groupby("winter")["ao"].mean()
    all_winter_mean = float(ao_clim_by_winter.mean())

    rows = []
    for o in onsets:
        w = int(winter_of(pd.DatetimeIndex([o]))[0])
        lo, hi = o + pd.Timedelta(days=WINDOW[0]), o + pd.Timedelta(days=WINDOW[1])
        win = ao[(ao.index >= lo) & (ao.index <= hi)]
        if len(win) < 8:
            continue
        doys = set(win["doy"])

        # conventional: DOY-matched days from winters with no event
        ctrl = nonssw[nonssw["doy"].isin(doys)]
        if len(ctrl) < 20:
            continue
        conventional = float(win["ao"].mean() - ctrl["ao"].mean())

        # corrected: same winter, days far from onset, DOY-adjusted using the
        # non-event seasonal cycle so the two comparisons see the same season
        own = ao[(ao["winter"] == w)]
        lag = (own.index - o).days.values
        base = own[np.abs(lag) > BASELINE_GAP]
        if len(base) < 20:
            continue
        clim = nonssw.groupby("doy")["ao"].mean()
        adj_win = win["ao"].mean() - clim.reindex(win["doy"]).mean()
        adj_base = base["ao"].mean() - clim.reindex(base["doy"]).mean()
        corrected = float(adj_win - adj_base)

        # ---- predictors ----
        pre = hf[(hf.index >= o - pd.Timedelta(days=45)) & (hf.index < o)]
        u10 = st.loc[(st.index >= o - pd.Timedelta(days=5))
                     & (st.index <= o + pd.Timedelta(days=5)), "uwnd_ms_10hPa"]
        z100 = st.loc[(st.index >= o) & (st.index <= o + pd.Timedelta(days=30)),
                      "hgt_m_100hPa"]
        rows.append({
            "onset": str(pd.Timestamp(o).date()), "winter": w,
            "conventional": round(conventional, 4),
            "corrected": round(corrected, 4),
            "design_sensitivity": round(conventional - corrected, 4),
            "precursor_vT": round(float(pre.mean()), 3) if len(pre) else np.nan,
            "vortex_u10": round(float(u10.mean()), 3) if len(u10) else np.nan,
            "polarcap_z100": round(float(z100.mean()), 2) if len(z100) else np.nan,
            "onset_doy": int(pd.Timestamp(o).dayofyear),
            "winter_AO_anom": round(float(ao_clim_by_winter.get(w, np.nan)
                                          - all_winter_mean), 4),
        })

    df = pd.DataFrame(rows).dropna()
    print(f"events usable: {len(df)} of {len(onsets)}")
    print(f"conventional mean {df['conventional'].mean():+.3f}, "
          f"corrected mean {df['corrected'].mean():+.3f}, "
          f"design sensitivity mean {df['design_sensitivity'].mean():+.3f}")

    preds = ["precursor_vT", "vortex_u10", "polarcap_z100",
             "onset_doy", "winter_AO_anom"]
    X = df[preds].values.astype(float)
    Xs = (X - X.mean(0)) / X.std(0)                 # standardise for comparability
    y = df["design_sensitivity"].values.astype(float)
    A = np.column_stack([np.ones(len(Xs)), Xs])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ beta
    # HC3-ish robust errors
    XtXi = np.linalg.pinv(A.T @ A)
    h = np.einsum("ij,jk,ik->i", A, XtXi, A)
    S = (A * (resid / (1 - np.clip(h, 0, .99)))[:, None])
    cov = XtXi @ (S.T @ S) @ XtXi
    se = np.sqrt(np.diag(cov))
    r2 = 1 - resid.var() / y.var()

    # permutation test on R^2 and each coefficient
    rng = np.random.default_rng(7)
    perm_r2, perm_b = [], []
    for _ in range(N_PERM):
        yp = rng.permutation(y)
        b, *_ = np.linalg.lstsq(A, yp, rcond=None)
        perm_r2.append(1 - (yp - A @ b).var() / yp.var())
        perm_b.append(b[1:])
    perm_b = np.array(perm_b)
    p_r2 = float((np.array(perm_r2) >= r2).mean())

    res = {"n_events": int(len(df)), "window": list(WINDOW),
           "mean_conventional": round(float(df["conventional"].mean()), 4),
           "mean_corrected": round(float(df["corrected"].mean()), 4),
           "mean_design_sensitivity": round(float(df["design_sensitivity"].mean()), 4),
           "r2": round(float(r2), 4), "r2_permutation_p": p_r2,
           "n_perm": N_PERM, "coefficients": {}, "per_event": rows}
    print(f"\nR2 = {r2:.3f}  (permutation P = {p_r2:.4f}, {N_PERM} shuffles)")
    print(f"{'predictor':16s} {'beta':>8s} {'robust SE':>10s} {'perm p':>8s}")
    print("-" * 48)
    for i, nm in enumerate(preds):
        b = float(beta[i + 1])
        pp = float((np.abs(perm_b[:, i]) >= abs(b)).mean())
        res["coefficients"][nm] = {"beta_standardised": round(b, 4),
                                   "robust_se": round(float(se[i + 1]), 4),
                                   "permutation_p": pp}
        print(f"{nm:16s} {b:+8.4f} {se[i+1]:10.4f} {pp:8.4f}"
              f"{'  *' if pp < 0.05 else ''}")

    # simple correlations, reported alongside so the regression is not the only view
    res["correlations_with_design_sensitivity"] = {
        nm: [round(float(stats.pearsonr(df[nm], y)[0]), 3),
             round(float(stats.pearsonr(df[nm], y)[1]), 4)] for nm in preds}
    print("\npairwise correlations with design sensitivity:")
    for nm, (r, p) in res["correlations_with_design_sensitivity"].items():
        print(f"  {nm:16s} r={r:+.3f}  p={p:.4f}")

    OUT.write_text(json.dumps(res, indent=2), encoding="utf8")
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    main()
