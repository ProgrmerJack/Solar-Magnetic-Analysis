#!/usr/bin/env python3
"""
consensus_monotonicity.py
=========================
Strengthens the marginal-event finding two ways:

  1. REPLICATION on NAO as well as AO. One index is an observation; two make a
     pattern.

  2. MONOTONICITY in detection count. The two-way split (>=4/6 vs 1-3/6) could be
     an artefact of where the line was drawn. If the pre-onset anomaly weakens
     steadily as more reanalyses agree an event occurred, that is far stronger
     evidence than any single threshold.

     Splitting 8 marginal events further would leave subsets of 3, which cannot
     support an event study. So monotonicity is tested per event instead: each
     event gets its own pre- and post-onset anomaly, measured against its OWN
     winter's baseline, and those are regressed on detection count. That uses all
     45 union events at full resolution.

PER-EVENT ANOMALIES
  pre  = mean index over -30..-16 d,  minus that winter's baseline (>75 d from
         onset), both day-of-year adjusted against the non-event seasonal cycle
  post = same for +15..+29 d

  A classification failure predicts: |pre| falls with detection count, |post|
  rises with it. Noise predicts neither.

Inference: Spearman rank correlation (detection count is ordinal and the
relationship need not be linear) plus a permutation test on the count labels.

Output: consensus_monotonicity.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "7_ensemble"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "06_simulation_validation"))
import calibrate_seasonality as C                   # noqa: E402

OUT = RESULTS / "consensus_monotonicity.json"
CAT = HERE.parents[0] / "02_event_catalogues" / "event_catalogue.csv"
IDX = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "daily_indices.parquet"
STRAT = ROOT / "data" / "processed" / "atmospheric" / "ncep_stratosphere.parquet"

SEASON = (11, 12, 1, 2, 3, 4)
PRE = (-30, -16)
POST = (15, 29)
BASELINE_GAP = 75
N_PERM = 10000


def winter_of(idx):
    return np.where(idx.month >= 11, idx.year + 1, idx.year)


def load_index(name):
    t = pd.read_parquet(IDX)
    s = t[t["index"] == name].set_index("date")["value"].sort_index()
    d = pd.DataFrame({"y": s.values}, index=s.index)
    d = d[np.isin(d.index.month, SEASON)]
    d["winter"] = winter_of(d.index)
    d["doy"] = d.index.dayofyear
    return d


def per_event(d, onsets, ssw_winters):
    """Pre/post anomalies vs own-winter baseline, DOY-adjusted."""
    nonssw = d[~d["winter"].isin(ssw_winters)]
    clim = nonssw.groupby("doy")["y"].mean()
    out = []
    for o in onsets:
        w = int(winter_of(pd.DatetimeIndex([o]))[0])
        own = d[d["winter"] == w]
        if len(own) < 60:
            out.append((np.nan, np.nan))
            continue
        lag = (own.index - o).days.values
        base = own[np.abs(lag) > BASELINE_GAP]
        if len(base) < 15:
            out.append((np.nan, np.nan))
            continue
        badj = base["y"].mean() - clim.reindex(base["doy"]).mean()
        vals = []
        for lo, hi in (PRE, POST):
            seg = own[(lag >= lo) & (lag <= hi)]
            if len(seg) < 6:
                vals.append(np.nan)
                continue
            vals.append(float(seg["y"].mean() - clim.reindex(seg["doy"]).mean() - badj))
        out.append(tuple(vals))
    return out


def main():
    cat = pd.read_csv(CAT)
    cat = cat[(~cat["final_warming_flag"]) & (cat["consensus_count"].fillna(0) >= 1)].copy()
    cat["n_detect"] = cat["consensus_count"].astype(float)
    cat["date"] = pd.to_datetime(cat["primary_central_date"])
    print(f"union events with a detection count: {len(cat)}")
    print("detection-count distribution:",
          cat["n_detect"].value_counts().sort_index().to_dict())

    st = pd.read_parquet(STRAT)
    if st.index.tz is not None:
        st.index = st.index.tz_convert("UTC").tz_localize(None)

    res = {"n_events": int(len(cat)), "pre_window": list(PRE),
           "post_window": list(POST), "n_perm": N_PERM, "indices": {}}

    for name in ("ao", "nao"):
        d = load_index(name)
        sub = cat[(cat["date"] >= d.index.min()) & (cat["date"] <= d.index.max())].copy()
        ssw_w = set(winter_of(pd.DatetimeIndex(sub["date"])))
        vals = per_event(d, pd.DatetimeIndex(sub["date"]), ssw_w)
        sub[f"pre"] = [v[0] for v in vals]
        sub[f"post"] = [v[1] for v in vals]
        s = sub.dropna(subset=["pre", "post"])
        print(f"\n=== {name.upper()}  ({len(s)} events with usable baselines) ===")

        entry = {"n_usable": int(len(s)), "by_count": {}, "monotonicity": {}}
        for k in sorted(s["n_detect"].unique()):
            g = s[s["n_detect"] == k]
            entry["by_count"][int(k)] = {
                "n": int(len(g)),
                "pre_mean": round(float(g["pre"].mean()), 3),
                "post_mean": round(float(g["post"].mean()), 3)}
            print(f"  detected by {int(k)}/6: n={len(g):2d}  "
                  f"pre {g['pre'].mean():+.3f}   post {g['post'].mean():+.3f}")

        for lab, col in (("pre", "pre"), ("post", "post")):
            rho, p_asym = stats.spearmanr(s["n_detect"], s[col])
            rng = np.random.default_rng(3)
            null = np.array([stats.spearmanr(rng.permutation(s["n_detect"].values),
                                             s[col].values)[0]
                             for _ in range(N_PERM)])
            p_perm = float((np.abs(null) >= abs(rho)).mean())
            entry["monotonicity"][lab] = {
                "spearman_rho": round(float(rho), 4),
                "p_asymptotic": float(p_asym),
                "p_permutation": p_perm}
            print(f"  Spearman(n_detect, {lab:4s}) rho={rho:+.3f}  "
                  f"perm P={p_perm:.4f}{'  *' if p_perm < 0.05 else ''}")

        # vortex wind vs detection count: is the classification story monotone too?
        u = []
        for o in s["date"]:
            w = st.loc[(st.index >= o - pd.Timedelta(days=5))
                       & (st.index <= o + pd.Timedelta(days=5)), "uwnd_ms_10hPa"]
            u.append(float(w.mean()) if len(w) else np.nan)
        s = s.assign(u10=u).dropna(subset=["u10"])
        rho_u, _ = stats.spearmanr(s["n_detect"], s["u10"])
        rng = np.random.default_rng(5)
        null_u = np.array([stats.spearmanr(rng.permutation(s["n_detect"].values),
                                           s["u10"].values)[0] for _ in range(N_PERM)])
        entry["u10_vs_count"] = {
            "spearman_rho": round(float(rho_u), 4),
            "p_permutation": float((np.abs(null_u) >= abs(rho_u)).mean()),
            "by_count": {int(k): round(float(g["u10"].mean()), 2)
                         for k, g in s.groupby("n_detect")}}
        print(f"  Spearman(n_detect, u10) rho={rho_u:+.3f}  "
              f"perm P={entry['u10_vs_count']['p_permutation']:.4f}")
        print(f"    u10 by count: {entry['u10_vs_count']['by_count']}")
        entry["per_event"] = s[["event_name", "primary_central_date", "n_detect",
                                "pre", "post", "u10"]].round(3).to_dict("records")
        res["indices"][name] = entry

    OUT.write_text(json.dumps(res, indent=2), encoding="utf8")
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    main()
