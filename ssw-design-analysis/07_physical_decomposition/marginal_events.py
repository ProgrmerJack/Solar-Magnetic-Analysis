#!/usr/bin/env python3
"""
marginal_events.py
==================
Localises the catalogue effect: WHICH events create the difference between the
union and consensus catalogues, and what is different about them?

BACKGROUND
  Catalogue choice moves the AO estimate 18-47x more than estimator choice.
  Under the union catalogue the pre-onset AO anomaly is significant; under the
  preregistered consensus catalogue it is not. Something about the events the
  union adds produces an apparent precursor.

FIRST HYPOTHESIS, TESTED AND REJECTED
  That reanalyses disagree on central dates, so mis-dated events leak their true
  response into the pre-onset window. Measured directly: the spread of central
  dates across the six reanalyses is **0 days for 23 of 37 events and <=2 days
  for 34**. Dating is not the explanation.

SECOND HYPOTHESIS, TESTED HERE
  The union adds events detected by only 1-3 reanalyses. If those marginal
  detections are weak or ambiguous vortex disturbances rather than clean
  warmings, they will carry a different -- and possibly mistimed -- surface
  signature, and including them changes the composite.

  This is a decomposition, not a screen: BOTH subsets are reported in full
  whatever they show, and the marginal set is small (n<=8), so its intervals
  will be wide and are quoted as such.

Output: marginal_events.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "7_ensemble"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "06_simulation_validation"))
import calibrate_seasonality as C                   # noqa: E402

OUT = RESULTS / "marginal_events.json"
CAT = HERE.parents[0] / "02_event_catalogues" / "event_catalogue.csv"
IDX = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "daily_indices.parquet"
STRAT = ROOT / "data" / "processed" / "atmospheric" / "ncep_stratosphere.parquet"
HF = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "heatflux_daily.parquet"

BINS = [(-60, -46), (-45, -31), (-30, -16), (-15, -1),
        (0, 14), (15, 29), (30, 44), (45, 60)]
LABELS = [f"{a:+d}..{b:+d}" for a, b in BINS]
BASELINE_GAP = 75
N_HARM = 3
N_BOOT = 1200
SEASON = (11, 12, 1, 2, 3, 4)
RE = ["ncep_ncar", "era40", "era_interim", "jra_55", "merra2", "era5"]


def winter_of(idx):
    return np.where(idx.month >= 11, idx.year + 1, idx.year)


def load_ao():
    t = pd.read_parquet(IDX)
    s = t[t["index"] == "ao"].set_index("date")["value"].sort_index()
    d = pd.DataFrame({"y": s.values}, index=s.index)
    d = d[np.isin(d.index.month, SEASON)]
    d["winter"] = winter_of(d.index)
    d["doy"] = d.index.dayofyear
    return d


def build(d, onsets):
    on = np.sort(np.array([np.datetime64(pd.Timestamp(o), "D") for o in onsets]))
    dd = d.index.values.astype("datetime64[D]")
    pos = np.searchsorted(on, dd)
    lo = np.clip(pos - 1, 0, len(on) - 1)
    hi = np.clip(pos, 0, len(on) - 1)
    dlo = (dd - on[lo]).astype(int)
    dhi = (dd - on[hi]).astype(int)
    lag = np.where(np.abs(dlo) <= np.abs(dhi), dlo, dhi)
    code = np.full(len(dd), -1)
    for i, (a, b) in enumerate(BINS):
        code[(lag >= a) & (lag <= b)] = i
    o = d.copy()
    o["bin"] = code
    o["is_baseline"] = np.abs(lag) > BASELINE_GAP
    return o[(o["bin"] >= 0) | o["is_baseline"]].copy()


def demean(A, codes):
    _, c = np.unique(codes, return_inverse=True)
    n = c.max() + 1
    cnt = np.bincount(c, minlength=n).astype(float)
    out = np.empty_like(A)
    for j in range(A.shape[1]):
        s = np.bincount(c, weights=A[:, j], minlength=n)
        out[:, j] = A[:, j] - (s / cnt)[c]
    return out


def fit(dd):
    y = dd["y"].values.astype(float)
    S = C.harmonics(dd["doy"].values, N_HARM)
    D = np.column_stack([(dd["bin"].values == i).astype(float)
                         for i in range(len(BINS))])
    Ad = demean(np.column_stack([y, D, S]), dd["winter"].values)
    try:
        b, *_ = np.linalg.lstsq(Ad[:, 1:], Ad[:, 0], rcond=None)
        return b[:len(BINS)]
    except Exception:
        return None


def run(d, onsets, seed=0):
    dd = build(d, onsets)
    est = fit(dd)
    if est is None:
        return None
    rng = np.random.default_rng(seed)
    winters = np.array(sorted(dd["winter"].unique()))
    idx = {w: np.flatnonzero((dd["winter"] == w).values) for w in winters}
    boot = []
    for _ in range(N_BOOT):
        pick = rng.choice(winters, len(winters), replace=True)
        rows = np.concatenate([idx[w] for w in pick])
        sub = dd.iloc[rows].copy()
        sub["winter"] = np.concatenate(
            [np.full(len(idx[w]), j) for j, w in enumerate(pick)])
        b = fit(sub)
        if b is not None and np.all(np.isfinite(b)):
            boot.append(b)
    boot = np.array(boot)
    return {lab: {"effect": round(float(est[i]), 4),
                  "CI95": [round(float(np.percentile(boot[:, i], 2.5)), 4),
                           round(float(np.percentile(boot[:, i], 97.5)), 4)],
                  "p_two_sided": float(2 * min((boot[:, i] >= 0).mean(),
                                               (boot[:, i] <= 0).mean()))}
            for i, lab in enumerate(LABELS)}


def main():
    cat = pd.read_csv(CAT)
    cat = cat[~cat["final_warming_flag"]]
    cat["n_detect"] = cat["consensus_count"]
    strong = cat[cat["n_detect"] >= 4]
    marginal = cat[(cat["n_detect"] >= 1) & (cat["n_detect"] <= 3)]
    print(f"strong (>=4/6 reanalyses): {len(strong)} events")
    print(f"marginal (1-3/6):          {len(marginal)} events")
    print("\nmarginal events:")
    print(marginal[["event_name", "primary_central_date", "n_detect"]]
          .to_string(index=False))

    # what is physically different about the marginal set?
    st = pd.read_parquet(STRAT)
    if st.index.tz is not None:
        st.index = st.index.tz_convert("UTC").tz_localize(None)
    hf = pd.read_parquet(HF).set_index("date")["vT_100hPa_45_75N"].sort_index()

    def props(dates):
        u, v = [], []
        for o in pd.to_datetime(dates):
            w = st.loc[(st.index >= o - pd.Timedelta(days=5))
                       & (st.index <= o + pd.Timedelta(days=5)), "uwnd_ms_10hPa"]
            p = hf[(hf.index >= o - pd.Timedelta(days=45)) & (hf.index < o)]
            if len(w):
                u.append(float(w.mean()))
            if len(p):
                v.append(float(p.mean()))
        return np.array(u), np.array(v)

    us, vs = props(strong["primary_central_date"])
    um, vm = props(marginal["primary_central_date"])
    from scipy import stats as sst
    phys = {
        "u10_strong_mean": round(float(us.mean()), 2) if len(us) else None,
        "u10_marginal_mean": round(float(um.mean()), 2) if len(um) else None,
        "u10_mannwhitney_p": float(sst.mannwhitneyu(us, um).pvalue) if len(um) else None,
        "vT_strong_mean": round(float(vs.mean()), 2) if len(vs) else None,
        "vT_marginal_mean": round(float(vm.mean()), 2) if len(vm) else None,
        "vT_mannwhitney_p": float(sst.mannwhitneyu(vs, vm).pvalue) if len(vm) else None,
    }
    print(f"\nvortex wind at onset  strong {phys['u10_strong_mean']} vs "
          f"marginal {phys['u10_marginal_mean']} m/s  "
          f"(P={phys['u10_mannwhitney_p']:.3f})")
    print(f"precursor v'T'        strong {phys['vT_strong_mean']} vs "
          f"marginal {phys['vT_marginal_mean']} K m/s  "
          f"(P={phys['vT_mannwhitney_p']:.3f})")

    ao = load_ao()
    res = {"n_strong": int(len(strong)), "n_marginal": int(len(marginal)),
           "marginal_events": marginal[["event_name", "primary_central_date",
                                        "n_detect"]].to_dict("records"),
           "physical_contrast": phys, "profiles": {}}

    for name, sub in (("strong_only", strong), ("marginal_only", marginal),
                      ("union_all", cat[cat["n_detect"] >= 1])):
        on = pd.DatetimeIndex(pd.to_datetime(sub["primary_central_date"]))
        on = on[(on >= ao.index.min()) & (on <= ao.index.max())]
        if len(on) < 3:
            continue
        # crc32, not hash(): see canonical_event_study.py
        prof = run(ao, on, seed=zlib.crc32(name.encode()) % 9999)
        if prof is None:
            continue
        res["profiles"][name] = {"n_events": len(on), "profile": prof}
        star = "".join("*" if prof[l]["p_two_sided"] < 0.05 else "." for l in LABELS)
        vals = " ".join(f"{prof[l]['effect']:+.2f}" for l in LABELS)
        print(f"  {name:14s} n={len(on):2d}  {vals}   {star}")

    OUT.write_text(json.dumps(res, indent=2), encoding="utf8")
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    main()
