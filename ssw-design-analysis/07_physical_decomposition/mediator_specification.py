#!/usr/bin/env python3
"""
mediator_specification.py
=========================
STRESS TEST on continuous_mediation.py before any claim is made from it.

That analysis found the stratosphere mediates only 14.2% of the
wave-driving-to-surface connection (CI -35% to +60%, p=0.50), with 86% direct.
Taken at face value that is the "common cause" reading of stratosphere-surface
coupling, which would be a significant claim.

WHY IT IS PROBABLY WRONG AS IT STANDS
  The mediator was the 10 hPa 60N zonal wind on a SINGLE DAY t. The vortex
  anomaly that actually couples to the surface persists for weeks and descends
  over ~10-20 days. A one-day snapshot is a noisy proxy for that sustained state.

  Measurement error in a MEDIATOR has a known and specific consequence: it
  attenuates the estimated indirect path toward zero and pushes the direct path
  away from zero. That is exactly the pattern reported. So the finding may be
  entirely an artefact of how M was defined, not a statement about the atmosphere.

WHAT THIS DOES
  Re-estimates the mediated share across mediator specifications that differ in
  how well they capture the sustained, descending anomaly:

    u10_day        10 hPa 60N wind, single day t          (the original)
    u10_0_10       same, mean over t..t+10
    u10_0_20       same, mean over t..t+20
    u10_0_30       same, mean over t..t+30
    z100_0_20      100 hPa polar-cap height, mean t..t+20  (lower-strat NAM)
    z10_0_20       10 hPa polar-cap height, mean t..t+20
    nam_composite  standardised mean of u10_0_20 and z100_0_20 -- two noisy
                   measures of one latent state, averaged to reduce error

  If the mediated share rises systematically as the mediator better captures the
  persistent anomaly, the original 14% was attenuation and the common-cause
  reading is NOT supported. If it stays near 14% across all of them, the finding
  survives its most obvious objection and becomes worth taking seriously.

  NOTE ON ORDERING: mediators averaged over t..t+k still precede the outcome
  window t+15..t+45 only partially for k>15. Specifications where the mediator
  window overlaps the outcome window are FLAGGED, because overlap induces
  mechanical correlation and inflates apparent mediation. Read those with
  suspicion; they are included to show the trend, not as estimates.

Output: mediator_specification.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "5_mechanism"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
sys.path.insert(0, str(HERE.parents[0] / "06_simulation_validation"))
import calibrate_seasonality as C                   # noqa: E402
import canonical_event_study as K                   # noqa: E402
import continuous_mediation as CM                   # noqa: E402

HF = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "heatflux_daily.parquet"
STRAT = ROOT / "data" / "processed" / "atmospheric" / "ncep_stratosphere.parquet"
X_LAG = (-15, -1)
Y_LAG = (15, 45)
N_BOOT = 1200
SEASON = (11, 12, 1, 2, 3)


def daily(s):
    return s.reindex(pd.date_range(s.index.min(), s.index.max(), freq="D"))


def build_all():
    ao = K.load_series("AO")["y"]
    hf = pd.read_parquet(HF).set_index("date")["vT_100hPa_45_75N"].sort_index()
    st = pd.read_parquet(STRAT)
    if st.index.tz is not None:
        st.index = st.index.tz_convert("UTC").tz_localize(None)

    u10 = daily(st["uwnd_ms_10hPa"].sort_index())
    z100 = daily(st["hgt_m_100hPa"].sort_index())
    z10 = daily(st["hgt_m_10hPa"].sort_index())
    aod = daily(ao)
    hfd = daily(hf)

    idx = ao.index.intersection(st.index).intersection(hf.index)
    idx = idx[np.isin(idx.month, SEASON)]
    d = pd.DataFrame(index=idx)
    d["winter"] = CM.winter_of(idx)
    d["doy"] = idx.dayofyear
    d["X"] = hfd.rolling(15, min_periods=10).mean().shift(1).reindex(idx).values
    fwd = aod[::-1].rolling(Y_LAG[1] - Y_LAG[0] + 1, min_periods=15).mean()[::-1]
    d["Y"] = fwd.shift(-Y_LAG[0]).reindex(idx).values

    def fwd_mean(s, k):
        r = s[::-1].rolling(k + 1, min_periods=max(3, k // 2))[::-1] if False else None
        return s[::-1].rolling(k + 1, min_periods=max(3, k // 2)).mean()[::-1].reindex(idx).values

    d["u10_day"] = u10.reindex(idx).values
    for k in (10, 20, 30):
        d[f"u10_0_{k}"] = fwd_mean(u10, k)
    # polar-cap height: sign-flip so that, like u10, HIGH = strong vortex
    d["z100_0_20"] = -fwd_mean(z100, 20)
    d["z10_0_20"] = -fwd_mean(z10, 20)
    return d


def med_share(d, mcol, seed=20260731):
    s = d[["winter", "doy", "X", "Y", mcol]].dropna().copy()
    for c in ("X", "Y", mcol):
        s[c] = (s[c] - s[c].mean()) / s[c].std()

    def f_tot(x):
        return CM.fit(x, ["X"])[0]

    def f_dir(x):
        return CM.fit(x, ["X", mcol])

    tot = float(f_tot(s))
    dm = f_dir(s)
    ind = tot - float(dm[0])
    rng = np.random.default_rng(seed)
    winters = np.array(sorted(s["winter"].unique()))
    idxw = {w: np.flatnonzero((s["winter"] == w).values) for w in winters}
    sh = []
    for _ in range(N_BOOT):
        pick = rng.choice(winters, len(winters), replace=True)
        sub = s.iloc[np.concatenate([idxw[w] for w in pick])].copy()
        sub["winter"] = np.concatenate(
            [np.full(len(idxw[w]), j) for j, w in enumerate(pick)])
        try:
            t = float(CM.fit(sub, ["X"])[0])
            dd = CM.fit(sub, ["X", mcol])
            if np.isfinite(t) and t != 0 and np.all(np.isfinite(dd)):
                sh.append((t - float(dd[0])) / t)
        except Exception:
            pass
    sh = np.array(sh)
    return {"total": round(tot, 4), "direct": round(float(dm[0]), 4),
            "indirect": round(ind, 4),
            "mediator_coef": round(float(dm[1]), 4),
            "mediated_share": round(ind / tot, 4) if tot else None,
            "share_CI95": [round(float(np.percentile(sh, 2.5)), 4),
                           round(float(np.percentile(sh, 97.5)), 4)],
            "n": int(len(s)), "corr_M_Y": round(float(s[mcol].corr(s["Y"])), 4),
            "corr_X_M": round(float(s[mcol].corr(s["X"])), 4)}


def main():
    d = build_all()
    specs = [("u10_day", False), ("u10_0_10", False), ("u10_0_20", True),
             ("u10_0_30", True), ("z100_0_20", True), ("z10_0_20", True)]
    # composite of two noisy measures of the same latent state
    z = d[["u10_0_20", "z100_0_20"]].copy()
    d["nam_composite"] = ((z["u10_0_20"] - z["u10_0_20"].mean()) / z["u10_0_20"].std()
                          + (z["z100_0_20"] - z["z100_0_20"].mean()) / z["z100_0_20"].std()) / 2
    specs.append(("nam_composite", True))

    print(f"{'mediator':15s} {'n':>6s} {'r(M,Y)':>8s} {'r(X,M)':>8s} {'total':>8s} "
          f"{'indirect':>9s} {'mediated share':>16s}  overlap")
    print("-" * 88)
    res = {"X_lag": list(X_LAG), "Y_lag": list(Y_LAG), "n_boot": N_BOOT, "specs": {}}
    for m, overlap in specs:
        r = med_share(d, m)
        r["mediator_window_overlaps_outcome"] = overlap
        res["specs"][m] = r
        ci = r["share_CI95"]
        print(f"{m:15s} {r['n']:6d} {r['corr_M_Y']:+8.3f} {r['corr_X_M']:+8.3f} "
              f"{r['total']:+8.3f} {r['indirect']:+9.3f} "
              f"{100*r['mediated_share']:7.1f}% [{100*ci[0]:+6.0f},{100*ci[1]:+6.0f}]"
              f"   {'YES' if overlap else 'no'}", flush=True)

    clean = [(m, res["specs"][m]) for m, o in specs if not o]
    print("\n  Non-overlapping specifications only (the interpretable ones):")
    for m, r in clean:
        print(f"    {m:15s} mediated share {100*r['mediated_share']:.1f}%")
    print("\n  If the share climbs steeply as the mediator captures more of the")
    print("  persistent anomaly, the original 14% was attenuation, not physics.")
    (RESULTS / "mediator_specification.json").write_text(json.dumps(res, indent=2),
                                                      encoding="utf8")
    print("\nSaved -> mediator_specification.json")


if __name__ == "__main__":
    main()
