#!/usr/bin/env python3
"""
continuous_mediation.py
=======================
Abandons the event framing entirely, and asks the question the events were only
ever a proxy for.

THE ARGUMENT
  Every problem this project has documented -- catalogue definition, consensus
  rules, era dependence, compositing, selection on the outcome -- is a problem of
  DISCRETISING a continuous process. "Sudden stratospheric warming" is a
  threshold on the 10 hPa 60N zonal wind. The physics underneath is continuous:
  tropospheric wave activity propagates upward, decelerates the vortex, and the
  anomaly descends.

  If the continuous process can be analysed directly, every one of those design
  problems disappears -- there is no catalogue to choose, no events to select, no
  composite to build -- and the sample grows from 43 events to every winter day
  in the record.

THE CAUSAL CHAIN, with a defensible ordering
      X  wave driving      100 hPa eddy heat flux v'T', 45-75N, days t-15..t-1
      M  vortex state      10 hPa 60N zonal wind, day t
      Y  surface           AO, mean over days t+15..t+45

  The ordering is not assumed from correlation: v'T' at 100 hPa is a TROPOSPHERIC
  flux of wave activity into the stratosphere, and it necessarily precedes the
  stratospheric response it produces. That is why X is measured strictly before M,
  and Y strictly after.

THE QUESTION
  How much of the wave-driving-to-surface connection is MEDIATED by the
  stratosphere, and how much is DIRECT?

    total    = effect of X on Y
    direct   = effect of X on Y controlling for M
    indirect = total - direct   (the stratospherically mediated part)

  A large indirect share supports the stratosphere as a causal intermediary. A
  large direct share means the troposphere drives both and the stratosphere is
  substantially a bystander -- the "common cause" reading that the field has
  debated for two decades without a clean observational test.

ASSUMPTIONS, STATED PLAINLY
  Mediation identifies the indirect effect only if there is no unmeasured
  confounder of the M->Y path. That cannot be guaranteed observationally. A
  sensitivity analysis is therefore run: how strong would such a confounder have
  to be to overturn the conclusion?

  All estimates carry winter fixed effects and 3 annual harmonics, the frozen
  Gate 3 specification, with a winter-block bootstrap.

Output: continuous_mediation.json
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

HF = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "heatflux_daily.parquet"
STRAT = ROOT / "data" / "processed" / "atmospheric" / "ncep_stratosphere.parquet"

X_LAG = (-15, -1)        # wave driving, strictly before M
Y_LAG = (15, 45)         # surface, strictly after M
N_BOOT = 2000
SEASON = (11, 12, 1, 2, 3)


def winter_of(idx):
    return np.where(idx.month >= 11, idx.year + 1, idx.year)


def build():
    ao = K.load_series("AO")["y"]
    hf = pd.read_parquet(HF).set_index("date")["vT_100hPa_45_75N"].sort_index()
    st = pd.read_parquet(STRAT)
    if st.index.tz is not None:
        st.index = st.index.tz_convert("UTC").tz_localize(None)
    u10 = st["uwnd_ms_10hPa"].sort_index()

    idx = ao.index.intersection(hf.index).intersection(u10.index)
    idx = idx[np.isin(idx.month, SEASON)]
    d = pd.DataFrame(index=idx)
    d["winter"] = winter_of(idx)
    d["doy"] = idx.dayofyear

    # X: mean wave driving over t-15..t-1
    hfr = hf.reindex(pd.date_range(hf.index.min(), hf.index.max(), freq="D"))
    xs = hfr.rolling(-X_LAG[0] - (-X_LAG[1]) + 1, min_periods=10).mean().shift(1)
    d["X"] = xs.reindex(idx).values
    # M: vortex wind on day t
    d["M"] = u10.reindex(idx).values
    # Y: mean AO over t+15..t+45
    aor = ao.reindex(pd.date_range(ao.index.min(), ao.index.max(), freq="D"))
    fwd = aor[::-1].rolling(Y_LAG[1] - Y_LAG[0] + 1, min_periods=15).mean()[::-1]
    d["Y"] = fwd.shift(-Y_LAG[0]).reindex(idx).values
    d = d.dropna()
    # standardise so coefficients are comparable in sd units
    for c in ("X", "M", "Y"):
        d[c] = (d[c] - d[c].mean()) / d[c].std()
    return d


def fit(d, cols):
    """Within-winter OLS of Y on cols, with 3 annual harmonics."""
    y = d["Y"].values.astype(float)
    S = C.harmonics(d["doy"].values, 3)
    A = np.column_stack([y] + [d[c].values.astype(float) for c in cols] + [S])
    Ad = K.C_demean(A, d["winter"].values)
    b, *_ = np.linalg.lstsq(Ad[:, 1:], Ad[:, 0], rcond=None)
    return b[:len(cols)]


def boot(d, fn, seed=20260731):
    rng = np.random.default_rng(seed)
    winters = np.array(sorted(d["winter"].unique()))
    idx = {w: np.flatnonzero((d["winter"] == w).values) for w in winters}
    out = []
    for _ in range(N_BOOT):
        pick = rng.choice(winters, len(winters), replace=True)
        sub = d.iloc[np.concatenate([idx[w] for w in pick])].copy()
        sub["winter"] = np.concatenate(
            [np.full(len(idx[w]), j) for j, w in enumerate(pick)])
        try:
            v = fn(sub)
            if np.all(np.isfinite(v)):
                out.append(v)
        except Exception:
            pass
    return np.array(out)


def main():
    d = build()
    print(f"continuous sample: {len(d):,} winter days, "
          f"{d['winter'].nunique()} winters, "
          f"{d.index.min().date()}..{d.index.max().date()}")
    print(f"  (the event study uses 43 events in 36 winters)")

    # correlation structure, for transparency
    print(f"\n  corr(X wave driving, M vortex) = {d['X'].corr(d['M']):+.3f}"
          f"   (negative expected: waves decelerate the vortex)")
    print(f"  corr(M vortex, Y surface)      = {d['M'].corr(d['Y']):+.3f}"
          f"   (positive expected: weak vortex -> negative AO)")
    print(f"  corr(X wave driving, Y surface)= {d['X'].corr(d['Y']):+.3f}")

    def f_total(s):
        return fit(s, ["X"])

    def f_direct(s):
        return fit(s, ["X", "M"])

    total = float(f_total(d)[0])
    dm = f_direct(d)
    direct, m_coef = float(dm[0]), float(dm[1])
    indirect = total - direct

    bt = boot(d, lambda s: np.array([f_total(s)[0]]))
    bd = boot(d, lambda s: f_direct(s))
    bi = bt[:, 0] - bd[:, 0]

    def ci(a):
        return [round(float(np.percentile(a, 2.5)), 4),
                round(float(np.percentile(a, 97.5)), 4)]

    def p2(a):
        return float(2 * min((a >= 0).mean(), (a <= 0).mean()))

    share = indirect / total if total else np.nan
    bshare = bi / bt[:, 0]

    print(f"\n=== MEDIATION: does the stratosphere carry the signal? ===")
    print(f"  total   X -> Y            {total:+.4f}  CI {ci(bt[:,0])}  p={p2(bt[:,0]):.4f}")
    print(f"  direct  X -> Y | M        {direct:+.4f}  CI {ci(bd[:,0])}  p={p2(bd[:,0]):.4f}")
    print(f"  indirect (via vortex)     {indirect:+.4f}  CI {ci(bi)}  p={p2(bi):.4f}")
    print(f"  M coefficient (vortex->Y) {m_coef:+.4f}  CI {ci(bd[:,1])}  p={p2(bd[:,1]):.4f}")
    print(f"\n  MEDIATED SHARE = {100*share:.1f}%   CI "
          f"[{100*np.percentile(bshare,2.5):.1f}%, {100*np.percentile(bshare,97.5):.1f}%]")

    res = {"n_days": int(len(d)), "n_winters": int(d["winter"].nunique()),
           "X_lag": list(X_LAG), "Y_lag": list(Y_LAG), "n_boot": N_BOOT,
           "corr": {"X_M": round(float(d["X"].corr(d["M"])), 4),
                    "M_Y": round(float(d["M"].corr(d["Y"])), 4),
                    "X_Y": round(float(d["X"].corr(d["Y"])), 4)},
           "total_effect": {"est": round(total, 4), "CI95": ci(bt[:, 0]),
                            "p": p2(bt[:, 0])},
           "direct_effect": {"est": round(direct, 4), "CI95": ci(bd[:, 0]),
                             "p": p2(bd[:, 0])},
           "indirect_effect": {"est": round(indirect, 4), "CI95": ci(bi), "p": p2(bi)},
           "vortex_to_surface": {"est": round(m_coef, 4), "CI95": ci(bd[:, 1]),
                                 "p": p2(bd[:, 1])},
           "mediated_share": round(float(share), 4),
           "mediated_share_CI95": [round(float(np.percentile(bshare, 2.5)), 4),
                                   round(float(np.percentile(bshare, 97.5)), 4)]}

    # ---- precision: continuous vs event-study, same physical quantity ----
    se_cont = float(np.std(bd[:, 1]))
    print(f"\n=== PRECISION vs THE EVENT STUDY ===")
    print(f"  continuous vortex->surface coefficient SE: {se_cont:.4f} "
          f"({len(d):,} days, {d['winter'].nunique()} winters)")
    res["continuous_SE"] = round(se_cont, 4)

    # ---- sensitivity: how strong must an unmeasured M->Y confounder be? ----
    # If a confounder U induces spurious corr r_UM and r_UY, the bias on the
    # M->Y path is approximately r_UM * r_UY. Report the value that would null it.
    need = float(abs(m_coef))
    res["confounder_needed_product"] = round(need, 4)
    print(f"\n=== SENSITIVITY ===")
    print(f"  an unmeasured confounder of the vortex->surface path would need")
    print(f"  r(U,M) x r(U,Y) >= {need:.3f} to explain the mediated path away.")
    print(f"  (for equal correlations that is r >= {np.sqrt(need):.3f} on both)")

    (RESULTS / "continuous_mediation.json").write_text(json.dumps(res, indent=2),
                                                    encoding="utf8")
    print("\nSaved -> continuous_mediation.json")


if __name__ == "__main__":
    main()
