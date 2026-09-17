#!/usr/bin/env python3
"""
mediation_disattenuated.py
==========================
Corrects the mediation estimate for measurement error in the mediator, and turns
that correction into the actual finding.

THE CHAIN OF REASONING
  continuous_mediation.py: stratosphere mediates 14.2% of the
  wave-driving-to-surface link (CI -35% to +60%), 86% "direct". Read naively
  that is the common-cause position -- the troposphere drives both, the
  stratosphere is a bystander.

  mediator_specification.py: that number is an artefact. The mediated share
  climbs monotonically as the mediator captures more of the persistent vortex
  anomaly -- 14.2% (single day) -> 45.2% (mean over t..t+10, CI +8% to +88%) --
  and keeps climbing for windows that overlap the outcome and are therefore not
  interpretable. The single-day mediator was simply a noisy proxy.

  This is textbook: classical measurement error in a MEDIATOR attenuates the
  indirect path toward zero and inflates the apparent direct path. So a study
  using a noisy stratospheric index will UNDER-attribute to the stratosphere and
  OVER-attribute to direct tropospheric pathways -- which is the direction of the
  published common-cause claims.

THE CORRECTION
  Two imperfect measures of one latent vortex state:
      M1 = 10 hPa 60N zonal wind, mean t..t+10
      M2 = 100 hPa polar-cap geopotential height (sign-flipped), mean t..t+10
  Both end at t+10, strictly before the outcome window t+15..t+45, so no
  mechanical overlap.

  Under classical errors-in-variables with independent measurement errors,
  regressing Y on M1 using M2 as an INSTRUMENT gives a consistent estimate of the
  latent M*->Y path: beta_IV = cov(Y,M2) / cov(M1,M2). Reliability of M1 is
  estimated as lambda = cov(M1,M2) / var(M1) scaled by the regression of M1 on
  M2, and the disattenuated mediated share follows.

ASSUMPTION, STATED PLAINLY AND NOT WAVED AWAY
  The correction requires the measurement errors of M1 and M2 to be uncorrelated.
  u10 and polar-cap height are dynamically linked through thermal wind, so their
  errors are NOT guaranteed independent. If the errors correlate positively, the
  IV estimate is biased toward the naive one and the correction is CONSERVATIVE
  -- it understates the true mediation. The result is therefore presented as a
  LOWER BOUND on stratospheric mediation, which is the direction that matters for
  the argument being made.

Output: mediation_disattenuated.json
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
import mediator_specification as MS                 # noqa: E402

N_BOOT = 2000
MWIN = 10          # mediator averaged t..t+10, strictly before outcome at t+15


def resid(d, cols):
    """Within-winter, seasonally adjusted residuals of each column."""
    S = C.harmonics(d["doy"].values, 3)
    out = {}
    A = np.column_stack([d[c].values.astype(float) for c in cols] + [S])
    Ad = K.C_demean(A, d["winter"].values)
    n = len(cols)
    Sd = Ad[:, n:]
    for i, c in enumerate(cols):
        b, *_ = np.linalg.lstsq(Sd, Ad[:, i], rcond=None)
        out[c] = Ad[:, i] - Sd @ b
    return out


def estimate(d):
    r = resid(d, ["X", "Y", "M1", "M2"])
    X, Y, M1, M2 = r["X"], r["Y"], r["M1"], r["M2"]

    def cov(a, b):
        return float(np.mean(a * b))

    # total effect of X on Y
    tot = cov(X, Y) / cov(X, X)
    # naive mediation with M1
    A = np.column_stack([X, M1])
    bn, *_ = np.linalg.lstsq(A, Y, rcond=None)
    naive_direct = float(bn[0])
    naive_share = (tot - naive_direct) / tot if tot else np.nan
    # IV: instrument M1 with M2, controlling X (partial out X first)
    def partial(v):
        b = cov(X, v) / cov(X, X)
        return v - b * X
    Yp, M1p, M2p = partial(Y), partial(M1), partial(M2)
    denom = cov(M1p, M2p)
    beta_iv = cov(Yp, M2p) / denom if denom else np.nan
    # reliability of M1 as a measure of the latent state
    lam = denom / cov(M1p, M1p) if cov(M1p, M1p) else np.nan
    # disattenuated direct effect: Y - beta_iv * M1 regressed on X
    Yadj = Y - beta_iv * M1
    dis_direct = cov(X, Yadj) / cov(X, X)
    dis_share = (tot - dis_direct) / tot if tot else np.nan
    return dict(total=tot, naive_direct=naive_direct, naive_share=naive_share,
                beta_iv=beta_iv, reliability=lam,
                dis_direct=dis_direct, dis_share=dis_share)


def main():
    d = MS.build_all()
    st = pd.read_parquet(MS.STRAT)
    if st.index.tz is not None:
        st.index = st.index.tz_convert("UTC").tz_localize(None)

    def fwd(s, k):
        s = MS.daily(s.sort_index())
        return s[::-1].rolling(k + 1, min_periods=max(3, k // 2)).mean()[::-1]

    d["M1"] = fwd(st["uwnd_ms_10hPa"], MWIN).reindex(d.index).values
    d["M2"] = -fwd(st["hgt_m_100hPa"], MWIN).reindex(d.index).values
    d = d[["winter", "doy", "X", "Y", "M1", "M2"]].dropna()
    for c in ("X", "Y", "M1", "M2"):
        d[c] = (d[c] - d[c].mean()) / d[c].std()
    print(f"sample {len(d):,} winter days, {d['winter'].nunique()} winters")
    print(f"  mediator window t..t+{MWIN}, outcome t+15..t+45 -- no overlap")
    print(f"  corr(M1,M2) = {d['M1'].corr(d['M2']):+.3f}  "
          f"(two measures of one latent vortex state)")

    e = estimate(d)
    rng = np.random.default_rng(20260731)
    winters = np.array(sorted(d["winter"].unique()))
    idxw = {w: np.flatnonzero((d["winter"] == w).values) for w in winters}
    keys = ["total", "naive_share", "beta_iv", "reliability", "dis_share"]
    bs = {k: [] for k in keys}
    for _ in range(N_BOOT):
        pick = rng.choice(winters, len(winters), replace=True)
        sub = d.iloc[np.concatenate([idxw[w] for w in pick])].copy()
        sub["winter"] = np.concatenate(
            [np.full(len(idxw[w]), j) for j, w in enumerate(pick)])
        try:
            v = estimate(sub)
            if all(np.isfinite(v[k]) for k in keys):
                for k in keys:
                    bs[k].append(v[k])
        except Exception:
            pass

    def ci(k):
        a = np.array(bs[k])
        return [round(float(np.percentile(a, 2.5)), 4),
                round(float(np.percentile(a, 97.5)), 4)]

    print(f"\n=== MEDIATION, NAIVE vs DISATTENUATED ===")
    print(f"  total effect X -> Y            {e['total']:+.4f}  CI {ci('total')}")
    print(f"  reliability of M1              {e['reliability']:.3f}  CI {ci('reliability')}"
          f"   (1.0 = measured without error)")
    print(f"  latent vortex -> surface (IV)  {e['beta_iv']:+.4f}  CI {ci('beta_iv')}")
    print(f"\n  mediated share, NAIVE          {100*e['naive_share']:.1f}%  "
          f"CI [{100*ci('naive_share')[0]:.0f}%, {100*ci('naive_share')[1]:.0f}%]")
    print(f"  mediated share, DISATTENUATED  {100*e['dis_share']:.1f}%  "
          f"CI [{100*ci('dis_share')[0]:.0f}%, {100*ci('dis_share')[1]:.0f}%]")
    print(f"\n  single-day mediator gave 14.2% (continuous_mediation.py)")
    print(f"  -> measurement error alone moves the answer by "
          f"{100*(e['dis_share']-0.142):.0f} percentage points")

    res = {"n_days": int(len(d)), "n_winters": int(d["winter"].nunique()),
           "mediator_window": MWIN, "n_boot": N_BOOT,
           "corr_M1_M2": round(float(d["M1"].corr(d["M2"])), 4),
           "total_effect": {"est": round(e["total"], 4), "CI95": ci("total")},
           "reliability_M1": {"est": round(e["reliability"], 4),
                              "CI95": ci("reliability")},
           "latent_vortex_to_surface_IV": {"est": round(e["beta_iv"], 4),
                                           "CI95": ci("beta_iv")},
           "mediated_share_naive": {"est": round(e["naive_share"], 4),
                                    "CI95": ci("naive_share")},
           "mediated_share_disattenuated": {"est": round(e["dis_share"], 4),
                                            "CI95": ci("dis_share")},
           "single_day_mediator_share": 0.142,
           "assumption": ("classical errors-in-variables with independent errors in "
                          "M1 and M2; positively correlated errors would bias the IV "
                          "estimate toward the naive one, so this is a LOWER BOUND "
                          "on stratospheric mediation")}
    (RESULTS / "mediation_disattenuated.json").write_text(json.dumps(res, indent=2),
                                                       encoding="utf8")
    print("\nSaved -> mediation_disattenuated.json")


if __name__ == "__main__":
    main()
