#!/usr/bin/env python3
"""
forecast_value_test.py
======================
DOES KNOWING THE EVENT IMPROVE A PROBABILISTIC FORECAST BEYOND THE SHIFT?

Plan approved 2026-09-25 (item C of the Nature Geoscience revision plan).

QUESTION
  For CMIP6 SSWs (Charlton-Polvani detection, 1,517 events, 20 members, 10
  models), compare two probabilistic forecasts of each event's surface response
  (annular mode, days +8..+52, in the member's own sigma units, NOT demeaned, so
  the common shift counts as forecast information):
    F0 SHIFT        Gaussian; mean = calendar regression (day-of-year sin/cos)
                    fitted on SSWs of the OTHER models; spread = its residual sd.
    F1 EVENT-AWARE  Gaussian; mean = ridge prediction from the stratospheric
                    predictors available at the stated lead plus calendar; spread
                    = out-of-fold residual sd within the training models.
  Tiers: P2 (information up to onset) primary, P1 (before onset) secondary.
  Score: CRPS on held-out events, leaving out one whole MODEL at a time (so
  CanESM5's ten members cannot dominate); CRPSS = 1 - CRPS(F1)/CRPS(F0).

PRE-SPECIFIED READOUTS
  (1) CRPSS at SSWs, with a model-cluster bootstrap interval: does event-aware
      information help at all?
  (2) SSW-specific value: CRPSS at SSWs minus CRPSS of the identical pipeline on
      size-matched event-free dates (same members, same counts and calendar days,
      pseudo-onsets > 135 d from every real onset); p = P(pseudo >= real).
  Falsifier of the thesis: (2) positive at p < 0.05 -- event-aware information is
  worth more at SSWs than on ordinary winter days.

Output: results/current/6_predictability/forecast_value_test.json
"""
import json
import multiprocessing as mp
import sys
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.linear_model import LinearRegression, RidgeCV
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "6_predictability"
sys.path.insert(0, str(HERE))
import within_model_check as W                       # noqa: E402

K_PSEUDO = 200
N_BOOT = 2000
WORKERS = 16
NAME = "forecast_value_test"
SEED = zlib.crc32(NAME.encode()) % (2 ** 32)
TIERS = {"P2 at-onset": 2, "P1 pre-onset": 1}


def crps_gauss(y, mu, sd):
    z = (y - mu) / sd
    return sd * (z * (2 * norm.cdf(z) - 1) + 2 * norm.pdf(z) - 1 / np.sqrt(np.pi))


def std_within(X, g):
    """Standardise each predictor within member (predictors only)."""
    Xo = X.astype(float).copy()
    for m in np.unique(g):
        i = g == m
        mu = np.nanmean(Xo[i], axis=0); sd = np.nanstd(Xo[i], axis=0)
        Xo[i] = (Xo[i] - mu) / np.where(sd > 1e-12, sd, 1.0)
    return Xo


def scores(X, cal, y, g, model):
    """Leave-one-model-out CRPS of F0 and F1 for every event (NaN rows dropped)."""
    ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
    X, cal, y, g, model = X[ok], cal[ok], y[ok], g[ok], model[ok]
    c0 = np.full(len(y), np.nan); c1 = np.full(len(y), np.nan)
    for mdl in np.unique(model):
        te, tr = model == mdl, model != mdl
        f0 = LinearRegression().fit(cal[tr], y[tr])
        mu0 = f0.predict(cal[te]); sd0 = np.std(y[tr] - f0.predict(cal[tr]), ddof=1)
        XA = np.column_stack([X, cal])
        f1 = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 4, 25)))
        # out-of-fold residual spread inside the training models
        oof = np.full(tr.sum(), np.nan); gi = g[tr]
        for a, b in GroupKFold(5).split(XA[tr], y[tr], gi):
            oof[b] = f1.fit(XA[tr][a], y[tr][a]).predict(XA[tr][b])
        sd1 = np.std(y[tr] - oof, ddof=1)
        mu1 = f1.fit(XA[tr], y[tr]).predict(XA[te])
        c0[te] = crps_gauss(y[te], mu0, sd0); c1[te] = crps_gauss(y[te], mu1, sd1)
    return c0, c1, model


def crpss(c0, c1):
    return float(1 - c1.sum() / c0.sum())


_S = {}


def _pseudo_task(k):
    """CRPSS of pseudo draw k (every member's k-th size-matched draw), per tier."""
    out = {}
    for tier, ci in _S["tier_ci"].items():
        X = np.concatenate([_S["M"][m][0][k][:, ci] for m in _S["members"]])
        cal = np.concatenate([_S["M"][m][0][k][:, _S["cal_ci"]] for m in _S["members"]])
        y = np.concatenate([_S["M"][m][1][k] for m in _S["members"]])
        g = np.concatenate([np.full(_S["M"][m][1].shape[1], m) for m in _S["members"]])
        mdl = np.array([s.split("_")[0] for s in g])
        c0, c1, _ = scores(std_within(X, g), cal, y, g, mdl)
        out[tier] = crpss(c0, c1)
    return out


def main():
    files = sorted(W.RAW.glob("*_zm.nc"))
    real_X, real_y, real_g, cols = [], [], [], None
    usable = []
    for f in files:
        c = W.prepare_member(f)
        if c is None:
            continue
        F = W.H.features_for(c["u"], c["plev"], c["lat"], c["idx"], c["on"])
        Y = W.C6.anom(c["am"], c["on"], c["cl"], W.OUT_WIN)
        real_X.append(F.values); real_y.append(Y); real_g.append(np.full(len(Y), f.stem))
        cols = list(F.columns); usable.append(f)
    X = np.concatenate(real_X); y = np.concatenate(real_y); g = np.concatenate(real_g)
    mdl = np.array([s.split("_")[0] for s in g])
    cal_ci = [cols.index("doy_sin"), cols.index("doy_cos")]
    tier_ci = {t: [cols.index(c) for c in W.P.tier_cols(cols, n)] for t, n in TIERS.items()}
    print(f"{len(y)} real events, {len(usable)} members, {len(np.unique(mdl))} models", flush=True)

    res = {"plan_approved": "2026-09-25", "n_events": int(np.isfinite(y).sum()),
           "n_members": len(usable), "n_models": int(len(np.unique(mdl))),
           "outcome_window_days": list(W.OUT_WIN), "k_pseudo": K_PSEUDO, "n_boot": N_BOOT,
           "seed": SEED, "cv": "leave one model out", "tiers": {}}

    # size-matched pseudo draws, each member's own seeded stream (as in J)
    W.K_DRAWS = K_PSEUDO
    # "fork": the pseudo-draw workers read _S, which Python 3.14's default
    # start method (no longer fork on Linux) would not pass on
    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(WORKERS, mp_context=ctx) as ex:
        M = {stem: (Fk, Yk) for stem, Fk, Yk, _ in ex.map(W.matched_draws, list(enumerate(map(str, usable))))}
    _S.update(M=M, members=[f.stem for f in usable], tier_ci=tier_ci, cal_ci=cal_ci)
    with ProcessPoolExecutor(WORKERS, mp_context=ctx) as ex:
        pseudo = list(ex.map(_pseudo_task, range(K_PSEUDO)))

    rng = np.random.default_rng(SEED)
    for tier, ci in tier_ci.items():
        c0, c1, m_te = scores(std_within(X[:, ci], g), X[:, cal_ci], y, g, mdl)
        real = crpss(c0, c1)
        um = np.unique(m_te); boot = []
        for _ in range(N_BOOT):
            pick = rng.choice(um, len(um), replace=True)
            idx = np.concatenate([np.flatnonzero(m_te == u) for u in pick])
            boot.append(crpss(c0[idx], c1[idx]))
        ps = np.array([p[tier] for p in pseudo])
        res["tiers"][tier] = {
            "n_features": len(ci),
            "CRPS_shift": round(float(np.mean(c0)), 5), "CRPS_event_aware": round(float(np.mean(c1)), 5),
            "CRPSS_ssw": round(real, 4),
            "CRPSS_ssw_CI95_model_bootstrap": [round(float(np.percentile(boot, 2.5)), 4),
                                               round(float(np.percentile(boot, 97.5)), 4)],
            "CRPSS_pseudo_mean": round(float(ps.mean()), 4),
            "CRPSS_pseudo_q025_q975": [round(float(np.percentile(ps, 2.5)), 4),
                                       round(float(np.percentile(ps, 97.5)), 4)],
            "ssw_specific_CRPSS": round(real - float(ps.mean()), 4),
            "p_pseudo_ge_ssw": round(float(np.mean(ps >= real)), 4)}
        r = res["tiers"][tier]
        print(f"{tier:14s} CRPSS at SSWs {real:+.4f} {r['CRPSS_ssw_CI95_model_bootstrap']} | "
              f"ordinary dates {r['CRPSS_pseudo_mean']:+.4f} {r['CRPSS_pseudo_q025_q975']} | "
              f"SSW-specific {r['ssw_specific_CRPSS']:+.4f}, p = {r['p_pseudo_ge_ssw']:.3f}", flush=True)
    (RESULTS / "forecast_value_test.json").write_text(json.dumps(res, indent=2), encoding="utf8",
                                                      newline="\n")
    print("Saved -> forecast_value_test.json")


if __name__ == "__main__":
    main()
