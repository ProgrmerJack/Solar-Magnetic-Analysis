#!/usr/bin/env python3
"""
stratifier_law_cmip6.py
=======================
THE POWER THE OBSERVATIONS DO NOT HAVE.

stratifier_inference.py established, on observations, that the between-group
surface contrast produced by splitting SSWs on a stratifier is predicted by
beta x dS (across 7 stratifiers under a true null: obs = 0.961 x pred + 0.002,
R^2 = 0.9945), and that the residual -- the part attributable to the event rather
than to the stratifier-surface regression -- is consistent with zero for every
stratifier.

That second statement is the weak one. With 29-43 events the null residual has an
SD of order 0.5-0.9 sigma, so "consistent with zero" bounds almost nothing. It is
an absence of evidence, and this project has already mistaken one of those for a
finding once (the AAO negative control) and nearly a second time (the bootstrap
"under-coverage" that was a sample-size confound).

The CMIP6 archive already pulled for this project carries 21 members across 11
models, ~165 winters each. That is of order 900 detected SSWs -- a factor ~30 more
than the observational catalogue, shrinking the residual CI by ~5.5x. It converts
"cannot distinguish from zero" into a real measurement of the residual.

WHAT IS TESTED, IDENTICALLY TO THE OBSERVATIONAL SCRIPT
  u10_pre     (-45,-16)  pre-onset vortex preconditioning     [not caused by the event]
  u10_post30  (0,30)     reversal depth                       [partly caused by the event]
  u10_post45  (0,45)     reversal persistence                 [partly caused by the event]
  AM_post     (8,52)     the outcome itself                   [the published criterion]

  For each: dY observed at real onsets, beta fitted ONLY on event-free pseudo-onset
  windows within the same member, dS, and residual = dY - beta x dS. The null
  residual distribution comes from pseudo-onsets drawn on event-cleaned model days,
  where the true differential effect is exactly zero by construction.

INTERPRETATION LIMIT, CARRIED OVER UNCHANGED
  For the post-onset stratospheric stratifiers the event partly CAUSES the
  stratifier, so beta x dS may contain a mediated causal path. Their residuals are
  an upper bound on contamination, not a causal effect. Only u10_pre (which the
  event cannot have caused) and AM_post (where the stratifier IS the outcome, so
  the whole term is selection) are causally interpretable.

  Models are not observations. What this establishes is whether the LAW holds in a
  system where the sample size permits testing it -- not what the real atmosphere's
  differential effect is.

Output: stratifier_law_cmip6.json
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
sys.path.insert(0, str(HERE.parents[0] / "07_physical_decomposition"))
import ensemble_precursor as EP                     # noqa: E402
import stratifier_bias_law as L                     # noqa: E402

RAW = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "raw" / "cmip6"
OUT_WIN = (8, 52)
INFLUENCE = (-60, 75)
SEASON = (11, 12, 1, 2, 3, 4)
STRATIFIERS = [
    ("u10_pre",    "pre",     (-45, -16)),
    ("u10_post30", "strat",   (0, 30)),
    ("u10_post45", "strat",   (0, 45)),
    ("AM_post",    "outcome", (8, 52)),
]
N_PSEUDO = 60          # pseudo-onset draws per member
SEED = 20260801


def influence_mask(index, onsets):
    dd = pd.DatetimeIndex(index).values.astype("datetime64[D]")
    m = np.zeros(len(dd), bool)
    for o in onsets:
        lag = (dd - np.datetime64(pd.Timestamp(o), "D")).astype(int)
        m |= (lag >= INFLUENCE[0]) & (lag <= INFLUENCE[1])
    return m


def zone_free_index(index, real_onsets):
    """Candidate pseudo-onset days whose OWN influence zone overlaps no real
    event's zone: |d - o| > INFLUENCE[1] - INFLUENCE[0] for every real onset o.

    Cleaning only the onset day (the previous practice) still let a pseudo-
    event's predictor, placebo and outcome windows run into a real SSW's
    response -- 6.3% of CMIP6 and 4.3% of ERA5 pseudo outcome windows, found by
    the methods audit 2026-09-25 -- which makes the null resemble real events and
    biases "event-specific" toward zero. A pseudo-event now gets the same
    exclusion zone a real one does; in practice it comes from an SSW-free winter.
    """
    dd = pd.DatetimeIndex(index).values.astype("datetime64[D]")
    ok = np.ones(len(dd), bool)
    sep = INFLUENCE[1] - INFLUENCE[0]
    for o in real_onsets:
        lag = (dd - np.datetime64(pd.Timestamp(o), "D")).astype(int)
        ok &= np.abs(lag) > sep
    return pd.DatetimeIndex(index)[ok]


def draw_pseudo(clean_index, doys, rng):
    by = {}
    for t in clean_index:
        by.setdefault(t.dayofyear, []).append(t)
    out = []
    for d in rng.permutation(doys):
        c = by.get(int(d))
        if c:
            out.append(c[rng.integers(len(c))])
    return pd.DatetimeIndex(out)


def win_mean(series, onsets, win, minfrac=3):
    out = []
    need = max(6, (win[1] - win[0]) // minfrac)
    for o in onsets:
        lo, hi = o + pd.Timedelta(days=win[0]), o + pd.Timedelta(days=win[1])
        v = series[(series.index >= lo) & (series.index <= hi)]
        out.append(float(v.mean()) if len(v) >= need else np.nan)
    return np.array(out)


def anom(series, onsets, clim, win):
    out = []
    need = max(6, (win[1] - win[0]) // 3)
    for o in onsets:
        lo, hi = o + pd.Timedelta(days=win[0]), o + pd.Timedelta(days=win[1])
        v = series[(series.index >= lo) & (series.index <= hi)]
        out.append(float(v.mean() - clim.reindex(v.index.dayofyear).mean())
                   if len(v) >= need else np.nan)
    return np.array(out)


def main():
    rng = np.random.default_rng(SEED)
    files = sorted(RAW.glob("*_zm.nc"))
    print(f"{len(files)} members on disk\n")

    real = {n: {"Y": [], "S": [], "mem": []} for n, *_ in STRATIFIERS}
    free = {n: {"Y": [], "S": []} for n, *_ in STRATIFIERS}
    nulldraw = {n: [] for n, *_ in STRATIFIERS}
    per_model = {}
    tot_ev = 0

    for f in files:
        model = f.name.split("_")[0]
        try:
            m = EP.load_member(f)
            full = m
        except Exception as e:
            print(f"  {f.name}: SKIP ({e})")
            continue
        m = m[np.isin(m.index.month, SEASON)].dropna()
        if len(m) < 2000:
            print(f"  {f.name}: SKIP (only {len(m)} days)")
            continue
        on = EP.detect_ssw(full["u10"].values, full.index)   # full daily series: CP07 needs contiguous days
        if len(on) < 15:
            print(f"  {f.name}: SKIP ({len(on)} events)")
            continue

        u10 = L.doy_standardise(m["u10"])
        am = m["am"]
        msk = influence_mask(am.index, on)
        clim = am[~msk].groupby(am[~msk].index.dayofyear).mean()
        clean = zone_free_index(am.index, on)
        doys = np.array([t.dayofyear for t in pd.DatetimeIndex(on)])

        Yr = anom(am, on, clim, OUT_WIN)
        amz = L.doy_standardise(am)
        for name, kind, win in STRATIFIERS:
            src = amz if name == "AM_post" else u10
            real[name]["Y"].append(Yr)
            real[name]["S"].append(win_mean(src, on, win))
            real[name]["mem"].append(np.full(len(on), f.name))

        # SAMPLE-SIZE MATCHING. The observed residual is computed on the POOLED
        # event set (~1900). If each null realisation used one member's ~90
        # pseudo-events the null would be sqrt(1900/90) ~ 4.6x too wide and every
        # p-value would be spuriously conservative -- the identical confound that
        # produced the phantom "bootstrap under-coverage" earlier in this project.
        # So draw k is stored PER MEMBER and pooled across members at index k,
        # giving null realisations of the same size as the observation.
        for k in range(N_PSEUDO):
            p = draw_pseudo(clean, doys, rng)
            if len(p) < 10:
                continue
            Yp = anom(am, p, clim, OUT_WIN)
            for name, kind, win in STRATIFIERS:
                src = amz if name == "AM_post" else u10
                Sp = win_mean(src, p, win)
                free[name]["Y"].append(Yp)
                free[name]["S"].append(Sp)
                nulldraw[name].append((k, Yp, Sp))

        per_model.setdefault(model, 0)
        per_model[model] += len(on)
        tot_ev += len(on)
        print(f"  {f.name:44s} {len(on):4d} events")

    print(f"\nTOTAL {tot_ev} events across {len(per_model)} models "
          f"(observations: 29-43)")

    res = {"n_events": int(tot_ev), "per_model": per_model,
           "n_pseudo_per_member": N_PSEUDO, "outcome_window": list(OUT_WIN),
           "stratifiers": {}}

    print(f"\n{'stratifier':12s} {'class':8s} {'dY_obs':>8s} {'dS':>7s} {'beta':>7s} "
          f"{'beta*dS':>8s} {'residual':>9s} {'null sd':>8s} {'p':>7s} {'%expl':>6s}")
    print("-" * 92)
    law_pred, law_obs = [], []
    for name, kind, win in STRATIFIERS:
        Y = np.concatenate(real[name]["Y"])
        Sv = np.concatenate(real[name]["S"])
        Yf = np.concatenate(free[name]["Y"])
        Sf = np.concatenate(free[name]["S"])
        ok = np.isfinite(Yf) & np.isfinite(Sf)
        beta = float(np.cov(Yf[ok], Sf[ok])[0, 1] / np.var(Sf[ok], ddof=1))

        dY, dS, n = L.contrast(Y, Sv)
        pb = beta * dS
        resid = dY - pb
        # pool each pseudo-draw index k across members -> null realisations of
        # the same size as the pooled observation
        byk = {}
        for k, yy, ss in nulldraw[name]:
            byk.setdefault(k, [[], []])
            byk[k][0].append(yy)
            byk[k][1].append(ss)
        nd = []
        for k, (ys, ss) in byk.items():
            d1, d2, nn = L.contrast(np.concatenate(ys), np.concatenate(ss))
            if np.isfinite(d1) and np.isfinite(d2):
                nd.append((d1, d2, nn))
        nd = np.array(nd)
        nres = nd[:, 0] - beta * nd[:, 1]
        p = float((np.abs(nres) >= abs(resid)).mean())
        expl = 100 * pb / dY if abs(dY) > 1e-9 else np.nan

        law_pred.append(float((beta * nd[:, 1]).mean()))
        law_obs.append(float(nd[:, 0].mean()))

        res["stratifiers"][name] = {
            "class": kind, "window": list(win), "n_events": int(n),
            "dY_observed": round(dY, 4), "dS": round(dS, 4),
            "beta": round(beta, 4), "predicted_bias": round(float(pb), 4),
            "residual": round(float(resid), 4),
            "null_residual_sd": round(float(nres.std()), 4),
            "null_realisation_size": int(np.median(nd[:, 2])),
            "n_null_realisations": int(len(nd)),
            "p_residual": p,
            "pct_explained_by_stratifier": round(float(expl), 1),
            "causally_interpretable": kind != "strat"}
        print(f"{name:12s} {kind:8s} {dY:+8.3f} {dS:+7.3f} {beta:+7.3f} {pb:+8.3f} "
              f"{resid:+9.3f} {nres.std():8.3f} {p:7.4f} {expl:5.0f}%")

    lp, lo_ = np.array(law_pred), np.array(law_obs)
    if len(lp) > 2:
        b = np.polyfit(lp, lo_, 1)
        r2 = float(np.corrcoef(lp, lo_)[0, 1] ** 2)
        res["law_across_stratifiers"] = {
            "slope": round(float(b[0]), 4), "intercept": round(float(b[1]), 4),
            "r2": round(r2, 4), "n_points": len(lp)}
        print(f"\nLAW ACROSS STRATIFIERS (null, {len(lp)} points): "
              f"obs = {b[0]:.3f} x pred {b[1]:+.3f}, R^2 = {r2:.4f}")
        print(f"  (observations gave slope 0.961, intercept +0.002, R^2 = 0.9945)")

    print("\n=== WHAT THE POWER BUYS ===")
    for name, kind, win in STRATIFIERS:
        r = res["stratifiers"][name]
        sd = r["null_residual_sd"]
        verdict = ("residual NONZERO" if r["p_residual"] < 0.05
                   else f"residual zero to +/-{1.96 * sd:.3f}")
        note = "" if kind != "strat" else "   (upper bound only, mediation)"
        print(f"  {name:12s} n={r['n_events']:4d}  {verdict}{note}")

    (RESULTS / "stratifier_law_cmip6.json").write_text(json.dumps(res, indent=2),
                                                    encoding="utf8")
    print("\nSaved -> stratifier_law_cmip6.json")


if __name__ == "__main__":
    main()
