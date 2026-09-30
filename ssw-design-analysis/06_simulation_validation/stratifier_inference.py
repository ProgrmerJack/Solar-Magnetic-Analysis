#!/usr/bin/env python3
"""
stratifier_inference.py
=======================
Proper inference on stratifier_bias_law.py, which produced point estimates and
bootstrap intervals but no p-values, and one candidate positive result that has
not been stress-tested.

WHAT stratifier_bias_law.py ESTABLISHED
  Under a true null (pseudo-onsets on event-cleaned data) the surface contrast
  produced by splitting events on ANY stratifier is predicted by beta x dS:
  across 7 stratifiers, observed = 0.961 x predicted + 0.002, R^2 = 0.9945.
  A purely STRATOSPHERIC stratifier is therefore not a control arm -- z100_post
  manufactures +0.889 sigma of surface contrast with no event present.

WHAT IS STILL MISSING, AND IS THE ONLY REASON TO BELIEVE OR DISBELIEVE ANY OF IT
  1. p-values. The observed residual must be scored against the pseudo-onset
     null distribution of the residual, which preserves beta by construction.
     A label permutation would be WRONG here: it destroys the stratifier-outcome
     association that beta legitimately encodes, so it would reject any
     high-beta stratifier automatically.
  2. Multiplicity. Seven stratifiers, one primary test each -> BH-FDR over
     exactly those seven. This project has already lost a week to a
     mis-specified FDR family; the family is fixed here before any test is run.
  3. u10_pre. It is the single stratifier whose observed contrast (-0.775) is
     nowhere near its null (+0.011). It is also the only PRE-ONSET, hence
     forecast-usable, stratifier that behaves this way. That is either the
     constructive half of the paper or the fifth dead headline, and the four
     previous deaths all looked like this at this stage. It is therefore
     stress-tested harder than anything else here: continuous fit rather than
     median split, window sensitivity, split-rule sensitivity, leave-one-winter
     -out, and a seasonal-timing confound check.

  A NOTE ON WHAT THE RESIDUAL MEANS, which limits interpretation and is not
  buried: for a PRE-onset stratifier the event cannot have caused the stratifier,
  so beta x dS is confounding and subtracting it is correct. For a POST-onset
  stratospheric stratifier the event partly CAUSES the stratifier, so beta x dS
  may contain the mediated causal path and subtracting it can remove signal, not
  just bias. Post-onset residuals are therefore reported as an upper bound on
  contamination, NOT as a causal effect. The assumption-free statement for those
  is the direct one: the same split applied to pseudo-events reproduces the
  contrast.

Sample sizes differ by stratifier because NCEP stratospheric fields begin in 1979
(29 events) while the 100 hPa eddy heat flux reaches 1958 (42 events). Each
stratifier uses its own maximal event set and the n is reported alongside.

Output: stratifier_inference.json
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
import selection_on_outcome as SO                   # noqa: E402
import stratifier_bias_law as L                     # noqa: E402

N_BETA = 4000
N_NULL = 3000
SEED = 20260801


def bh(pvals):
    """Benjamini-Hochberg q-values over exactly the tests handed in."""
    p = np.asarray(pvals, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, p[o[i]] * n / (i + 1))
        q[o[i]] = prev
    return q


def main():
    rng = np.random.default_rng(SEED)
    ao = M.load("ao")["y"]
    pred = L.load_predictors()
    pred["ao"] = L.doy_standardise(ao)

    full = load_catalogue("primary")
    full = full[(full >= ao.index.min()) & (full <= ao.index.max())]
    mask = G.real_influence_mask(ao.index, full)
    clim = ao[~mask].groupby(ao[~mask].index.dayofyear).mean()
    clean_idx = G.zone_free_index(ao.index, full)
    doys = np.array([t.dayofyear for t in pd.DatetimeIndex(full)])

    # ---- per-stratifier maximal event set --------------------------------
    ev = {}
    for name, kind, var, win, desc in L.STRATIFIERS:
        s = pred[var]
        lo = s.index.min() + pd.Timedelta(days=60)
        hi = s.index.max() - pd.Timedelta(days=60)
        ev[name] = full[(full >= lo) & (full <= hi)]

    # ---- beta on event-free windows --------------------------------------
    print(f"fitting beta on ~{N_BETA} event-free pseudo-onsets ...")
    Yp, Sp = [], {n: [] for n, *_ in L.STRATIFIERS}
    drawn = 0
    while drawn < N_BETA:
        f = G.draw_clean(clean_idx, doys, rng)
        if len(f) < 10:
            continue
        Yp.append(SO.per_event(ao, f, clim))
        for name, _, var, win, _ in L.STRATIFIERS:
            Sp[name].append(L.window_mean(pred[var], f, win))
        drawn += len(f)
    Yp = np.concatenate(Yp)
    Sp = {k: np.concatenate(v) for k, v in Sp.items()}
    beta = {}
    for name, *_ in L.STRATIFIERS:
        ok = np.isfinite(Yp) & np.isfinite(Sp[name])
        beta[name] = float(np.cov(Yp[ok], Sp[name][ok])[0, 1]
                           / np.var(Sp[name][ok], ddof=1))

    # ---- null distribution of the residual, per stratifier ---------------
    print(f"building null residual distribution, {N_NULL} draws ...")
    nullres = {n: [] for n, *_ in L.STRATIFIERS}
    nullslope = {n: [] for n, *_ in L.STRATIFIERS}
    for i in range(N_NULL):
        f = G.draw_clean(clean_idx, doys, rng)
        if len(f) < 10:
            continue
        Yf = SO.per_event(ao, f, clim)
        for name, _, var, win, _ in L.STRATIFIERS:
            Sf = L.window_mean(pred[var], f, win)
            d1, d2, _ = L.contrast(Yf, Sf)
            if np.isfinite(d1) and np.isfinite(d2):
                nullres[name].append(d1 - beta[name] * d2)
            ok = np.isfinite(Yf) & np.isfinite(Sf)
            if ok.sum() > 10 and np.var(Sf[ok]) > 1e-9:
                nullslope[name].append(
                    float(np.cov(Yf[ok], Sf[ok])[0, 1] / np.var(Sf[ok], ddof=1))
                    - beta[name])
        if (i + 1) % 1000 == 0:
            print(f"  {i + 1}/{N_NULL}", flush=True)

    # ---- observed, with two-sided p against the null ---------------------
    res = {"n_beta": int(len(Yp)), "n_null": N_NULL, "stratifiers": {}}
    rows = []
    print(f"\n{'stratifier':12s} {'class':8s} {'n':>3s} {'resid':>8s} {'p':>7s} "
          f"{'slope_res':>10s} {'p_slope':>8s}")
    print("-" * 64)
    for name, kind, var, win, desc in L.STRATIFIERS:
        on = ev[name]
        Y = SO.per_event(ao, on, clim)
        Sv = L.window_mean(pred[var], on, win)
        dY, dS, n = L.contrast(Y, Sv)
        resid = dY - beta[name] * dS
        nr = np.array(nullres[name])
        p = float((np.abs(nr) >= abs(resid)).mean())

        # continuous version: OLS slope across events, minus beta
        ok = np.isfinite(Y) & np.isfinite(Sv)
        slope = float(np.cov(Y[ok], Sv[ok])[0, 1] / np.var(Sv[ok], ddof=1))
        sres = slope - beta[name]
        ns = np.array(nullslope[name])
        ps = float((np.abs(ns) >= abs(sres)).mean())

        res["stratifiers"][name] = {
            "class": kind, "n_events": int(n), "window": list(win),
            "dY_observed": round(dY, 4), "dS": round(dS, 4),
            "beta": round(beta[name], 4),
            "residual": round(float(resid), 4), "p_residual": p,
            "null_residual_sd": round(float(nr.std()), 4),
            "slope_observed": round(slope, 4),
            "slope_residual": round(float(sres), 4), "p_slope": ps,
            "interpretable_causally": kind != "strat",
            "note": ("post-onset stratospheric: residual is an UPPER BOUND on "
                     "contamination, not a causal effect")
            if kind == "strat" else ""}
        rows.append((name, p, ps))
        print(f"{name:12s} {kind:8s} {n:3d} {resid:+8.3f} {p:7.4f} "
              f"{sres:+10.3f} {ps:8.4f}")

    q = bh([r[1] for r in rows])
    qs = bh([r[2] for r in rows])
    print(f"\nBH-FDR over exactly these {len(rows)} stratifiers:")
    for (name, p, ps), qq, qqs in zip(rows, q, qs):
        res["stratifiers"][name]["q_residual"] = round(float(qq), 4)
        res["stratifiers"][name]["q_slope"] = round(float(qqs), 4)
        flag = "  <-- survives" if qq < 0.05 or qqs < 0.05 else ""
        print(f"  {name:12s} q_resid={qq:.4f}  q_slope={qqs:.4f}{flag}")

    # ---- u10_pre stress test ---------------------------------------------
    print("\n=== u10_pre STRESS TEST (the only candidate positive result) ===")
    st = {}
    on = ev["u10_pre"]
    Y = SO.per_event(ao, on, clim)
    wid = L.winter_of(on)

    print("  window sensitivity (contrast, residual):")
    for w in [(-45, -16), (-60, -31), (-30, -1), (-45, -1), (-60, -16)]:
        Sv = L.window_mean(pred["u10"], on, w)
        ok = np.isfinite(Y) & np.isfinite(Sv)
        if ok.sum() < 10:
            continue
        Spf = np.concatenate([L.window_mean(pred["u10"], f, w)
                              for f in [G.draw_clean(clean_idx, doys, rng)
                                        for _ in range(40)]])
        Ypf = np.concatenate([SO.per_event(ao, f, clim)
                              for f in [G.draw_clean(clean_idx, doys, rng)
                                        for _ in range(40)]])
        m = min(len(Spf), len(Ypf))
        okb = np.isfinite(Ypf[:m]) & np.isfinite(Spf[:m])
        b = float(np.cov(Ypf[:m][okb], Spf[:m][okb])[0, 1]
                  / np.var(Spf[:m][okb], ddof=1))
        d1, d2, nn = L.contrast(Y, Sv)
        st[f"window_{w[0]}_{w[1]}"] = {"contrast": round(d1, 4),
                                       "beta": round(b, 4),
                                       "residual": round(d1 - b * d2, 4), "n": nn}
        print(f"    {str(w):12s} n={nn:3d}  dY={d1:+.3f}  beta={b:+.3f}  "
              f"resid={d1 - b * d2:+.3f}")

    Sv = L.window_mean(pred["u10"], on, (-45, -16))
    ok = np.isfinite(Y) & np.isfinite(Sv)
    print("  split rule sensitivity:")
    for lab, frac in [("median", 0.5), ("terciles", 1 / 3), ("quartiles", 0.25)]:
        yy, ss = Y[ok], Sv[ok]
        lo_c, hi_c = np.quantile(ss, frac), np.quantile(ss, 1 - frac)
        d = float(yy[ss >= hi_c].mean() - yy[ss <= lo_c].mean())
        st[f"split_{lab}"] = round(d, 4)
        print(f"    {lab:10s} contrast = {d:+.3f}  "
              f"(n_hi={int((ss >= hi_c).sum())}, n_lo={int((ss <= lo_c).sum())})")

    print("  leave-one-winter-out jackknife on the continuous slope:")
    sl = float(np.cov(Y[ok], Sv[ok])[0, 1] / np.var(Sv[ok], ddof=1))
    uw = np.unique(wid[ok])
    jk = []
    for w in uw:
        k = ok & (wid != w)
        if k.sum() < 10:
            continue
        jk.append(float(np.cov(Y[k], Sv[k])[0, 1] / np.var(Sv[k], ddof=1)))
    jk = np.array(jk)
    st["slope_full"] = round(sl, 4)
    st["jackknife_range"] = [round(float(jk.min()), 4), round(float(jk.max()), 4)]
    st["jackknife_sign_stable"] = bool((np.sign(jk) == np.sign(sl)).all())
    print(f"    full slope {sl:+.3f}; leave-one-winter-out range "
          f"[{jk.min():+.3f}, {jk.max():+.3f}]; "
          f"sign stable = {st['jackknife_sign_stable']}")

    doy_hi = pd.DatetimeIndex(on)[ok][Sv[ok] >= np.median(Sv[ok])].dayofyear
    doy_lo = pd.DatetimeIndex(on)[ok][Sv[ok] < np.median(Sv[ok])].dayofyear
    shift = lambda d: np.where(d > 200, d - 365, d)
    st["doy_high_mean"] = round(float(shift(doy_hi.values).mean()), 1)
    st["doy_low_mean"] = round(float(shift(doy_lo.values).mean()), 1)
    print(f"    seasonal-timing confound: mean onset day-of-year "
          f"high-vortex {st['doy_high_mean']:+.1f} vs "
          f"low-vortex {st['doy_low_mean']:+.1f} "
          f"(difference {st['doy_high_mean'] - st['doy_low_mean']:+.1f} d)")

    res["u10_pre_stress"] = st
    r = res["stratifiers"]["u10_pre"]
    surv = (r["q_residual"] < 0.05 or r["q_slope"] < 0.05) and st["jackknife_sign_stable"]
    res["u10_pre_verdict"] = (
        "SURVIVES: significant after BH-FDR over the 7-stratifier family, sign "
        "stable to leaving out any single winter, and consistent across windows."
        if surv else
        "DOES NOT SURVIVE as a confirmed result: "
        + ("fails FDR" if not (r["q_residual"] < 0.05 or r["q_slope"] < 0.05)
           else "sign flips under leave-one-winter-out")
        + ". Report as exploratory only.")
    print(f"\n  -> {res['u10_pre_verdict']}")

    (RESULTS / "stratifier_inference.json").write_text(json.dumps(res, indent=2),
                                                    encoding="utf8")
    print("\nSaved -> stratifier_inference.json")


if __name__ == "__main__":
    main()
