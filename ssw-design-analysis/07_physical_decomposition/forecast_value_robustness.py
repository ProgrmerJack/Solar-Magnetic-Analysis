#!/usr/bin/env python3
"""
forecast_value_robustness.py
============================
DOES THE FORECAST-VALUE RESULT (Q) SURVIVE THE REFEREE'S OBJECTIONS?

Plan approved 2026-09-29. Result Q (forecast_value_test.py) found that stratospheric
information adds <1% CRPS skill beyond the shifted distribution after CMIP6 SSWs,
less than on ordinary winter days. Objections, each answered here with the same
pipeline (Gaussian SHIFT vs EVENT-AWARE forecast, leave-one-MODEL-out CRPS):
  O1 post-onset stratosphere left out            -> tier P3 (adds days 0..30)
  O2 ordinary days are SSW-free winters, enriched in strong vortices, so the
     predictors span a wider range there        -> a WEAK-VORTEX null: zone-free
     dates whose member u(10 hPa, 60N) lies in that member's lowest 15% of
     zone-free November-March days, one per real onset within +-30 calendar days;
     and the across-event s.d. of the predictors at SSWs vs both nulls
  O3 ten of twenty members are CanESM5           -> CRPSS with models weighted
     equally, and with CanESM5 excluded
  O4 "event information" is only zonal-mean wind -> predictor set PLUS: the base
     set, wave-driving proxies (zonal-wind tendency at 10/50/100 hPa, 60N and cap;
     10-100 hPa shear; 75N-45N geometry, as in predictability_wave_driving.py)
     and the pre-onset SURFACE annular mode (a tropospheric precursor), all in
     the same windows and tiers.
  Written before any of these was computed. Sanity: the base set must reproduce
  Q's P2 and P1 CRPSS at SSWs.
  CORRECTION 2026-09-29 (review): the surface-precursor windows were first taken
  from the November-April series, so onsets before about 4 December had NaN
  precursors and 204 of 1,517 events were silently dropped from the PLUS set.
  The precursor now comes from the full-year annular-mode series with its own
  day-of-year climatology (days outside every event's -60..+75 influence zone),
  so both sets score the same events; n per set is recorded.

READOUTS for each set (base, plus) x tier (P1, P2, P3): CRPSS at SSWs (pooled,
model-bootstrap 95% interval; model-balanced; without CanESM5); ordinary-date and
weak-vortex nulls (200 draws each): mean, 2.5-97.5%, p = P(null >= SSW).

Output: results/current/6_predictability/forecast_value_robustness.json
"""
import json
import multiprocessing as mp
import sys
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "6_predictability"
sys.path.insert(0, str(HERE))
import forecast_value_test as FV                     # noqa: E402
import predictability_ceiling as P                   # noqa: E402
import predictability_wave_driving as WD             # noqa: E402
import within_model_check as W                       # noqa: E402

NAME = "forecast_value_robustness"
SEED = zlib.crc32(NAME.encode()) % (2 ** 32)
K = 200
N_BOOT = 2000
WORKERS = 16
WEAK_Q = 0.15
CAL_WIN = 30
TIERS = {"P1 pre-onset": 1, "P2 at-onset": 2, "P3 + post-onset stratosphere": 3}
_S = {}


def extra_features(c, on):
    """Wave-driving proxies and the pre-onset surface annular mode, windowed like P."""
    lat, plev, u, idx = c["lat"], c["plev"], c["u"], c["idx"]
    i60, i45, i75 = (int(np.argmin(np.abs(lat - x))) for x in (60.0, 45.0, 75.0))
    cap = lat >= 65.0; wcap = np.cos(np.deg2rad(lat[cap]))
    feat, lev = {}, {}
    for L in P.LEVELS_PA:
        j = int(np.argmin(np.abs(plev - L)))
        if not np.isclose(plev[j], L, rtol=0.1):
            continue
        s60 = pd.Series(u[:, j, i60], index=idx)
        scap = pd.Series(np.average(u[:, j, cap], weights=wcap, axis=1), index=idx)
        lev[L] = (s60, j)
        for src, tag in ((s60, "60N"), (scap, "cap")):
            zt = WD.doy_std(src.diff())
            for wn, win, _ in P.WINDOWS:
                feat[f"dudt_{P.LEVEL_LAB[L]}_{tag}_{wn}"] = P._wm(zt, on, win)
    if 1000.0 in lev and 10000.0 in lev:
        sh = WD.doy_std(lev[1000.0][0] - lev[10000.0][0])
        for wn, win, _ in P.WINDOWS:
            feat[f"shear10_100_{wn}"] = P._wm(sh, on, win)
    if 1000.0 in lev:
        j = lev[1000.0][1]
        gm = WD.doy_std(pd.Series(u[:, j, i75] - u[:, j, i45], index=idx))
        for wn, win, _ in P.WINDOWS:
            feat[f"geom75_45_{wn}"] = P._wm(gm, on, win)
    am = c["am_full"]
    ama = am - c["cl_full"].reindex(am.index.dayofyear).values
    for wn, win, t in P.WINDOWS:
        if t <= 2:                                   # a PRECURSOR: before/at onset only
            feat[f"amsurf_{wn}"] = P._wm(ama, on, win)
    return pd.DataFrame(feat)


def with_full_am(c, f):
    """Add the full-year annular mode and its climatology (outside influence zones)."""
    m = W.EP.load_member(Path(f))
    am = m["am"].dropna()
    msk = W.C6.influence_mask(am.index, c["on"])
    c["am_full"] = am
    c["cl_full"] = am[~msk].groupby(am[~msk].index.dayofyear).mean()
    return c


def build(c, on):
    base = W.H.features_for(c["u"], c["plev"], c["lat"], c["idx"], on)
    y = W.C6.anom(c["am"], on, c["cl"], W.OUT_WIN)
    ext = extra_features(c, on)
    return base, pd.concat([base.reset_index(drop=True), ext.reset_index(drop=True)], axis=1), y


def u10_60(c):
    j = int(np.argmin(np.abs(c["plev"] - 1000.0))); i = int(np.argmin(np.abs(c["lat"] - 60.0)))
    return pd.Series(c["u"][:, j, i], index=c["idx"])


def member_draws(task):
    """K ordinary and K weak-vortex pseudo draws for one member, own seeded stream."""
    i, f = task
    c = with_full_am(W.prepare_member(Path(f)), f)
    rng = np.random.default_rng(np.random.SeedSequence(SEED).spawn(i + 1)[i])
    n = len(c["on"])
    u = u10_60(c)
    cln = c["cln"][np.isin(c["cln"].month, (11, 12, 1, 2, 3))]
    uc = u.reindex(cln)
    weak = cln[(uc < np.nanquantile(uc, WEAK_Q)).values]
    wdoy = np.array([t.dayofyear for t in weak])
    out = {"ordinary": [], "weak": []}
    for _ in range(K):
        for kind in ("ordinary", "weak"):
            if kind == "ordinary":
                p_on = W.C6.draw_pseudo(c["cln"], c["dy"], rng)
            else:
                sel = []
                for d in rng.permutation(c["dy"]):
                    dd = np.abs(wdoy - d); dd = np.minimum(dd, 365 - dd)
                    ok = np.flatnonzero(dd <= CAL_WIN)
                    if len(ok):
                        sel.append(weak[ok[rng.integers(len(ok))]])
                p_on = pd.DatetimeIndex(sel)
            base, plus, y = build(c, p_on)
            pad = n - len(y)
            if pad > 0:
                base = pd.concat([base, pd.DataFrame(np.nan, index=range(pad), columns=base.columns)],
                                 ignore_index=True)
                plus = pd.concat([plus, pd.DataFrame(np.nan, index=range(pad), columns=plus.columns)],
                                 ignore_index=True)
                y = np.concatenate([y, np.full(pad, np.nan)])
            out[kind].append((base.values[:n], plus.values[:n], np.asarray(y)[:n]))
    return Path(f).stem, out, list(build(c, c["on"][:3])[1].columns)


def score_set(X, cal, y, g, balanced=False, drop=None):
    mdl = np.array([s.split("_")[0] for s in g])
    keep = np.ones(len(y), bool) if drop is None else mdl != drop
    c0, c1, m_te = FV.scores(FV.std_within(X[keep], g[keep]), cal[keep], y[keep], g[keep], mdl[keep])
    ok = np.isfinite(c0) & np.isfinite(c1)
    if not balanced:
        return FV.crpss(c0[ok], c1[ok]), c0, c1, m_te
    w = np.zeros(len(c0))
    for m in np.unique(m_te):
        i = (m_te == m) & ok
        w[i] = 1.0 / i.sum()
    return float(1 - np.sum(w * np.nan_to_num(c1)) / np.sum(w * np.nan_to_num(c0))), c0, c1, m_te


def _null_task(args):
    kind, k = args
    out = {}
    for sname, si in (("base", 0), ("plus", 1)):
        for tier, cols in _S["tier_cols"][sname].items():
            X = np.concatenate([_S["D"][m][kind][k][si][:, cols] for m in _S["members"]])
            cal = np.concatenate([_S["D"][m][kind][k][si][:, _S["cal"][sname]] for m in _S["members"]])
            y = np.concatenate([_S["D"][m][kind][k][2] for m in _S["members"]])
            g = np.concatenate([np.full(len(_S["D"][m][kind][k][2]), m) for m in _S["members"]])
            out[(sname, tier, "pooled")] = score_set(X, cal, y, g)[0]
            out[(sname, tier, "balanced")] = score_set(X, cal, y, g, balanced=True)[0]
    return kind, k, out


def main():
    files = sorted(W.RAW.glob("*_zm.nc"))
    mem, Xb, Xp, Y, G = [], [], [], [], []
    for f in files:
        c = W.prepare_member(f)
        if c is None:
            continue
        c = with_full_am(c, f)
        base, plus, y = build(c, c["on"])
        Xb.append(base.values); Xp.append(plus.values); Y.append(y); G.append(np.full(len(y), f.stem))
        mem.append(f); cb, cp = list(base.columns), list(plus.columns)
    Xb, Xp, y, g = np.concatenate(Xb), np.concatenate(Xp), np.concatenate(Y), np.concatenate(G)
    cal = {"base": [cb.index("doy_sin"), cb.index("doy_cos")], "plus": [cp.index("doy_sin"), cp.index("doy_cos")]}
    tier_cols = {s: {t: [cols.index(x) for x in P.tier_cols(cols, n)] for t, n in TIERS.items()}
                 for s, cols in (("base", cb), ("plus", cp))}
    XX = {"base": Xb, "plus": Xp}
    print(f"{np.isfinite(y).sum()} events, {len(mem)} members; base {len(cb)} / plus {len(cp)} columns", flush=True)

    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(WORKERS, mp_context=ctx) as ex:
        D = {stem: out for stem, out, _ in ex.map(member_draws, list(enumerate(map(str, mem))))}
    _S.update(D=D, members=[f.stem for f in mem], tier_cols=tier_cols, cal=cal)
    tasks = [(kind, k) for kind in ("ordinary", "weak") for k in range(K)]
    with ProcessPoolExecutor(WORKERS, mp_context=ctx) as ex:
        nulls = list(ex.map(_null_task, tasks))

    rng = np.random.default_rng(SEED)
    res = {"plan_approved": "2026-09-29", "seed": SEED, "k_draws": K, "n_boot": N_BOOT,
           "weak_quantile": WEAK_Q, "n_events": int(np.isfinite(y).sum()), "n_members": len(mem),
           "columns": {"base": len(cb), "plus": len(cp)}, "results": {}}
    # predictor spread at SSWs vs nulls (raw columns, tier P2, base set)
    ci2 = tier_cols["base"]["P2 at-onset"]
    sd_ssw = np.nanstd(Xb[:, ci2], axis=0)
    spread = {}
    for kind in ("ordinary", "weak"):
        sds = np.nanmean([np.nanstd(np.concatenate([D[m][kind][k][0] for m in _S["members"]])[:, ci2], axis=0)
                          for k in range(0, K, 10)], axis=0)
        spread[kind] = round(float(np.nanmedian(sd_ssw / sds)), 4)
    res["predictor_sd_ratio_ssw_over_null_median"] = spread
    for sname in ("base", "plus"):
        for tier, cols in tier_cols[sname].items():
            X = XX[sname][:, cols]; C = XX[sname][:, cal[sname]]
            real, c0, c1, m_te = score_set(X, C, y, g)
            ok = np.isfinite(c0) & np.isfinite(c1)
            um = np.unique(m_te[ok]); boot = []
            for _ in range(N_BOOT):
                idx = np.concatenate([np.flatnonzero((m_te == u) & ok) for u in rng.choice(um, len(um))])
                boot.append(FV.crpss(c0[idx], c1[idx]))
            bal = score_set(X, C, y, g, balanced=True)[0]
            nocan = score_set(X, C, y, g, drop="CanESM5")[0]
            r = {"n_events_scored": int(np.isfinite(c0).sum()), "CRPSS_ssw": round(real, 4),
                 "CRPSS_ssw_CI95_model_bootstrap": [round(float(q), 4) for q in np.percentile(boot, [2.5, 97.5])],
                 "CRPSS_ssw_model_balanced": round(bal, 4), "CRPSS_ssw_without_CanESM5": round(nocan, 4)}
            for kind in ("ordinary", "weak"):
                for agg in ("pooled", "balanced"):
                    v = np.array([o[(sname, tier, agg)] for kd, k, o in nulls if kd == kind])
                    ref = real if agg == "pooled" else bal
                    r[f"null_{kind}_{agg}"] = {"mean": round(float(v.mean()), 4),
                                               "q025_q975": [round(float(q), 4) for q in np.percentile(v, [2.5, 97.5])],
                                               "p_null_ge_ssw": round(float(np.mean(v >= ref)), 4)}
            res["results"][f"{sname} | {tier}"] = r
            print(f"{sname:4s} {tier:30s} SSW {real:+.4f} {r['CRPSS_ssw_CI95_model_bootstrap']} bal {bal:+.4f} "
                  f"noCan {nocan:+.4f} | ordinary {r['null_ordinary_pooled']['mean']:+.4f} "
                  f"p={r['null_ordinary_pooled']['p_null_ge_ssw']} | weak {r['null_weak_pooled']['mean']:+.4f} "
                  f"p={r['null_weak_pooled']['p_null_ge_ssw']}", flush=True)
    print("predictor sd ratio SSW/null (median, P2 base):", spread)
    (RESULTS / "forecast_value_robustness.json").write_text(json.dumps(res, indent=2), encoding="utf8",
                                                            newline="\n")
    print("Saved -> forecast_value_robustness.json")


if __name__ == "__main__":
    main()
