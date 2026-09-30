#!/usr/bin/env python3
"""
obs_model_compatibility.py
==========================
IS THE OBSERVED RECORD A TYPICAL 43-EVENT DRAW FROM THE CMIP6 SSW POPULATION, AND
DOES THE OBSERVATIONAL PIPELINE RECOVER THE TRUTH WHEN A MODEL PLAYS REALITY?

Plan approved 2026-09-30. Written before it was run.

DATA (as is_downward_propagation_a_class.py)
  Observations: CPC AO daily, primary catalogue (43 events), per-event mean over
  days +8..+52, anomalies from the event-free day-of-year climatology; 300
  day-of-year-matched event-free sets (gate3_clean_null).
  CMIP6: each member's annular-mode index, Charlton-Polvani events (1,517), the
  same outcome; per member 20 event-free sets (stratifier_law_cmip6.draw_pseudo).
  The observed index is an EOF of 1000 hPa height and the CMIP6 one of zonal-mean
  sea-level pressure; both are standardised; this is stated as a limit.

STATISTICS (the observational pipeline, computed identically on any event set)
  T1 pure-shift KS statistic (events against the event-free pool plus the shift)
  T2 variance ratio, events / event-free pool
  T3 implied shift, mean(events) - mean(pool)
  T4 downward rate: conditions 1-2 of the criterion on the same series over days
     8-52 (mean negative, more than half of days negative)
  T5 share of the DW-NDW contrast (same outcome) reproduced by event-free sets

A. COMPATIBILITY. 10,000 draws of 43 events without replacement from the pooled
   CMIP6 events (with the event-free sets of the members drawn) and, for each
   model with at least 43 events, 2,000 draws from that model. For each statistic,
   the two-sided percentile of the observed value; combined by the minimum of the
   five two-sided p values, calibrated on the same draws. Reading: p < 0.05 means
   the observed record is not a typical CMIP6 draw on that statistic. Passing
   shows compatibility, which is weaker than exchangeability, and is described so.

B. PERFECT MODEL. Each model in turn is "reality": 1,000 draws of 43 of its events,
   the observational pipeline applied, and the answers compared with that model's
   full-population values:
   PM1 size of the pure-shift KS test at 5% (rejection rate);
   PM2 how often the observed downward rate of the draw falls outside the 95%
       interval of shifted event-free sets (the two-thirds test's error rate);
   PM3 coverage of the model's full-population variance ratio by the 43-event
       bootstrap 95% interval (2,000 resamples).
   Pass, fixed now: PM1 and PM2 rates within [0.01, 0.10]; PM3 coverage >= 0.90.

IMPLEMENTATION NOTES (fixed before the run)
  T4 uses the daily anomaly series (full record / full year, so March windows are
  complete; an event is classified when 90% of its window days exist). T5's
  event-free contrast is the contrast within the whole event-free pool of the
  data set (its population value), the same definition for observations and every
  draw; it therefore differs slightly from the published-criterion audit's
  set-average. For pooled CMIP6 draws the reference pool is all CMIP6 event-free
  events; for per-model draws, that model's. PM2 uses 200 shifted event-free sets
  of the draw's size from the model's pool, shifted by the draw's own shift.
  T3 is expressed in units of the data set's own event-free s.d. (the observed AO
  and the CMIP6 annular mode are standardised differently, so raw differences of
  means are not comparable); T1, T2, T4 and T5 are scale-free.

Output: results/current/6_predictability/obs_model_compatibility.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "6_predictability"
PHYS = ROOT / "ssw-design-analysis" / "07_physical_decomposition"
AUD = ROOT / "ssw-design-analysis" / "08_literature_audit"
for d in (HERE, PHYS, AUD, HERE.parents[0] / "02_event_catalogues", HERE.parents[0] / "05_corrected_estimators"):
    sys.path.insert(0, str(d))
from build_catalogue import load_catalogue          # noqa: E402
import gate3_clean_null as G                        # noqa: E402
import multi_index_event_study as M                 # noqa: E402
import selection_on_outcome as SO                   # noqa: E402
import ensemble_precursor as EP                     # noqa: E402
import stratifier_law_cmip6 as C6                   # noqa: E402
import criterion_sweep as CS                        # noqa: E402  (window matrices)

NAME = "obs_model_compatibility"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
RAW = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "raw" / "cmip6"
OUT_WIN = (8, 52)
SEASON = (11, 12, 1, 2, 3, 4)
N_OBS_SETS = 300
N_MEMBER_SETS = 20
N_POOLED, N_PER_MODEL, N_PM = 10000, 2000, 1000
N_PM_SETS, N_PM_BOOT = 200, 2000
K = 43
STATS = ("T1_ks", "T2_var_ratio", "T3_shift", "T4_rate_c12", "T5_share_reproduced")


def c12(Mw):
    """Conditions 1-2 over days 8-52 on daily anomalies (tau 0, f 0.5)."""
    lab, ok = CS.classify(Mw, None, 52, 0.0, 0.5, False)
    return lab, ok


def contrast(y, lab, ok):
    u = ok & np.isfinite(y)
    if (lab & u).sum() < 3 or ((~lab) & u).sum() < 3:
        return np.nan
    return float(y[lab & u].mean() - y[(~lab) & u].mean())


class Pool:
    """An event-free reference population: outcomes, labels, pool contrast."""
    def __init__(self, y, Mw):
        ok = np.isfinite(y)
        self.y, self.M = y[ok], Mw[ok]
        self.lab, self.ok = c12(self.M)
        self.contrast = contrast(self.y, self.lab, self.ok)
        self.mean, self.var = float(self.y.mean()), float(self.y.var(ddof=1))


def statistics(y, Mw, pool):
    ok = np.isfinite(y)
    y, Mw = y[ok], Mw[ok]
    sh = float(y.mean() - pool.mean)
    lab, okl = c12(Mw)
    c = contrast(y, lab, okl)
    return {"T1_ks": float(stats.ks_2samp(y, pool.y + sh).statistic),
            "T2_var_ratio": float(y.var(ddof=1) / pool.var),
            "T3_shift": sh / np.sqrt(pool.var),
            "T4_rate_c12": float(lab[okl].mean()),
            "T5_share_reproduced": float(pool.contrast / c) if np.isfinite(c) and c != 0 else np.nan}


def two_sided(draws, v):
    d = np.asarray(draws, float); d = d[np.isfinite(d)]
    return float(min(1.0, 2 * min(np.mean(d <= v), np.mean(d >= v))))


def load_obs(rng):
    ao = M.load("ao")["y"]
    real = load_catalogue("primary")
    real = real[(real >= ao.index.min()) & (real <= ao.index.max())]
    mask = G.real_influence_mask(ao.index, real)
    clim = ao[~mask].groupby(ao[~mask].index.dayofyear).mean()
    an = ao - clim.reindex(ao.index.dayofyear).values
    clean = G.zone_free_index(ao.index, real)
    by = G.build_by_doy(clean)
    doys = np.array([t.dayofyear for t in pd.DatetimeIndex(real)])
    ys, ms = [], []
    for _ in range(N_OBS_SETS):
        f = G.draw_clean(clean, doys, rng, by)
        if len(f) >= 10:
            ys.append(SO.per_event(ao, f, clim)); ms.append(CS.window_matrix(an, f))
    return (SO.per_event(ao, real, clim), CS.window_matrix(an, real),
            Pool(np.concatenate(ys), np.vstack(ms)))


def load_cmip6(rng):
    ev_y, ev_M, ev_mod, ps_y, ps_M, ps_mod = [], [], [], [], [], []
    for f in sorted(RAW.glob("*_zm.nc")):
        try:
            full = EP.load_member(f)
        except Exception:
            continue
        m = full[np.isin(full.index.month, SEASON)].dropna()
        if len(m) < 2000:
            continue
        on = EP.detect_ssw(full["u10"].values, full.index)
        if len(on) < 15:
            continue
        am = m["am"]
        msk = C6.influence_mask(am.index, on)
        cl = am[~msk].groupby(am[~msk].index.dayofyear).mean()
        cln = C6.zone_free_index(am.index, on)
        dy = np.array([t.dayofyear for t in pd.DatetimeIndex(on)])
        an = full["am"] - cl.reindex(full["am"].index.dayofyear).values
        mod = f.stem.split("_")[0]
        ev_y.append(C6.anom(am, on, cl, OUT_WIN)); ev_M.append(CS.window_matrix(an, on))
        ev_mod += [mod] * len(on)
        for _ in range(N_MEMBER_SETS):
            p = C6.draw_pseudo(cln, dy, rng)
            if len(p) >= 10:
                ps_y.append(C6.anom(am, p, cl, OUT_WIN)); ps_M.append(CS.window_matrix(an, p))
                ps_mod += [mod] * len(p)
    return (np.concatenate(ev_y), np.vstack(ev_M), np.array(ev_mod),
            np.concatenate(ps_y), np.vstack(ps_M), np.array(ps_mod))


def main():
    rng = np.random.default_rng(SEED)
    yo, Mo, pool_o = load_obs(rng)
    obs = statistics(yo, Mo, pool_o)
    print(f"observations: {int(np.isfinite(yo).sum())} events; {obs}", flush=True)
    ey, eM, emod, py, pM, pmod = load_cmip6(rng)
    fin = np.isfinite(ey)
    ey, eM, emod = ey[fin], eM[fin], emod[fin]
    pool_all = Pool(py, pM)
    full = statistics(ey, eM, pool_all)
    print(f"CMIP6: {len(ey)} events; full population {full}", flush=True)
    res = {"plan_approved": "2026-09-30", "seed": SEED, "k_events_per_draw": K,
           "observations": obs, "cmip6_full_population": full, "n_cmip6_events": int(len(ey)),
           "compatibility": {}, "perfect_model": {}}

    def compat(label, idx_pool, pool, n_draws):
        draws = {s: [] for s in STATS}
        for _ in range(n_draws):
            i = rng.choice(idx_pool, K, replace=False)
            st = statistics(ey[i], eM[i], pool)
            for s in STATS:
                draws[s].append(st[s])
        draws = {s: np.array(v) for s, v in draws.items()}
        ps = {s: two_sided(draws[s], obs[s]) for s in STATS}
        # calibrate the minimum p on the draws themselves
        mins = np.array([min(two_sided(draws[s], draws[s][b]) for s in STATS)
                         for b in range(min(n_draws, 2000))])
        pmin = min(ps.values())
        out = {"n_draws": n_draws, "p_two_sided": {s: round(v, 4) for s, v in ps.items()},
               "draw_q025_q975": {s: [round(float(q), 4) for q in np.nanpercentile(draws[s], [2.5, 97.5])]
                                  for s in STATS},
               "min_p": round(pmin, 4), "p_combined": round(float(np.mean(mins <= pmin)), 4)}
        print(f"  {label}: {out['p_two_sided']} combined {out['p_combined']}", flush=True)
        return out

    res["compatibility"]["pooled"] = compat("pooled", np.arange(len(ey)), pool_all, N_POOLED)
    models = sorted(set(emod))
    for mdl in models:
        idx = np.flatnonzero(emod == mdl)
        if len(idx) < K:
            res["compatibility"][mdl] = {"n_events": int(len(idx)), "skipped": "fewer than 43 events"}
            continue
        pm = Pool(py[pmod == mdl], pM[pmod == mdl])
        res["compatibility"][mdl] = {"n_events": int(len(idx)), **compat(mdl, idx, pm, N_PER_MODEL)}

    # ---------------- perfect model
    for mdl in models:
        idx = np.flatnonzero(emod == mdl)
        if len(idx) < K:
            continue
        pm = Pool(py[pmod == mdl], pM[pmod == mdl])
        truth = statistics(ey[idx], eM[idx], pm)
        ks_rej, pm2_out, pm3_cov = 0, 0, 0
        for _ in range(N_PM):
            i = rng.choice(idx, K, replace=False)
            y, Mw = ey[i], eM[i]
            sh = float(y.mean() - pm.mean)
            ks_rej += stats.ks_2samp(y, pm.y + sh).pvalue < 0.05
            lab, okl = c12(Mw)
            r_obs = lab[okl].mean()
            rates = []
            for _ in range(N_PM_SETS):
                j = rng.integers(0, len(pm.y), K)
                lb, ok = c12(pm.M[j] + sh)
                rates.append(lb[ok].mean())
            lo, hi = np.percentile(rates, [2.5, 97.5])
            pm2_out += not (lo <= r_obs <= hi)
            b = rng.integers(0, K, (N_PM_BOOT, K))
            vr = y[b].var(axis=1, ddof=1) / pm.var
            lo, hi = np.percentile(vr, [2.5, 97.5])
            pm3_cov += lo <= truth["T2_var_ratio"] <= hi
        r = {"n_events": int(len(idx)), "full_population": {k: round(v, 4) for k, v in truth.items()},
             "PM1_ks_rejection_rate": round(ks_rej / N_PM, 4),
             "PM2_two_thirds_error_rate": round(pm2_out / N_PM, 4),
             "PM3_var_ratio_coverage": round(pm3_cov / N_PM, 4)}
        r["pass"] = bool(0.01 <= r["PM1_ks_rejection_rate"] <= 0.10 and 0.01 <= r["PM2_two_thirds_error_rate"] <= 0.10
                         and r["PM3_var_ratio_coverage"] >= 0.90)
        res["perfect_model"][mdl] = r
        print(f"  perfect model {mdl} (n={len(idx)}): {r}", flush=True)
    (RESULTS / "obs_model_compatibility.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print("Saved -> obs_model_compatibility.json")


if __name__ == "__main__":
    main()
