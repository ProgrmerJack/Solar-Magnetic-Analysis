#!/usr/bin/env python3
"""
criterion_sweep.py
==================
DOES THE "DOWNWARD" RATE AND CONTRAST TRACK ONE SHIFTED POPULATION FOR EVERY
REASONABLE VERSION OF THE CRITERION, OR ONLY FOR THE ONE WE CHOSE?

Plan approved 2026-09-30. Written before the sweep was run.

QUESTION
  era5_recompute_and_two_thirds.py showed, for the Karpechko et al. (2017)
  criterion as published, that event-free dates displaced by the measured surface
  shift reproduce the observed downward rate, and that event-free dates reproduce
  97% of the class contrast. A referee can answer that the criterion was badly
  tuned. If the label reflected two kinds of event, some versions of the criterion
  would separate SSWs from the shifted population (a plateau: a rate insensitive
  to the threshold, or a contrast the null cannot reach). If it is a threshold on
  one shifted population, every version should track the shifted null smoothly.

GRID (fixed now; 108 specifications in ERA5)
  classification window days 8..W, W in {25, 35, 52}
  level of conditions 1-2: NAM 1000 hPa or 850 hPa (ERA5, era5_nam_daily.parquet)
  condition 1 threshold: window mean below tau, tau in {0, -0.25, -0.5} (sigma)
  condition 2: share of days negative above f, f in {0.5, 0.6, 0.7}
  condition 3 (150 hPa NAM negative on more than 70% of days): off / on
  The published criterion is (52, 1000, 0, 0.5, on).

EVENTS AND NULL (as era5_recompute_and_two_thirds.py)
  The 39 primary-catalogue events inside the ERA5 NAM record; event-free dates
  from gate3_clean_null (more than 135 days from every catalogued onset), drawn
  day-of-year matched, 400 sets. Shift per specification: mean of the window-mean
  NAM at events minus that at event-free dates (each level measured on its own;
  with condition 3 on, the 150 hPa field is shifted by its own measured shift).
  Outcome for the contrast: the per-event 1000 hPa NAM of selection_on_outcome
  (its outcome window), as in the published-criterion audit.

STATISTICS PER SPECIFICATION
  rate at SSWs; rate of shifted event-free sets (mean, 2.5-97.5%);
  z = (rate_SSW - mean) / s.d. over shifted sets;
  contrast DW - NDW at SSWs and the share reproduced by unshifted event-free sets.

TESTS
  S1 (primary, conditions 1-2 only, 54 specifications): max |z| over
     specifications, referred to the same statistic computed for each shifted
     event-free set against the others (the joint null, correlation between
     specifications kept); p = P(null max|z| >= observed).
  S2 (primary): mean z over the 54 (a systematic excess), same reference.
  S3 (secondary): S1 and S2 with condition 3 on (54 specifications); the
     stratospheric shift is not a surface shift, so a failure here speaks to the
     stratospheric condition, not the surface classification.
  Descriptive: the share of the contrast reproduced, per specification.
  Coverage rule (implementation): an event is classified when at least 90% of its
  window days exist (the published audit's "at least 20 days" would drop every
  18-day window at W = 25); ERA5 NAM is complete, so this binds only in CMIP6.
  CMIP6 anomalies: each member's annular mode minus its event-free day-of-year
  mean, divided by the s.d. of those daily anomalies on event-free days, so tau is
  in sigma as in ERA5; the full-year series is used so March windows are complete.
  CMIP6 (secondary): conditions 1-2 on each member's annular mode, W x tau x f
     (27 specifications), 1,517 events, 20 pseudo sets per member; S1 and S2.
  Reading, fixed now: S1 or S2 at p < 0.05 in ERA5 means some version of the
  criterion separates SSWs from one shifted population, and the paper says so.

POST-HOC CALIBRATION (added 2026-09-30 after a code review; labelled so in the
JSON): the registered S1/S2 reference scores each event-free set against the
others at the SSWs' measured shift, but the SSWs' own window mean equals that
shift by construction, so the observed z is squeezed toward 0 and the test cannot
reject (reviewer's simulation: 0/200 rejections under the null). The calibrated
reference replays the whole procedure with each event-free set b playing the SSWs
(displaced by the measured shift, then its own shift re-estimated against the
other sets, exactly as for the real events); p = share of those replays with a
statistic at least as extreme. Power: the same replay with a planted two-class
structure (two thirds of the displaced set shifted by a further -0.6 sigma, one
third by +1.2 sigma, keeping the mean shift), detection at calibrated p < 0.05.

Output: results/current/9_literature/criterion_sweep.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "9_literature"
SIM = ROOT / "ssw-design-analysis" / "06_simulation_validation"
PHYS = ROOT / "ssw-design-analysis" / "07_physical_decomposition"
for d in (SIM, PHYS, HERE.parents[0] / "02_event_catalogues", HERE.parents[0] / "05_corrected_estimators"):
    sys.path.insert(0, str(d))
from build_catalogue import load_catalogue          # noqa: E402
import gate3_clean_null as G                        # noqa: E402
import selection_on_outcome as SO                   # noqa: E402

NAME = "criterion_sweep"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
ERA5 = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "era5_nam_daily.parquet"
RAW = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "raw" / "cmip6"
LO, HI = 8, 52
WS, TAUS, FS = (25, 35, 52), (0.0, -0.25, -0.5), (0.5, 0.6, 0.7)
LEVELS = ("nam_1000", "nam_850")
N_SETS = 400
N_SETS_CMIP6 = 20


def window_matrix(series, onsets):
    """Rows: onsets; columns: days LO..HI after onset (NaN outside the record)."""
    idx = pd.DatetimeIndex(series.index)
    pos = idx.get_indexer(pd.DatetimeIndex(onsets).normalize())
    v = series.values.astype(float)
    M = np.full((len(pos), HI - LO + 1), np.nan)
    for j, p in enumerate(pos):
        if p < 0:
            continue
        a, b = p + LO, p + HI + 1
        seg = v[max(a, 0):min(b, len(v))]
        M[j, max(0, -a):max(0, -a) + len(seg)] = seg
    return M


def classify(Ms, Mt, W, tau, f, c3, ds=0.0, dt=0.0):
    """Conditions 1-2 (and 3) on window days LO..W. Returns (label, usable)."""
    k = W - LO + 1
    a = Ms[:, :k] + ds
    n = np.isfinite(a).sum(1)
    ok = n >= 0.9 * k
    with np.errstate(invalid="ignore", divide="ignore"):
        c1 = np.nanmean(np.where(np.isfinite(a), a, np.nan), axis=1) < tau
        c2 = (a < 0).sum(1) / np.maximum(n, 1) > f
        lab = c1 & c2
        if c3:
            t = Mt[:, :k] + dt
            nt = np.isfinite(t).sum(1)
            c3v = np.where(nt >= 0.9 * k, (t < 0).sum(1) / np.maximum(nt, 1) > 0.7, True)
            lab = lab & c3v
    return lab, ok


def wmean(M, W):
    k = W - LO + 1
    with np.errstate(invalid="ignore"):
        return np.nanmean(M[:, :k], axis=1)


def specs(levels, c3s):
    return [(W, L, tau, f, c3) for W in WS for L in levels for tau in TAUS for f in FS for c3 in c3s]


def joint_test(obs, sets):
    """obs: rates per spec; sets: (n_sets, n_spec). max|z| and mean z against the
    same statistics of each set referred to the others (leave-one-out moments)."""
    n = len(sets)
    mu, sd = sets.mean(0), sets.std(0, ddof=1)
    z_obs = (obs - mu) / np.where(sd > 0, sd, np.nan)
    zb = []
    for b in range(n):
        o = np.delete(sets, b, axis=0)
        m, s = o.mean(0), o.std(0, ddof=1)
        zb.append((sets[b] - m) / np.where(s > 0, s, np.nan))
    zb = np.array(zb)
    mx_o, mx_b = np.nanmax(np.abs(z_obs)), np.nanmax(np.abs(zb), axis=1)
    mz_o, mz_b = np.nanmean(z_obs), np.nanmean(zb, axis=1)
    return {"n_specs": int(len(obs)), "n_sets": int(n),
            "max_abs_z": round(float(mx_o), 3), "p_max_abs_z": round(float(np.mean(mx_b >= mx_o)), 4),
            "mean_z": round(float(mz_o), 3), "p_mean_z_excess": round(float(np.mean(mz_b >= mz_o)), 4),
            "p_mean_z_deficit": round(float(np.mean(mz_b <= mz_o)), 4)}, z_obs


def era5_arm(rng):
    e = pd.read_parquet(ERA5)
    allev = load_catalogue("primary")
    real = allev[(allev >= e.index.min() + pd.Timedelta(days=70))
                 & (allev <= e.index.max() - pd.Timedelta(days=70))]
    mask = G.real_influence_mask(e.index, allev)
    clean = G.zone_free_index(e.index, allev)
    by = G.build_by_doy(clean)
    doys = np.array([t.dayofyear for t in pd.DatetimeIndex(real)])
    clim = e["nam_1000"][~mask].groupby(e["nam_1000"][~mask].index.dayofyear).mean()
    sets = [G.draw_clean(clean, doys, rng, by) for _ in range(N_SETS)]
    Mr = {c: window_matrix(e[c], real) for c in (*LEVELS, "nam_150")}
    Mp = [{c: window_matrix(e[c], s) for c in (*LEVELS, "nam_150")} for s in sets]
    Yr = SO.per_event(e["nam_1000"], real, clim)
    Yp = [SO.per_event(e["nam_1000"], s, clim) for s in sets]
    # measured shifts per (window, level): events minus pooled event-free dates
    shift = {}
    for W in WS:
        for c in (*LEVELS, "nam_150"):
            pool = np.concatenate([wmean(m[c], W) for m in Mp])
            shift[(W, c)] = float(np.nanmean(wmean(Mr[c], W)) - np.nanmean(pool))
    out, rows = {}, []
    sp = specs(LEVELS, (False, True))
    obs = np.zeros(len(sp)); sh = np.zeros((N_SETS, len(sp)))
    for j, (W, L, tau, f, c3) in enumerate(sp):
        lab, ok = classify(Mr[L], Mr["nam_150"], W, tau, f, c3)
        obs[j] = lab[ok].mean()
        cons = []
        for b, m in enumerate(Mp):
            lb, okb = classify(m[L], m["nam_150"], W, tau, f, c3, shift[(W, L)],
                               shift[(W, "nam_150")] if c3 else 0.0)
            sh[b, j] = lb[okb].mean()
            l0, ok0 = classify(m[L], m["nam_150"], W, tau, f, c3)
            y = Yp[b]; u = ok0 & np.isfinite(y)
            if (l0 & u).sum() >= 3 and ((~l0) & u).sum() >= 3:
                cons.append(y[l0 & u].mean() - y[(~l0) & u].mean())
        u = ok & np.isfinite(Yr)
        c_real = float(Yr[lab & u].mean() - Yr[(~lab) & u].mean()) if (lab & u).sum() >= 3 and ((~lab) & u).sum() >= 3 else np.nan
        rows.append({"W": W, "level": L, "tau": tau, "f": f, "cond3": c3,
                     "n_events": int(ok.sum()), "rate_ssw": round(float(obs[j]), 4),
                     "rate_shifted_mean": round(float(sh[:, j].mean()), 4),
                     "rate_shifted_q025_q975": [round(float(q), 4) for q in np.percentile(sh[:, j], [2.5, 97.5])],
                     "shift_surface": round(shift[(W, L)], 4),
                     "contrast_ssw": None if not np.isfinite(c_real) else round(c_real, 4),
                     "contrast_eventfree_mean": round(float(np.mean(cons)), 4) if cons else None,
                     "share_reproduced_pct": round(float(100 * np.mean(cons) / c_real), 1)
                     if cons and np.isfinite(c_real) and abs(c_real) > 1e-9 else None})
    prim = np.array([not r["cond3"] for r in rows])
    s1, z1 = joint_test(obs[prim], sh[:, prim])
    s3, z3 = joint_test(obs[~prim], sh[:, ~prim])
    for r, z in zip([r for r in rows if not r["cond3"]], z1):
        r["z"] = round(float(z), 3)
    for r, z in zip([r for r in rows if r["cond3"]], z3):
        r["z"] = round(float(z), 3)
    # ---- post-hoc calibration (see docstring)
    rcal = np.random.default_rng(zlib.crc32(b"criterion_sweep|calibration"))
    N_CAL = 200

    def replay(b, planted):
        """Set b plays the SSWs: displaced by the measured shift (optionally with
        a planted two-class split), own shift re-estimated against the others."""
        m = Mp[b]
        cls = rcal.uniform(size=m["nam_1000"].shape[0]) < 2 / 3
        extra = np.where(cls, -0.6, 1.2)[:, None] if planted else 0.0
        others = [k for k in range(N_SETS) if k != b][:N_CAL]
        stats_ = {}
        for part, sel in (("c12", prim), ("c3", ~prim)):
            obs_b, sh_b = [], []
            for j in np.flatnonzero(sel):
                W, L, tau, f, c3 = sp[j]
                Ms = m[L] + shift[(W, L)] + extra
                Mt = m["nam_150"] + (shift[(W, "nam_150")] if c3 else 0.0)
                own = float(np.nanmean(wmean(Ms, W)) - np.nanmean(np.concatenate([wmean(Mp[k][L], W) for k in others])))
                own_t = float(np.nanmean(wmean(Mt, W)) - np.nanmean(np.concatenate([wmean(Mp[k]["nam_150"], W) for k in others])))
                lab, ok = classify(Ms, Mt, W, tau, f, c3)
                obs_b.append(lab[ok].mean())
                col = []
                for k in others:
                    lk, okk = classify(Mp[k][L], Mp[k]["nam_150"], W, tau, f, c3, own, own_t if c3 else 0.0)
                    col.append(lk[okk].mean())
                sh_b.append(col)
            obs_b, sh_b = np.array(obs_b), np.array(sh_b).T
            mu, sd = sh_b.mean(0), sh_b.std(0, ddof=1)
            z = (obs_b - mu) / np.where(sd > 0, sd, np.nan)
            stats_[part] = (float(np.nanmax(np.abs(z))), float(np.nanmean(z)),
                            float(np.mean((obs_b < np.percentile(sh_b, 2.5, axis=0))
                                          | (obs_b > np.percentile(sh_b, 97.5, axis=0)))))
        return stats_

    null = [replay(b, False) for b in range(N_CAL)]
    plant = [replay(b, True) for b in range(N_CAL)]
    calib = {"label": "post hoc, added after code review", "n_replays": N_CAL}
    for part, s_obs in (("c12", s1), ("c3", s3)):
        mx = np.array([r[part][0] for r in null]); mz = np.array([r[part][1] for r in null])
        fo = np.array([r[part][2] for r in null])
        pmx = lambda v: float(np.mean(mx >= v))
        pmz = lambda v: float(np.mean(np.abs(mz - np.mean(mz)) >= abs(v - np.mean(mz))))
        calib[part] = {"p_max_abs_z_calibrated": round(pmx(s_obs["max_abs_z"]), 4),
                       "p_mean_z_calibrated_two_sided": round(pmz(s_obs["mean_z"]), 4),
                       "null_max_abs_z_q95": round(float(np.percentile(mx, 95)), 3),
                       "null_share_outside_95_mean": round(float(fo.mean()), 4),
                       "power_planted_classes_max_abs_z": round(float(np.mean([pmx(r[part][0]) < 0.05 for r in plant])), 3),
                       "power_planted_classes_mean_z": round(float(np.mean([pmz(r[part][1]) < 0.05 for r in plant])), 3),
                       "planted_share_outside_95_mean": round(float(np.mean([r[part][2] for r in plant])), 4)}
    pub = [r for r in rows if (r["W"], r["level"], r["tau"], r["f"], r["cond3"]) == (52, "nam_1000", 0.0, 0.5, True)][0]
    return {"n_events": int(len(real)), "n_sets": N_SETS, "S1_S2_conditions_1_2": s1,
            "calibration_posthoc": calib,
            "S3_with_condition_3": s3, "published_spec": pub,
            "share_reproduced_range_pct": [min(r["share_reproduced_pct"] for r in rows if r["share_reproduced_pct"] is not None),
                                           max(r["share_reproduced_pct"] for r in rows if r["share_reproduced_pct"] is not None)],
            "specs": rows}


def cmip6_arm(rng):
    import ensemble_precursor as EP
    import stratifier_law_cmip6 as C6
    Mr, Mp = [], [[] for _ in range(N_SETS_CMIP6)]
    n_mem = 0
    for fpath in sorted(RAW.glob("*_zm.nc")):
        try:
            full = EP.load_member(fpath)
        except Exception:
            continue
        if len(full[np.isin(full.index.month, (11, 12, 1, 2, 3, 4))].dropna()) < 2000:
            continue
        on = EP.detect_ssw(full["u10"].values, full.index)
        if len(on) < 15:
            continue
        am = full["am"]
        msk = C6.influence_mask(am.index, on)
        cl = am[~msk].groupby(am[~msk].index.dayofyear).mean()
        an = am - cl.reindex(am.index.dayofyear).values
        an = an / np.nanstd(an[~msk].values)
        cln = C6.zone_free_index(am.index, on)
        cln = cln[np.isin(pd.DatetimeIndex(cln).month, (11, 12, 1, 2, 3, 4))]
        dy = np.array([t.dayofyear for t in pd.DatetimeIndex(on)])
        Mr.append(window_matrix(an, on))
        for k in range(N_SETS_CMIP6):
            Mp[k].append(window_matrix(an, C6.draw_pseudo(cln, dy, rng)))
        n_mem += 1
    Mr = np.vstack(Mr); Mp = [np.vstack(m) for m in Mp]
    shift = {W: float(np.nanmean(wmean(Mr, W)) - np.nanmean(np.concatenate([wmean(m, W) for m in Mp])))
             for W in WS}
    sp = [(W, tau, f) for W in WS for tau in TAUS for f in FS]
    obs = np.zeros(len(sp)); sh = np.zeros((N_SETS_CMIP6, len(sp))); rows = []
    for j, (W, tau, f) in enumerate(sp):
        lab, ok = classify(Mr, None, W, tau, f, False)
        obs[j] = lab[ok].mean()
        for k, m in enumerate(Mp):
            lb, okb = classify(m, None, W, tau, f, False, shift[W])
            sh[k, j] = lb[okb].mean()
        rows.append({"W": W, "tau": tau, "f": f, "n_events": int(ok.sum()),
                     "rate_ssw": round(float(obs[j]), 4),
                     "rate_shifted_mean": round(float(sh[:, j].mean()), 4),
                     "rate_shifted_min_max": [round(float(sh[:, j].min()), 4), round(float(sh[:, j].max()), 4)],
                     "shift": round(shift[W], 4)})
    s, z = joint_test(obs, sh)
    for r, zz in zip(rows, z):
        r["z"] = round(float(zz), 3)
    return {"n_members": n_mem, "n_events": int(len(Mr)), "n_sets": N_SETS_CMIP6,
            "S1_S2": s, "specs": rows}


def main():
    rng = np.random.default_rng(SEED)
    res = {"plan_approved": "2026-09-30", "seed": SEED, "grid": {"W": WS, "tau": TAUS, "f": FS,
           "levels": LEVELS, "cond3": [False, True]}, "window_start": LO,
           "inputs": [str(ERA5.relative_to(ROOT)), "03_data_ingestion/raw/cmip6/*_zm.nc"]}
    res["era5"] = era5_arm(rng)
    e = res["era5"]
    print(f"ERA5 {e['n_events']} events, {e['n_sets']} shifted event-free sets")
    print(f"  S1/S2 conditions 1-2: {e['S1_S2_conditions_1_2']}")
    print(f"  S3 with condition 3:  {e['S3_with_condition_3']}")
    print(f"  published spec: {e['published_spec']}")
    print(f"  POST-HOC calibration: {e['calibration_posthoc']}")
    print(f"  share of contrast reproduced, range over specs: {e['share_reproduced_range_pct']}")
    for r in e["specs"]:
        print(f"   W{r['W']} {r['level']} tau{r['tau']:+.2f} f{r['f']} c3={int(r['cond3'])}: "
              f"SSW {r['rate_ssw']:.3f} shifted {r['rate_shifted_mean']:.3f} "
              f"{r['rate_shifted_q025_q975']} z {r['z']:+.2f} share {r['share_reproduced_pct']}", flush=True)
    res["cmip6"] = cmip6_arm(rng)
    c = res["cmip6"]
    print(f"CMIP6 {c['n_events']} events, {c['n_members']} members: {c['S1_S2']}")
    for r in c["specs"]:
        print(f"   W{r['W']} tau{r['tau']:+.2f} f{r['f']}: SSW {r['rate_ssw']:.3f} "
              f"shifted {r['rate_shifted_mean']:.3f} z {r['z']:+.2f}", flush=True)
    (RESULTS / "criterion_sweep.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print("Saved -> criterion_sweep.json")


if __name__ == "__main__":
    main()
