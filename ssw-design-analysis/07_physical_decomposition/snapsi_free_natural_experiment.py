#!/usr/bin/env python3
"""
snapsi_free_natural_experiment.py
=================================
A CONFIRMATORY TEST, ON DATA NEVER ANALYSED, OF THE LAW THE INTERVENTIONS
SUGGESTED: THE DOWNWARD-NON-DOWNWARD CONTRAST IS WHAT ONE POPULATION, SHIFTED AND
RESCALED, GIVES -- FORCING MOVES IT ONLY THROUGH LOCATION AND SPREAD.

Registered 2026-10-03, committed BEFORE any member of SNAPSI's `free` arm was
downloaded (03_data_ingestion/acquire_snapsi_surface.py --free for psl;
acquire_snapsi_u10.py for u at 10 hPa, 60N), and before the pooled law test (L)
or the regional-class test (G6 of snapsi_contrast_generality) was computed.

WHY
  In every intervention set analysed so far (SNAPSI zonal-mean and full nudging,
  the SH warming, the strongly forced pairs, January 2021), departures from an
  unchanged contrast tracked the change in spread (r 0.95-0.99) -- found AFTER the
  tests. A law found after the fact must be tested on data that played no part in
  finding it. SNAPSI's `free` ensembles (no nudging) start from the same initial
  states as the nudged ones; some members develop the observed SSW and some do
  not. That is a natural experiment in up to ten models -- one of them
  (NOAA-GFDL) used nowhere else in this project.

DATA  SNAPSI `free`, NH initialisations s20180125, s20180208, s20181213, s20190108;
  every centre with free psl and ua except NRL (excluded throughout the project).
  Outcome: polar-cap (60-90N) mean sea-level pressure NAM proxy, exactly as result L
  (window post-onset days +8..+25 after the observed central date, 6-hourly).
  Stratosphere: daily-mean zonal-mean u at 10 hPa, 60N (latitude row nearest 60N).
MEMBERS (classified on the STRATOSPHERE only, never on the outcome)
  SSW member      westerly on lead day 0, and a westerly-to-easterly change of the
                  daily mean within central date -10 .. +10 days
  non-SSW member  westerly on every day from lead 0 to central date +10
  others dropped. An ensemble enters if it has >= 5 of each ("mixed").
  Standardisation: by the non-SSW members' window means (mean 0, s.d. 1).
  DW: conditions 1-2 of the Karpechko criterion, as result L.
TESTS (unit: mixed ensemble; intervals from 10,000 resamples of centres, and,
  if fewer than 6 centres are mixed, centre means are reported instead of an
  interval -- FAILURES 2026-10-03)
  F1 rate: DW share, SSW vs non-SSW members; shift of the mean.
  F2 PRIMARY, the law: residual = C_SSW - C_LS, where C_LS is the contrast of the
     non-SSW members' 6-hourly series mapped x -> m_s + k (x - m_n) (k = ratio of
     s.d. of window means, m = means), classified identically.
     READING, fixed now: interval inside +-0.25 sigma -> "the law holds on
     independent data"; interval excluding 0 with |mean| > 0.25 -> "the law fails:
     the contrast carries information beyond location and spread"; otherwise
     inconclusive.
  F3 the narrower claim: paired C_SSW - C_non; back-shift Delta_b (SSW members
     minus the shift, classified) - C_non; and Delta_b against the law's
     prediction C_non (k - 1): slope and residual.
  F3 READING, added 2026-10-03 before any free-arm result existed (theory note
     theory_spread_law.md: for one location-scale population the contrast's
     elasticity with respect to the spread is exactly 1; for two response types
     separated by Delta whose mixing weight the SSW changes it is about 0.60
     (Delta = 2), 0.37 (Delta = 3), 0.16 (Delta = 4)): the slope of Delta_b on the
     law's prediction C_non (k - 1), with a 95% interval from resampling centres
     (ensembles, if fewer than 6 centres). Interval entirely above 0.60 -> "two
     response types separated by 2 s.d. or more are excluded"; interval entirely
     below 0.60 -> "consistent with forced mixing of separated types"; otherwise
     inconclusive. Whether the interval contains 1 is reported separately.
     (Numbers now reproduced in the repository by threshold_law_theory.py: 0.58,
     0.35, 0.14. Corrected 2026-10-03, still before any free-arm result: the first version read
     an interval excluding 1 as evidence against one population, but mixtures
     predict slopes BELOW 1; the trigger was the post-hoc slope 1.14 [1.04, 1.24] (a
     quick 4,000-resample bootstrap; the stored value is 1.14 [1.05, 1.25])
     on the 36 SNAPSI zonal-mean pairs, data already seen, disclosed here.)
  F4 one population or a mixture: within each ensemble, SSW and non-SSW window
     means each standardised by their own mean and s.d.; pooled residual quantile
     difference at q = 0.05..0.95, and the Anderson-Darling k-sample statistic with
     2,000 within-ensemble label permutations. A mixture of two response types
     whose weight the SSW changes would change the SHAPE: a quantile difference
     beyond +-0.25 sigma, or permutation p < 0.05, refutes one population.
  F5 balance: mean NAM proxy over lead days 0..(central date - 1) (before onset),
     SSW minus non-SSW members, in the same units. |F5| > 0.25 sigma means members
     that develop the SSW already differed, and F1-F3 are then also reported net
     of it (by adjusting each group's mean); F2 and F4 are scale-free in the
     location and are reported as they are.
L POOLED LAW TEST (registered here, computed after F): the F2 residual of every
  intervention set -- SNAPSI zonal-mean (location-scale residual of the 25 pairs,
  snapsi_contrast_generality G3 NH_25), SH (G3), full nudging (recomputed the same
  way on the nudged-full pairs that form a contrast), January 2021 (each pair, from
  davis2021_intervention_test arrays), and F -- averaged with equal weight per SET,
  with a 95% interval from resampling within sets. Same reading as F2.
  Secondary: the same for the translation-only residual (no rescaling), which the
  law predicts to be biased exactly where spread changes.

PRE-RUN CHECK (--selftest, before any data): 12 synthetic ensembles of 25 + 25
  members; one widened population gives F2 +0.008 (s.d. 0.26) and passes F4
  (AD p 0.985); a 50/50 two-type mixture fails F4 (AD p 0.0085, max dq 0.35).

Output: results/current/8_experiment/snapsi_free_natural_experiment.json
"""
import json
import sys
import warnings
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "8_experiment"
ING = HERE.parents[0] / "03_data_ingestion"
FREE_RED = ING / "_snapsi_reduced_free"
U10 = ING / "_snapsi_u10_free"
sys.path.insert(0, str(HERE))
import snapsi_selection_test as L                   # noqa: E402

warnings.filterwarnings("ignore", message=".*(variant|p-value (capped|floored)).*")
NAME = "snapsi_free_natural_experiment"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
N_BOOT, N_PERM, MIN_EACH, TOL, MIN_CENTRES = 10000, 2000, 5, 0.25, 6
NH_INITS = list(L.ONSET)
QS = np.round(np.arange(0.05, 0.951, 0.05), 2)


# ---------------------------------------------------------------- data
def load_psl(centre, init):
    fs = sorted(FREE_RED.glob(f"{centre}_free_{init}_*.parquet"))
    return pd.concat([pd.read_parquet(f) for f in fs]) if fs else None


def load_u10(centre, init):
    fs = sorted(U10.glob(f"u10_{centre}_free_{init}_*.parquet"))
    return pd.concat([pd.read_parquet(f) for f in fs]) if fs else None


def classify_members(u, init):
    """member -> 'ssw' | 'non' | None from daily-mean u(10 hPa, 60N)."""
    off = L.post_onset_offset(init)                  # central date, lead days
    d = u.assign(day=np.floor(u.lead_days).astype(int)).groupby(["member", "day"])["u10_60N"].mean().unstack()
    out = {}
    for m, row in d.iterrows():
        row = row.dropna()
        if 0 not in row.index or row.loc[0] <= 0 or row.index.max() < off + 10:
            out[m] = None
            continue
        win = row.loc[max(1, off - 10): off + 10]
        prev = row.shift(1).loc[win.index]
        rev = bool(((win < 0) & (prev > 0)).any())
        west = bool((row.loc[0: off + 10] > 0).all())
        out[m] = "ssw" if rev else ("non" if west else None)
    return out


def member_series(df, init, base, sd):
    off = L.post_onset_offset(init)
    lo, hi = L.WINDOW[0] + off, L.WINDOW[1] + off
    w = df[(df.lead_days >= lo) & (df.lead_days <= hi)]
    return {m: -(g.psl_cap.values - base) / sd for m, g in w.groupby("member")}


def pre_onset(df, init, base, sd):
    off = L.post_onset_offset(init)
    w = df[(df.lead_days >= 0) & (df.lead_days < off)]
    return {m: float(np.mean(-(g.psl_cap.values - base) / sd)) for m, g in w.groupby("member")}


# ---------------------------------------------------------------- statistics
def classify(series):
    means = np.array([s.mean() for s in series])
    return means, (means < 0) & np.array([(s < 0).mean() > 0.5 for s in series])


def contrast(series):
    means, dw = classify(series)
    if dw.sum() < 3 or (~dw).sum() < 3:
        return np.nan
    return float(means[dw].mean() - means[~dw].mean())


def ensemble_stats(S, N):
    ms, dws = classify(S); mn, dwn = classify(N)
    m_s, m_n = float(ms.mean()), float(mn.mean())
    k = float(ms.std(ddof=1) / mn.std(ddof=1))
    s = m_s - m_n
    out = {"n_ssw": len(S), "n_non": len(N), "dw_rate_ssw": float(dws.mean()), "dw_rate_non": float(dwn.mean()),
           "shift": s, "sd_ratio": k, "C_ssw": contrast(S), "C_non": contrast(N),
           "C_LS": contrast([m_s + k * (x - m_n) for x in N]), "C_backshift": contrast([x - s for x in S]),
           "C_shift_null": contrast([x + s for x in N])}
    out["F2_resid_law"] = out["C_ssw"] - out["C_LS"]
    out["paired_diff"] = out["C_ssw"] - out["C_non"]
    out["Delta_b"] = out["C_backshift"] - out["C_non"]
    out["law_pred_Delta_b"] = out["C_non"] * (k - 1)
    zs = (ms - m_s) / ms.std(ddof=1); zn = (mn - m_n) / mn.std(ddof=1)
    out["_z"] = (zs, zn)
    return out


def summarise(df, col, rng):
    d = df.dropna(subset=[col])
    if d.empty:
        return None
    cs = d.centre.unique()
    res = {"mean": round(float(d[col].mean()), 4), "n_ensembles": int(len(d)), "n_centres": int(len(cs)),
           "centre_means": {c: round(float(g[col].mean()), 4) for c, g in d.groupby("centre")}}
    if len(cs) >= MIN_CENTRES:
        by = {c: d.loc[d.centre == c, col].values for c in cs}
        bs = [np.mean(np.concatenate([by[c] for c in rng.choice(cs, len(cs))])) for _ in range(N_BOOT)]
        res["ci95"] = [round(float(q), 4) for q in np.quantile(bs, [0.025, 0.975])]
    else:
        res["ci95"] = None
        res["note"] = f"fewer than {MIN_CENTRES} centres: centre means, no interval"
    return res


def reading(r):
    if r is None or r.get("ci95") is None:
        return "no interval (too few centres)"
    lo, hi = r["ci95"]
    if -TOL <= lo and hi <= TOL:
        return "the law holds on independent data"
    if (lo > 0 or hi < 0) and abs(r["mean"]) > TOL:
        return "the law fails"
    return "inconclusive"


def shape_test(zs_list, rng):
    """F4: pooled residual quantiles and permutation AD over ensembles."""
    from scipy.stats import anderson_ksamp
    zs = np.concatenate([a for a, _ in zs_list]); zn = np.concatenate([b for _, b in zs_list])
    dq = np.quantile(zs, QS) - np.quantile(zn, QS)
    ad0 = anderson_ksamp([zs, zn]).statistic
    exceed = 0
    for _ in range(N_PERM):
        a, b = [], []
        for x, y in zs_list:
            z = np.concatenate([x, y]); p = rng.permutation(len(z))
            a.append(z[p[:len(x)]]); b.append(z[p[len(x):]])
        exceed += anderson_ksamp([np.concatenate(a), np.concatenate(b)]).statistic >= ad0
    return {"quantile_diff": dict(zip([str(q) for q in QS], np.round(dq, 4).tolist())),
            "max_abs_quantile_diff": round(float(np.max(np.abs(dq))), 4),
            "AD_stat": round(float(ad0), 4), "AD_perm_p": round((exceed + 1) / (N_PERM + 1), 4),
            "reading": ("one population (shape unchanged)"
                        if np.max(np.abs(dq)) <= TOL and (exceed + 1) / (N_PERM + 1) >= 0.05
                        else "shape differs: refutes one population")}


def davis_ls_residuals():
    """January 2021 pairs: nudged contrast minus the location-scale null of the
    no-SSW arm (daily NAM proxy), and the translation-only residual."""
    import davis2021_intervention_test as D
    d = pd.read_parquet(D.DATA)
    out = []
    for name, (a_s, a_c, ref) in D.PAIRS.items():
        Zs, _ = D.arm_arrays(d, a_s, ref); Zc, _ = D.arm_arrays(d, a_c, ref)
        cm = Zc.mean(axis=1); base, sd = float(cm.mean()), float(cm.std(ddof=1))
        Ns = -(Zs.values - base) / sd; Nc = -(Zc.values - base) / sd
        ms, mc = Ns.mean(1), Nc.mean(1)
        k = ms.std(ddof=1) / mc.std(ddof=1)
        cs = contrast(list(Ns))
        out.append({"unit": name, "resid_LS": cs - contrast(list(ms.mean() + k * (Nc - mc.mean()))),
                    "resid_shift": cs - contrast(list(Nc + ms.mean() - mc.mean()))})
    return pd.DataFrame(out)


def pooled_law(f_rows, rng):
    """L: equal weight per intervention set, bootstrap within sets."""
    gen = json.loads((RESULTS / "snapsi_contrast_generality.json").read_text())
    nh = pd.DataFrame(gen["pairs_NH"]); nh = nh[~nh.strong]
    sh = pd.DataFrame(gen["G3_third_event_SH"]["sh_per_centre"])
    full = pd.DataFrame(gen["G5_nudged_full"]["control"]["pairs"]) if isinstance(gen.get("G5_nudged_full"), dict) else None
    sets = {"SNAPSI_zonal_NH_25": (nh.resid_locscale.values, nh.resid_shift.values),
            "SNAPSI_SH": (sh.resid_locscale.values, sh.resid_shift.values)}
    if full is not None:
        sets["SNAPSI_full_3D"] = (full.resid_locscale.values, full.resid_shift.values)
    dv = davis_ls_residuals()
    sets["CESM2_Jan2021"] = (dv.resid_LS.values, dv.resid_shift.values)
    if f_rows is not None and len(f_rows):
        sets["SNAPSI_free_natural"] = (f_rows.F2_resid_law.values,
                                       (f_rows.C_ssw - contrast_shift_only(f_rows)).values)
    out = {"sets": {}}
    for which, j in (("location_scale", 0), ("shift_only", 1)):
        clean = {k: v[j][np.isfinite(v[j])] for k, v in sets.items()}
        clean = {k: v for k, v in clean.items() if len(v)}
        pt = float(np.mean([v.mean() for v in clean.values()]))
        bs = [np.mean([rng.choice(v, len(v)).mean() for v in clean.values()]) for _ in range(N_BOOT)]
        r = {"mean": round(pt, 4), "ci95": [round(float(q), 4) for q in np.quantile(bs, [0.025, 0.975])],
             "n_sets": len(clean), "per_set": {k: [round(float(v.mean()), 4), int(len(v))] for k, v in clean.items()}}
        r["reading"] = reading(r)
        out[which] = r
    # POST HOC (added 2026-10-04 after the run, for Fig. 2b): a 95% interval for each
    # set's mean residual, from its own random stream; the pooled numbers are unchanged
    rp = np.random.default_rng(SEED + 70)
    for which, j in (("location_scale", 0), ("shift_only", 1)):
        out[which]["per_set_ci95"] = {}
        for k, v in sets.items():
            x = v[j][np.isfinite(v[j])]
            if len(x) >= 2:
                bs = [rp.choice(x, len(x)).mean() for _ in range(N_BOOT)]
                out[which]["per_set_ci95"][k] = [round(float(q), 4) for q in np.quantile(bs, [0.025, 0.975])]
    return out


def contrast_shift_only(f_rows):
    return f_rows.get("C_shift_null", pd.Series(np.nan, index=f_rows.index))


def selftest():
    """Pre-run check (registered): synthetic one-population ensembles (widened)
    must give F2 ~ 0 and pass F4; a 50/50 two-type mixture must fail F4."""
    rng = np.random.default_rng(1)
    def ens(n, mu, sd, mix=False):
        out = []
        for _ in range(n):
            m = (mu + sd * rng.normal()) if not mix else \
                (mu - 1.5 + 0.5 * rng.normal() if rng.random() < 0.5 else mu + 1.0 + 0.5 * rng.normal())
            out.append(m + rng.normal(0, 1.2, 72))
        return out
    res, zl, zm = [], [], []
    for _ in range(12):
        st = ensemble_stats(ens(25, -0.8, 1.3), ens(25, 0, 1.0)); zl.append(st.pop("_z")); res.append(st["F2_resid_law"])
        st2 = ensemble_stats(ens(25, -0.8, 1.0, mix=True), ens(25, 0, 1.0)); zm.append(st2.pop("_z"))
    a, b = shape_test(zl, rng), shape_test(zm, rng)
    print(f"selftest: F2 mean {np.nanmean(res):+.3f} (s.d. {np.nanstd(res):.2f}); one population AD p "
          f"{a['AD_perm_p']}, max dq {a['max_abs_quantile_diff']}; mixture AD p {b['AD_perm_p']}, "
          f"max dq {b['max_abs_quantile_diff']}")
    return 0 if abs(np.nanmean(res)) < 0.1 and a["AD_perm_p"] >= 0.05 and b["AD_perm_p"] < 0.05 else 1


def main():
    if "--selftest" in sys.argv:
        return selftest()
    rng = np.random.default_rng(SEED)
    centres = sorted({p.name.split("_")[0] for p in FREE_RED.glob("*.parquet")} - {"NRL"})
    if not centres:
        sys.exit("no free-arm data: run acquire_snapsi_surface.py --free and acquire_snapsi_u10.py")
    rows, zl, bal = [], [], []
    for c in centres:
        for init in NH_INITS:
            p, u = load_psl(c, init), load_u10(c, init)
            if p is None or u is None or not L.spans_window(p, init):
                continue
            lab = classify_members(u, init)
            non = [m for m, v in lab.items() if v == "non"]; ssw = [m for m, v in lab.items() if v == "ssw"]
            row = {"centre": c, "init": init, "n_members": len(lab), "n_ssw_all": len(ssw), "n_non_all": len(non)}
            pm = p[p.member.isin(non)]
            if len(ssw) >= MIN_EACH and len(non) >= MIN_EACH:
                wm = L.window_means(pm, init)
                base, sd = float(wm.mean()), float(wm.std(ddof=1))
                Sd = member_series(p[p.member.isin(ssw)], init, base, sd)
                Nd = member_series(pm, init, base, sd)
                st = ensemble_stats(list(Sd.values()), list(Nd.values()))
                zl.append(st.pop("_z"))
                pre = pre_onset(p, init, base, sd)
                st["F5_pre_onset_diff"] = float(np.mean([pre[m] for m in ssw if m in pre]) -
                                                np.mean([pre[m] for m in non if m in pre]))
                row.update(st); row["mixed"] = True
            else:
                row["mixed"] = False
            rows.append(row)
    df = pd.DataFrame(rows)
    mx = df[df.mixed] if "mixed" in df else df.iloc[0:0]
    res = {"registered": "2026-10-03, committed before any free-arm member was downloaded",
           "seed": SEED, "n_boot": N_BOOT, "n_perm": N_PERM, "tolerance_sigma": TOL,
           "n_ensembles": int(len(df)), "n_mixed": int(len(mx)), "centres_mixed": sorted(mx.centre.unique().tolist())}
    for col in ("F2_resid_law", "paired_diff", "Delta_b", "shift", "dw_rate_ssw", "dw_rate_non",
                "sd_ratio", "F5_pre_onset_diff"):
        res[col] = summarise(mx, col, rng)
    res["F2_reading"] = reading(res["F2_resid_law"])
    d3 = mx.dropna(subset=["Delta_b", "law_pred_Delta_b"])
    if len(d3) > 2:
        b = np.polyfit(d3.law_pred_Delta_b, d3.Delta_b, 1)
        r3 = np.random.default_rng(SEED + 3)
        units = d3.centre.unique() if d3.centre.nunique() >= MIN_CENTRES else d3.index.values
        key = "centre" if d3.centre.nunique() >= MIN_CENTRES else None
        bs = []
        for _ in range(N_BOOT):
            pick = r3.choice(units, len(units))
            dd = pd.concat([d3[d3.centre == u] for u in pick]) if key else d3.loc[pick]
            if dd.law_pred_Delta_b.nunique() > 2:
                bs.append(np.polyfit(dd.law_pred_Delta_b, dd.Delta_b, 1)[0])
        ci = [round(float(q), 4) for q in np.quantile(bs, [0.025, 0.975])]
        res["F3_Delta_b_vs_law"] = {"slope": round(float(b[0]), 4), "slope_ci95": ci, "resampled": key or "ensembles",
                                    "intercept": round(float(b[1]), 4),
                                    "r": round(float(np.corrcoef(d3.law_pred_Delta_b, d3.Delta_b)[0, 1]), 4),
                                    "mean_residual": round(float((d3.Delta_b - d3.law_pred_Delta_b).mean()), 4),
                                    "contains_1": bool(ci[0] <= 1 <= ci[1]),
                                    "reading": ("two response types separated by 2 s.d. or more are excluded"
                                                if ci[0] > 0.60 else
                                                "consistent with forced mixing of separated types"
                                                if ci[1] < 0.60 else "inconclusive")}
    if zl:
        res["F4_shape"] = shape_test(zl, rng)
    # L is computed with its own random stream, after F
    res["L_pooled_law"] = pooled_law(mx if len(mx) else None, np.random.default_rng(SEED + 7))
    res["ensembles"] = df.round(4).to_dict(orient="records")
    for k in ("n_ensembles", "n_mixed", "centres_mixed", "F2_resid_law", "F2_reading", "Delta_b",
              "F5_pre_onset_diff", "F3_Delta_b_vs_law", "F4_shape", "L_pooled_law"):
        print(k, json.dumps(res.get(k))[:400], flush=True)
    (RESULTS / f"{NAME}.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf8", newline="\n")
    print(f"Saved -> {NAME}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
