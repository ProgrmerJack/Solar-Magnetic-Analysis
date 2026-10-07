#!/usr/bin/env python3
"""
snapsi_contrast_generality.py
=============================
DOES AN IMPOSED SSW LEAVE THE CLASS STRUCTURE UNCHANGED WHERE THE FORCING IS
STRONGEST, IN A THIRD EVENT, AND FOR A VARIABLE NOT USED TO CLASSIFY?

Registered 2026-10-02 (revision 5), committed before any of these quantities was
computed. Answers three objections to result L (snapsi_selection_test):
  (i)   its paired contrast uses 25 of 36 NH pairs; the other 11 nudged ensembles
        are 96-100% DW and cannot form a contrast, so invariance might hold only
        where the forcing is weaker;
  (ii)  the causal core is two NH events;
  (iii) the observational test of a variable not used to classify
        (criterion_regional_contrast) has power 0.30.
Survey: snapsi_contrast_diagnosis D2 shifts the control forward on the 23 pairs
where both classes exist; nothing back-shifts the nudged members, uses the 11
excluded pairs, transforms the SH control by location and scale, or applies the
matched shifted null to SNAPSI temperature. Hence this script.

UNITS AND DATA exactly as result L: window-mean polar-cap NAM proxy over post-onset
days +8..+25 (6-hourly), standardised by the pair's control ensemble; DW =
Karpechko conditions 1-2 (window mean < 0 and > 50% of 6-hourly values < 0);
a contrast needs >= 3 members per class. Imposed effect s = mean nudged minus mean
control window mean. Intervals: 10,000 resamples of centres (pairs nested in
centres), as L. Tolerance 0.25 sigma, as fixed for result Z1 (residual quantiles).

G1 PRIMARY -- the 11 strongly forced pairs (back-shift test)
  Every 6-hourly value of every nudged member minus s, classified identically;
  Delta_b = contrast(back-shifted nudged) - contrast(control). If the imposed SSW
  acts on the class structure only by translating the distribution, Delta_b ~ 0
  in every pair, including those where the unshifted nudged ensemble has no NDW
  members. Reported for the 11 pairs (primary), the 25 and all 36.
  Reading, fixed now: 95% interval of the 11-pair mean inside [-0.25, +0.25]
  -> "class structure unchanged where the forcing is strongest"; interval
  excludes 0 AND |mean| > 0.25 -> "the forcing changes the class structure there;
  the contrast claim must be restricted to the 25 pairs"; otherwise inconclusive.
  Sensitivity (would G1 see it?): nudged members' deviations from the ensemble-
  mean series inflated by k = sqrt(1.25) and sqrt(1.5) (forced variance of 0.25
  and 0.5 sigma^2) before back-shifting; detected if the interval excludes 0.
  Also: pooled variance ratio nudged/control in the 11 pairs, members re-centred
  on their ensemble means, 10,000 resamples of members within ensembles.
G2 dose-response: slope of Delta_b on s over 36 pairs, and of the paired
  difference C_nudged - C_control on s over the 25; centre bootstrap. A slope
  interval that covers 0 = no evidence that the contrast moves with forcing
  strength; the fitted Delta at the mean s of the 11 pairs is reported.
G3 THIRD EVENT (SH minor warming, s20190829, 8 centres)
  (a) shift-only null: control + s, classified; residual C_nudged - C_null.
  (b) location-scale null: control 6-hourly values x -> m_n + k (x - m_c), with
      m and k = sd_n / sd_c from the window means, classified; residual likewise.
  (c) Delta_b as in G1.
  Same three for the 25 NH pairs for comparison. Reading, fixed now: if (b)'s
  interval lies inside +-0.25 while (a)'s excludes 0, the SH contrast change is
  the widening of one population -- forcing changes the contrast there only
  through the spread, which is itself a forced change, and the paper must say so.
  Caveat stated now: (b) matches the first two moments by construction, so it
  tests only shape beyond them and condition 2.
G4 REGIONAL TEMPERATURE, high-power analogue of criterion_regional_contrast
  T = northern-Eurasian (50-65N, 10-130E) window-mean tas, post-onset days 8-24,
  standardised by the pair's control (also reported in K).
  (a) PRIMARY: observed C_T = mean T(DW) - mean T(NDW) among nudged members;
      matched shifted null C_T0 = the same among control members whose NAM
      series is displaced by s and classified identically, each keeping its own
      temperature; residual C_T - C_T0 over pairs where both have >= 3 per
      class. Reading, fixed now: interval inside +-0.25 -> "the regional class
      contrast is what selection plus the SSW-free circulation-temperature
      relation produce"; interval excludes 0 AND |mean| > 0.25 -> "the label
      marks regional information beyond the shifted null"; otherwise
      inconclusive. Power: the same statistic with a planted extra cooling of
      0.5 sigma (and of the observational test's 1.32 K, converted per pair) on
      nudged DW members.
  (b) intervention: C_T(nudged) - C_T(control), paired, the 25 pairs; and with
      back-shifted nudged members, all 36 pairs.
  Secondary: the other three regions of snapsi_regional_test.
G5 nudged-full (the full 3-D stratosphere, vortex geometry included), registered
  now, run only once its members are cached: result L's paired contrast and G1's
  Delta_b for nudged-full against control (control-full for UKMO, both reported)
  -- ECMWF and UKMO at the four NH initialisations, Meteo-France s20190108.
  NOAA-GFDL has no control arm and is not used.

G6 OTHER IMPACT CLASSES (registered 2026-10-03, committed before it was computed;
  own random stream): is the law a property of outcome-defined classes in general,
  not of the SSW label? 24 Northern Hemisphere regions (latitude bands 30-50N and
  50-70N x twelve 30-degree longitude sectors), window-mean 2 m temperature over
  post-onset days 8-24, standardised by each pair's control. Class: "cold outcome"
  = below the control ensemble's lower tercile (a fixed threshold, like the label's
  zero). In each of the 36 NH pairs and each region, where both arms have >= 3
  members per class: contrast C = mean(cold) - mean(not cold); residual of the
  SSW-imposed contrast against (a) the location-scale null of the control (the
  law) and (b) the shift-only null; and the change of the cold-class rate.
  Summary over all region-pairs (bootstrap over centres, regions kept together).
  READING, fixed now: law residual inside +-0.25 -> "the law holds for regional
  impact classes"; excluding 0 with |mean| > 0.25 -> "it fails"; else
  inconclusive. Secondary: the share of region-pairs whose |law residual| is
  below 0.25, against that of the shift-only null.

G7 WHY IS THE ZONAL-PAIR ELASTICITY 1.14, NOT 1? (registered 2026-10-04 before it
  was computed; own random stream). The law's elasticity of exactly 1 holds for a
  pure threshold on one location-scale variable. The downward label adds condition
  2 (more than half the 6-hourly values negative), and forcing may scale the
  spread of member means without scaling each member's within-window noise.
  (a) The same slope (back-shift difference on C_control (k - 1), 36 pairs, centre
      bootstrap) with condition 1 only (window mean < 0).
  (b) Calibration under the exact law: in each pair, synthetic forced members built
      from resampled control members, Y_i(t) = s + mbar + k (m_i - mbar) + (X_i(t) - m_i)
      (member means scaled by the pair's observed k, within-window noise kept), and
      a variant scaling the whole series, Y_i(t) = s + mbar + k (X_i(t) - mbar);
      the slope recomputed over the 36 pairs, 200 replicates of each.
  READING, fixed now: if the observed 1.14 lies inside the 95% range of the
  means-only calibration, the excess over 1 is what the two-condition label gives
  under one population; if (a) contains 1, the excess comes from condition 2. If
  1.14 lies outside both calibrations, the excess is unexplained and is reported.

G8 SHAPE TEST WITH POWER, FOR SSWs (registered 2026-10-04 before it was computed;
  own random stream). The natural experiment's shape test is underpowered against
  asymmetric two-type mixtures (threshold_law_theory shape_power: 0-0.47 for types
  1.5-2 s.d. apart). The 36 SNAPSI zonal-nudging pairs (~50 + ~50 members) are a
  design as large as the human-forcing one. Within each pair both arms are
  standardised by their own mean and s.d. of window means; pooled residual quantiles
  (q 0.05..0.95) and the Anderson-Darling k-sample statistic, 2,000 within-pair
  permutations. Power for this design (36 units of 50 + 50; reference weight 0.3,
  forced 0.6; Delta in {1.5, 2, 3}, width ratio r in {0.5, 1, 2}; 100 replicates of
  200 permutations). READING, fixed now: if p >= 0.05 and max |dq| <= 0.25, every
  (Delta, r) with power >= 0.8 is excluded for the imposed SSWs.

G9 DOES THE TIMING OF THE FORCED RESPONSE EXPLAIN THE 1.14? (registered 2026-10-07
  before it was computed; own random stream; run with --g9). G7's exact-law calibration
  moves every member by ONE shift over the whole window (1.05 [0.96, 1.13]); the
  observed two-condition slope 1.1434 lies just outside. The label's condition 2 counts
  negative 6-hourly values, so a forced shift that varies within the window acts on it
  differently from a constant one. Per pair, the ensemble-mean series of each arm and the
  across-member s.d. series, smoothed by a centred 5-day (20-step) running mean:
  (c2) timing only: Y_i(t) = s(t) + mbar + k (m_i - mbar) + (X_i(t) - m_i), with
       s(t) = smoothed nudged minus control mean series, k the pair's s.d. ratio of
       member means (G7's means-only construction with s(t) for s);
  (c1) time-varying marginals: Y_i(t) = mn(t) + (sn(t) / sc(t)) (X_i(t) - mc(t))
       (forced mean and spread matched at every time; control temporal dependence kept);
  each from resampled control members, slope over the 36 pairs, 200 replicates.
  HELD OUT (secondary): s(t) and k fitted on a random half of the nudged members, the
  slope computed from the other half against control, c2 calibrated in the same sizes;
  10 splits; reported: share of splits with the held-out slope inside its range.
  DIAGNOSTIC: longest run of consecutive negative values per member, nudged vs c1
  synthetic (mean difference, pooled over pairs) -- order-sensitive, fitted by neither.
  READING, fixed now: 1.1434 inside the c2 95% range -> "the excess is the timing of the
  forced response"; else inside c1 -> "explained by time-varying marginals (timing and
  spread)"; else "unexplained by the marginals: a change of temporal dependence".

Output: results/current/8_experiment/snapsi_contrast_generality.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "8_experiment"
sys.path.insert(0, str(HERE))
import snapsi_selection_test as L                   # noqa: E402
import snapsi_regional_test as RT                   # noqa: E402

NAME = "snapsi_contrast_generality"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
N_BOOT = 10000
MIN_CLASS = 3
TOL = 0.25
NH_INITS = list(L.ONSET)
SH_INIT = "s20190829"
FULL = {"ECMWF": NH_INITS, "UKMO": NH_INITS, "Meteo-France": ["s20190108"]}
REG = "NEURASIA"
# nudged-full / control-full live in their own cache (acquire_snapsi_surface.py --full);
# location only, set before G5 was run -- the registered design is unchanged
FULL_RED = HERE.parents[0] / "03_data_ingestion" / "_snapsi_reduced_full"


def load(centre, arm, init):
    if arm.endswith("-full"):
        fs = sorted(FULL_RED.glob(f"{centre}_{arm}_{init}_*.parquet"))
        return pd.concat([pd.read_parquet(f) for f in fs]) if fs else None
    return L.load(centre, arm, init)


# ---------------------------------------------------------------- member data
def member_nam(centre, init, arms=("nudged", "control"), ref="control"):
    """{arm: {member: standardised 6-hourly NAM proxy over the window}}."""
    r = load(centre, ref, init)
    if r is None or not L.spans_window(r, init):
        return None
    rm = L.window_means(r, init)
    base, sd = float(rm.mean()), float(rm.std(ddof=1))
    off = L.post_onset_offset(init)
    lo, hi = L.WINDOW[0] + off, L.WINDOW[1] + off
    out = {}
    for arm in arms:
        d = load(centre, arm, init)
        if d is None or not L.spans_window(d, init):
            return None
        w = d[(d["lead_days"] >= lo) & (d["lead_days"] <= hi)]
        out[arm] = {m: -(g["psl_cap"].values - base) / sd for m, g in w.groupby("member")}
    return out


def classify(series):
    means = np.array([s.mean() for s in series])
    dw = (means < 0) & np.array([(s < 0).mean() > 0.5 for s in series])
    return means, dw


def contrast(series, values=None):
    """DW - NDW contrast of `values` (default: the window means themselves)."""
    means, dw = classify(series)
    v = means if values is None else np.asarray(values)
    if dw.sum() < MIN_CLASS or (~dw).sum() < MIN_CLASS:
        return np.nan
    return float(v[dw].mean() - v[~dw].mean())


def inflate(series, k):
    """Members' deviations from the ensemble-mean series scaled by k."""
    n = min(len(s) for s in series)
    X = np.array([s[:n] for s in series])
    return list(X.mean(0) + k * (X - X.mean(0)))


# ---------------------------------------------------------------- statistics
def cboot(df, col, rng, stat=np.mean):
    df = df.dropna(subset=[col])
    if df.empty:
        return None
    cs = df["centre"].unique()
    by = {c: df.loc[df.centre == c, col].values for c in cs}
    bs = np.array([stat(np.concatenate([by[c] for c in rng.choice(cs, len(cs))])) for _ in range(N_BOOT)])
    return {"mean": round(float(stat(df[col].values)), 4),
            "ci95": [round(float(q), 4) for q in np.quantile(bs, [0.025, 0.975])],
            "n_pairs": int(len(df)), "n_centres": int(len(cs))}


def slope_boot(df, x, y, rng):
    df = df.dropna(subset=[x, y])
    cs = df["centre"].unique()
    by = {c: df[df.centre == c] for c in cs}
    b0 = float(np.polyfit(df[x], df[y], 1)[0])
    bs = []
    for _ in range(N_BOOT):
        d = pd.concat([by[c] for c in rng.choice(cs, len(cs))])
        if d[x].nunique() > 2:
            bs.append(np.polyfit(d[x], d[y], 1)[0])
    a = np.polyfit(df[x], df[y], 1)
    return {"slope": round(b0, 4), "ci95": [round(float(q), 4) for q in np.quantile(bs, [0.025, 0.975])],
            "intercept": round(float(a[1]), 4), "n_pairs": int(len(df))}


def reading(ci, mean, tol=TOL):
    if tol * -1 <= ci[0] and ci[1] <= tol:
        return "within tolerance"
    if (ci[0] > 0 or ci[1] < 0) and abs(mean) > tol:
        return "departure"
    return "inconclusive"


def pooled_var_ratio(ens, rng):
    vr = float(np.var(np.concatenate([e[0] for e in ens]), ddof=1) /
               np.var(np.concatenate([e[1] for e in ens]), ddof=1))
    bs = []
    for _ in range(N_BOOT):
        a = np.concatenate([rng.choice(e[0], len(e[0])) for e in ens])
        b = np.concatenate([rng.choice(e[1], len(e[1])) for e in ens])
        bs.append(np.var(a, ddof=1) / np.var(b, ddof=1))
    return {"ratio": round(vr, 4), "ci95": [round(float(q), 4) for q in np.quantile(bs, [0.025, 0.975])],
            "n_ensembles": len(ens)}


# ---------------------------------------------------------------- per pair
def pair_row(centre, init, nm, arm_n="nudged", arm_c="control"):
    nud, con = list(nm[arm_n].values()), list(nm[arm_c].values())
    mn, dwn = classify(nud)
    mc, dwc = classify(con)
    s = float(mn.mean() - mc.mean())
    k = float(mn.std(ddof=1) / mc.std(ddof=1))
    m_n, m_c = float(mn.mean()), float(mc.mean())
    row = {"centre": centre, "init": init, "s": s, "sd_ratio": k,
           "dw_rate_nudged": float(dwn.mean()), "dw_rate_control": float(dwc.mean()),
           "C_nudged": contrast(nud), "C_control": contrast(con),
           "C_backshift": contrast([x - s for x in nud]),
           "C_null_shift": contrast([x + s for x in con]),
           "C_null_locscale": contrast([m_n + k * (x - m_c) for x in con])}
    row["paired_diff"] = row["C_nudged"] - row["C_control"]
    row["Delta_b"] = row["C_backshift"] - row["C_control"]
    row["resid_shift"] = row["C_nudged"] - row["C_null_shift"]
    row["resid_locscale"] = row["C_nudged"] - row["C_null_locscale"]
    for lab, kk in (("k1.25", np.sqrt(1.25)), ("k1.5", np.sqrt(1.5))):
        infl = inflate(nud, kk)
        si = float(np.mean([x.mean() for x in infl]) - m_c)
        row[f"Delta_b_{lab}"] = contrast([x - si for x in infl]) - row["C_control"]
    return row, (mn - m_n, mc - m_c)


def regional_row(centre, init, nm):
    """G4 quantities for one NH pair (T standardised by the control)."""
    tn, tc = RT.region_series(centre, "nudged", init), RT.region_series(centre, "control", init)
    if tn.empty or tc.empty:
        return None
    nud = {m: v for m, v in nm["nudged"].items() if m in tn.index}
    con = {m: v for m, v in nm["control"].items() if m in tc.index}
    if len(nud) < 10 or len(con) < 10:
        return None
    mu, sdK = float(tc.loc[list(con), REG].mean()), float(tc.loc[list(con), REG].std(ddof=1))
    Tn = (tn.loc[list(nud), REG].values - mu) / sdK
    Tc = (tc.loc[list(con), REG].values - mu) / sdK
    sn, sc = list(nud.values()), list(con.values())
    mn, dwn = classify(sn)
    mc, _ = classify(sc)
    s = float(mn.mean() - mc.mean())
    row = {"centre": centre, "init": init, "sd_K": sdK, "s": s,
           "CT_nudged": contrast(sn, Tn), "CT_control": contrast(sc, Tc),
           "CT_null": contrast([x + s for x in sc], Tc),
           "CT_backshift": contrast([x - s for x in sn], Tn)}
    row["resid"] = row["CT_nudged"] - row["CT_null"]
    row["resid_plant_0.5sigma"] = contrast(sn, Tn - 0.5 * dwn) - row["CT_null"]
    row["resid_plant_1.32K"] = contrast(sn, Tn - (1.3192 / sdK) * dwn) - row["CT_null"]
    row["paired_diff"] = row["CT_nudged"] - row["CT_control"]
    row["Delta_b"] = row["CT_backshift"] - row["CT_control"]
    for r in RT.REGIONS:
        if r == REG:
            continue
        m2, s2 = float(tc.loc[list(con), r].mean()), float(tc.loc[list(con), r].std(ddof=1))
        Tn2 = (tn.loc[list(nud), r].values - m2) / s2
        Tc2 = (tc.loc[list(con), r].values - m2) / s2
        row[f"resid_{r}"] = contrast(sn, Tn2) - contrast([x + s for x in sc], Tc2)
    return row


# ---------------------------------------------------------------- G7
def contrast_c1(series):
    means = np.array([x.mean() for x in series]); dw = means < 0
    if dw.sum() < MIN_CLASS or (~dw).sum() < MIN_CLASS:
        return np.nan
    return float(means[dw].mean() - means[~dw].mean())


def g7_units(nud, con, cfun):
    mn = np.array([x.mean() for x in nud]); mc = np.array([x.mean() for x in con])
    s = float(mn.mean() - mc.mean()); k = float(mn.std(ddof=1) / mc.std(ddof=1))
    cc = cfun(con)
    return cfun([x - s for x in nud]) - cc, cc * (k - 1)


def slope(xs, ys):
    xs, ys = np.asarray(xs), np.asarray(ys); ok = np.isfinite(xs) & np.isfinite(ys)
    return float(np.polyfit(xs[ok], ys[ok], 1)[0]) if ok.sum() > 3 else np.nan


def elasticity_diagnosis(usable, strong, rng, n_rep=200):
    pairs = []
    for c in usable:
        for init in NH_INITS:
            nm = member_nam(c, init)
            if nm is not None:
                pairs.append((c, list(nm["nudged"].values()), list(nm["control"].values())))
    out = {}
    # (a) condition 1 only, with the centre bootstrap
    rows = [(c, *g7_units(n, k, contrast_c1)) for c, n, k in pairs]
    d = pd.DataFrame(rows, columns=["centre", "y", "x"]).dropna()
    cs = d.centre.unique(); bs = []
    for _ in range(N_BOOT):
        dd = pd.concat([d[d.centre == cc] for cc in rng.choice(cs, len(cs))])
        bs.append(slope(dd.x, dd.y))
    out["a_condition1_only"] = {"slope": round(slope(d.x, d.y), 4),
                                "ci95": [round(float(q), 4) for q in np.nanquantile(bs, [0.025, 0.975])],
                                "n_pairs": int(len(d))}
    # (b) calibration under the exact law
    obs = [(n, k) for _, n, k in pairs]
    ks = [np.std([x.mean() for x in n], ddof=1) / np.std([x.mean() for x in k], ddof=1) for n, k in obs]
    ss = [np.mean([x.mean() for x in n]) - np.mean([x.mean() for x in k]) for n, k in obs]
    for variant in ("means_only", "whole_series"):
        sl = []
        for _ in range(n_rep):
            xs, ys = [], []
            for (nud, con), k, s in zip(obs, ks, ss):
                n0 = min(len(x) for x in con)
                C = [con[i][:n0] for i in rng.integers(len(con), size=len(con))]
                F = [con[i][:n0] for i in rng.integers(len(con), size=len(nud))]
                mbar = np.mean([x.mean() for x in con])
                if variant == "means_only":
                    Y = [s + mbar + k * (x.mean() - mbar) + (x - x.mean()) for x in F]
                else:
                    Y = [s + mbar + k * (x - mbar) for x in F]
                y_, x_ = g7_units(Y, C, contrast)
                ys.append(y_); xs.append(x_)
            sl.append(slope(xs, ys))
        sl = np.asarray(sl, float)
        out[f"b_calibration_{variant}"] = {"mean_slope": round(float(np.nanmean(sl)), 4),
                                           "range95": [round(float(q), 4) for q in np.nanquantile(sl, [0.025, 0.975])],
                                           "n_rep": n_rep}
    obs_slope = 1.1434
    r = out["b_calibration_means_only"]["range95"]
    out["reading"] = ("the excess over 1 is what the two-condition label gives under one population"
                      if r[0] <= obs_slope <= r[1] else
                      "the excess comes from condition 2" if out["a_condition1_only"]["ci95"][0] <= 1 <= out["a_condition1_only"]["ci95"][1]
                      else "the excess is unexplained")
    out["observed_two_condition_slope"] = obs_slope
    print("G7:", out, flush=True)
    return out


# ---------------------------------------------------------------- G8
def _ad_perm(units, rng, n_perm):
    from scipy.stats import anderson_ksamp
    import warnings
    warnings.filterwarnings("ignore")
    A = np.concatenate([a for a, _ in units]); B = np.concatenate([b for _, b in units])
    ad0 = anderson_ksamp([A, B]).statistic; ex = 0
    for _ in range(n_perm):
        a_, b_ = [], []
        for x, y in units:
            z = np.concatenate([x, y]); pp = rng.permutation(len(z)); a_.append(z[pp[:len(x)]]); b_.append(z[pp[len(x):]])
        ex += anderson_ksamp([np.concatenate(a_), np.concatenate(b_)]).statistic >= ad0
    return ad0, (ex + 1) / (n_perm + 1)


def shape_with_power(usable, rng):
    units = []
    for c in usable:
        for init in NH_INITS:
            nm = member_nam(c, init)
            if nm is None:
                continue
            mn = np.array([x.mean() for x in nm["nudged"].values()]); mc = np.array([x.mean() for x in nm["control"].values()])
            units.append(((mn - mn.mean()) / mn.std(ddof=1), (mc - mc.mean()) / mc.std(ddof=1)))
    qs = np.round(np.arange(0.05, 0.951, 0.05), 2)
    A = np.concatenate([a for a, _ in units]); B = np.concatenate([b for _, b in units])
    dq = np.quantile(A, qs) - np.quantile(B, qs)
    ad0, pval = _ad_perm(units, rng, 2000)
    out = {"n_pairs": len(units), "AD_stat": round(float(ad0), 4), "AD_perm_p": round(float(pval), 4),
           "max_abs_quantile_diff": round(float(np.max(np.abs(dq))), 4),
           "quantile_diff": dict(zip([str(q) for q in qs], np.round(dq, 4).tolist()))}
    sizes = [(len(a), len(b)) for a, b in units]
    def draw(n, w, D, r):
        z = rng.random(n) < w
        return np.where(z, D + r * rng.normal(size=n), rng.normal(size=n))
    def std(x):
        return (x - x.mean()) / x.std(ddof=1)
    power = {}
    for D in (1.5, 2.0, 3.0):
        for r in (0.5, 1.0, 2.0):
            rej = []
            for _ in range(100):
                u = [(std(draw(nf, 0.6, D, r)), std(draw(nr, 0.3, D, r))) for nf, nr in sizes]
                rej.append(_ad_perm(u, rng, 200)[1] < 0.05)
            power[f"D{D}_r{r}"] = round(float(np.mean(rej)), 3)
    out["power"] = power
    ok = pval >= 0.05 and np.max(np.abs(dq)) <= TOL
    out["no_shape_change"] = bool(ok)
    out["excluded_types"] = sorted(k for k, v in power.items() if ok and v >= 0.8)
    print("G8:", {k: out[k] for k in ("n_pairs", "AD_perm_p", "max_abs_quantile_diff", "power", "excluded_types")}, flush=True)
    return out


# ---------------------------------------------------------------- G8x (POST HOC, 2026-10-07)
def _g8x_job(args):
    """Monte Carlo distribution of the G8 Anderson-Darling statistic under one two-type
    mixture, in the observed design (pair sizes); returns P(AD_sim <= AD_obs)."""
    D, r, wc, wn, sizes, ad_obs, seed, n_sim = args
    from scipy.stats import anderson_ksamp
    import warnings
    warnings.filterwarnings("ignore")
    rng = np.random.default_rng(seed)
    def draw(n, w):
        z = rng.random(n) < w
        return np.where(z, D + r * rng.normal(size=n), rng.normal(size=n))
    def std(x):
        return (x - x.mean()) / x.std(ddof=1)
    le = 0
    for _ in range(n_sim):
        A = np.concatenate([std(draw(nn, wn)) for nn, _ in sizes]); B = np.concatenate([std(draw(nc, wc)) for _, nc in sizes])
        le += anderson_ksamp([A, B]).statistic <= ad_obs
    sd_c = np.sqrt((1 - wc) + wc * r * r + wc * (1 - wc) * D * D)
    return D, r, wc, wn, (le + 1) / (n_sim + 1), (wn - wc) * D / sd_c


def g8_inversion(usable, n_sim=400):
    """POST HOC (2026-10-07, after review): G8 excluded types by 'no rejection + power
    >= 0.8 at one weight change', which is not a confidence-set exclusion, and with two
    arms arbitrary two-type mixing cannot be rejected at all. Here: test inversion within
    a stated family -- two Gaussian types D apart (units of the first type's s.d.), width
    ratio r, control weight wc, nudged weight wn, common to all pairs. A family member is
    REJECTED at 5% if the observed AD statistic lies below the 5th percentile of its
    distribution in the observed design (Monte Carlo, n_sim draws, within-arm
    standardisation exactly as G8). Reported: the surviving set, and, against each
    member's implied mean shift (control s.d.), whether mixtures able to produce the
    imposed SSW shift survive."""
    from scipy.stats import anderson_ksamp
    units = []; shifts = []
    for c in usable:
        for init in NH_INITS:
            nm = member_nam(c, init)
            if nm is None:
                continue
            mn = np.array([x.mean() for x in nm["nudged"].values()]); mc = np.array([x.mean() for x in nm["control"].values()])
            units.append(((mn - mn.mean()) / mn.std(ddof=1), (mc - mc.mean()) / mc.std(ddof=1)))
            shifts.append((mn.mean() - mc.mean()) / mc.std(ddof=1))
    A = np.concatenate([a for a, _ in units]); B = np.concatenate([b for _, b in units])
    ad_obs = float(anderson_ksamp([A, B]).statistic)
    sizes = [(len(a), len(b)) for a, b in units]
    ws = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
    grid = [(D, r, wc, wn) for D in (1.0, 1.5, 2.0, 2.5, 3.0, 4.0) for r in (0.5, 1.0, 2.0) for wc in ws for wn in ws if wc != wn]
    jobs = [(D, r, wc, wn, sizes, ad_obs, SEED + 81 + i, n_sim) for i, (D, r, wc, wn) in enumerate(grid)]
    import multiprocessing as mp
    with mp.get_context("fork").Pool(8) as pool:
        outs = pool.map(_g8x_job, jobs)
    obs_shift = float(np.mean(np.abs(shifts)))
    rows = [{"D": D, "r": r, "wc": wc, "wn": wn, "p": round(pv, 4), "implied_shift": round(float(sh), 3), "rejected": bool(pv < 0.05)}
            for D, r, wc, wn, pv, sh in outs]
    surv = [x for x in rows if not x["rejected"]]
    big = [x for x in rows if abs(x["implied_shift"]) >= 0.5 * obs_shift]
    out = {"note": "POST HOC test inversion (2026-10-07); family: two Gaussian types, common weights; replaces the "
                   "'power >= 0.8' exclusion of G8", "AD_obs": round(ad_obs, 4), "n_pairs": len(units), "n_sim": n_sim,
           "observed_mean_abs_shift_sd": round(obs_shift, 3), "n_family": len(rows), "n_rejected": sum(x["rejected"] for x in rows),
           "n_able_to_produce_half_the_shift": len(big), "n_of_those_rejected": sum(x["rejected"] for x in big),
           "survivors_able_to_produce_half_the_shift": [x for x in big if not x["rejected"]],
           "max_abs_implied_shift_among_survivors": round(max((abs(x["implied_shift"]) for x in surv), default=0.0), 3),
           "by_D_r": {f"D{D}_r{r}": {"n": sum(1 for x in rows if x["D"] == D and x["r"] == r),
                                     "rejected": sum(1 for x in rows if x["D"] == D and x["r"] == r and x["rejected"]),
                                     "max_abs_dw_surviving": round(max((abs(x["wn"] - x["wc"]) for x in surv if x["D"] == D and x["r"] == r), default=0.0), 2)}
                      for D in (1.0, 1.5, 2.0, 2.5, 3.0, 4.0) for r in (0.5, 1.0, 2.0)},
           "grid": rows}
    print("G8x:", {k: v for k, v in out.items() if k not in ("grid", "by_D_r")}, flush=True)
    print("G8x by D,r:", out["by_D_r"], flush=True)
    return out


# ---------------------------------------------------------------- G9 (registered 2026-10-07)
def _smooth(x, w=20):
    k = np.ones(w) / w
    pad = np.concatenate([np.full(w // 2, x[0]), x, np.full(w - w // 2 - 1, x[-1])])
    return np.convolve(pad, k, mode="valid")


def _arm_stats(series):
    n0 = min(len(x) for x in series)
    X = np.array([x[:n0] for x in series])
    return X, _smooth(X.mean(0)), _smooth(X.std(0, ddof=1))


def _longest_neg(x):
    best = run = 0
    for v in x:
        run = run + 1 if v < 0 else 0; best = max(best, run)
    return best


def timing_diagnosis(usable, rng, n_rep=200, n_split=10):
    pairs = []
    for c in usable:
        for init in NH_INITS:
            nm = member_nam(c, init)
            if nm is not None:
                pairs.append((list(nm["nudged"].values()), list(nm["control"].values())))
    obs = slope(*zip(*[g7_units(n, k, contrast)[::-1] for n, k in pairs]))
    def synth(con_sel, n_out, Xn_fit, variant):
        Xc, mc, sc = _arm_stats(con_sel)
        Xn, mn, sn = Xn_fit
        n0 = min(Xc.shape[1], len(mn))
        F = Xc[rng.integers(len(Xc), size=n_out), :n0]
        if variant == "c2":
            mbar = Xc.mean(); k = Xn.mean(1).std(ddof=1) / Xc.mean(1).std(ddof=1)
            st = (mn - mc)[:n0]
            return [st + mbar + k * (f.mean() - mbar) + (f - f.mean()) for f in F]
        return [mn[:n0] + (sn[:n0] / sc[:n0]) * (f - mc[:n0]) for f in F]
    out = {"observed_slope_36": round(obs, 4)}
    for variant in ("c2", "c1"):
        sl, runs = [], []
        for _ in range(n_rep):
            xs, ys = [], []
            for nud, con in pairs:
                C = [con[i] for i in rng.integers(len(con), size=len(con))]
                Y = synth(con, len(nud), _arm_stats(nud), variant)
                y_, x_ = g7_units(Y, C, contrast); ys.append(y_); xs.append(x_)
                if variant == "c1" and len(runs) < 36 * 20:
                    runs.append(np.mean([_longest_neg(v) for v in nud]) - np.mean([_longest_neg(v) for v in Y]))
            sl.append(slope(xs, ys))
        sl = np.asarray(sl, float)
        r = [round(float(q), 4) for q in np.nanquantile(sl, [0.025, 0.975])]
        out[variant] = {"mean_slope": round(float(np.nanmean(sl)), 4), "range95": r, "n_rep": n_rep,
                        "observed_inside": bool(r[0] <= obs <= r[1])}
        if variant == "c1":
            out["longest_negative_run_nudged_minus_c1_steps"] = {"mean": round(float(np.mean(runs)), 3),
                                                                 "ci95_over_draws": [round(float(q), 3) for q in np.quantile(runs, [0.025, 0.975])]}
    # held out: fit on half of the nudged members, score the other half
    inside = []; held = []
    for _ in range(n_split):
        split = [rng.permutation(len(n)) for n, _ in pairs]
        A = [[n[i] for i in sp[: len(n) // 2]] for (n, _), sp in zip(pairs, split)]
        Bh = [[n[i] for i in sp[len(n) // 2:]] for (n, _), sp in zip(pairs, split)]
        ob = slope(*zip(*[g7_units(b, k, contrast)[::-1] for b, (_, k) in zip(Bh, pairs)]))
        sl = []
        for _ in range(n_rep // 4):
            xs, ys = [], []
            for a, b, (_, con) in zip(A, Bh, pairs):
                C = [con[i] for i in rng.integers(len(con), size=len(con))]
                Y = synth(con, len(b), _arm_stats(a), "c2")
                y_, x_ = g7_units(Y, C, contrast); ys.append(y_); xs.append(x_)
            sl.append(slope(xs, ys))
        r = np.nanquantile(sl, [0.025, 0.975]); held.append(round(ob, 4)); inside.append(bool(r[0] <= ob <= r[1]))
    out["held_out_c2"] = {"slopes_heldout_half": held, "share_inside": round(float(np.mean(inside)), 3), "n_split": n_split}
    out["reading"] = ("the excess is the timing of the forced response" if out["c2"]["observed_inside"] else
                      "explained by time-varying marginals (timing and spread)" if out["c1"]["observed_inside"] else
                      "unexplained by the marginals: a change of temporal dependence")
    print("G9:", out, flush=True)
    return out


# ---------------------------------------------------------------- G6
G6_REGIONS = {f"{la0}-{la0 + 20}N_{lo0}-{lo0 + 30}E": (la0, la0 + 20, lo0, lo0 + 30)
              for la0 in (30, 50) for lo0 in range(0, 360, 30)}


def g6_series(centre, exp, init):
    """member -> window-mean tas (K) for every G6 region (days 8-24 after onset)."""
    off = L.post_onset_offset(init)
    lo, hi = RT.TAS_DAYS[0] + off, RT.TAS_DAYS[1] + off
    out = {}
    for f in sorted(RT.TAS.glob(f"{centre}_{exp}_{init}_*.npz")):
        with np.load(f) as z:
            day, tas, n, lat, lon = z["lead_day"], z["tas"], z["n_steps"], z["lat"], z["lon"]
        sel = (day >= lo) & (day <= hi)
        if sel.sum() != hi - lo + 1 or (n[sel] != 4).any():
            continue
        field = tas[sel].mean(axis=0)
        row = {}
        for r, (a, b, c, d) in G6_REGIONS.items():
            mi = (lat >= a) & (lat < b); mj = (lon >= c) & (lon < d)
            w = np.cos(np.deg2rad(lat[mi]))[:, None] * np.ones(mj.sum())[None, :]
            row[r] = float((field[np.ix_(mi, mj)] * w).sum() / w.sum())
        out[f.stem.split("_", 3)[3]] = row
    return pd.DataFrame.from_dict(out, orient="index")


def class_contrast(x, thr):
    cold = x < thr
    if cold.sum() < MIN_CLASS or (~cold).sum() < MIN_CLASS:
        return np.nan
    return float(x[cold].mean() - x[~cold].mean())


def regional_classes(usable, rng):
    rows = []
    for c in usable:
        for init in NH_INITS:
            tn, tc = g6_series(c, "nudged", init), g6_series(c, "control", init)
            if len(tn) < 10 or len(tc) < 10:
                continue
            for r in G6_REGIONS:
                mu, sd = tc[r].mean(), tc[r].std(ddof=1)
                xn, xc = ((tn[r] - mu) / sd).values, ((tc[r] - mu) / sd).values
                thr = float(np.quantile(xc, 1 / 3))
                m_n, m_c, k = xn.mean(), xc.mean(), xn.std(ddof=1) / xc.std(ddof=1)
                cn = class_contrast(xn, thr)
                rows.append({"centre": c, "init": init, "region": r, "rate_nudged": float((xn < thr).mean()),
                             "rate_control": float((xc < thr).mean()), "shift": float(m_n - m_c), "sd_ratio": float(k),
                             "resid_law": cn - class_contrast(m_n + k * (xc - m_c), thr),
                             "resid_shift": cn - class_contrast(xc + (m_n - m_c), thr)})
    df = pd.DataFrame(rows)
    out = {"n_region_pairs": int(len(df)), "regions": list(G6_REGIONS)}
    for col in ("resid_law", "resid_shift", "shift", "sd_ratio"):
        d = df.dropna(subset=[col]); cs = d.centre.unique()
        by = {cc: d.loc[d.centre == cc, col].values for cc in cs}
        bs = [np.mean(np.concatenate([by[cc] for cc in rng.choice(cs, len(cs))])) for _ in range(N_BOOT)]
        out[col] = {"mean": round(float(d[col].mean()), 4), "ci95": [round(float(q), 4) for q in np.quantile(bs, [0.025, 0.975])],
                    "n": int(len(d)), "n_centres": int(len(cs))}
    df["rate_change"] = df.rate_nudged - df.rate_control
    out["rate_change"] = {"mean": round(float(df.rate_change.mean()), 4),
                          "abs_mean": round(float(df.rate_change.abs().mean()), 4)}
    ok = df.dropna(subset=["resid_law", "resid_shift"])
    out["share_abs_resid_lt_tol"] = {"law": round(float((ok.resid_law.abs() < TOL).mean()), 4),
                                     "shift_only": round(float((ok.resid_shift.abs() < TOL).mean()), 4), "n": int(len(ok))}
    r = out["resid_law"]
    out["reading"] = ("the law holds for regional impact classes" if -TOL <= r["ci95"][0] and r["ci95"][1] <= TOL
                      else "it fails" if (r["ci95"][0] > 0 or r["ci95"][1] < 0) and abs(r["mean"]) > TOL
                      else "inconclusive")
    out["per_region_mean_resid_law"] = {k: round(float(v), 4) for k, v in df.groupby("region").resid_law.mean().items()}
    print("G6:", {k: out[k] for k in ("n_region_pairs", "resid_law", "resid_shift", "rate_change",
                                      "share_abs_resid_lt_tol", "reading")}, flush=True)
    return out


# ---------------------------------------------------------------- self-test
def selftest():
    """Synthetic pairs: a pure translation must give Delta_b ~ 0 even when the
    shifted ensemble has no NDW members; an added forced spread must be seen."""
    rng = np.random.default_rng(3)
    ok0 = ok1 = 0
    means = []
    for rep in range(100):
        d0, d1 = [], []
        for c in range(9):
            def members(n, sd_offset=1.0):          # member offset + 6-hourly noise
                return [rng.normal(0, sd_offset) + rng.normal(0, 1.2, 72) for _ in range(n)]
            con = members(50)
            s = rng.uniform(-3.5, -2.5)                 # strong enough to leave no NDW
            nud = [x + s for x in members(50)]
            nud_w = [x + s for x in members(50, np.sqrt(2.0))]   # forced variance added
            cc = contrast(con)
            d0.append(contrast([x - s for x in nud]) - cc)
            si = np.mean([x.mean() for x in nud_w])
            d1.append(contrast([x - si for x in nud_w]) - cc)
        means.append(np.nanmean(d0))
        ok0 += abs(means[-1]) < TOL
        ok1 += np.nanmean(d1) < -TOL
    bias = float(np.mean(means))
    print(f"selftest: translation |Delta_b| < {TOL} in {ok0}/100 (bias {bias:+.3f}); "
          f"doubled variance gives Delta_b < -{TOL} in {ok1}/100")
    return ok0 >= 95 and ok1 >= 95 and abs(bias) < 0.03


# ---------------------------------------------------------------- main
def main():
    if "--selftest" in sys.argv:
        return 0 if selftest() else 1
    if "--g9" in sys.argv:                      # registered 2026-10-07; writes only its own block
        res = json.loads((RESULTS / f"{NAME}.json").read_text())
        centres = sorted({p.name.split("_")[0] for p in L.RED.glob("*.parquet")})
        res["G9_timing_diagnosis"] = timing_diagnosis([c for c in centres if not L.corruption_guard(c)[0]],
                                                      np.random.default_rng(SEED + 9))
        (RESULTS / f"{NAME}.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
        print(f"Saved -> {NAME}.json (G9)")
        return 0
    if "--g8x" in sys.argv:                     # POST HOC block only; registered results untouched
        res = json.loads((RESULTS / f"{NAME}.json").read_text())
        centres = sorted({p.name.split("_")[0] for p in L.RED.glob("*.parquet")})
        res["G8x_post_hoc_inversion"] = g8_inversion([c for c in centres if not L.corruption_guard(c)[0]])
        (RESULTS / f"{NAME}.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
        print(f"Saved -> {NAME}.json (G8x)")
        return 0
    rng = np.random.default_rng(SEED)
    sel = json.loads((RESULTS / "snapsi_selection_test.json").read_text())
    strong = {(e["centre"], e["init"]) for e in sel["all_ensembles"]
              if e["hemisphere"] == "NH" and e["arm"] == "nudged" and not e["contrast_estimable"]}
    centres = sorted({p.name.split("_")[0] for p in L.RED.glob("*.parquet")})
    usable = [c for c in centres if not L.corruption_guard(c)[0]]

    rows, sh_rows, reg_rows, ens_strong = [], [], [], []
    for c in usable:
        for init in NH_INITS + [SH_INIT]:
            nm = member_nam(c, init)
            if nm is None:
                continue
            row, (dn, dc) = pair_row(c, init, nm)
            if init == SH_INIT:
                sh_rows.append(row)
                continue
            row["strong"] = (c, init) in strong
            rows.append(row)
            if row["strong"]:
                ens_strong.append((dn, dc))
            rr = regional_row(c, init, nm)
            if rr is not None:
                reg_rows.append(rr)
    df, sh, rg = pd.DataFrame(rows), pd.DataFrame(sh_rows), pd.DataFrame(reg_rows)
    assert len(df) == 36 and df.strong.sum() == 11, (len(df), df.strong.sum())

    res = {"registered": "2026-10-02, committed before the run", "seed": SEED, "n_boot": N_BOOT,
           "tolerance_sigma": TOL, "window_nam": list(L.WINDOW), "window_tas": list(RT.TAS_DAYS),
           "inputs": ["ssw-design-analysis/03_data_ingestion/_snapsi_reduced",
                      "ssw-design-analysis/03_data_ingestion/_snapsi_tas",
                      "results/current/8_experiment/snapsi_selection_test.json"],
           "strong_pairs": sorted(f"{c}|{i}" for c, i in strong)}
    g1 = {}
    for lab, d in (("strong_11", df[df.strong]), ("other_25", df[~df.strong]), ("all_36", df)):
        g1[lab] = {k: cboot(d, k, rng) for k in ("Delta_b", "Delta_b_k1.25", "Delta_b_k1.5", "s", "sd_ratio")}
        g1[lab]["dw_rate_nudged_mean"] = round(float(d.dw_rate_nudged.mean()), 4)
    p = g1["strong_11"]["Delta_b"]
    g1["reading_primary"] = reading(p["ci95"], p["mean"])
    g1["variance_ratio_strong_11"] = pooled_var_ratio(ens_strong, rng)
    g1["composition_strong_11"] = {"by_centre": df[df.strong].centre.value_counts().to_dict(),
                                   "by_init": df[df.strong].init.value_counts().to_dict()}
    res["G1_strongly_forced"] = g1

    res["G2_dose_response"] = {"Delta_b_on_s_36": slope_boot(df, "s", "Delta_b", rng),
                               "paired_diff_on_s_25": slope_boot(df[~df.strong], "s", "paired_diff", rng),
                               "mean_s_strong_11": round(float(df[df.strong].s.mean()), 4),
                               "mean_s_other_25": round(float(df[~df.strong].s.mean()), 4)}
    g2 = res["G2_dose_response"]["paired_diff_on_s_25"]
    res["G2_dose_response"]["paired_diff_fitted_at_strong_mean_s"] = round(
        g2["intercept"] + g2["slope"] * res["G2_dose_response"]["mean_s_strong_11"], 4)

    g3 = {"SH": {k: cboot(sh, k, rng) for k in ("C_nudged", "C_control", "resid_shift", "resid_locscale",
                                                "Delta_b", "s", "sd_ratio")},
          "NH_25": {k: cboot(df[~df.strong], k, rng) for k in ("resid_shift", "resid_locscale", "Delta_b")}}
    a, b = g3["SH"]["resid_shift"], g3["SH"]["resid_locscale"]
    g3["reading"] = ("SH contrast change is the widening of one population"
                     if reading(b["ci95"], b["mean"]) == "within tolerance" and (a["ci95"][0] > 0 or a["ci95"][1] < 0)
                     else f"shift-only {reading(a['ci95'], a['mean'])}; location-scale {reading(b['ci95'], b['mean'])}")
    g3["sh_per_centre"] = sh.round(4).to_dict(orient="records")
    res["G3_third_event_SH"] = g3

    g4 = {"n_pairs_with_tas": int(len(rg)), "mean_control_sd_K": round(float(rg.sd_K.mean()), 4)}
    for k in ("resid", "resid_plant_0.5sigma", "resid_plant_1.32K", "CT_nudged", "CT_null", "CT_control",
              "paired_diff", "Delta_b") + tuple(f"resid_{r}" for r in RT.REGIONS if r != REG):
        g4[k] = cboot(rg, k, rng)
    g4["Delta_b_strong_11"] = cboot(rg[[(c, i) in strong for c, i in zip(rg.centre, rg.init)]], "Delta_b", rng)
    g4["reading_primary"] = reading(g4["resid"]["ci95"], g4["resid"]["mean"])
    g4["planted_detected"] = {k: bool(g4[k]["ci95"][1] < 0) for k in ("resid_plant_0.5sigma", "resid_plant_1.32K")}
    g4["resid_K_approx"] = round(g4["resid"]["mean"] * g4["mean_control_sd_K"], 3)
    res["G4_regional_temperature"] = g4

    g5 = {}
    # own stream (set before G5 first ran), so running G5 leaves every other interval unchanged
    rng5 = np.random.default_rng(SEED + 5)
    for c, inits in FULL.items():
        for init in inits:
            for ref in (["control-full", "control"] if c == "UKMO" else ["control"]):
                nm = member_nam(c, init, arms=("nudged-full", ref), ref=ref)
                if nm is None:
                    continue
                row, _ = pair_row(c, init, nm, "nudged-full", ref)
                g5.setdefault(ref if ref == "control-full" else "control", []).append(row)
    res["G5_nudged_full"] = ({k: {"pairs": pd.DataFrame(v).round(4).to_dict(orient="records"),
                                  **{s: cboot(pd.DataFrame(v), s, rng5) for s in ("paired_diff", "Delta_b", "s")}}
                              for k, v in g5.items()}
                             if g5 else "not run: nudged-full members not cached (CEDA token)")
    # POST HOC (added 2026-10-02 after the registered run, not part of the registration):
    # G1's per-pair Delta_b tracks the change in spread; a contrast proportional to
    # the spread predicts Delta_b = C_control * (sd_ratio - 1). Residual after it:
    df["spread_pred"] = df.C_control * (df.sd_ratio - 1)
    df["Delta_b_net_spread"] = df.Delta_b - df.spread_pred
    res["G1x_post_hoc_spread"] = {
        "note": "post hoc, after the registered run",
        **{lab: {"r_Delta_b_vs_spread_pred": round(float(np.corrcoef(d.Delta_b, d.spread_pred)[0, 1]), 4),
                 "Delta_b_net_spread": cboot(d, "Delta_b_net_spread", rng)}
           for lab, d in (("strong_11", df[df.strong]), ("other_25", df[~df.strong]), ("all_36", df))}}
    # elasticity slope (Delta_b on the spread prediction), post hoc, own random stream
    r11 = np.random.default_rng(SEED + 11)
    b0 = float(np.polyfit(df.spread_pred, df.Delta_b, 1)[0]); bs = []
    cs_ = df.centre.unique()
    for _ in range(N_BOOT):
        dd = pd.concat([df[df.centre == c] for c in r11.choice(cs_, len(cs_))])
        bs.append(np.polyfit(dd.spread_pred, dd.Delta_b, 1)[0])
    res["G1x_post_hoc_spread"]["elasticity_slope_36"] = {
        "slope": round(b0, 4), "ci95": [round(float(q), 4) for q in np.quantile(bs, [0.025, 0.975])]}
    print("G1x post hoc:", res["G1x_post_hoc_spread"])
    # G6 (registered 2026-10-03): regional cold-outcome classes, own random stream
    res["G6_regional_classes"] = regional_classes(usable, np.random.default_rng(SEED + 6))
    # POST HOC (added 2026-10-03 after G5 ran, not part of the registration): per-centre
    # means (with 2-3 centres a centre bootstrap gives only the extreme centre means) and
    # the spread check of G1x for the nudged-full pairs. No random numbers are drawn.
    if g5:
        g5x = {"note": "post hoc, after the registered run"}
        for k, v in g5.items():
            d5 = pd.DataFrame(v)
            d5["spread_pred"] = d5.C_control * (d5.sd_ratio - 1)
            g5x[k] = {"centre_means": {c: {m: round(float(g[m].mean()), 4) for m in ("paired_diff", "Delta_b")}
                                       for c, g in d5.groupby("centre")},
                      "r_Delta_b_vs_spread_pred": round(float(np.corrcoef(d5.Delta_b, d5.spread_pred)[0, 1]), 4)
                      if len(d5) > 2 else None,
                      "Delta_b_net_spread_mean": round(float((d5.Delta_b - d5.spread_pred).mean()), 4)}
        res["G5x_post_hoc"] = g5x
        print("G5x post hoc:", g5x)
    res["G7_elasticity_diagnosis"] = elasticity_diagnosis(usable, strong, np.random.default_rng(SEED + 7))
    res["G8_shape_with_power"] = shape_with_power(usable, np.random.default_rng(SEED + 8))
    res["pairs_NH"] = df.round(4).to_dict(orient="records")
    res["pairs_regional"] = rg.round(4).to_dict(orient="records")

    for lab in ("strong_11", "other_25", "all_36"):
        q = g1[lab]
        print(f"G1 {lab:9s} Delta_b {q['Delta_b']['mean']:+.3f} {q['Delta_b']['ci95']}  "
              f"k1.25 {q['Delta_b_k1.25']['mean']:+.3f} {q['Delta_b_k1.25']['ci95']}  "
              f"k1.5 {q['Delta_b_k1.5']['mean']:+.3f} {q['Delta_b_k1.5']['ci95']}  s {q['s']['mean']:+.2f}")
    print("G1 reading:", g1["reading_primary"], "| var ratio strong:", g1["variance_ratio_strong_11"])
    print("G2:", res["G2_dose_response"])
    print("G3 SH:", {k: (v["mean"], v["ci95"]) for k, v in g3["SH"].items()}, "|", g3["reading"])
    print("G3 NH:", {k: (v["mean"], v["ci95"]) for k, v in g3["NH_25"].items()})
    print("G4:", {k: (v["mean"], v["ci95"], v["n_pairs"]) for k, v in g4.items() if isinstance(v, dict) and "mean" in v})
    print("G4 reading:", g4["reading_primary"], g4["planted_detected"], "resid ~K", g4["resid_K_approx"])
    print("G5:", "run" if g5 else res["G5_nudged_full"])
    (RESULTS / f"{NAME}.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print(f"Saved -> {NAME}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
