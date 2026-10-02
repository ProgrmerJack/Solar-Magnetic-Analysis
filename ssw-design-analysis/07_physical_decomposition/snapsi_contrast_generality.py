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


# ---------------------------------------------------------------- member data
def member_nam(centre, init, arms=("nudged", "control"), ref="control"):
    """{arm: {member: standardised 6-hourly NAM proxy over the window}}."""
    r = L.load(centre, ref, init)
    if r is None or not L.spans_window(r, init):
        return None
    rm = L.window_means(r, init)
    base, sd = float(rm.mean()), float(rm.std(ddof=1))
    off = L.post_onset_offset(init)
    lo, hi = L.WINDOW[0] + off, L.WINDOW[1] + off
    out = {}
    for arm in arms:
        d = L.load(centre, arm, init)
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
    for c, inits in FULL.items():
        for init in inits:
            for ref in (["control-full", "control"] if c == "UKMO" else ["control"]):
                nm = member_nam(c, init, arms=("nudged-full", ref), ref=ref)
                if nm is None:
                    continue
                row, _ = pair_row(c, init, nm, "nudged-full", ref)
                g5.setdefault(ref if ref == "control-full" else "control", []).append(row)
    res["G5_nudged_full"] = ({k: {"pairs": pd.DataFrame(v).round(4).to_dict(orient="records"),
                                  **{s: cboot(pd.DataFrame(v), s, rng) for s in ("paired_diff", "Delta_b", "s")}}
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
    print("G1x post hoc:", res["G1x_post_hoc_spread"])
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
