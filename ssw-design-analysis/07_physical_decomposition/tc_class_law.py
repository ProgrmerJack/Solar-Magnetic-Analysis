#!/usr/bin/env python3
"""
tc_class_law.py
===============
DOES THE THRESHOLD-CLASS LAW HOLD OUTSIDE THE STRATOSPHERE, AND CAN ITS TEST TELL A
THRESHOLD ON ONE POPULATION FROM GENUINE TYPES? TROPICAL-CYCLONE CLASSES.

Registered 2026-10-08, committed (with its validation gate) BEFORE any IBTrACS or ERSST
value was fetched (03_data_ingestion/acquire_tc_data.py).

CLASSES (outcome-defined, from the literature)
  RI     rapid intensification: the storm's largest 24-h over-water increase of maximum
         sustained wind >= 30 kt (Kaplan & DeMaria 2003, ~95th percentile of 24-h changes)
  MAJOR  lifetime maximum intensity (LMI) >= 96 kt (Saffir-Simpson category 3+)
  Literature positions: intensification rates show "no evidence of significant departures
  from exponential distributions" (Kowch & Emanuel 2015) -> one population expected;
  the LMI distribution is bimodal, "two types of storms" (RI vs non-RI; Lee et al. 2016)
  -> a genuine two-type mixture is a live possibility. The test is run on both.
DATA  IBTrACS v04r01 (all basins), seasons 1982-2024, USA agency winds (USA_WIND, 1-min,
  kt; the agency field available in every basin); storms reaching >= 34 kt. Per storm:
  Y_RI = largest V(t+24 h) - V(t) over synoptic times (00/06/12/18 UTC) with both points
  over water (DIST2LAND > 0); Y_LMI = maximum USA_WIND. Genesis = first record >= 34 kt;
  SST_g = NOAA ERSST v5 monthly mean (2 deg; the 1-deg OISST v2 monthly file ends in 2023)
  at the genesis month and the nearest cell with a value.
FORCING LEVELS  terciles of SST_g over all storms (cool L, middle M, warm H). Observational:
  levels differ also in basin, latitude and season; this is a test of how the outcome
  distribution changes with the thermodynamic environment, not an intervention.
TESTS (unit: storm; uncertainty: 1,999 bootstrap resamples of whole seasons, applied
  jointly to all levels, because storms of a season share an environment)
  C1 PRIMARY  fixed two-type mixing, raw-moment collinearity T across L, M, H (as
     mixture_moment_test, years pooled within a level), for Y_RI and for Y_LMI; Holm over
     the two outcomes. Category: "incompatible with fixed two-type mixing" (Holm p < 0.05);
     "compatible" (95% interval of rho inside [-0.25, +0.25]; rho = 0 under mixing,
     1 under pure translation); otherwise "unresolved".
  C2 SECONDARY  the class law from L to H: in units of L, H's contrast (class mean minus
     rest) against L's values mapped by H's location and scale and classified at the same
     threshold; residual read against +-0.25 sigma as damip_class_law ("holds"/"fails"/
     "inconclusive"), for RI (on Y_RI) and MAJOR (on Y_LMI).
  C3 GROUND TRUTH (descriptive; the RI label is observable per storm): LMI mean and s.d.
     of RI and of non-RI storms in each level, and Kolmogorov-Smirnov p across levels
     within each type. Fixed types would keep both within-type distributions while the RI
     share changes. AGREEMENT (fixed now): C1 on Y_LMI is "right" if it says incompatible
     and either type's KS p < 0.05 (Holm over the two types), or says compatible and
     neither does; "wrong" in the two opposite cases; otherwise "not decidable".
READING of the pair (fixed now): "the test discriminates" if one outcome is incompatible
  with fixed mixing and the other is compatible; "one population in both" if both are
  incompatible; "mixing in both" if both compatible; otherwise "not resolved".
VALIDATION GATE (--validate; before any data): 1,000 simulations of 43 seasons x 80 storms
  (season effects s.d. 0.3), levels by tercile of a storm covariate that raises the weight
  of a second type (gamma components, separations as LMI types); false-rejection rate of
  the season-bootstrap T test at nominal 0.05: one-sided 95% Clopper-Pearson upper bound
  must be < 0.07. Power reported against one shifted population. A failed gate stops the
  analysis and is reported.
  GATE HISTORY (before any data, before this file was first committed): v1 used the
  centred season-bootstrap p-value of mixture_moment_test; it FAILED (false rejection
  0.058, upper bound 0.0717). With 43 clusters that p-value is liberal. v2, used here:
  Wald statistic T / bootstrap s.d., referred to Student t with (seasons - 1) d.f.;
  the v1 p-value is stored as p_T_centred_v1 and not used.

Output: results/current/8_experiment/tc_class_law.json
Run:    tc_class_law.py --validate | tc_class_law.py
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta, ks_2samp
from scipy.stats import t as t_dist

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "8_experiment"
OUT = RESULTS / "tc_class_law.json"
ING = HERE.parents[0] / "03_data_ingestion"
TRACKS = ING / "tc_storms.parquet"            # written by acquire_tc_data.py
NAME = "tc_class_law"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
B, N_SIM, GATE, TOL = 1999, 1000, 0.07, 0.25
RI_KT, MAJOR_KT = 30.0, 96.0


# ------------------------------------------------------------------ statistics
def level_moments(y, lev, yr, w):
    """weighted mean and second raw moment per level; w = season weights (bootstrap counts)."""
    ww = w[yr]
    out = []
    for k in range(3):
        m = (lev == k) & np.isfinite(y)
        s = ww[m].sum()
        out.append(((ww[m] * y[m]).sum() / s, (ww[m] * y[m] ** 2).sum() / s))
    return out


def T_of(mom):
    (m0, a0), (m1, a1), (m2, a2) = mom
    return (a1 - a0) * (m2 - m0) - (a2 - a0) * (m1 - m0)


def moment_test(y, lev, yr, n_years, rng):
    """y standardised by level 0; season bootstrap; returns T, p, rho and interval."""
    ok = np.isfinite(y)
    m0, s0 = y[ok & (lev == 0)].mean(), y[ok & (lev == 0)].std(ddof=1)
    z = (y - m0) / s0
    one = level_moments(z, lev, yr, np.ones(n_years))
    t0 = T_of(one)
    mu = [m for m, _ in one]
    trans = -(mu[1] - mu[0]) * (mu[2] - mu[0]) * (mu[2] - mu[1])
    tb = np.empty(B); rb = np.empty(B)
    for b in range(B):
        w = np.bincount(rng.integers(n_years, size=n_years), minlength=n_years).astype(float)
        tb[b] = T_of(level_moments(z, lev, yr, w))
        rb[b] = tb[b] / trans if abs(trans) > 1e-12 else np.nan
    p_c = float((1 + np.sum(np.abs(tb - t0) >= abs(t0))) / (B + 1))
    p = float(2 * t_dist.sf(abs(t0) / tb.std(ddof=1), n_years - 1))
    return {"T": round(float(t0), 5), "p_T": round(p, 4), "p_T_centred_v1": round(p_c, 4), "rho": round(float(t0 / trans), 4) if abs(trans) > 1e-12 else None,
            "rho_ci95": [round(float(x), 4) for x in np.nanquantile(rb, [0.025, 0.975])] if abs(trans) > 1e-12 else None,
            "level_means_sd_units": [round(float(x), 4) for x in mu]}


def contrast(x, thr):
    up = x >= thr
    if up.sum() < 3 or (~up).sum() < 3:
        return np.nan
    return float(x[up].mean() - x[~up].mean())


def law_residual(y, lev, yr, n_years, thr, rng):
    def stat(w):
        ww = w[yr]
        L = np.repeat(y[(lev == 0)], ww[lev == 0].astype(int)); H = np.repeat(y[(lev == 2)], ww[lev == 2].astype(int))
        L, H = L[np.isfinite(L)], H[np.isfinite(H)]
        mL, sL = L.mean(), L.std(ddof=1)
        l, h = (L - mL) / sL, (H - mL) / sL; c = (thr - mL) / sL
        k = h.std(ddof=1) / l.std(ddof=1)
        return contrast(h, c) - contrast(h.mean() + k * (l - l.mean()), c), float((H >= thr).mean()), float((L >= thr).mean())
    r0, rh, rl = stat(np.ones(n_years))
    bs = [stat(np.bincount(rng.integers(n_years, size=n_years), minlength=n_years).astype(float))[0] for _ in range(B)]
    lo, hi = np.nanquantile(bs, [0.025, 0.975])
    read = ("holds" if -TOL <= lo and hi <= TOL else "fails" if (lo > 0 or hi < 0) and abs(r0) > TOL else "inconclusive")
    return {"residual_sd_units": round(float(r0), 4), "ci95": [round(float(lo), 4), round(float(hi), 4)], "reading": read,
            "class_rate_L": round(rl, 4), "class_rate_H": round(rh, 4)}


def holm(ps):
    ps = np.asarray(ps); o = np.argsort(ps); adj = np.empty_like(ps); run = 0.0
    for r, i in enumerate(o):
        run = max(run, min(1.0, (len(ps) - r) * ps[i])); adj[i] = run
    return adj


# ------------------------------------------------------------------ validation
def simulate(rng, mixing=True, n_years=43, per=80):
    yr = np.repeat(np.arange(n_years), per)
    cov = rng.normal(size=yr.size) + 0.5 * rng.normal(size=n_years)[yr]     # environment, season-correlated
    lev = np.digitize(cov, np.quantile(cov, [1 / 3, 2 / 3]))
    seff = 0.3 * rng.normal(size=n_years)[yr]
    if mixing:
        p = 1 / (1 + np.exp(-(cov * 1.2 - 0.5)))                            # weight of the second type rises with cov
        typ = rng.random(yr.size) < p
        y = np.where(typ, 3.0 + rng.gamma(4, 0.35, yr.size), rng.gamma(3, 0.4, yr.size)) + seff
    else:
        y = 1.2 * cov + rng.gamma(3, 0.5, yr.size) + seff
    return y, lev, yr, n_years


def validate():
    rng = np.random.default_rng(SEED)
    rej_mix = sum(moment_test(*simulate(rng, True), rng=rng)["p_T"] < 0.05 for _ in range(N_SIM))
    rej_one = sum(moment_test(*simulate(rng, False), rng=rng)["p_T"] < 0.05 for _ in range(N_SIM // 5))
    ub = float(beta.ppf(0.95, rej_mix + 1, N_SIM - rej_mix))
    val = {"n_sim": N_SIM, "B": B, "false_rejection_rate": rej_mix / N_SIM, "upper95": round(ub, 4), "gate": GATE,
           "passes": ub < GATE, "power_vs_one_population": round(rej_one / (N_SIM // 5), 4)}
    print(json.dumps(val, indent=1))
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    res.setdefault("validation_history", []).append({"v1_centred_p": {"false_rejection_rate": 0.058, "upper95": 0.0717, "passes": False}}) if not res.get("validation_history") else None
    res["validation"] = val; res["method_fixed"] = "2026-10-08, before any data"
    OUT.write_text(json.dumps(res, indent=1), encoding="utf8", newline="\n")
    return 0


# ------------------------------------------------------------------ main
def main():
    if "--validate" in sys.argv:
        return validate()
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    if not res.get("validation", {}).get("passes"):
        sys.exit("validation gate not passed (run --validate)")
    d = pd.read_parquet(TRACKS)
    d = d[(d.season >= 1982) & (d.season <= 2024) & np.isfinite(d.sst_g)].reset_index(drop=True)
    cuts = np.quantile(d.sst_g, [1 / 3, 2 / 3])
    lev = np.digitize(d.sst_g.values, cuts)
    seasons = np.sort(d.season.unique()); yr = np.searchsorted(seasons, d.season.values); ny = len(seasons)
    rng = np.random.default_rng(SEED + 1)
    res.update({"registered": "2026-10-08, committed before any data", "seed": SEED, "B": B, "n_storms": int(len(d)),
                "n_seasons": int(ny), "sst_tercile_cuts_C": [round(float(x), 3) for x in cuts],
                "levels": {k: {"n": int((lev == i).sum()), "ri_share": round(float((d.y_ri.values[lev == i] >= RI_KT).mean()), 4),
                               "major_share": round(float((d.y_lmi.values[lev == i] >= MAJOR_KT).mean()), 4),
                               "mean_lmi": round(float(d.y_lmi.values[lev == i].mean()), 2)}
                           for i, k in enumerate(("L", "M", "H"))}})
    c1 = {o: moment_test(d[o].values.astype(float), lev, yr, ny, rng) for o in ("y_ri", "y_lmi")}
    adj = holm([c1[o]["p_T"] for o in c1])
    for (o, c), a in zip(c1.items(), adj):
        ci = c["rho_ci95"]
        c["p_T_holm"] = round(float(a), 4)
        c["category"] = ("incompatible with fixed two-type mixing" if a < 0.05 else
                         "compatible with fixed two-type mixing" if ci and -TOL <= ci[0] and ci[1] <= TOL else "unresolved")
    res["C1_mixing"] = c1
    res["C2_law"] = {"RI": law_residual(d.y_ri.values.astype(float), lev, yr, ny, RI_KT, rng),
                     "MAJOR": law_residual(d.y_lmi.values.astype(float), lev, yr, ny, MAJOR_KT, rng)}
    ri = d.y_ri.values >= RI_KT; g = {}
    for typ, m in (("RI", ri), ("nonRI", ~ri & np.isfinite(d.y_ri.values))):
        per = [d.y_lmi.values[m & (lev == i)] for i in range(3)]
        g[typ] = {"lmi_mean": [round(float(x.mean()), 2) for x in per], "lmi_sd": [round(float(x.std(ddof=1)), 2) for x in per],
                  "n": [int(len(x)) for x in per], "ks_p_L_vs_H": round(float(ks_2samp(per[0], per[2]).pvalue), 6)}
    ks_holm = holm([g["RI"]["ks_p_L_vs_H"], g["nonRI"]["ks_p_L_vs_H"]])
    differ = bool((ks_holm < 0.05).any()); cat = c1["y_lmi"]["category"]
    g["ks_p_holm"] = [round(float(x), 6) for x in ks_holm]
    g["agreement_C1_lmi"] = ("right" if ("incompatible" in cat and differ) or ("compatible with" in cat and not differ) else
                             "wrong" if ("incompatible" in cat and not differ) or ("compatible with" in cat and differ) else
                             "not decidable")
    res["C3_ground_truth"] = g
    cats = [c1["y_ri"]["category"], c1["y_lmi"]["category"]]
    inc = ["incompatible" in c for c in cats]; com = ["compatible with" in c for c in cats]
    res["reading_pair"] = ("the test discriminates" if (inc[0] and com[1]) or (inc[1] and com[0]) else
                           "one population in both" if all(inc) else "mixing in both" if all(com) else "not resolved")
    print(json.dumps({k: res[k] for k in ("levels", "C1_mixing", "C2_law", "C3_ground_truth", "reading_pair")}, indent=1))
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf8", newline="\n")
    print(f"Saved -> {OUT.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
