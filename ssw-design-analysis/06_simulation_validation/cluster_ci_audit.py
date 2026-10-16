#!/usr/bin/env python3
"""
cluster_ci_audit.py
===================
ARE THE HEADLINE SNAPSI INTERVALS, BUILT BY RESAMPLING EIGHT OR NINE MODELLING CENTRES,
TOO NARROW? AND DOES ANY READING CHANGE UNDER AN INTERVAL THAT IS VALID FOR SO FEW CLUSTERS?

Registered 2026-10-09 (plan approved by the user), committed BEFORE any interval below was
recomputed. Trigger: the validation gate of class_forcing_content P1 found the percentile
centre bootstrap liberal with nine centres (false exclusion 8.5%, FAILURES 2026-10-09).

STATISTICS (per-unit values already stored; the unit is a centre x initialisation pair or,
  for the natural experiment, a mixed ensemble; the estimator is the mean over units)
  S1 paired contrast difference, nudged - control, 25 NH pairs     (generality pairs_NH paired_diff, pairs forming both contrasts)
  S2 shift-null residual C_nudged - C_null_shift, NH                 (pairs_NH resid_shift)
  S3 location-scale residual, NH                                     (pairs_NH resid_locscale)
  S4 back-shift difference, the 11 strongly forced pairs             (pairs_NH Delta_b, strong)
  S5 northern-Eurasian temperature class residual                    (pairs_regional resid)
  S6 natural-experiment law residual, mixed free ensembles           (free ensembles F2_resid_law)
  S7 Southern Hemisphere location-scale residual                     (G3 sh_per_centre resid_locscale)
  Each published value is reproduced first (mean must match to 1e-3), otherwise the audit
  of that statistic stops and the mismatch is reported.
METHODS
  A  percentile bootstrap resampling centres (10,000), as published
  B  cluster-robust t interval: mean over units, variance from centre sums with the
     small-sample factor C/(C-1) (CR1), Student t with C - 1 degrees of freedom
VALIDATION (per statistic, in its own layout: units per centre as observed)
  2,000 synthetic data sets with a true mean of 0: centre effects N(0, s_b^2) plus unit noise,
  s_b and s_w from the one-way random-effects ANOVA of the observed values (s_b^2 floored
  at 0), noise Gaussian and, separately, Student t with 4 d.f. scaled to s_w; coverage of
  the 95% interval by A and B. A method is valid for a statistic if, in both noise cases,
  the one-sided 95% Clopper-Pearson upper bound of non-coverage is below 0.07.
REPORT  for each statistic: both intervals, validity, and its reading under the published
  rule (tolerance +-0.25 sigma: inside -> holds; excluding 0 with |mean| > 0.25 -> fails;
  else inconclusive) with the valid method (B if both are valid); a change of reading is
  flagged. Nothing in the paper is changed by this script.
POST HOC (2026-10-16, after review; labelled in the JSON): the cluster-robust interval of the
  planted-cooling check of S5 (pairs_regional resid_plant_1.32K and resid_plant_0.5sigma).
REGISTERED AMENDMENT R10 (2026-10-16, plan approved by the user; committed before any interval below was
  computed; run with --r10, which leaves every value above unchanged). Disclosed: before registration a
  literature note rescaled the stored S5 interval to a 90% level (about [-0.24, +0.07]); nothing else here
  had been computed. Same seven statistics, same data.
  C  sign-flip randomization interval: test of mean = theta by flipping the signs of the centre sums of
     (x - theta), all 2^C patterns, inverted over a grid of theta (step se/200, +-6 se)
  D  restricted wild cluster bootstrap-t, Rademacher weights, all 2^C patterns, studentised by CR1,
     inverted over the same grid
  E  CR1 with KMA and UKMO as one cluster (both run GloSea6-family systems), t with C' - 1 d.f.
  F  CR1 90% interval (two one-sided tests at 0.05): "equivalent" if inside +-0.25 (registered secondary;
     the registered primary reading stays the 95% one)
  VALIDATION as above (2,000 sets per noise case, gate: one-sided 95% Clopper-Pearson upper bound of
  non-coverage < 0.07 for C, D, E; < 0.12 for the 90% interval F), in three layouts: as observed
  ("indep"), and with KMA and UKMO given one shared centre effect ("shared"), for B and E as well.
  READING per statistic and method under the published rule; "robust" if every valid method gives the
  reading of B; S5's equivalence reading is reported as secondary whatever it is.
Output: results/current/3_calibration/cluster_ci_audit.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
from scipy.stats import beta, t as tdist

ROOT = Path(__file__).resolve().parents[2]
R8 = ROOT / "results" / "current" / "8_experiment"
OUT = ROOT / "results" / "current" / "3_calibration" / "cluster_ci_audit.json"
NAME = "cluster_ci_audit"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
N_BOOT, N_SIM, TOL, GATE = 10000, 2000, 0.25, 0.07


def load():
    g = json.loads((R8 / "snapsi_contrast_generality.json").read_text())
    f = json.loads((R8 / "snapsi_free_natural_experiment.json").read_text())
    nh = g["pairs_NH"]; rg = g["pairs_regional"]
    def col(rows, key, cond=lambda r: True):
        v = [(r["centre"], r[key]) for r in rows if cond(r) and r.get(key) is not None and np.isfinite(r[key])]
        return [c for c, _ in v], np.array([x for _, x in v], float)
    sh = g["G3_third_event_SH"]["sh_per_centre"]
    S = {"S1_paired_diff": (col(nh, "paired_diff"), None),
         "S2_resid_shift_NH": (col(nh, "resid_shift"), g["G3_third_event_SH"]["NH_25"]["resid_shift"]),
         "S3_resid_locscale_NH": (col(nh, "resid_locscale"), g["G3_third_event_SH"]["NH_25"]["resid_locscale"]),
         "S4_Delta_b_strong11": (col(nh, "Delta_b", lambda r: r.get("strong")), g["G1_strongly_forced"]["strong_11"]["Delta_b"]),
         "S5_regional_resid": (col(rg, "resid"), g["G4_regional_temperature"]["resid"]),
         "S6_free_law_resid": (col([e for e in f["ensembles"] if e.get("mixed")], "F2_resid_law"), f["F2_resid_law"]),
         "S7_SH_resid_locscale": (col(sh, "resid_locscale"), g["G3_third_event_SH"]["SH"]["resid_locscale"])}
    sel = json.loads((R8 / "snapsi_selection_test.json").read_text())["paired_NH"]
    S["S1_paired_diff"] = (S["S1_paired_diff"][0], {"mean": sel["paired_difference"], "ci95": sel["paired_difference_CI95_centre_bootstrap"]})
    return S


def boot_ci(c, x, rng, n=N_BOOT):
    """percentile bootstrap of the mean over units, resampling centres (vectorised)."""
    cs = np.unique(c)
    sums = np.array([x[c == k].sum() for k in cs]); cnt = np.array([np.sum(c == k) for k in cs])
    idx = rng.integers(len(cs), size=(n, len(cs)))
    bs = sums[idx].sum(1) / cnt[idx].sum(1)
    return np.quantile(bs, [0.025, 0.975])


def cr1_ci(c, x):
    cs = np.unique(c); C = len(cs); m = x.mean(); N = len(x)
    s = np.array([np.sum(x[c == k] - m) for k in cs])
    se = np.sqrt(C / (C - 1) * np.sum(s ** 2)) / N
    q = tdist.ppf(0.975, C - 1)
    return np.array([m - q * se, m + q * se])


def anova(c, x):
    cs = np.unique(c); n = np.array([np.sum(c == k) for k in cs]); C = len(cs); N = len(x)
    means = np.array([x[c == k].mean() for k in cs])
    ssw = sum(np.sum((x[c == k] - x[c == k].mean()) ** 2) for k in cs); msw = ssw / max(N - C, 1)
    ssb = np.sum(n * (means - x.mean()) ** 2); msb = ssb / (C - 1)
    n0 = (N - np.sum(n ** 2) / N) / (C - 1)
    return max((msb - msw) / n0, 0.0), msw


def reading(m, ci):
    lo, hi = ci
    if -TOL <= lo and hi <= TOL:
        return "holds"
    if (lo > 0 or hi < 0) and abs(m) > TOL:
        return "fails"
    return "inconclusive"


def main():
    rng = np.random.default_rng(SEED)
    res = {"registered": "2026-10-09, committed before any interval was recomputed", "seed": SEED, "n_boot": N_BOOT,
           "n_sim": N_SIM, "statistics": {}}
    ub = lambda k, n: float(beta.ppf(0.95, k + 1, n - k)) if k < n else 1.0
    for name, ((cent, x), pub) in load().items():
        c = np.array(cent); m = float(x.mean())
        r = {"n_units": int(len(x)), "n_centres": int(len(np.unique(c))), "mean": round(m, 4),
             "published": pub}
        if pub and abs(pub["mean"] - m) > 1e-3:
            r["status"] = "STOPPED: published mean not reproduced"
            res["statistics"][name] = r; print(name, r, flush=True); continue
        a = boot_ci(c, x, rng); b = cr1_ci(c, x)
        sb2, sw2 = anova(c, x)
        cov = {}
        for noise in ("gauss", "t4"):
            miss = {"A": 0, "B": 0}
            cs = np.unique(c)
            for _ in range(N_SIM):
                eff = {k: rng.normal(0, np.sqrt(sb2)) for k in cs}
                e = rng.normal(size=len(x)) if noise == "gauss" else rng.standard_t(4, len(x)) / np.sqrt(2.0)
                y = np.array([eff[k] for k in c]) + np.sqrt(sw2) * e
                ca = boot_ci(c, y, rng, 999); cb = cr1_ci(c, y)
                miss["A"] += not (ca[0] <= 0 <= ca[1]); miss["B"] += not (cb[0] <= 0 <= cb[1])
            cov[noise] = {k: {"non_coverage": v / N_SIM, "upper95": round(ub(v, N_SIM), 4)} for k, v in miss.items()}
        valid = {k: all(cov[n][k]["upper95"] < GATE for n in cov) for k in ("A", "B")}
        use = "B" if valid["B"] else ("A" if valid["A"] else None)
        ci_use = {"A": a, "B": b}.get(use)
        r.update({"ci_A_centre_bootstrap": [round(float(v), 4) for v in a], "ci_B_cluster_robust_t": [round(float(v), 4) for v in b],
                  "s_between": round(float(np.sqrt(sb2)), 4), "s_within": round(float(np.sqrt(sw2)), 4),
                  "validation": cov, "valid": valid, "method_used": use,
                  "reading_published_interval": reading(m, pub["ci95"]) if pub else None,
                  "reading_valid_interval": reading(m, ci_use) if use else "no valid method"})
        r["reading_changed"] = bool(pub and use and r["reading_published_interval"] != r["reading_valid_interval"])
        res["statistics"][name] = r
        print(name, json.dumps({k: r[k] for k in ("mean", "ci_A_centre_bootstrap", "ci_B_cluster_robust_t", "valid",
                                                    "reading_published_interval", "reading_valid_interval")}), flush=True)
    g = json.loads((R8 / "snapsi_contrast_generality.json").read_text())
    ph = {}
    for key in ("resid_plant_1.32K", "resid_plant_0.5sigma"):
        v = [(r["centre"], r[key]) for r in g["pairs_regional"] if r.get(key) is not None and np.isfinite(r[key])]
        c = np.array([a for a, _ in v]); x = np.array([b for _, b in v])
        ph[key] = {"mean": round(float(x.mean()), 4), "ci_B_cluster_robust_t": [round(float(q), 4) for q in cr1_ci(c, x)],
                   "n_units": int(len(x)), "n_centres": int(len(np.unique(c)))}
    res["post_hoc_planted_S5"] = ph
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf8", newline="\n")
    print(f"Saved -> {OUT.relative_to(ROOT)}")
    return 0


MERGE = {"KMA": "KMA+UKMO", "UKMO": "KMA+UKMO"}


def _grid(c, x):
    b = cr1_ci(c, x); m = x.mean(); se = (b[1] - m) / tdist.ppf(0.975, len(np.unique(c)) - 1)
    return m, se, m + se * np.arange(-1200, 1201) / 200.0


def _patterns(C):
    return np.array([[1.0 if (i >> k) & 1 else -1.0 for k in range(C)] for i in range(2 ** C)])


def _invert(theta, p):
    keep = theta[p >= 0.05]
    return np.array([keep.min(), keep.max()]) if keep.size else np.array([np.nan, np.nan])


def signflip_ci(c, x):
    cs = np.unique(c); N = len(x); m, se, th = _grid(c, x)
    S = np.array([x[c == k].sum() for k in cs]); n = np.array([np.sum(c == k) for k in cs])
    s = S[None, :] - n[None, :] * th[:, None]                    # grid x C
    T = np.abs(m - th); Ts = np.abs(s @ _patterns(len(cs)).T) / N  # grid x 2^C
    return _invert(th, np.mean(Ts >= T[:, None] - 1e-12, axis=1))


def wcb_ci(c, x):
    cs = np.unique(c); C = len(cs); N = len(x); m, se, th = _grid(c, x)
    S = np.array([x[c == k].sum() for k in cs]); n = np.array([np.sum(c == k) for k in cs])
    t_obs = np.abs(m - th) / se
    s = S[None, :] - n[None, :] * th[:, None]; E = _patterns(C)
    es = s[:, None, :] * E[None, :, :]                           # grid x 2^C x C
    d = es.sum(2) / N
    sest = np.sqrt(C / (C - 1) * np.sum((es - n[None, None, :] * d[:, :, None]) ** 2, axis=2)) / N
    ts = np.abs(d) / np.where(sest > 0, sest, np.inf)
    return _invert(th, np.mean(ts >= t_obs[:, None] - 1e-12, axis=1))


def merged_ci(c, x):
    return cr1_ci(np.array([MERGE.get(k, k) for k in c]), x)


def cr1_ci90(c, x):
    cs = np.unique(c); C = len(cs); m = x.mean(); N = len(x)
    s = np.array([np.sum(x[c == k] - m) for k in cs])
    se = np.sqrt(C / (C - 1) * np.sum(s ** 2)) / N; q = tdist.ppf(0.95, C - 1)
    return np.array([m - q * se, m + q * se])


def equivalence(ci):
    return "equivalent" if -TOL <= ci[0] and ci[1] <= TOL else "not shown equivalent"


def r10():
    rng = np.random.default_rng(SEED + 10)
    res = json.loads(OUT.read_text())
    out = {"registered": "2026-10-16 (amendment R10), committed before any of these intervals was computed",
           "disclosed": "a literature note had rescaled the stored S5 interval to about 90% before registration",
           "n_sim": N_SIM, "statistics": {}}
    ub = lambda k, n: float(beta.ppf(0.95, k + 1, n - k)) if k < n else 1.0
    meth = {"B": cr1_ci, "C": signflip_ci, "D": wcb_ci, "E": merged_ci, "F": cr1_ci90}
    for name, ((cent, x), pub) in load().items():
        c = np.array(cent); m = float(x.mean()); cs = np.unique(c)
        if pub and abs(pub["mean"] - m) > 1e-3:
            continue
        cis = {k: f(c, x) for k, f in meth.items()}
        sb2, sw2 = anova(c, x)
        val = {}
        for lay in ("indep", "shared"):
            for noise in ("gauss", "t4"):
                miss = {k: 0 for k in meth}
                for _ in range(N_SIM):
                    eff = {k: rng.normal(0, np.sqrt(sb2)) for k in cs}
                    if lay == "shared" and "KMA" in eff and "UKMO" in eff:
                        eff["UKMO"] = eff["KMA"]
                    e = rng.normal(size=len(x)) if noise == "gauss" else rng.standard_t(4, len(x)) / np.sqrt(2.0)
                    y = np.array([eff[k] for k in c]) + np.sqrt(sw2) * e
                    for k, f in meth.items():
                        ci = f(c, y); miss[k] += not (ci[0] <= 0 <= ci[1])
                val[f"{lay}_{noise}"] = {k: {"non_coverage": v / N_SIM, "upper95": round(ub(v, N_SIM), 4)} for k, v in miss.items()}
        valid = {k: all(val[s][k]["upper95"] < (0.12 if k == "F" else GATE) for s in val) for k in meth}
        rd = {k: reading(m, cis[k]) for k in "BCDE" if valid[k]}
        r = {"mean": round(m, 4), "n_units": int(len(x)), "n_centres": int(len(cs)),
             "ci": {k: [round(float(v), 4) for v in cis[k]] for k in meth}, "validation": val, "valid": valid,
             "readings_valid_methods": rd, "robust": bool(valid["B"] and all(v == rd["B"] for v in rd.values())),
             "equivalence_90_secondary": equivalence(cis["F"]) if valid["F"] else "F not valid"}
        out["statistics"][name] = r
        print(name, json.dumps({k: r[k] for k in ("mean", "ci", "valid", "readings_valid_methods", "robust", "equivalence_90_secondary")}), flush=True)
    res["r10_few_cluster_robustness"] = out
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf8", newline="\n")
    print(f"Saved -> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(r10() if "--r10" in sys.argv else main())
