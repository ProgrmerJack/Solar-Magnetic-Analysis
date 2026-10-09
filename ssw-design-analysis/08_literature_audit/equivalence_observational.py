#!/usr/bin/env python3
"""
equivalence_observational.py
============================
ARE THE INCONCLUSIVE OBSERVATIONAL CLASS-CONTRAST TESTS EVIDENCE OF NO SSW-SPECIFIC
DIFFERENCE, OR ONLY INSENSITIVE? Equivalence tests (TOST; Lakens 2017) and Bayes factors
(Dienes 2014) on the matched-null tests of criterion_regional_contrast (ERA5 temperature)
and criterion_precip_contrast (E-OBS rainfall), 39 SSWs.

POST HOC relative to those tests (their p values were seen). Bounds, priors and readings
below were fixed and committed on 2026-10-09 BEFORE any equivalence or Bayes-factor value
was computed. Inputs: the two producers' null replicates (03_data_ingestion/_null_replicates/,
written by a re-run that reproduced their JSONs exactly) and their JSONs.

MODEL  The producers' power analysis plants an extra effect delta on downward events of a
  null replicate, which moves its class contrast by exactly delta. So under effect delta the
  observed contrast C_obs is distributed as C_null + delta, where C_null are the 4,000 null
  replicates (selection plus the ordinary circulation-outcome relation, held-out winters).
TOST  bounds +-D, D = 0.25 s.d. of the outcome on November-March event-free dates (the
  tolerance used throughout the paper). p_lo = P(C_null - D >= C_obs), p_hi = P(C_null + D <=
  C_obs); p_TOST = max. Reported also: the 90% interval of delta by inversion,
  [C_obs - q95(C_null), C_obs - q05(C_null)], and the smallest D that would be declared
  equivalent (max |end|), in outcome units and in s.d.
BAYES FACTOR  H0: delta = 0; H1: delta half-normal in the direction of the alternative the
  power analysis planted (temperature: downward colder, delta < 0; rainfall: downward
  wetter, delta > 0; wet-spell length: longer), SD tau = 0.5 event-free s.d. (the planted
  effect; Dienes's default of SD = the expected effect). Likelihood: Gaussian kernel density of
  C_null evaluated at C_obs - delta. BF01 = L(0) / integral L(delta) p1(delta) d delta.
  Sensitivity: tau = 0.25 and 1.0 s.d.; a Gaussian approximation of C_null.
READINGS (fixed now)  TOST p < 0.05: "equivalent within +-0.25 s.d."; otherwise "equivalence
  not shown". BF01 > 3: "evidence of no SSW-specific difference"; BF01 < 1/3: "evidence of a
  difference"; otherwise "insensitive".
TARGETS  primary: northern-Eurasian temperature days 8-24 (cross-validated null, surface
  criterion) and Iberian rainfall days 1-25; secondary (descriptive, no correction): the other
  regions and windows of both producers, and Iberian wet-spell length.

Output: results/current/9_literature/equivalence_observational.json
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import gaussian_kde, norm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RES = ROOT / "results" / "current" / "9_literature"
NULLS = HERE.parents[0] / "03_data_ingestion" / "_null_replicates"
OUT = RES / "equivalence_observational.json"
D_SD, TAU_SD, TAU_SENS = 0.25, 0.5, (0.25, 1.0)
PRIMARY = ("T|NEURASIA_d8_24", "PR|IBERIA_d1_25")


def bf01(c_obs, tau, sign, dens):
    """BF01 for delta = 0 against delta = sign * |N(0, tau)| (half-normal); dens = density of C_null."""
    a = np.linspace(0.0, 8 * tau, 4001)                       # |delta|
    m1 = np.trapezoid(dens(c_obs - sign * a) * 2 * norm.pdf(a, 0, tau), a)
    return float(dens(np.array([c_obs]))[0] / m1)


def analyse(cn, c_obs, sd, sign):
    D, tau = D_SD * sd, TAU_SD * sd
    p_lo = float(np.mean(cn - D >= c_obs)); p_hi = float(np.mean(cn + D <= c_obs))
    q05, q95 = np.quantile(cn, [0.05, 0.95])
    ci90 = [float(c_obs - q95), float(c_obs - q05)]
    smallest = max(abs(ci90[0]), abs(ci90[1]))
    kde = gaussian_kde(cn)
    gau = norm(cn.mean(), cn.std(ddof=1)).pdf                                # Gaussian approximation
    b = bf01(c_obs, tau, sign, kde)
    sens = {f"tau_{t}sd": round(bf01(c_obs, t * sd, sign, kde), 3) for t in TAU_SENS}
    sens["gaussian_tau_0.5sd"] = round(bf01(c_obs, tau, sign, gau), 3)
    p_tost = max(p_lo, p_hi)
    return {"C_obs": round(c_obs, 4), "C_null_mean": round(float(cn.mean()), 4), "delta_hat": round(float(c_obs - cn.mean()), 4),
            "event_free_sd": round(sd, 4), "D": round(D, 4), "p_TOST": round(p_tost, 4), "delta_ci90": [round(x, 4) for x in ci90],
            "smallest_equivalence_bound": round(smallest, 4), "smallest_bound_sd": round(smallest / sd, 3),
            "tau": round(tau, 4), "direction": "+" if sign > 0 else "-", "BF01": round(b, 3), "BF01_sensitivity": sens,
            "TOST_reading": "equivalent within +-0.25 s.d." if p_tost < 0.05 else "equivalence not shown",
            "BF_reading": ("evidence of no SSW-specific difference" if b > 3 else "evidence of a difference" if b < 1 / 3 else "insensitive"),
            "n_null": int(len(cn))}


def main():
    rt = json.loads((RES / "criterion_regional_contrast.json").read_text())["primary_surface"]["cv"]
    rp = json.loads((RES / "criterion_precip_contrast.json").read_text())["contrast"]
    zt = np.load(NULLS / "criterion_regional_contrast.npz"); zp = np.load(NULLS / "criterion_precip_contrast.npz")
    out = {"label": "post hoc relative to the matched-null tests; bounds, priors and readings committed before computation",
           "D_sd": D_SD, "tau_sd": TAU_SD, "primary": list(PRIMARY), "results": {}}
    for k, rec in rt.items():
        cn = zt[f"primary_surface|cv|{k}"]
        assert abs(cn.mean() - rec["C_null_mean_K"]) < 1e-3
        out["results"][f"T|{k}"] = analyse(cn, rec["C_obs_K"], rec["event_free_sd_K"], -1)
    for k, rec in rp.items():
        if "event_free_sd" not in rec:
            continue
        cn = zp[k]
        assert abs(cn.mean() - rec["C_null_mean"]) < 1e-3
        out["results"][f"PR|{k}"] = analyse(cn, rec["C_obs"], rec["event_free_sd"], +1 if not k.startswith("SCAND") else -1)
        if k.startswith("SCAND"):
            # added after the run: the docstring registers "downward wetter" for all rainfall, the code (also
            # committed before the run) used "drier" for Scandinavia; both directions are reported
            out["results"][f"PR|{k}"]["BF01_docstring_direction_wetter"] = round(bf01(rec["C_obs"], TAU_SD * rec["event_free_sd"], +1, gaussian_kde(cn)), 3)
    for k, v in out["results"].items():
        print(k, {x: v[x] for x in ("delta_hat", "p_TOST", "smallest_bound_sd", "BF01", "TOST_reading", "BF_reading")})
    OUT.write_text(json.dumps(out, indent=1), encoding="utf8", newline="\n")
    print(f"Saved -> {OUT.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
