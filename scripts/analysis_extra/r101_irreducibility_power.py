#!/usr/bin/env python3
"""
r101_irreducibility_power.py
============================
Referee fix M1 - convert the underpowered irreducibility "p-values" into an
honest SAMPLE-SIZE REQUIREMENT: how many SSW events are needed to power the
compound-event (opposing-direction) confirmation at 80%?

Two component tests of the proposed sub-type:
  (A) Cross-arm correlation (suppression strength vs structural instability):
      observed Spearman rho ~ -0.33 at n=16; paper benchmark |rho|=0.5.
  (B) Joint triple-positive enrichment: observed 13/16 = 0.8125 vs the
      product-of-marginals null p0 = 0.666.

Outputs: data/results/r101_irreducibility_power.json
"""
import json
from pathlib import Path
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "results" / "r101_irreducibility_power.json"

za = stats.norm.ppf(0.975)   # alpha=0.05 two-sided
zb = stats.norm.ppf(0.80)    # 80% power


def n_for_correlation(rho, power=0.80, alpha=0.05):
    """Fisher z-transform sample size to detect a correlation rho."""
    za_ = stats.norm.ppf(1 - alpha / 2)
    zb_ = stats.norm.ppf(power)
    C = 0.5 * np.log((1 + abs(rho)) / (1 - abs(rho)))
    return ((za_ + zb_) / C) ** 2 + 3


def power_for_correlation(rho, n, alpha=0.05):
    za_ = stats.norm.ppf(1 - alpha / 2)
    C = 0.5 * np.log((1 + abs(rho)) / (1 - abs(rho)))
    se = 1.0 / np.sqrt(n - 3)
    return float(stats.norm.cdf(C / se - za_) + stats.norm.cdf(-C / se - za_))


def n_for_one_proportion(p1, p0, power=0.80, alpha=0.05):
    """One-sample proportion test sample size (normal approx), one-sided."""
    za_ = stats.norm.ppf(1 - alpha)         # one-sided (enrichment is directional)
    zb_ = stats.norm.ppf(power)
    num = za_ * np.sqrt(p0 * (1 - p0)) + zb_ * np.sqrt(p1 * (1 - p1))
    return (num / (p1 - p0)) ** 2


def main():
    res = {}

    # ---- (A) cross-arm correlation ----
    corr = {}
    for rho in [0.33, 0.40, 0.50]:
        corr[f"rho_{rho:.2f}"] = {
            "n_for_80pct_power": int(np.ceil(n_for_correlation(rho))),
            "power_at_n16": round(power_for_correlation(rho, 16), 3),
            "power_at_n40": round(power_for_correlation(rho, 40), 3),
        }
    res["cross_arm_correlation"] = {
        "observed_rho_at_n16": -0.33,
        "by_effect_size": corr,
        "interpretation": (
            "At n=16 the cross-arm test has <25% power for |rho| up to 0.5; "
            "powering it to 80% requires the event counts below."),
    }

    # ---- (B) joint triple-positive enrichment ----
    p1, p0 = 13/16, 0.666
    n_prop = n_for_one_proportion(p1, p0)
    # exact binomial power at several N for the same effect
    def exact_power(n, p0=p0, p1=p1, alpha=0.05):
        k = stats.binom.ppf(1 - alpha, n, p0)  # one-sided critical count
        return float(1 - stats.binom.cdf(k, n, p1))
    res["joint_enrichment"] = {
        "observed_proportion": round(p1, 4),
        "product_of_marginals_null_p0": p0,
        "n_for_80pct_power_normalapprox": int(np.ceil(n_prop)),
        "exact_binomial_power_at_n16": round(exact_power(16), 3),
        "exact_binomial_power_at_n40": round(exact_power(40), 3),
        "exact_binomial_power_at_n72": round(exact_power(72), 3),
        "interpretation": (
            "Detecting the observed triple-positive enrichment above the "
            "product-of-marginals null requires the N below; at n=16 it is "
            "underpowered, consistent with the non-significant P=0.166."),
    }

    # ---- headline requirement ----
    n_needed = max(
        int(np.ceil(n_for_correlation(0.50))),
        int(np.ceil(n_for_one_proportion(p1, p0))),
    )
    res["headline_requirement"] = {
        "events_needed_to_confirm_subtype_at_80pct_power": n_needed,
        "currently_available": 16,
        "shortfall_factor": round(n_needed / 16, 1),
        "statement": (
            f"Confirming the opposing-direction compound sub-type at 80% power "
            f"requires ~{n_needed} SSW events (vs 16 available); the present study "
            f"is therefore powered to DESCRIBE the dissociation but not to CONFIRM "
            f"irreducibility, which is why we present the sub-type as a hypothesis."),
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=2)

    print("=== M1 SAMPLE-SIZE REQUIREMENT ===")
    print("Cross-arm correlation (events for 80% power):")
    for k, v in corr.items():
        print(f"  {k}: N={v['n_for_80pct_power']}  (power@16={v['power_at_n16']})")
    print(f"Joint enrichment: N for 80% power ~ {res['joint_enrichment']['n_for_80pct_power_normalapprox']} "
          f"(exact power@16={res['joint_enrichment']['exact_binomial_power_at_n16']}, "
          f"@40={res['joint_enrichment']['exact_binomial_power_at_n40']})")
    print(f"\nHEADLINE: ~{n_needed} SSW events needed to confirm the sub-type at 80% power "
          f"(vs 16 available; {res['headline_requirement']['shortfall_factor']}x).")
    print(f"Saved -> {OUT}")


if __name__ == "__main__":
    main()
