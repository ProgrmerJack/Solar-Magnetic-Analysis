"""
GEE Power Simulation for SSW-Avalanche Manuscript (R60)

Estimates the statistical power of GEE with winter-level exchangeable
correlation to detect the observed SSW effect (IRR=0.70) given n_eff=16
SSW events across 21 winters.

Key result: Power is only 26%, and the observed P=0.18 is the MEDIAN
expected p-value under H1 — the GEE null is a power artifact, not
evidence against the association.

Output: data/results/r60_gee_power_simulation.json
"""

import numpy as np
from scipy import stats
import json

np.random.seed(42)

true_irr = 0.70
true_beta = np.log(true_irr)
n_events = 16
alpha = 0.05
n_sim = 50000

# Standard errors from actual GEE output
naive_se = 0.21
gee_se = 0.27

# Analytical power
z_crit = stats.norm.ppf(1 - alpha / 2)
ncp = true_beta / gee_se
power_analytical = stats.norm.cdf(ncp - z_crit) + stats.norm.cdf(-ncp - z_crit)

# Parametric bootstrap
power_count = 0
p_values = []
for _ in range(n_sim):
    beta_hat = np.random.normal(true_beta, gee_se)
    z_stat = beta_hat / gee_se
    p_val = 2 * (1 - stats.norm.cdf(abs(z_stat)))
    p_values.append(p_val)
    if p_val < alpha:
        power_count += 1

sim_power = power_count / n_sim
p_values = np.array(p_values)

percentile_of_observed = np.mean(p_values <= 0.18) * 100

# Power requirements
z_80 = stats.norm.ppf(0.80)
beta_needed = -(z_crit + z_80) * gee_se
irr_needed = np.exp(beta_needed)
se_needed = abs(true_beta) / (z_crit + z_80)
n_needed = n_events * (gee_se / se_needed) ** 2

results = {
    "true_irr": true_irr,
    "true_beta": round(true_beta, 4),
    "gee_se": gee_se,
    "naive_se": naive_se,
    "design_effect": round((gee_se / naive_se) ** 2, 2),
    "analytical_power": round(power_analytical, 4),
    "simulated_power": round(sim_power, 4),
    "observed_p": 0.18,
    "percentile_of_observed_p_under_H1": round(percentile_of_observed, 1),
    "irr_needed_for_80pct_power": round(irr_needed, 3),
    "n_events_needed_for_80pct_power": int(np.ceil(n_needed)),
    "median_p_under_H1": round(float(np.median(p_values)), 3),
    "mean_p_under_H1": round(float(np.mean(p_values)), 3),
    "pct_significant_under_H1": round(sim_power * 100, 1),
}

with open("data/results/r60_gee_power_simulation.json", "w") as f:
    json.dump(results, f, indent=2)

print(json.dumps(results, indent=2))
