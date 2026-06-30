"""
Brown's method for combining dependent p-values from multi-system evidence.

Addresses reviewer concern that standard Fisher combination assumes independence
between geographic measurement systems that share the same SSW event calendar.

Brown's method (Brown 1975; Poole et al. 2016) adjusts the Fisher chi-squared
test statistic variance for known inter-system correlations using the 
Kost-McDermott covariance approximation.

Results reported in main text §"Convergent evidence synthesis."
"""
import numpy as np
from scipy import stats
import json

# Three geographically independent measurement systems
streams = {
    "Swiss SLF sign-randomisation": {"p": 0.0006, "type": "output_suppression"},
    "US SNOTEL SWE loading": {"p": 2.5e-7, "type": "input_loading"},
    "French S2M/Crocus danger": {"p": 0.011, "type": "output_suppression"},
}

p_values = [s["p"] for s in streams.values()]
k = len(p_values)

# Standard Fisher combination (assumes independence)
T_fisher = -2 * sum(np.log(p) for p in p_values)
df_fisher = 2 * k
p_fisher = stats.chi2.sf(T_fisher, df_fisher)

# Brown's method with Kost-McDermott covariance approximation
# cov(-2*ln(p_i), -2*ln(p_j)) ≈ 3.25*rho + 0.75*rho^2
# where rho is the correlation between the underlying test statistics

def browns_method(p_values, rho_matrix):
    """Compute Brown's method p-value given correlations between test statistics."""
    k = len(p_values)
    T = -2 * sum(np.log(p) for p in p_values)
    E_T = 2 * k  # Expected value under H0

    # Compute variance adjustment
    total_cov = 0.0
    for i in range(k):
        for j in range(i + 1, k):
            rho = rho_matrix[i][j]
            cov_ij = 3.25 * rho + 0.75 * rho**2
            total_cov += cov_ij

    var_T = 4 * k + 2 * total_cov  # Under independence: 4k; Brown: 4k + 2*sum(cov)

    # Scaled chi-squared approximation
    c = var_T / (2 * E_T)
    f = 2 * E_T**2 / var_T

    p_brown = stats.chi2.sf(T / c, f)
    return {
        "T": float(T),
        "variance": float(var_T),
        "effective_df": float(f),
        "scale_factor": float(c),
        "p_value": float(p_brown),
    }


# Conservative scenario: rho=0.5 (Swiss-French), rho=0.2 (cross-Atlantic)
rho_conservative = [[0, 0.5, 0.2], [0.5, 0, 0.2], [0.2, 0.2, 0]]
result_conservative = browns_method(p_values, rho_conservative)

# Ultra-conservative: rho=0.7 (Swiss-French), rho=0.4 (cross-Atlantic)
rho_ultra = [[0, 0.7, 0.4], [0.7, 0, 0.4], [0.4, 0.4, 0]]
result_ultra = browns_method(p_values, rho_ultra)

# Two-stream (Swiss + French only, same hypothesis direction)
p_two = [0.0006, 0.011]
rho_two_conservative = [[0, 0.5], [0.5, 0]]
rho_two_ultra = [[0, 0.7], [0.7, 0]]
result_two_cons = browns_method(p_two, rho_two_conservative)
result_two_ultra = browns_method(p_two, rho_two_ultra)

results = {
    "method": "Brown's method for dependent p-value combination",
    "reference": "Brown (1975) Biometrics 31:987-992; Poole et al. (2016) Bioinformatics 32:i430-i436",
    "streams": {name: {"p_value": s["p"], "type": s["type"]} for name, s in streams.items()},
    "fisher_independent": {
        "T": float(T_fisher),
        "df": df_fisher,
        "p_value": float(p_fisher),
        "note": "Assumes full independence; likely too liberal",
    },
    "browns_conservative": {
        "correlations": {"Swiss-French": 0.5, "Swiss-US": 0.2, "French-US": 0.2},
        "justification": "European streams share blocking; cross-Atlantic correlation is weak",
        **result_conservative,
    },
    "browns_ultra_conservative": {
        "correlations": {"Swiss-French": 0.7, "Swiss-US": 0.4, "French-US": 0.4},
        "justification": "Extreme dependence well beyond published SSW teleconnection correlations",
        **result_ultra,
    },
    "two_stream_conservative": {
        "streams": "Swiss + French only (same direction: suppression)",
        "correlation": 0.5,
        **result_two_cons,
    },
    "two_stream_ultra_conservative": {
        "streams": "Swiss + French only (same direction: suppression)",
        "correlation": 0.7,
        **result_two_ultra,
    },
    "conclusion": (
        "Even under ultra-conservative dependence (rho=0.7 European, 0.4 cross-Atlantic), "
        "the null hypothesis of no systematic snow/avalanche response to SSW weather is "
        f"rejected at P = {result_ultra['p_value']:.2e}. The two-stream same-direction "
        f"sensitivity (Swiss + French only, rho=0.7) gives P = {result_two_ultra['p_value']:.2e}."
    ),
}

# Print summary
print("=" * 70)
print("BROWN'S METHOD: DEPENDENT P-VALUE COMBINATION")
print("=" * 70)
print(f"\nStreams:")
for name, s in streams.items():
    print(f"  {name}: P = {s['p']:.2e} ({s['type']})")
print(f"\nFisher (independent): P = {p_fisher:.2e} (T={T_fisher:.2f}, df={df_fisher})")
print(f"\nBrown's (conservative, rho=0.5/0.2): P = {result_conservative['p_value']:.2e} (eff.df={result_conservative['effective_df']:.2f})")
print(f"Brown's (ultra-conservative, rho=0.7/0.4): P = {result_ultra['p_value']:.2e} (eff.df={result_ultra['effective_df']:.2f})")
print(f"\nTwo-stream Swiss+French (rho=0.5): P = {result_two_cons['p_value']:.2e}")
print(f"Two-stream Swiss+French (rho=0.7): P = {result_two_ultra['p_value']:.2e}")
print(f"\nConclusion: All scenarios reject H0 at conventional alpha=0.05")

# Save results
output_path = "../../data/results/r86_browns_method.json"
with open(output_path, "w") as f:
    json.dump(results, f, indent=2)
print(f"\nResults saved to {output_path}")
