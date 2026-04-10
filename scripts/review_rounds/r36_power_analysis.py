"""
Formal power analysis for multi-country replication.
Shows: (1) What n is needed for 80% power to detect Swiss effect size
       (2) Whether Norway/Utah are consistent with Swiss effect size given their n
       (3) Prediction for future SSW events
"""
import numpy as np
from scipy import stats

print("=" * 60)
print("FORMAL POWER ANALYSIS")
print("=" * 60)

# Swiss parameters
d_swiss = -1.06  # Cohen's d on log(RR)
mean_logRR_swiss = -1.14  # mean log(RR)
sd_logRR_swiss = abs(mean_logRR_swiss / d_swiss)  # ~1.075

print(f"\nSwiss effect: d = {d_swiss:.2f}, mean log(RR) = {mean_logRR_swiss:.2f}")
print(f"SD of log(RR): {sd_logRR_swiss:.3f}")

# Part 1: Sample size needed for sign test (14/16 = 87.5% decrease)
print("\n--- Part 1: Sample size for sign test ---")
p_decrease = 14/16  # observed probability of decrease
for n in range(3, 30):
    # Power = P(k >= ceil(n/2) + 1 | p = p_decrease) for one-sided sign test
    # Actually, P(rejecting H0 of 50%) at alpha=0.05
    from scipy.stats import binomtest, binom
    
    # Critical value: smallest k such that P(X >= k | p=0.5) <= 0.05
    for k_crit in range(n+1):
        if 1 - binom.cdf(k_crit - 1, n, 0.5) <= 0.05:
            break
    
    # Power = P(X >= k_crit | p = p_decrease)
    power = 1 - binom.cdf(k_crit - 1, n, p_decrease)
    if n in [4, 5, 8, 10, 12, 15, 16, 20, 25]:
        print(f"  n={n:2d}: critical k={k_crit}, power = {power:.3f}" + (" ***" if power >= 0.8 else ""))

# Part 2: Power for t-test on log(RR) 
print("\n--- Part 2: Sample size for t-test on log(RR) ---")
for n in [4, 5, 6, 8, 10, 12, 15, 16, 20]:
    # Power of one-sample t-test
    # Under H1: t ~ noncentral t with ncp = d * sqrt(n)
    ncp = abs(d_swiss) * np.sqrt(n)
    t_crit = stats.t.ppf(0.975, df=n-1)
    power = 1 - stats.nct.cdf(t_crit, df=n-1, nc=ncp)
    print(f"  n={n:2d}: ncp={ncp:.2f}, power = {power:.3f}" + (" ***" if power >= 0.8 else ""))

# Part 3: Are Norway/Utah consistent with Swiss?
print("\n--- Part 3: Consistency of Norway/Utah with Swiss ---")

# Norway: 4 events, all 4 decreased, danger-level metric
# If true p_decrease = 0.875 (Swiss), P(4/4 decrease) = 0.875^4 = 0.586
p_swiss = 0.875
p_4_of_4 = p_swiss ** 4
print(f"If p_decrease = {p_swiss}: P(4/4 decrease) = {p_4_of_4:.3f}")
print(f"  Norway (4/4): consistent (P > 0.5)")

# Utah: 4 events, 3 decrease (1 unknown direction exactly)
p_3_of_4 = binom.pmf(3, 4, p_swiss) + binom.pmf(4, 4, p_swiss)
print(f"If p_decrease = {p_swiss}: P(>=3/4 decrease) = {p_3_of_4:.3f}")
print(f"  Utah (>=3/4): consistent (P > 0.5)")

# EAWS: 2 events, both show gradient
# Overall consistency: 32/36 country-event pairs show decrease
k_total, n_total = 32, 36
# Expected under Swiss p_decrease
expected = n_total * p_swiss
print(f"\nCombined: {k_total}/{n_total} decrease (expected {expected:.0f}/{n_total} under Swiss rate)")
_, p_binom = binomtest(k_total, n_total, p_swiss).statistic, binomtest(k_total, n_total, p_swiss).pvalue
print(f"  Binomial test vs Swiss rate: P = {binomtest(k_total, n_total, p_swiss).pvalue:.3f}")
print(f"  Result: {'consistent' if binomtest(k_total, n_total, p_swiss).pvalue > 0.05 else 'inconsistent'}")

# Part 4: Minimum detectable effect size for Norway/Utah
print("\n--- Part 4: Minimum detectable effect (80% power) ---")
for n in [4, 8, 16, 32]:
    # For sign test
    for k_crit in range(n+1):
        if 1 - binom.cdf(k_crit - 1, n, 0.5) <= 0.05:
            break
    
    # Find p such that power = 0.8
    for p_test in np.arange(0.5, 1.0, 0.001):
        power = 1 - binom.cdf(k_crit - 1, n, p_test)
        if power >= 0.8:
            break
    print(f"  n={n:2d}: need p_decrease >= {p_test:.3f} for 80% power (Swiss: {p_swiss})")

# Part 5: Prediction for next 5 SSW events
print("\n--- Part 5: Prediction for next 5 SSW events ---")
print(f"Under Swiss effect (p_decrease = {p_swiss}):")
for k in range(6):
    prob = binom.pmf(k, 5, 1-p_swiss)  # k non-decreasing
    prob_at_least = 1 - binom.cdf(k-1, 5, 1-p_swiss) if k > 0 else 1
    print(f"  P(exactly {k} exceptions out of 5) = {prob:.3f}")
print(f"  P(all 5 decrease) = {p_swiss**5:.3f}")
print(f"  P(4+ decrease) = {sum(binom.pmf(k, 5, p_swiss) for k in [4,5]):.3f}")

print("\n--- Part 6: Effect size estimation precision ---")
# How many events needed for 80% CI width <= 0.5 on d?
for n in [8, 16, 32, 50, 64]:
    se_d = np.sqrt(1/n + d_swiss**2 / (2*n))
    ci_width = 2 * 1.96 * se_d
    print(f"  n={n:2d}: SE(d) = {se_d:.3f}, 95% CI width = {ci_width:.2f}")

print("\nDone.")
