"""
Recompute BH-FDR correction for 17 primary tests (Norway MW excluded).

Norway MW P<10^-6 was removed from the primary FDR family because it
conflates 4 SSW events into a single day-level comparison, violating
the independence assumption for valid BH correction. The event-level
sign test (4/4, P=0.125) and region-pair sign test (20/23, P=0.0005)
remain in the family.

Output: EDT7 FDR values for paper/main.tex
"""
import numpy as np

raw_p = np.array([
    3.6e-5,   # Combined sign (22/24)
    0.0005,   # Swiss permutation (10,000)
    0.0005,   # Norway pairs sign (20/23)
    0.0006,   # Swiss t-test (log RR) — exact 0.000639
    0.001,    # Swiss Wilcoxon (log RR) — exact 0.001007
    0.004,    # Swiss sign test (14/16) — exact 0.004181
    0.004,    # Phase: onset (14/16)
    0.006,    # Rutschblock MW (2-sided)
    0.021,    # Phase: pre-SSW (13/16)
    0.077,    # Phase: post (12/16)
    0.077,    # Phase: late (12/16)
    0.125,    # Utah sign (4/4)
    0.32,     # ERA5 T2m (event-level)
    0.79,     # AO mediation
    0.86,     # Sintering model
    0.92,     # NAO mediation
    0.999,    # ERA5 snowfall
])

names = [
    "Combined sign (22/24)",
    "Swiss permutation (10,000)",
    "Norway pairs sign (20/23)",
    "Swiss t-test (log RR)",
    "Swiss Wilcoxon (log RR)",
    "Swiss sign test (14/16)",
    "Phase: onset (14/16)",
    "Rutschblock MW (2-sided)",
    "Phase: pre-SSW (13/16)",
    "Phase: post (12/16)",
    "Phase: late (12/16)",
    "Utah sign (4/4)",
    "ERA5 T2m (event-level)",
    "AO mediation",
    "Sintering model",
    "NAO mediation",
    "ERA5 snowfall",
]

m = len(raw_p)
print(f"Number of tests in primary FDR family: {m}")
print(f"(Norway MW excluded from primary family)\n")

# BH step-up
fdr = np.zeros(m)
fdr[-1] = min(1.0, raw_p[-1])
for i in range(m - 2, -1, -1):
    fdr[i] = min(fdr[i + 1], m / (i + 1) * raw_p[i])

print(f"{'Test':<35s}  {'Raw P':>12s}  {'FDR P':>8s}  {'Status':<10s}")
print("-" * 75)
for name, p, f in zip(names, raw_p, fdr):
    sig = "Significant" if f < 0.05 else "Not sig."
    print(f"{name:<35s}  {p:>12.6f}  {f:>8.4f}  {sig}")

# How many significant after FDR?
n_sig = sum(1 for f in fdr if f < 0.05)
print(f"\nSignificant tests at FDR < 0.05: {n_sig}/{m}")
