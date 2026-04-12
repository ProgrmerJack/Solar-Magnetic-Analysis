"""
forest_plot_cross_regional.py — Forest plot of cross-regional SSW-avalanche effect sizes
==========================================================================================
Generates a forest plot showing country-level effect sizes with 95% CIs,
making geographic heterogeneity (including Canada's opposite sign) immediately visible.

Output: data/figures/forest_plot_cross_regional.pdf
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Data: country/dataset, effect size (log RR or standardised d), 95% CI lower, upper, metric type
data = [
    # Country, estimate, CI_lo, CI_hi, label, n_events
    ("Switzerland (occurrence)", np.log(0.32), np.log(0.20), np.log(0.54), "log(RR)", 16),
    ("Utah (occurrence)", np.log(0.34), np.log(0.08), np.log(1.40), "log(RR)", 4),
    ("Norway (danger Δ)", -0.67, -1.35, 0.01, "Cohen's d", 4),
    ("France–N. Alps (danger Δ)", -0.15, -0.40, 0.10, "Δ danger", 56),
    ("France–all Alps (danger Δ)", -0.08, -0.25, 0.09, "Δ danger", 88),
    ("Canada (danger Δ)", 0.37, -0.10, 0.84, "Δ danger", 2),
]

fig, ax = plt.subplots(figsize=(8, 5))

y_positions = list(range(len(data) - 1, -1, -1))
colors = ['#2166ac', '#2166ac', '#4393c3', '#4393c3', '#92c5de', '#d6604d']

for i, (name, est, lo, hi, metric, n) in enumerate(data):
    y = y_positions[i]
    color = colors[i]
    ax.plot([lo, hi], [y, y], color=color, linewidth=2, solid_capstyle='round')
    ax.plot(est, y, 'o', color=color, markersize=8, zorder=5)
    ax.text(-3.0, y, f"{name} (n={n})", va='center', ha='left', fontsize=9)
    sign = '+' if est > 0 else ''
    ax.text(2.0, y, f"{sign}{est:.2f} [{lo:.2f}, {hi:.2f}]", va='center', ha='left', fontsize=8, family='monospace')

ax.axvline(0, color='grey', linestyle='--', linewidth=0.8, zorder=0)
ax.set_xlim(-3.2, 3.5)
ax.set_ylim(-0.8, len(data) - 0.2)
ax.set_yticks([])
ax.set_xlabel('Effect size (log RR or standardised difference)', fontsize=10)
ax.set_title('Cross-regional SSW–avalanche effect sizes', fontsize=11, fontweight='bold')

ax.annotate('← Suppression', xy=(-1.5, -0.6), fontsize=8, color='#2166ac', ha='center')
ax.annotate('Enhancement →', xy=(1.5, -0.6), fontsize=8, color='#d6604d', ha='center')

plt.tight_layout()
out_path = "data/figures/forest_plot_cross_regional.pdf"
plt.savefig(out_path, dpi=300, bbox_inches='tight')
print(f"Saved: {out_path}")
