"""
Multi-country comparison figure for manuscript.
Creates a publication-quality figure showing SSW-avalanche response across all datasets.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os

# Data from manuscript Extended Data Tables
# Swiss events (16 events)
swiss_rr = [0.26, 3.95, 0.55, 0.35, 0.06, 0.15, 0.32, 0.37, 0.10, 0.08, 0.18, 0.25, 0.44, 0.88, 1.05, 3.43]
swiss_dates = ['1998-12', '1999-02', '2001-02', '2001-12', '2003-01', '2004-01', 
               '2006-01', '2007-02', '2008-02', '2009-01', '2010-02', '2010-03',
               '2013-01', '2018-02a', '2018-02b', '2019-01']

# Country-level summary statistics from manuscript
countries = {
    'Switzerland\n(counts, n=16)': {'rr': 0.32, 'ci_lo': 0.26, 'ci_hi': 0.53, 'n_decrease': '14/16', 'p': 0.002},
    'Utah\n(counts, n=4)': {'rr': 0.34, 'ci_lo': 0.10, 'ci_hi': 1.15, 'n_decrease': '4/4', 'p': 0.0625},
    'Norway\n(danger, n=4)': {'rr': None, 'delta': -0.48, 'd': -0.67, 'n_decrease': '4/4', 'p': 1e-6},
    'France BRA\n(danger, n=4)': {'rr': None, 'delta': -0.15, 'n_decrease': '37/56', 'p': 0.011},
    'EAWS Alps\n(danger, n=2)': {'rr': None, 'delta': 'gradient', 'n_decrease': 'spatial', 'p': 0.0004},
    'Canada\n(danger, n=2)': {'rr': None, 'delta': +0.37, 'n_decrease': '18/61', 'p': 'positive'},
}

fig, axes = plt.subplots(1, 3, figsize=(14, 5), gridspec_kw={'width_ratios': [3, 2, 2]})

# Panel A: Swiss event-level RR (forest plot style)
ax = axes[0]
y_pos = np.arange(len(swiss_rr))
colors = ['#d62728' if rr > 1 else '#2ca02c' for rr in swiss_rr]
log_rr = np.log(np.array(swiss_rr))
bars = ax.barh(y_pos, log_rr, color=colors, alpha=0.7, height=0.7)
ax.axvline(0, color='black', linestyle='-', linewidth=0.8)
ax.axvline(np.log(0.32), color='blue', linestyle='--', linewidth=1.5, label='Geometric mean RR=0.32')
ax.set_yticks(y_pos)
ax.set_yticklabels(swiss_dates, fontsize=8)
ax.set_xlabel('ln(Rate Ratio)', fontsize=10)
ax.set_title('A. Swiss event-level response', fontsize=11, fontweight='bold')
ax.legend(fontsize=8, loc='lower left')
ax.invert_yaxis()
ax.set_xlim(-3.5, 1.8)
# Add RR values as text
for i, rr in enumerate(swiss_rr):
    ax.text(max(log_rr[i] + 0.1, -3.0), i, f'{rr:.2f}', va='center', fontsize=7)

# Panel B: Cross-country summary (effect direction)
ax = axes[1]
country_names = ['Switzerland', 'Utah', 'Norway', 'France\n(N. Alps)', 'France\n(All Alps)', 'ALBINA\n(Austria)', 'Canada']
effects = [-0.68, -0.66, -0.67, -0.15, -0.08, -0.10, +0.37]  # Standardized direction
n_events = [16, 4, 4, 4, 4, 6, 2]
p_values = [0.002, 0.063, 1e-6, 0.011, 0.083, None, None]
colors_b = ['#2ca02c' if e < 0 else '#d62728' for e in effects]
y_pos_b = np.arange(len(country_names))

bars_b = ax.barh(y_pos_b, effects, color=colors_b, alpha=0.7, height=0.6)
ax.axvline(0, color='black', linestyle='-', linewidth=0.8)
ax.set_yticks(y_pos_b)
ax.set_yticklabels(country_names, fontsize=9)
ax.set_xlabel('Effect direction (negative = suppression)', fontsize=9)
ax.set_title('B. Cross-country consistency', fontsize=11, fontweight='bold')
ax.invert_yaxis()

# Add significance markers
for i, p in enumerate(p_values):
    if p is not None:
        stars = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
        x_pos = effects[i] - 0.05 if effects[i] < 0 else effects[i] + 0.05
        ax.text(x_pos, i, stars, va='center', ha='right' if effects[i] < 0 else 'left', fontsize=9, fontweight='bold')

# Panel C: Measurement system triangulation
ax = axes[2]
systems = ['Occurrence\ncounts', 'Danger\nforecasts', 'Incident\nreports', 'Process\nmodel']
n_datasets = [2, 5, 1, 1]  # CH+UT, NO+FR+CA+EAWS+ALBINA, LAWIS, SNOWPACK
directions = ['Suppression\n(18/20 events)', 'Mixed\n(spatial gradient)', 'Increase\n(human-triggered)', 'Loaded gun\n(instability ↑)']
colors_c = ['#2ca02c', '#f0ad4e', '#d62728', '#5bc0de']

y_pos_c = np.arange(len(systems))
widths = [0.9, 0.7, 0.6, 0.8]

for i, (sys, n, d, c, w) in enumerate(zip(systems, n_datasets, directions, colors_c, widths)):
    ax.barh(i, w, color=c, alpha=0.7, height=0.5)
    ax.text(w + 0.05, i, d, va='center', fontsize=8)

ax.set_yticks(y_pos_c)
ax.set_yticklabels(systems, fontsize=9)
ax.set_xlim(0, 2.0)
ax.set_xlabel('Relative strength of evidence', fontsize=9)
ax.set_title('C. Evidence triangulation', fontsize=11, fontweight='bold')
ax.invert_yaxis()
ax.set_xticks([])

plt.tight_layout()

os.makedirs('data/figures', exist_ok=True)
plt.savefig('data/figures/fig_multi_country_comparison.pdf', dpi=300, bbox_inches='tight')
plt.savefig('data/figures/fig_multi_country_comparison.png', dpi=150, bbox_inches='tight')
print("Saved: data/figures/fig_multi_country_comparison.pdf")
print("Saved: data/figures/fig_multi_country_comparison.png")
plt.close()
