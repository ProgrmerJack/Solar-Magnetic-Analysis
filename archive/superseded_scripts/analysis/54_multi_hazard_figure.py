#!/usr/bin/env python3
"""Script 54: Multi-Hazard Regime Shift Figure

Creates a publication-quality figure showing how SSW-driven regime
redistribution simultaneously shifts multiple hazard-relevant weather
variables. This is the key figure for Impact/Scope arguments.
"""
import json, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from scipy import stats
from pathlib import Path

warnings.filterwarnings('ignore')

ROOT = Path(__file__).resolve().parents[2]

# Load cross-hazard results
with open(ROOT / 'data/results/51_cross_hazard_regime.json') as f:
    xhaz = json.load(f)

# Create figure
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# ============ Panel A: Multi-hazard regime shift bar chart ============
ax = axes[0, 0]

hazards = [
    ('Temperature\n(cold-wave risk)', xhaz['cold_wave']['difference_K'], 'K', xhaz['cold_wave']['p_value']),
    ('Rain fraction\n(flood/landslide)', (xhaz['rain_on_snow']['ratio'] - 1) * 100, '%', 0.1),
    ('Snow depth\n(spring flood)', xhaz['snow_loading']['change_pct'], '%', xhaz['snow_loading']['p_value']),
    ('Wind speed\n(infrastructure)', xhaz['wind_hazard']['change_pct'], '%', xhaz['wind_hazard']['p_value']),
    ('Cold stagnation\n(air quality)', (xhaz['air_quality']['ratio'] - 1) * 100, '%', 0.1),
]

names = [h[0] for h in hazards]
values = [h[1] for h in hazards]
p_vals = [h[3] for h in hazards]

colors = ['#2166ac' if v < 0 else '#b2182b' for v in values]
bars = ax.barh(range(len(names)), values, color=colors, edgecolor='black', linewidth=0.5, alpha=0.8)

# Add significance markers
for i, (val, p) in enumerate(zip(values, p_vals)):
    marker = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    x_pos = val + (0.3 if val >= 0 else -0.3)
    ax.text(x_pos, i, marker, ha='left' if val >= 0 else 'right', va='center', fontsize=9, fontweight='bold')

ax.set_yticks(range(len(names)))
ax.set_yticklabels(names, fontsize=9)
ax.axvline(0, color='black', linewidth=0.8)
ax.set_xlabel('Change during SSW windows', fontsize=10)
ax.set_title('a  Multi-hazard weather regime shift', fontsize=11, fontweight='bold', loc='left')
ax.invert_yaxis()

# ============ Panel B: Event-level consistency ============
ax = axes[0, 1]

events = xhaz['multi_hazard_summary']['event_details']
event_dates = [e['event_date'] for e in events]
temp_anoms = [e['temp_anom_K'] for e in events]
precip_ratios = [e['precip_ratio'] for e in events]

# Scatter: temperature anomaly vs precipitation ratio
ax.scatter(temp_anoms, precip_ratios, s=80, c='steelblue', edgecolors='black', 
          linewidths=0.5, alpha=0.8, zorder=5)

# Add reference lines
ax.axhline(1, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
ax.axvline(0, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)

# Quadrant labels
ax.text(-2.5, 0.6, 'Cold & Dry\n(loaded-gun)', fontsize=8, ha='center', color='#2166ac', fontweight='bold')
ax.text(2, 1.3, 'Warm & Wet\n(trigger-rich)', fontsize=8, ha='center', color='#b2182b', fontweight='bold')
ax.text(-2.5, 1.3, 'Cold & Wet\n(snow loading)', fontsize=8, ha='center', color='gray')
ax.text(2, 0.6, 'Warm & Dry\n(fire risk)', fontsize=8, ha='center', color='gray')

n_cold_dry = sum(1 for t, p in zip(temp_anoms, precip_ratios) if t < 0 and p < 1)
ax.text(0.98, 0.02, f'{n_cold_dry}/{len(events)} events in\ncold-dry quadrant',
       transform=ax.transAxes, fontsize=8, ha='right', va='bottom',
       bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

ax.set_xlabel('Temperature anomaly (K)', fontsize=10)
ax.set_ylabel('Precipitation ratio (SSW/control)', fontsize=10)
ax.set_title('b  Event-level hazard-vector consistency', fontsize=11, fontweight='bold', loc='left')

# ============ Panel C: Loaded-gun mechanism schematic ============
ax = axes[1, 0]
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.set_aspect('equal')
ax.axis('off')

# Draw mechanism boxes
boxes = [
    (1.5, 8.5, 'SSW Event\n(Stratospheric\nwarming)', '#fee090', 'bold'),
    (1.5, 6.0, 'Regime\nRedistribution\n(cold-dry ↑)', '#abd9e9', 'bold'),
    (5.0, 3.5, 'Triggers ↓\n(warming −57%\nrain-on-snow −29%\nsolar −17.5%)', '#d73027', 'normal'),
    (1.5, 3.5, 'Loading ↑\n(snowfall +5%\nwind +4.4%\nPWL +9.4%)', '#4575b4', 'normal'),
    (3.2, 1.0, 'LOADED GUN\nNatural release ↓ (RR=0.32)\nHuman-triggered risk ↑', '#fee090', 'bold'),
]

for x, y, text, color, weight in boxes:
    bbox = dict(boxstyle='round,pad=0.4', facecolor=color, alpha=0.8, edgecolor='black')
    ax.text(x, y, text, ha='center', va='center', fontsize=8, fontweight=weight, bbox=bbox)

# Arrows
arrow_style = dict(arrowstyle='->', color='black', linewidth=1.5)
ax.annotate('', xy=(1.5, 7.2), xytext=(1.5, 7.8), arrowprops=arrow_style)
ax.annotate('', xy=(1.5, 4.5), xytext=(1.5, 5.2), arrowprops=arrow_style)
ax.annotate('', xy=(5.0, 4.5), xytext=(3.5, 5.2), arrowprops=arrow_style)
ax.annotate('', xy=(3.2, 2.0), xytext=(2.0, 2.8), arrowprops=arrow_style)
ax.annotate('', xy=(3.2, 2.0), xytext=(4.5, 2.8), arrowprops=arrow_style)

# Generalisation arrow
ax.annotate('Generalises to:\n• Permafrost thaw\n• River ice jams\n• Wildfire weather',
           xy=(8.0, 1.0), fontsize=7, ha='center', va='center',
           bbox=dict(boxstyle='round,pad=0.3', facecolor='#f7f7f7', alpha=0.9, edgecolor='gray'))
ax.annotate('', xy=(7.0, 1.0), xytext=(5.5, 1.0), 
           arrowprops=dict(arrowstyle='->', color='gray', linewidth=1, linestyle='--'))

ax.set_title('c  Loaded-gun mechanism: decoupled loading–trigger', fontsize=11, fontweight='bold', loc='left')

# ============ Panel D: Simultaneous hazard shifts histogram ============
ax = axes[1, 1]

sim_shifts = xhaz['multi_hazard_summary']
events_data = sim_shifts['event_details']

# Calculate shifts per event
shifts_per_event = []
for e in events_data:
    n = 0
    if e['temp_anom_K'] < -0.5:
        n += 1
    if e['precip_ratio'] < 0.85:
        n += 1
    if e['wind_change_pct'] > 5:
        n += 1
    if e['snow_depth_change_pct'] > 2:
        n += 1
    shifts_per_event.append(n)

bins = range(0, 5)
counts = [shifts_per_event.count(i) for i in bins]
colors_hist = ['#f7f7f7', '#d1e5f0', '#67a9cf', '#2166ac']
ax.bar(bins, counts, color=colors_hist[:len(bins)], edgecolor='black', linewidth=0.5, width=0.7)

for i, c in enumerate(counts):
    if c > 0:
        ax.text(i, c + 0.2, str(c), ha='center', fontsize=10, fontweight='bold')

ax.set_xlabel('Number of simultaneous hazard-domain shifts', fontsize=10)
ax.set_ylabel('Number of SSW events', fontsize=10)
ax.set_xticks(range(4))
ax.set_title('d  Multi-hazard simultaneity per SSW event', fontsize=11, fontweight='bold', loc='left')

mean_shifts = np.mean(shifts_per_event)
frac_multi = np.mean([s >= 2 for s in shifts_per_event])
ax.text(0.98, 0.98, f'Mean shifts: {mean_shifts:.1f}\n≥2 shifts: {frac_multi*100:.0f}% of events',
       transform=ax.transAxes, fontsize=9, ha='right', va='top',
       bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

plt.tight_layout()

# Save
out_pdf = ROOT / 'data/figures/fig_multi_hazard_regime.pdf'
out_png = ROOT / 'data/figures/fig_multi_hazard_regime.png'
fig.savefig(out_pdf, dpi=300, bbox_inches='tight')
fig.savefig(out_png, dpi=150, bbox_inches='tight')
plt.close()

print(f"Saved to {out_pdf}")
print(f"Saved to {out_png}")
