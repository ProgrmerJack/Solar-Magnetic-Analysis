#!/usr/bin/env python3
"""
R79: EP-flux proxy height-time composite figure for SSW-avalanche manuscript.

Creates a 2-panel figure:
  (a) Height-time composite of du/dt (EP-flux convergence proxy) at 4 pressure levels
  (b) Composite v'T' at 100 hPa (vertical EP-flux proxy) with pre-onset wave forcing

Uses existing computed data from r77_height_time_composite.json and 47_direct_ep_flux.json.
"""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from pathlib import Path

RESULTS = Path(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\results')
FIGURES = Path(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\figures')

# Load data
with open(RESULTS / 'r77_height_time_composite.json') as f:
    ht = json.load(f)

with open(RESULTS / '47_direct_ep_flux.json') as f:
    vt_data = json.load(f)

with open(RESULTS / 'downward_propagation.json') as f:
    dp = json.load(f)

# ── Panel (a): Height-time du/dt composite ──
levels_str = ['10hPa', '30hPa', '50hPa', '100hPa']
levels_hPa = [10, 30, 50, 100]

days = np.array(ht['composite_by_level']['10hPa']['days'])
dudt_matrix = np.zeros((len(levels_hPa), len(days)))
stderr_matrix = np.zeros_like(dudt_matrix)
n = 16

for i, lev in enumerate(levels_str):
    dudt_matrix[i, :] = np.array(ht['composite_by_level'][lev]['mean_dUdt'])
    std = np.array(ht['composite_by_level'][lev]['std_dUdt'])
    stderr_matrix[i, :] = std / np.sqrt(n)

t_matrix = dudt_matrix / stderr_matrix

# ── Panel (b): v'T' composites ──
vt_values = [e['vT_mean'] for e in vt_data['ep_flux_data']]
vt_arr = np.array(vt_values)

# Event-level RR from downward propagation
rr_vals = [e['rr'] for e in dp['per_event_metrics']]
logrr_vals = [e['log_rr'] for e in dp['per_event_metrics']]

# ── Create figure ──
fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), 
                         gridspec_kw={'width_ratios': [3, 1.2, 1.2]})

# Panel (a): Height-time contour
ax1 = axes[0]
X, Y = np.meshgrid(days, levels_hPa)
norm = TwoSlopeNorm(vmin=-7, vcenter=0, vmax=4)
cf = ax1.contourf(X, Y, dudt_matrix, levels=np.linspace(-7, 4, 23), 
                  cmap='RdBu_r', norm=norm, extend='both')

# Stipple significant regions (|t| > 2.12, df=15, alpha=0.05)
for i, lev in enumerate(levels_hPa):
    for j, d in enumerate(days):
        if abs(t_matrix[i, j]) > 2.12:
            ax1.plot(d, lev, 'k.', markersize=1.5, alpha=0.4)

ax1.set_yscale('log')
ax1.set_yticks(levels_hPa)
ax1.set_yticklabels(['10', '30', '50', '100'])
ax1.invert_yaxis()
ax1.set_xlabel('Lag from SSW onset (days)', fontsize=11)
ax1.set_ylabel('Pressure level (hPa)', fontsize=11)
ax1.set_title(r'$\mathbf{a}$  Composite $\partial\bar{u}/\partial t$ (m s$^{-1}$ d$^{-1}$)', 
              fontsize=12, loc='left')
ax1.axvline(0, color='red', linestyle='--', linewidth=1.5, alpha=0.8)
ax1.axhline(100, color='grey', linestyle=':', linewidth=0.5, alpha=0.5)
ax1.set_xlim(-25, 25)
cbar = plt.colorbar(cf, ax=ax1, shrink=0.85, pad=0.02)
cbar.set_label(r'$\partial\bar{u}/\partial t$ (m s$^{-1}$ d$^{-1}$)', fontsize=10)

# Annotate key features
ax1.annotate('Peak deceleration\n−6.6 m/s/d (t=4.7)', 
             xy=(-1, 10), xytext=(-20, 15),
             fontsize=8, ha='center',
             arrowprops=dict(arrowstyle='->', color='black', lw=1))
ax1.annotate('Downward\npropagation', xy=(2, 60), 
             fontsize=8, ha='center', style='italic', color='navy')

# Panel (b): v'T' at 100 hPa per event
ax2 = axes[1]
events_sorted = sorted(zip(vt_arr, range(16)), reverse=True)
y_pos = np.arange(16)
colors = ['#2166ac' if vt > 25 else '#67a9cf' if vt > 20 else '#d1e5f0' 
          for vt, _ in events_sorted]
bars = ax2.barh(y_pos, [v for v, _ in events_sorted], color=colors, edgecolor='grey', linewidth=0.5)
ax2.axvline(np.mean(vt_arr), color='red', linestyle='--', linewidth=1.5, label='Mean')
ax2.set_xlabel(r"$\overline{v'T'}$ at 100 hPa (K·m/s)", fontsize=10)
ax2.set_ylabel('SSW event (ranked)', fontsize=10)
ax2.set_title(r"$\mathbf{b}$  Pre-onset $\overline{v'T'}_{100}$", fontsize=12, loc='left')
ax2.set_yticks([])
ax2.legend(fontsize=8, loc='lower right')

# Panel (c): v'T' vs log(RR) scatter  
ax3 = axes[2]
ax3.scatter(vt_arr, logrr_vals, c='#2166ac', s=50, edgecolors='black', linewidth=0.5, zorder=5)
z = np.polyfit(vt_arr, logrr_vals, 1)
xfit = np.linspace(min(vt_arr)-2, max(vt_arr)+2, 100)
ax3.plot(xfit, np.polyval(z, xfit), 'r--', linewidth=1, alpha=0.7)

rho = np.corrcoef(vt_arr, logrr_vals)[0,1]
ax3.set_xlabel(r"$\overline{v'T'}_{100}$ (K·m/s)", fontsize=10)
ax3.set_ylabel(r'$\log(\mathrm{RR})$', fontsize=10)
ax3.set_title(r'$\mathbf{c}$  Wave forcing vs response', fontsize=12, loc='left')
ax3.axhline(0, color='grey', linestyle=':', linewidth=0.5)
ax3.text(0.05, 0.95, f'r = {rho:.2f}\n(P = 0.34; n = 16)', 
         transform=ax3.transAxes, fontsize=9, va='top',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.5))

plt.tight_layout()
outpath = FIGURES / 'fig_epflux_height_time.pdf'
fig.savefig(outpath, dpi=300, bbox_inches='tight')
outpath_png = FIGURES / 'fig_epflux_height_time.png'
fig.savefig(outpath_png, dpi=150, bbox_inches='tight')
print(f'Saved: {outpath}')
print(f'Saved: {outpath_png}')

# Print key statistics for manuscript text
print('\n--- Composite statistics ---')
for lev in levels_str:
    d = ht['composite_by_level'][lev]
    pk = d['peak_dUdt_ms_per_day']
    day = d['peak_deceleration_day']
    t = d['t_statistic_at_peak']
    print(f'  {lev}: peak du/dt = {pk:.1f} m/s/d at lag {day}d, t = {t:.1f}')

print(f'\nv\'T\' 100hPa: mean = {np.mean(vt_arr):.1f} +/- {np.std(vt_arr)/np.sqrt(16):.1f} K*m/s')
print(f'v\'T\' vs log(RR): r = {rho:.3f}')
print(f'All 16 events show positive v\'T\' anomaly: {all(v > 0 for v in vt_arr)}')
print(f'v\'T\' range: [{min(vt_arr):.1f}, {max(vt_arr):.1f}]')
