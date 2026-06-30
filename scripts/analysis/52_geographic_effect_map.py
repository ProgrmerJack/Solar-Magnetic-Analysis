#!/usr/bin/env python3
"""Script 52: Geographic Effect Map

Creates a geographic map showing the spatial pattern of SSW-avalanche
effects across European Alpine regions and international validation sites.
"""
import json, warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from pathlib import Path

warnings.filterwarnings('ignore')

ROOT = Path(__file__).resolve().parents[2]

# Regional effect data from manuscript and results
# Format: name, lat, lon, effect_size (RR or d), p_value, n_events, metric_type
regions = [
    # Switzerland (primary)
    {'name': 'Switzerland\n(Primary)', 'lat': 46.8, 'lon': 8.2, 
     'effect': 0.32, 'p': 0.001, 'n': 16, 'metric': 'RR', 'category': 'primary'},
    
    # EAWS regions (from SI tables)
    {'name': 'W. Switzerland', 'lat': 46.3, 'lon': 7.0,
     'effect': -0.85, 'p': 0.01, 'n': 2, 'metric': 'd', 'category': 'eaws'},
    {'name': 'E. Switzerland', 'lat': 46.8, 'lon': 9.5,
     'effect': -0.45, 'p': 0.05, 'n': 2, 'metric': 'd', 'category': 'eaws'},
    {'name': 'W. Austria', 'lat': 47.1, 'lon': 11.0,
     'effect': -0.30, 'p': 0.10, 'n': 2, 'metric': 'd', 'category': 'eaws'},
    {'name': 'Bavaria', 'lat': 47.4, 'lon': 11.8,
     'effect': -0.20, 'p': 0.15, 'n': 2, 'metric': 'd', 'category': 'eaws'},
    
    # International validation
    {'name': 'Norway', 'lat': 61.5, 'lon': 8.5,
     'effect': 0.45, 'p': 0.06, 'n': 4, 'metric': 'RR', 'category': 'validation'},
    {'name': 'French\nN. Alps', 'lat': 45.2, 'lon': 6.5,
     'effect': 0.52, 'p': 0.09, 'n': 4, 'metric': 'RR', 'category': 'validation'},
    {'name': 'Utah\n(SNOTEL)', 'lat': 40.5, 'lon': -111.5,
     'effect': 0.48, 'p': 0.08, 'n': 4, 'metric': 'RR', 'category': 'validation'},
    
    # Canada (structured boundary)
    {'name': 'Canada\n(Québec)', 'lat': 48.5, 'lon': -71.0,
     'effect': 0.55, 'p': 0.10, 'n': 2, 'metric': 'RR', 'category': 'boundary_neg'},
    {'name': 'Canada\n(BC/AB)', 'lat': 51.5, 'lon': -116.5,
     'effect': 1.45, 'p': 0.15, 'n': 2, 'metric': 'RR', 'category': 'boundary_pos'},
]

# Create figure with two panels: Europe detail + hemispheric overview
fig, (ax_eur, ax_hem) = plt.subplots(1, 2, figsize=(14, 7),
    gridspec_kw={'width_ratios': [1.3, 1]})

# ============ PANEL A: European Alps Detail ============
ax = ax_eur

# Simple geographic background
ax.set_xlim(4, 16)
ax.set_ylim(44.5, 48.5)
ax.set_aspect(1.4)

# Add Alpine arc approximation
alpine_lons = [6.0, 7.5, 8.5, 10.0, 11.5, 13.0, 14.5, 15.5]
alpine_lats = [44.5, 46.0, 46.5, 47.0, 47.2, 47.1, 46.8, 46.5]
ax.fill_between(alpine_lons, [y-0.5 for y in alpine_lats], alpine_lats, 
                alpha=0.1, color='gray', label='Alpine arc')
ax.plot(alpine_lons, alpine_lats, 'k-', alpha=0.3, linewidth=1)

# Country boundaries (simplified)
# Switzerland
ch_lons = [6.0, 7.5, 8.5, 10.5, 10.5, 9.5, 8.0, 6.0]
ch_lats = [46.0, 45.8, 46.0, 46.5, 47.8, 47.6, 47.8, 47.3]
ax.plot(ch_lons + [ch_lons[0]], ch_lats + [ch_lats[0]], 'k-', alpha=0.2, linewidth=0.5)

# Plot European regions
for r in regions:
    if r['category'] not in ('primary', 'eaws', 'validation') or r['lon'] < 4 or r['lon'] > 16:
        continue
    if r['lat'] < 44.5 or r['lat'] > 48.5:
        continue
    
    # Size based on sample size
    size = max(80, r['n'] * 40)
    
    # Color based on effect direction and significance
    if r['metric'] == 'RR':
        suppressed = r['effect'] < 1
    else:
        suppressed = r['effect'] < 0
    
    if suppressed:
        color = '#2166ac' if r['p'] < 0.05 else '#67a9cf'
    else:
        color = '#b2182b' if r['p'] < 0.05 else '#ef8a62'
    
    edge = 'black' if r['p'] < 0.05 else 'gray'
    lw = 2 if r['category'] == 'primary' else 1
    
    ax.scatter(r['lon'], r['lat'], s=size, c=color, edgecolors=edge, 
              linewidths=lw, zorder=5, alpha=0.85)
    
    # Label
    offset_y = 0.25 if r['lat'] > 46.5 else -0.35
    fontsize = 7 if r['category'] != 'primary' else 9
    fontweight = 'bold' if r['category'] == 'primary' else 'normal'
    
    label = r['name']
    if r['metric'] == 'RR':
        label += f"\nRR={r['effect']:.2f}"
    else:
        label += f"\nd={r['effect']:.2f}"
    
    ax.annotate(label, (r['lon'], r['lat']), 
               xytext=(0, offset_y*60), textcoords='offset points',
               fontsize=fontsize, fontweight=fontweight, ha='center', va='center',
               bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='none'))

ax.set_xlabel('Longitude (°E)', fontsize=10)
ax.set_ylabel('Latitude (°N)', fontsize=10)
ax.set_title('a  European Alps: Spatial Pattern of SSW–Avalanche Response', 
             fontsize=11, fontweight='bold', loc='left')
ax.grid(True, alpha=0.2)

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#2166ac', 
           markeredgecolor='black', markersize=10, label='Suppressed (P<0.05)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#67a9cf',
           markeredgecolor='gray', markersize=10, label='Suppressed (P≥0.05)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#ef8a62',
           markeredgecolor='gray', markersize=10, label='Elevated (P≥0.05)'),
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=8, framealpha=0.9)

# ============ PANEL B: Hemispheric Overview ============
ax2 = ax_hem
ax2.set_xlim(-130, 20)
ax2.set_ylim(35, 70)
ax2.set_aspect(1.4)

# Plot all regions including North America
for r in regions:
    size = max(100, r['n'] * 30)
    
    if r['metric'] == 'RR':
        suppressed = r['effect'] < 1
    else:
        suppressed = r['effect'] < 0
    
    if suppressed:
        color = '#2166ac' if r['p'] < 0.05 else '#67a9cf'
    else:
        color = '#b2182b' if r['p'] < 0.05 else '#ef8a62'
    
    edge = 'black' if r['p'] < 0.05 else 'gray'
    lw = 2 if r['category'] == 'primary' else 1
    
    ax2.scatter(r['lon'], r['lat'], s=size, c=color, edgecolors=edge,
               linewidths=lw, zorder=5, alpha=0.85)
    
    label = r['name'].replace('\n', ' ')
    if r['metric'] == 'RR':
        label += f" RR={r['effect']:.2f}"
    else:
        label += f" d={r['effect']:.2f}"
    
    ax2.annotate(label, (r['lon'], r['lat']),
                xytext=(10, 5), textcoords='offset points',
                fontsize=7, ha='left',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='none'))

# Add continental outlines (very simplified)
ax2.axhline(y=45, color='gray', alpha=0.1, linewidth=0.5)
ax2.axvline(x=-80, color='gray', alpha=0.1, linewidth=0.5)

# Climate zone annotations
ax2.annotate('Continental\nsnow climate', xy=(-115, 53), fontsize=8, 
            color='darkred', ha='center', style='italic', alpha=0.6)
ax2.annotate('Maritime/\ntransitional', xy=(-70, 50), fontsize=8,
            color='darkblue', ha='center', style='italic', alpha=0.6)
ax2.annotate('Alpine\ncontinental', xy=(10, 47.5), fontsize=8,
            color='darkblue', ha='center', style='italic', alpha=0.6)

ax2.set_xlabel('Longitude (°)', fontsize=10)
ax2.set_ylabel('Latitude (°N)', fontsize=10)
ax2.set_title('b  Hemispheric Context: Climate-Zone Dependence', 
             fontsize=11, fontweight='bold', loc='left')
ax2.grid(True, alpha=0.2)

# Add inset text explaining the pattern
textbox = ("Continental snow climates (BC/AB):\n"
           "  SSW → deeper cold → MORE avalanches\n"
           "Maritime/transitional climates:\n"
           "  SSW → trigger suppression → FEWER avalanches")
props = dict(boxstyle='round', facecolor='lightyellow', alpha=0.9, edgecolor='gray')
ax2.text(0.02, 0.02, textbox, transform=ax2.transAxes, fontsize=7,
        verticalalignment='bottom', bbox=props)

plt.tight_layout()

# Save
out_pdf = ROOT / 'data/figures/fig_geographic_effect_map.pdf'
out_png = ROOT / 'data/figures/fig_geographic_effect_map.png'
fig.savefig(out_pdf, dpi=300, bbox_inches='tight')
fig.savefig(out_png, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved to {out_pdf}")
print(f"Saved to {out_png}")
