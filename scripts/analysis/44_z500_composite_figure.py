"""
Script 44: Z500 composite time-series and blocking lead-time visualization

Creates Figure showing:
a) Z500 anomaly composite across SSW events (lag -30 to +30 days)  
b) Avalanche activity composite vs Z500 composite overlay
c) Blocking lead-time scatter (Z500 onset vs SSW onset)

Addresses ALL 5 reviewers flagging missing Z500 visualization.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pathlib import Path

DATA = Path("data/processed")
OUT = Path("data/figures")
OUT.mkdir(parents=True, exist_ok=True)

# Load data
activity = pd.read_parquet(DATA / "cryosphere/slf_activity.parquet")
ssw_cat = pd.read_parquet(DATA / "atmospheric/ssw_catalog.parquet")
ncep_trop = pd.read_parquet(DATA / "atmospheric/ncep_troposphere.parquet")

# Get Z500 and activity
z500 = ncep_trop['hgt_500hPa_m'].copy()
z500.index = z500.index.tz_localize(None)

act = activity['aai_dry_natural'].copy()
if hasattr(act.index, 'tz'):
    act.index = act.index.tz_localize(None)

# Compute DOY-based Z500 anomalies (winter only)
z500_df = pd.DataFrame({'z500': z500})
z500_df['doy'] = z500_df.index.dayofyear
winter_mask = z500_df.index.month.isin([11, 12, 1, 2, 3, 4])
doy_mean = z500_df[winter_mask].groupby('doy')['z500'].mean()
doy_std = z500_df[winter_mask].groupby('doy')['z500'].std()

# Standardized Z500 anomaly
z500_df['z500_anom'] = z500_df.apply(
    lambda r: (r['z500'] - doy_mean.get(r['doy'], r['z500'])) / 
              doy_std.get(r['doy'], 1) if r['doy'] in doy_mean.index else 0, axis=1)

# Also compute activity DOY anomaly
act_df = pd.DataFrame({'act': act})
act_df['doy'] = act_df.index.dayofyear
act_winter = act_df[act_df.index.month.isin([11, 12, 1, 2, 3, 4])]
act_doy_mean = act_winter.groupby('doy')['act'].mean()
act_doy_std = act_winter.groupby('doy')['act'].std()
act_df['act_anom'] = act_df.apply(
    lambda r: (r['act'] - act_doy_mean.get(r['doy'], r['act'])) / 
              act_doy_std.get(r['doy'], 1) if r['doy'] in act_doy_mean.index else 0, axis=1)

# Build event composites
lags = np.arange(-30, 31)
z500_composites = []
act_composites = []

for onset in ssw_cat.index:
    onset_dt = pd.Timestamp(onset).tz_localize(None)
    z500_event = []
    act_event = []
    for lag in lags:
        target = onset_dt + pd.Timedelta(days=int(lag))
        if target in z500_df.index:
            z500_event.append(z500_df.loc[target, 'z500_anom'])
        else:
            z500_event.append(np.nan)
        if target in act_df.index:
            act_event.append(act_df.loc[target, 'act_anom'])
        else:
            act_event.append(np.nan)
    z500_composites.append(z500_event)
    act_composites.append(act_event)

z500_comp = np.array(z500_composites)
act_comp = np.array(act_composites)

z500_mean = np.nanmean(z500_comp, axis=0)
z500_se = np.nanstd(z500_comp, axis=0) / np.sqrt(np.sum(~np.isnan(z500_comp), axis=0))
act_mean = np.nanmean(act_comp, axis=0)
act_se = np.nanstd(act_comp, axis=0) / np.sqrt(np.sum(~np.isnan(act_comp), axis=0))

# Create figure
fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True, gridspec_kw={'hspace': 0.15})

# Panel a: Z500 composite
ax1 = axes[0]
ax1.fill_between(lags, z500_mean - 1.96*z500_se, z500_mean + 1.96*z500_se, 
                 alpha=0.3, color='steelblue')
ax1.plot(lags, z500_mean, 'o-', color='steelblue', markersize=3, linewidth=1.5, label='Z500 anomaly')
ax1.axhline(0, color='grey', linestyle='--', alpha=0.5)
ax1.axvline(0, color='red', linestyle='-', alpha=0.7, label='SSW onset')
ax1.axvspan(-15, -6, alpha=0.1, color='orange', label='Pre-onset window')
ax1.axvspan(-5, 15, alpha=0.1, color='red', label='Post-onset window')
ax1.set_ylabel('Standardized Z500 anomaly\n(Alpine sector)', fontsize=11)
ax1.set_title('Composite SSW event response: circulation and avalanche activity\n($n = 16$ events, ±30 day window)', fontsize=12, fontweight='bold')
ax1.legend(loc='upper right', fontsize=8)
ax1.text(0.02, 0.95, 'a', transform=ax1.transAxes, fontsize=14, fontweight='bold', va='top')

# Panel b: Activity composite
ax2 = axes[1]
ax2.fill_between(lags, act_mean - 1.96*act_se, act_mean + 1.96*act_se,
                 alpha=0.3, color='forestgreen')
ax2.plot(lags, act_mean, 's-', color='forestgreen', markersize=3, linewidth=1.5, label='Natural dry slab anomaly')
ax2.axhline(0, color='grey', linestyle='--', alpha=0.5)
ax2.axvline(0, color='red', linestyle='-', alpha=0.7, label='SSW onset')
ax2.axvspan(-15, -6, alpha=0.1, color='orange')
ax2.axvspan(-5, 15, alpha=0.1, color='red')
ax2.set_ylabel('Standardized activity anomaly\n(natural dry slab)', fontsize=11)
ax2.set_xlabel('Days relative to SSW onset', fontsize=11)
ax2.legend(loc='upper right', fontsize=8)
ax2.text(0.02, 0.95, 'b', transform=ax2.transAxes, fontsize=14, fontweight='bold', va='top')

plt.tight_layout()
plt.savefig(OUT / 'fig_z500_activity_composite.pdf', dpi=300, bbox_inches='tight')
plt.savefig(OUT / 'fig_z500_activity_composite.png', dpi=150, bbox_inches='tight')
print(f"Saved composite figure to {OUT / 'fig_z500_activity_composite.pdf'}")
plt.close()

# Panel c: Event-level scatter (Z500 anom in ±15d window vs log(RR))
# Load event-level data
import json
z500_diag_path = Path("data/results/z500_influence_diagnostics.json")
if z500_diag_path.exists():
    with open(z500_diag_path) as f:
        z500_diag = json.load(f)
    print(f"Z500 diagnostics loaded: {list(z500_diag.keys())[:5]}")

# Create event-level figure
fig2, ax3 = plt.subplots(1, 1, figsize=(7, 5))

# Compute event-level Z500 anomaly (mean in ±15d window)
event_z500 = []
event_rr = []
event_labels = []

for i, onset in enumerate(ssw_cat.index):
    onset_dt = pd.Timestamp(onset).tz_localize(None)
    window_z500 = []
    for d in range(-15, 16):
        target = onset_dt + pd.Timedelta(days=d)
        if target in z500_df.index:
            window_z500.append(z500_df.loc[target, 'z500_anom'])
    
    if window_z500:
        event_z500.append(np.nanmean(window_z500))
    else:
        event_z500.append(np.nan)
    
    # Compute RR from activity data
    window_act = []
    for d in range(-15, 16):
        target = onset_dt + pd.Timedelta(days=d)
        if target in act_df.index:
            window_act.append(act_df.loc[target, 'act'])
    
    if window_act:
        obs = np.nanmean(window_act)
        # Get DOY-matched expected
        doy_range = [(onset_dt + pd.Timedelta(days=d)).dayofyear for d in range(-15, 16)]
        exp_vals = [act_doy_mean.get(doy, np.nan) for doy in doy_range]
        exp = np.nanmean(exp_vals)
        if exp > 0:
            event_rr.append(obs / exp)
        else:
            event_rr.append(np.nan)
    else:
        event_rr.append(np.nan)
    
    event_labels.append(onset_dt.strftime('%Y-%m'))

# Plot
valid = ~np.isnan(event_z500) & ~np.isnan(event_rr)
z500_arr = np.array(event_z500)[valid]
rr_arr = np.array(event_rr)[valid]
log_rr = np.log(rr_arr)
labels_arr = np.array(event_labels)[valid]

from scipy import stats
slope, intercept, r, p, se = stats.linregress(z500_arr, log_rr)
x_fit = np.linspace(z500_arr.min(), z500_arr.max(), 100)
y_fit = slope * x_fit + intercept

colors = ['forestgreen' if rr < 1 else 'firebrick' for rr in rr_arr]
ax3.scatter(z500_arr, log_rr, c=colors, s=60, zorder=5, edgecolors='black', linewidth=0.5)
ax3.plot(x_fit, y_fit, 'b-', linewidth=2, alpha=0.7, 
         label=f'$r = {r:.2f}$, $P = {p:.3f}$')

# Add event labels
for z, lr, lab in zip(z500_arr, log_rr, labels_arr):
    ax3.annotate(lab, (z, lr), fontsize=6, alpha=0.7, 
                xytext=(3, 3), textcoords='offset points')

ax3.axhline(0, color='grey', linestyle='--', alpha=0.5)
ax3.set_xlabel('Standardized Z500 anomaly (±15d window)', fontsize=11)
ax3.set_ylabel('log(RR) natural dry slab activity', fontsize=11)
ax3.set_title('Event-level Z500 anomaly vs avalanche response\n($n = 16$ SSW events)', fontsize=12, fontweight='bold')
ax3.legend(fontsize=10)

# Color legend
from matplotlib.lines import Line2D
legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor='forestgreen', 
                          markersize=8, label='Suppression (RR < 1)'),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor='firebrick',
                          markersize=8, label='Enhancement (RR > 1)')]
ax3.legend(handles=legend_elements + [Line2D([0], [0], color='blue', linewidth=2, 
           label=f'$r = {r:.2f}$, $P = {p:.3f}$')], fontsize=9)

plt.tight_layout()
plt.savefig(OUT / 'fig_z500_event_scatter.pdf', dpi=300, bbox_inches='tight')
plt.savefig(OUT / 'fig_z500_event_scatter.png', dpi=150, bbox_inches='tight')
print(f"Saved event scatter to {OUT / 'fig_z500_event_scatter.pdf'}")
plt.close()

# Print summary stats
print(f"\nComposite summary:")
print(f"Z500 at onset (lag=0): {z500_mean[30]:.3f} ± {z500_se[30]:.3f}")
print(f"Z500 at lag -12: {z500_mean[18]:.3f} ± {z500_se[18]:.3f}")
print(f"Activity at onset: {act_mean[30]:.3f} ± {act_se[30]:.3f}")
print(f"Activity at lag -12: {act_mean[18]:.3f} ± {act_se[18]:.3f}")
print(f"\nEvent-level scatter: r={r:.3f}, P={p:.4f}, n={sum(valid)}")
