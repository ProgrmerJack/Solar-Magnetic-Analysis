"""
Script 45: Loaded-gun divergence figure

Shows the dual-metric signature:
- Natural dry slab activity DECREASES during SSW windows
- Fatal accidents INCREASE during SSW windows

This is the central mechanistic prediction and the most compelling visualization
for the loaded-gun concept.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from scipy import stats

DATA = Path("data/processed")
OUT = Path("data/figures")
OUT.mkdir(parents=True, exist_ok=True)

# Load all datasets
activity = pd.read_parquet(DATA / "cryosphere/slf_activity.parquet")
accidents = pd.read_parquet(DATA / "cryosphere/slf_accidents.parquet")
ssw_cat = pd.read_parquet(DATA / "atmospheric/ssw_catalog.parquet")

# Prepare activity data
act = activity['aai_dry_natural'].copy()
if hasattr(act.index, 'tz'):
    act.index = act.index.tz_localize(None)
act_df = pd.DataFrame({'act': act})
act_df['doy'] = act_df.index.dayofyear

# DOY means for activity
act_doy_mean = act_df.groupby('doy')['act'].mean()

# Prepare accident data
acc = accidents.reset_index()
acc['date'] = pd.to_datetime(acc['date']).dt.tz_localize(None)

# Daily accident counts
acc_daily = acc.groupby('date').agg(
    n_acc=('number_dead', 'count'),
    n_dead=('number_dead', 'sum')
).reindex(pd.date_range('1970-01-01', '2025-04-30'), fill_value=0)
acc_daily.index.name = 'date'
acc_daily['doy'] = acc_daily.index.dayofyear

# DOY means for accidents (winter only)
winter_acc = acc_daily[acc_daily.index.month.isin([11,12,1,2,3,4])]
acc_doy_mean = winter_acc.groupby('doy')['n_acc'].mean()

# Compute event-level metrics for each SSW event
events = []
for onset in ssw_cat.index:
    onset_dt = pd.Timestamp(onset).tz_localize(None)
    
    # Activity RR (±15 day window, DOY-matched)
    obs_act = []
    exp_act = []
    for d in range(-15, 16):
        target = onset_dt + pd.Timedelta(days=d)
        if target in act_df.index:
            obs_act.append(act_df.loc[target, 'act'])
            doy = target.dayofyear
            if doy in act_doy_mean.index:
                exp_act.append(act_doy_mean[doy])
    
    if obs_act and exp_act and np.mean(exp_act) > 0:
        act_rr = np.mean(obs_act) / np.mean(exp_act)
    else:
        act_rr = np.nan
    
    # Accident RR (±15 day window, DOY-matched)
    obs_acc = []
    exp_acc = []
    for d in range(-15, 16):
        target = onset_dt + pd.Timedelta(days=d)
        if target in acc_daily.index:
            obs_acc.append(acc_daily.loc[target, 'n_acc'])
            doy = target.dayofyear
            if doy in acc_doy_mean.index:
                exp_acc.append(acc_doy_mean[doy])
    
    if obs_acc and exp_acc and np.mean(exp_acc) > 0:
        acc_rr = np.mean(obs_acc) / np.mean(exp_acc)
    else:
        acc_rr = np.nan
    
    events.append({
        'onset': onset_dt,
        'year': onset_dt.year,
        'act_rr': act_rr,
        'acc_rr': acc_rr,
        'has_activity': not np.isnan(act_rr),
        'label': onset_dt.strftime('%Y-%m')
    })

df = pd.DataFrame(events)

# Print summary
print("=== LOADED-GUN DIVERGENCE ===\n")
print(f"Events with activity data: {df['has_activity'].sum()}")
print(f"Events with accident data: {df['acc_rr'].notna().sum()}")

# Activity RR stats (where available)
act_valid = df[df['has_activity']]
print(f"\nNatural activity (n={len(act_valid)} events):")
print(f"  Mean RR: {act_valid['act_rr'].mean():.3f}")
print(f"  Events with RR < 1: {(act_valid['act_rr'] < 1).sum()}/{len(act_valid)}")

# Accident RR stats
acc_valid = df[df['acc_rr'].notna()]
print(f"\nFatal accidents (n={len(acc_valid)} events):")
print(f"  Mean RR: {acc_valid['acc_rr'].mean():.3f}")  
print(f"  Events with RR > 1: {(acc_valid['acc_rr'] > 1).sum()}/{len(acc_valid)}")

# Events with BOTH metrics
both = df[df['has_activity'] & df['acc_rr'].notna()]
divergent = both[(both['act_rr'] < 1) & (both['acc_rr'] > 1)]
print(f"\nEvents with both metrics: {len(both)}")
print(f"Events showing divergence (act↓, acc↑): {len(divergent)}/{len(both)}")
print(f"Divergence fraction: {len(divergent)/len(both)*100:.0f}%")

# Create figure
fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={'width_ratios': [2, 1]})

# Panel a: Event-level bar chart
ax1 = axes[0]
# Sort by onset date
df_sorted = df.sort_values('onset')

# Plot all events that have accident data
x = np.arange(len(df_sorted))
width = 0.35

# Activity bars (only where available)
act_vals = df_sorted['act_rr'].values
acc_vals = df_sorted['acc_rr'].values

has_act = ~np.isnan(act_vals)
has_acc = ~np.isnan(acc_vals)

# Activity bars
for i, (has, val) in enumerate(zip(has_act, act_vals)):
    if has:
        color = 'forestgreen' if val < 1 else 'lightcoral'
        ax1.bar(i - width/2, np.log2(val), width, color=color, edgecolor='black', linewidth=0.5, alpha=0.8)

# Accident bars
for i, (has, val) in enumerate(zip(has_acc, acc_vals)):
    if has:
        color = 'firebrick' if val > 1 else 'lightblue'
        ax1.bar(i + width/2, np.log2(val), width, color=color, edgecolor='black', linewidth=0.5, alpha=0.8)

ax1.axhline(0, color='black', linewidth=1)
ax1.set_xticks(x)
ax1.set_xticklabels(df_sorted['label'], rotation=90, fontsize=7)
ax1.set_ylabel('$\\log_2$(Rate Ratio)', fontsize=11)
ax1.set_title('The Loaded-Gun Divergence:\nNatural activity ↓ while fatal accidents ↑', fontsize=12, fontweight='bold')

# Custom legend
act_patch = mpatches.Patch(facecolor='forestgreen', edgecolor='black', label='Natural activity (RR < 1)')
act_up = mpatches.Patch(facecolor='lightcoral', edgecolor='black', label='Natural activity (RR > 1)')
acc_patch = mpatches.Patch(facecolor='firebrick', edgecolor='black', label='Fatal accidents (RR > 1)')
acc_down = mpatches.Patch(facecolor='lightblue', edgecolor='black', label='Fatal accidents (RR < 1)')
ax1.legend(handles=[act_patch, act_up, acc_patch, acc_down], fontsize=8, loc='lower left')
ax1.text(0.02, 0.98, 'a', transform=ax1.transAxes, fontsize=14, fontweight='bold', va='top')

# Panel b: Summary violin/box plot
ax2 = axes[1]

# Prepare data for box plots
act_log_rr = np.log2(act_valid['act_rr'].values)
acc_log_rr = np.log2(acc_valid['acc_rr'].dropna().values)

bp = ax2.boxplot([act_log_rr, acc_log_rr], 
                  positions=[1, 2], widths=0.6,
                  patch_artist=True,
                  medianprops=dict(color='black', linewidth=2))

bp['boxes'][0].set_facecolor('forestgreen')
bp['boxes'][0].set_alpha(0.6)
bp['boxes'][1].set_facecolor('firebrick')
bp['boxes'][1].set_alpha(0.6)

# Add individual points
np.random.seed(42)
jitter1 = np.random.normal(0, 0.05, len(act_log_rr))
jitter2 = np.random.normal(0, 0.05, len(acc_log_rr))
ax2.scatter(1 + jitter1, act_log_rr, color='forestgreen', alpha=0.5, s=20, zorder=5)
ax2.scatter(2 + jitter2, acc_log_rr, color='firebrick', alpha=0.5, s=20, zorder=5)

ax2.axhline(0, color='grey', linestyle='--', alpha=0.5)
ax2.set_xticks([1, 2])
ax2.set_xticklabels(['Natural\nactivity\n($n$=16)', 'Fatal\naccidents\n($n$=29)'], fontsize=10)
ax2.set_ylabel('$\\log_2$(Rate Ratio)', fontsize=11)
ax2.set_title('Summary', fontsize=12, fontweight='bold')

# Add statistical annotations
act_median = np.median(act_log_rr)
acc_median = np.median(acc_log_rr)
ax2.annotate(f'Median RR={2**act_median:.2f}', xy=(1, act_median), xytext=(1.3, act_median-0.5),
            fontsize=8, ha='left', arrowprops=dict(arrowstyle='->', color='forestgreen'))
ax2.annotate(f'Median RR={2**acc_median:.2f}', xy=(2, acc_median), xytext=(2.3, acc_median+0.3),
            fontsize=8, ha='left', arrowprops=dict(arrowstyle='->', color='firebrick'))

ax2.text(0.02, 0.98, 'b', transform=ax2.transAxes, fontsize=14, fontweight='bold', va='top')

plt.tight_layout()
plt.savefig(OUT / 'fig_loaded_gun_divergence.pdf', dpi=300, bbox_inches='tight')
plt.savefig(OUT / 'fig_loaded_gun_divergence.png', dpi=150, bbox_inches='tight')
print(f"\nSaved figure to {OUT / 'fig_loaded_gun_divergence.pdf'}")
plt.close()

# Compute correlation between the two metrics for events with both
if len(both) >= 5:
    r, p = stats.pearsonr(both['act_rr'], both['acc_rr'])
    print(f"\nCorrelation (act_rr vs acc_rr) for n={len(both)} events: r={r:.3f}, P={p:.3f}")
