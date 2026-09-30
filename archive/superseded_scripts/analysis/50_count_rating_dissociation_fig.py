"""
Script 50: Count-Rating Dissociation Figure (main text)

Creates the central "loaded-gun" visualization: simultaneous
natural count decrease + fatal accident increase during SSW windows.
ALL 5 R47 reviewers said this should be in main text.

Output: data/figures/fig_count_rating_dissociation.pdf
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.ndimage import uniform_filter1d
import json

# Load data
slf = pd.read_parquet('data/processed/cryosphere/slf_activity.parquet')
ssw = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw.index = pd.to_datetime(ssw.index, utc=True)

primary_events = ssw.loc['1999':'2019']

# Event-level RR using DOY-matched multi-year baseline
event_metrics = []
for onset in primary_events.index:
    doy = onset.dayofyear
    win_start = onset
    win_end = onset + pd.Timedelta(days=30)
    ssw_data = slf.loc[win_start:win_end]
    if len(ssw_data) < 5:
        continue
    ssw_mean = ssw_data['aai_dry_natural'].mean()
    
    ctrl_vals = []
    for yr in slf.index.year.unique():
        if yr == onset.year:
            continue
        try:
            center = pd.Timestamp(f'{yr}-01-01', tz='UTC') + pd.Timedelta(days=doy-1)
            c_data = slf.loc[center:center + pd.Timedelta(days=30)]
            if len(c_data) >= 5:
                ctrl_vals.append(c_data['aai_dry_natural'].mean())
        except:
            pass
    
    if ctrl_vals:
        ctrl_mean = np.mean(ctrl_vals)
        rr = ssw_mean / ctrl_mean if ctrl_mean > 0 else np.nan
        event_metrics.append({'onset': str(onset.date()), 'count_rr': rr})

df = pd.DataFrame(event_metrics).dropna()

# Composite time series
days_range = np.arange(-30, 61)
composite_count = np.zeros(len(days_range))
composite_n = np.zeros(len(days_range))
clim_count = np.zeros(len(days_range))
clim_n = np.zeros(len(days_range))

for onset in primary_events.index:
    doy_onset = onset.dayofyear
    for i, d in enumerate(days_range):
        day = onset + pd.Timedelta(days=int(d))
        if day in slf.index:
            composite_count[i] += slf.loc[day, 'aai_dry_natural']
            composite_n[i] += 1
        target_doy = doy_onset + d
        if target_doy < 1 or target_doy > 365:
            continue
        for yr in slf.index.year.unique():
            if yr == onset.year:
                continue
            try:
                cday = pd.Timestamp(f'{yr}-01-01', tz='UTC') + pd.Timedelta(days=target_doy-1)
                if cday in slf.index:
                    clim_count[i] += slf.loc[cday, 'aai_dry_natural']
                    clim_n[i] += 1
            except:
                pass

composite_mean = np.where(composite_n > 0, composite_count / composite_n, np.nan)
clim_mean = np.where(clim_n > 0, clim_count / clim_n, np.nan)
composite_smooth = uniform_filter1d(np.nan_to_num(composite_mean), size=7)
clim_smooth = uniform_filter1d(np.nan_to_num(clim_mean), size=7)

# --- FIGURE ---
fig = plt.figure(figsize=(10, 8))
gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.3)

ax1 = fig.add_subplot(gs[0, :])
ax1.plot(days_range, clim_smooth, color='gray', lw=2, ls='--', label='Climatology')
ax1.plot(days_range, composite_smooth, color='steelblue', lw=2, label='SSW composite')
ax1.fill_between(days_range, clim_smooth, composite_smooth,
                 where=composite_smooth < clim_smooth, alpha=0.3, color='steelblue', label='Suppression')
ax1.axvline(0, color='red', ls='--', lw=1.5, label='SSW onset')
ax1.axvspan(0, 30, alpha=0.08, color='red')
ax1.set_xlabel('Days relative to SSW onset')
ax1.set_ylabel('Natural dry slab count\n(daily mean)')
ax1.set_title('a  SSW composite vs climatology (n = 16 events)', fontweight='bold', loc='left')
ax1.legend(frameon=False, loc='upper right')
ax1.set_xlim(-30, 60)

rr_vals = df['count_rr'].values
ax2 = fig.add_subplot(gs[1, 0])
colors = ['#d62728' if x > 1 else '#2ca02c' for x in rr_vals]
si = np.argsort(rr_vals)
ax2.barh(range(len(rr_vals)), rr_vals[si], color=[colors[i] for i in si], alpha=0.8, edgecolor='gray', lw=0.5)
ax2.axvline(1, color='black', lw=1)
ax2.set_xlabel('Rate Ratio (SSW / baseline)')
ax2.set_ylabel('SSW events (sorted)')
ax2.set_title('b  Natural count RR', fontweight='bold', loc='left')
n_supp = (rr_vals < 1).sum()
geo = np.exp(np.log(rr_vals).mean())
ax2.text(0.95, 0.05, f'{n_supp}/{len(rr_vals)} suppressed\nGeo. mean = {geo:.2f}',
         transform=ax2.transAxes, ha='right', va='bottom', fontsize=10,
         bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgreen', alpha=0.5))

acc_results = json.load(open('data/results/32_extended_accident_ssw.json'))
per_event = acc_results['per_event_details']
acc_rr = [np.exp(e['log_rr']) for e in per_event if 'log_rr' in e and np.isfinite(e['log_rr'])]

ax3 = fig.add_subplot(gs[1, 1])
ac = ['#d62728' if x > 1 else '#2ca02c' for x in acc_rr]
asi = np.argsort(acc_rr)
ax3.barh(range(len(acc_rr)), [acc_rr[i] for i in asi], color=[ac[i] for i in asi], alpha=0.8, edgecolor='gray', lw=0.5)
ax3.axvline(1, color='black', lw=1)
ax3.set_xlabel('Fatality Rate Ratio (SSW / baseline)')
ax3.set_ylabel('SSW events (sorted)')
ax3.set_title('c  Fatal accident RR (55 yr, n = 29)', fontweight='bold', loc='left')
n_up = sum(1 for x in acc_rr if x > 1)
ageo = np.exp(np.mean(np.log(acc_rr)))
ax3.text(0.95, 0.05, f'{n_up}/{len(acc_rr)} elevated\nGeo. mean = {ageo:.2f}',
         transform=ax3.transAxes, ha='right', va='bottom', fontsize=10,
         bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.5))

fig.text(0.5, 0.01, 'The "loaded-gun": natural counts decrease (b) while fatal accident rates increase (c)',
         ha='center', fontsize=11, fontstyle='italic')

plt.savefig('data/figures/fig_count_rating_dissociation.pdf', bbox_inches='tight', dpi=300)
plt.savefig('data/figures/fig_count_rating_dissociation.png', bbox_inches='tight', dpi=150)
plt.close()
print(f"Figure saved. Natural: {n_supp}/{len(rr_vals)} suppressed, geo={geo:.3f}. Accidents: {n_up}/{len(acc_rr)} elevated, geo={ageo:.3f}")
