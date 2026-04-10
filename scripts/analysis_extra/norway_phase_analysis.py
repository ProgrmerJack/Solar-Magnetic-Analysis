"""
Norwegian NVE phase decomposition analysis.
Analyze pre-SSW vs post-SSW danger level patterns.
Also: daily lag analysis for mechanism timing.
"""
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import binom
import json

# Load Norwegian NVE data
nve = pd.read_csv('data/cryosphere/norway_nve/nve_ssw_analysis.csv')
print(f'NVE data: {len(nve)} rows')
print(f'Columns: {nve.columns.tolist()}')
print(f'SSW dates: {nve["ssw_date"].unique()}')
print(f'Window types: {nve["window_type"].unique()}')
print()

# SSW events in Norwegian data
ssw_dates = nve['ssw_date'].unique()
print(f'SSW events: {len(ssw_dates)}')

# Phase analysis: pre vs post for each event
for sd in ssw_dates:
    if sd == 'control':
        continue
    event = nve[nve['ssw_date'] == sd]
    pre = event[event['window_type'] == 'ssw']
    
    # Parse dates to compute pre/post
    event_dates = pd.to_datetime(event['date'])
    ssw_dt = pd.to_datetime(sd)
    
    pre_mask = event_dates < ssw_dt
    post_mask = event_dates >= ssw_dt
    
    pre_data = event[pre_mask.values]
    post_data = event[post_mask.values]
    
    pre_danger = pre_data['danger_level'].mean() if len(pre_data) > 0 else np.nan
    post_danger = post_data['danger_level'].mean() if len(post_data) > 0 else np.nan
    
    print(f'  {sd}: PRE danger={pre_danger:.2f} (n={len(pre_data)}) '
          f'POST danger={post_danger:.2f} (n={len(post_data)})')

# Overall pre vs post
all_ssw = nve[nve['window_type'] == 'ssw']
all_ctrl = nve[nve['window_type'] == 'ctrl']

print(f'\n=== OVERALL ===')
print(f'SSW window mean danger: {all_ssw["danger_level"].mean():.3f} (n={len(all_ssw)})')
print(f'Control mean danger:    {all_ctrl["danger_level"].mean():.3f} (n={len(all_ctrl)})')

# Pre-SSW vs Post-SSW across all events
all_dates = pd.to_datetime(all_ssw['date'])
all_ssw_dates = pd.to_datetime(all_ssw['ssw_date'])
pre_mask = all_dates < all_ssw_dates
post_mask = all_dates >= all_ssw_dates

pre_danger = all_ssw.loc[pre_mask.values, 'danger_level']
post_danger = all_ssw.loc[post_mask.values, 'danger_level']

print(f'\nPRE-SSW danger: {pre_danger.mean():.3f} (n={len(pre_danger)})')
print(f'POST-SSW danger: {post_danger.mean():.3f} (n={len(post_danger)})')
u_stat, p_val = stats.mannwhitneyu(pre_danger, post_danger, alternative='two-sided')
print(f'MW test (PRE vs POST): U={u_stat:.0f}, P={p_val:.4f}')

# Pre-SSW vs control
u_pre, p_pre = stats.mannwhitneyu(pre_danger, all_ctrl['danger_level'], alternative='less')
print(f'PRE vs Control: U={u_pre:.0f}, P={p_pre:.4f}')

# Post-SSW vs control
u_post, p_post = stats.mannwhitneyu(post_danger, all_ctrl['danger_level'], alternative='less')
print(f'POST vs Control: U={u_post:.0f}, P={p_post:.4f}')

# Daily lag analysis: compute mean danger level by day relative to SSW onset
print(f'\n=== DAILY LAG ANALYSIS ===')
lag_data = []
for _, row in all_ssw.iterrows():
    date = pd.to_datetime(row['date'])
    ssw_dt = pd.to_datetime(row['ssw_date'])
    lag = (date - ssw_dt).days
    lag_data.append({'lag': lag, 'danger': row['danger_level']})

lag_df = pd.DataFrame(lag_data)
print('Day-by-day danger levels relative to SSW onset:')
for lag in range(-15, 16):
    sub = lag_df[lag_df['lag'] == lag]
    if len(sub) > 0:
        mean_d = sub['danger'].mean()
        n = len(sub)
        print(f'  Day {lag:+3d}: danger={mean_d:.2f} (n={n})')

# Control daily analysis
ctrl_lag = []
for _, row in all_ctrl.iterrows():
    date = pd.to_datetime(row['date'])
    ssw_dt = pd.to_datetime(row['ssw_date'])
    lag = (date - ssw_dt).days
    ctrl_lag.append({'lag': lag, 'danger': row['danger_level']})

ctrl_df = pd.DataFrame(ctrl_lag)
ctrl_mean = ctrl_df.groupby('lag')['danger'].mean()
print(f'\nControl mean danger by lag: {ctrl_mean.mean():.3f} (range: {ctrl_mean.min():.2f}-{ctrl_mean.max():.2f})')

# Save results
output = {
    'overall': {
        'ssw_mean_danger': float(all_ssw['danger_level'].mean()),
        'ctrl_mean_danger': float(all_ctrl['danger_level'].mean()),
        'pre_mean_danger': float(pre_danger.mean()),
        'post_mean_danger': float(post_danger.mean()),
        'pre_vs_ctrl_p': float(p_pre),
        'post_vs_ctrl_p': float(p_post),
        'pre_vs_post_p': float(p_val),
    },
    'daily_lag': {str(lag): float(lag_df[lag_df['lag'] == lag]['danger'].mean()) 
                  for lag in range(-15, 16) if len(lag_df[lag_df['lag'] == lag]) > 0}
}

with open('data/results/norway_phase_analysis.json', 'w') as f:
    json.dump(output, f, indent=2)
print('\nSaved to data/results/norway_phase_analysis.json')
