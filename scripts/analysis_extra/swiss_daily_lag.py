"""
Swiss daily lag analysis - compute avalanche rate by day relative to SSW onset.
Compare with Norwegian pattern to show cross-country consistency.
"""
import pandas as pd
import numpy as np
from scipy import stats
import json

panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
ssw = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw.index = ssw.index.tz_localize(None)
aval = 'dry_natural_size_1234'

# Only winter SSW events that overlap panel
panel_range = panel.dropna(subset=[aval])
ssw_in_range = ssw[(ssw.index >= panel_range.index.min()) & (ssw.index <= panel_range.index.max())]
print(f'SSW events in Swiss data: {len(ssw_in_range)}')

# Daily lag analysis
lag_data = []
for sd in ssw_in_range.index:
    for lag in range(-30, 31):
        date = sd + pd.Timedelta(days=lag)
        if date in panel_range.index:
            val = panel_range.loc[date, aval]
            if not np.isnan(val):
                lag_data.append({'lag': lag, 'aval': val, 'ssw_date': str(sd.date())})

lag_df = pd.DataFrame(lag_data)
print(f'Total lag observations: {len(lag_df)}')

# Compute control rate (seasonal mean for same DOYs across all non-SSW years)
# Simpler: use the full winter mean as baseline
winter = panel_range.copy()
winter['month'] = winter.index.month
winter = winter[(winter['month'] >= 11) | (winter['month'] <= 4)]
baseline_rate = winter[aval].mean()
print(f'Baseline winter rate: {baseline_rate:.3f}')

# Daily lag profile
print('\n=== SWISS DAILY LAG PROFILE ===')
print(f'{"Day":>5} {"Rate":>8} {"n":>4} {"RR":>8} {"vs_base_P":>10}')
daily_rates = {}
for lag in range(-30, 31):
    sub = lag_df[lag_df['lag'] == lag]
    if len(sub) > 0:
        rate = sub['aval'].mean()
        n = len(sub)
        rr = rate / baseline_rate if baseline_rate > 0 else np.nan
        # One-sample test vs baseline
        if n > 1 and sub['aval'].std() > 0:
            t, p = stats.ttest_1samp(sub['aval'], baseline_rate)
            p_str = f'{p:.4f}'
        else:
            p_str = 'n/a'
        daily_rates[lag] = {'rate': rate, 'n': n, 'rr': rr}
        marker = '*' if isinstance(p_str, str) and p_str != 'n/a' and float(p_str) < 0.05 else ''
        print(f'{lag:+5d} {rate:8.3f} {n:4d} {rr:8.3f} {p_str:>10} {marker}')

# Phase summary
pre_data = lag_df[(lag_df['lag'] >= -15) & (lag_df['lag'] < 0)]
post_data = lag_df[(lag_df['lag'] >= 0) & (lag_df['lag'] <= 15)]
late_data = lag_df[(lag_df['lag'] > 15) & (lag_df['lag'] <= 30)]

print(f'\n=== PHASE SUMMARY ===')
print(f'PRE [-15,0):  rate={pre_data["aval"].mean():.3f}, n={len(pre_data)}, '
      f'RR={pre_data["aval"].mean()/baseline_rate:.3f}')
print(f'POST [0,+15]: rate={post_data["aval"].mean():.3f}, n={len(post_data)}, '
      f'RR={post_data["aval"].mean()/baseline_rate:.3f}')
print(f'LATE (+15,+30]: rate={late_data["aval"].mean():.3f}, n={len(late_data)}, '
      f'RR={late_data["aval"].mean()/baseline_rate:.3f}')

# Pre vs Post comparison
u, p_pp = stats.mannwhitneyu(pre_data['aval'], post_data['aval'], alternative='two-sided')
print(f'PRE vs POST: MW P={p_pp:.4f}')

# Both vs baseline
u_pre, p_pre = stats.mannwhitneyu(pre_data['aval'], winter[aval].dropna(), alternative='less')
u_post, p_post = stats.mannwhitneyu(post_data['aval'], winter[aval].dropna(), alternative='less')
print(f'PRE vs baseline: MW P={p_pre:.4f}')
print(f'POST vs baseline: MW P={p_post:.4f}')

# Event-by-event pre vs post
print(f'\n=== EVENT-BY-EVENT PRE vs POST ===')
event_results = []
for sd in ssw_in_range.index:
    pre = lag_df[(lag_df['ssw_date'] == str(sd.date())) & (lag_df['lag'] >= -15) & (lag_df['lag'] < 0)]
    post = lag_df[(lag_df['ssw_date'] == str(sd.date())) & (lag_df['lag'] >= 0) & (lag_df['lag'] <= 15)]
    if len(pre) > 0 and len(post) > 0:
        pre_rate = pre['aval'].mean()
        post_rate = post['aval'].mean()
        event_results.append({
            'date': str(sd.date()),
            'pre_rate': pre_rate,
            'post_rate': post_rate,
            'pre_lower': pre_rate < post_rate,
        })
        marker = 'PRE<POST' if pre_rate < post_rate else 'POST<PRE'
        print(f'  {sd.date()}: PRE={pre_rate:.2f} POST={post_rate:.2f} [{marker}]')

n_pre_lower = sum(1 for e in event_results if e['pre_lower'])
n_total = len(event_results)
print(f'\n  PRE lower in {n_pre_lower}/{n_total} events')
p_sign = float(binom.sf(n_pre_lower - 1, n_total, 0.5)) if n_pre_lower > 0 else 1.0
from scipy.stats import binom
p_sign = float(binom.sf(n_pre_lower - 1, n_total, 0.5))
print(f'  Sign test: P={p_sign:.4f}')

# Save
output = {
    'baseline_rate': float(baseline_rate),
    'n_events': int(len(ssw_in_range)),
    'phases': {
        'pre': {'rate': float(pre_data['aval'].mean()), 'n': int(len(pre_data)),
                'rr': float(pre_data['aval'].mean()/baseline_rate), 'p_vs_baseline': float(p_pre)},
        'post': {'rate': float(post_data['aval'].mean()), 'n': int(len(post_data)),
                 'rr': float(post_data['aval'].mean()/baseline_rate), 'p_vs_baseline': float(p_post)},
        'late': {'rate': float(late_data['aval'].mean()), 'n': int(len(late_data)),
                 'rr': float(late_data['aval'].mean()/baseline_rate)},
    },
    'pre_vs_post_p': float(p_pp),
    'daily_rates': {str(k): v for k, v in daily_rates.items()},
}

with open('data/results/swiss_daily_lag.json', 'w') as f:
    json.dump(output, f, indent=2, default=str)
print('\nSaved to data/results/swiss_daily_lag.json')
