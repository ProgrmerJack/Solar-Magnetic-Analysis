"""
R52 Pre-Event Snowpack Comparison: SSW vs Non-SSW Years

Tests whether antecedent snowpack conditions (snow depth, SWE, PWL prevalence,
stability indices) differ between SSW and non-SSW years in the 30-to-16 day
pre-onset window. If they do NOT differ, this rules out pre-event snowpack
as a confounder — the trigger-suppression effect operates on the meteorological
forcing window alone, strengthening the mechanistic argument.

Output: data/results/r52_pre_event_snowpack.json
"""
import pandas as pd
import numpy as np
from scipy import stats
import json, os

# Load data
snowpack = pd.read_csv('data/cryosphere/envidat/weather_snowpack_danger.csv',
                       parse_dates=['datum'])
ssw_cat = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')

# SSW onset dates within SNOWPACK period (2004-2013 based on manuscript)
ssw_dates = ssw_cat.index.tz_localize(None)
snowpack_range = (snowpack['datum'].min(), snowpack['datum'].max())
print(f"SNOWPACK date range: {snowpack_range[0]} to {snowpack_range[1]}")

ssw_in_range = ssw_dates[(ssw_dates >= snowpack_range[0]) & 
                          (ssw_dates <= snowpack_range[1])]
print(f"SSW events in SNOWPACK range: {len(ssw_in_range)}")
for d in ssw_in_range:
    print(f"  {d.strftime('%Y-%m-%d')}")

# Variables to compare
vars_to_test = {
    'HS_mod': 'Snow depth (cm)',
    'SWE': 'Snow water equivalent (kg/m²)',
    'HN24': 'New snow 24h (cm)',
    'pwl_100': 'PWL prevalence (top 100cm)',
    'ssi_pwl': 'Stability index (PWL)',
    'sk38_pwl': 'Skier stability index',
    'Sn': 'Natural stability',
    'Sd': 'Deformation stability',
    'hoar_size': 'Surface hoar size',
}

# Pre-event window: -30 to -16 days before SSW onset
# This is BEFORE the SSW window, so it captures pre-existing snowpack
pre_window = (-30, -16)

results = {}

# For each SSW event, get station-day data in pre-event window
ssw_pre_records = []
for onset in ssw_in_range:
    start = onset + pd.Timedelta(days=pre_window[0])
    end = onset + pd.Timedelta(days=pre_window[1])
    mask = (snowpack['datum'] >= start) & (snowpack['datum'] <= end)
    subset = snowpack.loc[mask].copy()
    subset['ssw_event'] = onset.strftime('%Y-%m-%d')
    subset['is_ssw'] = True
    ssw_pre_records.append(subset)

ssw_pre = pd.concat(ssw_pre_records, ignore_index=True) if ssw_pre_records else pd.DataFrame()
print(f"\nSSW pre-event station-days: {len(ssw_pre)}")

# For control: same DOY windows from non-SSW winters
# Get DOYs for each SSW pre-event window
ctrl_records = []
all_winters = snowpack['datum'].dt.year.unique()
ssw_years = set(ssw_in_range.year)

for onset in ssw_in_range:
    doy_start = (onset + pd.Timedelta(days=pre_window[0])).dayofyear
    doy_end = (onset + pd.Timedelta(days=pre_window[1])).dayofyear
    
    for year in all_winters:
        if year in ssw_years:
            continue
        # Match by DOY
        mask = (snowpack['datum'].dt.year == year) & \
               (snowpack['datum'].dt.dayofyear >= doy_start) & \
               (snowpack['datum'].dt.dayofyear <= doy_end)
        subset = snowpack.loc[mask].copy()
        if len(subset) > 0:
            subset['ssw_event'] = 'control'
            subset['is_ssw'] = False
            ctrl_records.append(subset)

ctrl_pre = pd.concat(ctrl_records, ignore_index=True) if ctrl_records else pd.DataFrame()
print(f"Control pre-event station-days: {len(ctrl_pre)}")

# Compare each variable
print(f"\n{'Variable':<30} {'SSW mean':>10} {'Ctrl mean':>10} {'Diff%':>8} {'P-value':>10} {'Cohen d':>8}")
print("-" * 80)

for var, label in vars_to_test.items():
    if var not in snowpack.columns:
        print(f"{label:<30} {'N/A':>10}")
        continue
    
    ssw_vals = ssw_pre[var].dropna()
    ctrl_vals = ctrl_pre[var].dropna()
    
    if len(ssw_vals) < 10 or len(ctrl_vals) < 10:
        print(f"{label:<30} insufficient data")
        continue
    
    ssw_mean = ssw_vals.mean()
    ctrl_mean = ctrl_vals.mean()
    diff_pct = ((ssw_mean - ctrl_mean) / ctrl_mean * 100) if ctrl_mean != 0 else 0
    
    # Mann-Whitney U test (robust to non-normality)
    u_stat, p_val = stats.mannwhitneyu(ssw_vals, ctrl_vals, alternative='two-sided')
    
    # Cohen's d
    pooled_std = np.sqrt((ssw_vals.std()**2 + ctrl_vals.std()**2) / 2)
    d = (ssw_mean - ctrl_mean) / pooled_std if pooled_std > 0 else 0
    
    print(f"{label:<30} {ssw_mean:10.2f} {ctrl_mean:10.2f} {diff_pct:7.1f}% {p_val:10.4f} {d:8.3f}")
    
    results[var] = {
        'label': label,
        'ssw_mean': round(float(ssw_mean), 3),
        'ctrl_mean': round(float(ctrl_mean), 3),
        'diff_pct': round(float(diff_pct), 1),
        'mann_whitney_p': round(float(p_val), 4),
        'cohen_d': round(float(d), 3),
        'n_ssw': int(len(ssw_vals)),
        'n_ctrl': int(len(ctrl_vals)),
    }

# Event-level aggregation: mean pre-event snowpack per SSW event vs per control winter
print("\n\n=== EVENT-LEVEL COMPARISON ===")
print("Mean pre-event snowpack depth per SSW event vs per control winter:")

ssw_event_means = ssw_pre.groupby('ssw_event')['HS_mod'].mean()
# For controls, group by year
ctrl_pre_copy = ctrl_pre.copy()
ctrl_pre_copy['year'] = ctrl_pre_copy['datum'].dt.year
ctrl_year_means = ctrl_pre_copy.groupby('year')['HS_mod'].mean()

print(f"\nSSW event-level snow depths (n={len(ssw_event_means)}):")
for ev, val in ssw_event_means.items():
    print(f"  {ev}: {val:.1f} cm")
print(f"  Mean: {ssw_event_means.mean():.1f} cm")

print(f"\nControl winter-level snow depths (n={len(ctrl_year_means)}):")
print(f"  Mean: {ctrl_year_means.mean():.1f} cm")
print(f"  Range: {ctrl_year_means.min():.1f} - {ctrl_year_means.max():.1f} cm")

# Event-level test
t_stat, p_event = stats.mannwhitneyu(ssw_event_means.values, ctrl_year_means.values, 
                                      alternative='two-sided')
d_event = (ssw_event_means.mean() - ctrl_year_means.mean()) / \
          np.sqrt((ssw_event_means.std()**2 + ctrl_year_means.std()**2) / 2)
print(f"\nEvent-level Mann-Whitney P = {p_event:.4f}, Cohen's d = {d_event:.3f}")

results['event_level_HS'] = {
    'ssw_mean': round(float(ssw_event_means.mean()), 1),
    'ctrl_mean': round(float(ctrl_year_means.mean()), 1),
    'n_ssw_events': int(len(ssw_event_means)),
    'n_ctrl_winters': int(len(ctrl_year_means)),
    'mann_whitney_p': round(float(p_event), 4),
    'cohen_d': round(float(d_event), 3),
}

# SWE event-level
if 'SWE' in snowpack.columns:
    ssw_swe = ssw_pre.groupby('ssw_event')['SWE'].mean()
    ctrl_swe = ctrl_pre_copy.groupby('year')['SWE'].mean()
    t_swe, p_swe = stats.mannwhitneyu(ssw_swe.values, ctrl_swe.values, alternative='two-sided')
    d_swe = (ssw_swe.mean() - ctrl_swe.mean()) / \
            np.sqrt((ssw_swe.std()**2 + ctrl_swe.std()**2) / 2)
    print(f"\nSWE event-level: SSW={ssw_swe.mean():.1f}, Ctrl={ctrl_swe.mean():.1f}")
    print(f"  Mann-Whitney P = {p_swe:.4f}, Cohen's d = {d_swe:.3f}")
    
    results['event_level_SWE'] = {
        'ssw_mean': round(float(ssw_swe.mean()), 1),
        'ctrl_mean': round(float(ctrl_swe.mean()), 1),
        'mann_whitney_p': round(float(p_swe), 4),
        'cohen_d': round(float(d_swe), 3),
    }

# PWL event-level
if 'pwl_100' in snowpack.columns:
    ssw_pwl = ssw_pre.groupby('ssw_event')['pwl_100'].mean()
    ctrl_pwl = ctrl_pre_copy.groupby('year')['pwl_100'].mean()
    t_pwl, p_pwl = stats.mannwhitneyu(ssw_pwl.values, ctrl_pwl.values, alternative='two-sided')
    d_pwl = (ssw_pwl.mean() - ctrl_pwl.mean()) / \
            np.sqrt((ssw_pwl.std()**2 + ctrl_pwl.std()**2) / 2)
    print(f"\nPWL event-level: SSW={ssw_pwl.mean():.3f}, Ctrl={ctrl_pwl.mean():.3f}")
    print(f"  Mann-Whitney P = {p_pwl:.4f}, Cohen's d = {d_pwl:.3f}")
    
    results['event_level_PWL'] = {
        'ssw_mean': round(float(ssw_pwl.mean()), 3),
        'ctrl_mean': round(float(ctrl_pwl.mean()), 3),
        'mann_whitney_p': round(float(p_pwl), 4),
        'cohen_d': round(float(d_pwl), 3),
    }

# Summary
n_nonsig = sum(1 for v in results.values() if isinstance(v, dict) and 'mann_whitney_p' in v and v['mann_whitney_p'] > 0.05)
n_total = sum(1 for v in results.values() if isinstance(v, dict) and 'mann_whitney_p' in v)
print(f"\n=== SUMMARY ===")
print(f"{n_nonsig}/{n_total} comparisons non-significant (P > 0.05)")
print("This rules out pre-event snowpack as a confounder if most are non-significant.")

# Save results
os.makedirs('data/results', exist_ok=True)
with open('data/results/r52_pre_event_snowpack.json', 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nResults saved to data/results/r52_pre_event_snowpack.json")
