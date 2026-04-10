"""R28 Formal Causal Mediation + Effective Sample Size Analysis"""
import pandas as pd
import numpy as np
import json
from scipy import stats

# Load panel data
panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
panel = panel.reset_index()
if 'time' in panel.columns:
    panel = panel.rename(columns={'time': 'date'})
panel['date'] = pd.to_datetime(panel['date']).dt.tz_localize(None)

# Load ERA5 for regime classification
era5 = pd.read_parquet('data/processed/era5_swiss_alps_extended.parquet')
era5 = era5.reset_index()
era5['date'] = pd.to_datetime(era5['date']).dt.tz_localize(None)

# Merge ERA5 into panel
panel = panel.merge(era5[['date', 't2m_K', 'tp_mm']], on='date', how='left')

# SSW catalog
ssw = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw = ssw.reset_index()
if 'onset_date' in ssw.columns:
    ssw['onset_date'] = pd.to_datetime(ssw['onset_date']).dt.tz_localize(None)

winter = panel[panel['is_winter'] == 1].copy()
winter['doy'] = winter['date'].dt.dayofyear

# Define SSW windows
ssw_windows = []
for _, row in ssw.iterrows():
    onset = row['onset_date']
    start = onset - pd.Timedelta(days=15)
    end = onset + pd.Timedelta(days=15)
    ssw_windows.append((start, end, onset))

winter['in_ssw'] = False
winter['ssw_event'] = None
for start, end, onset in ssw_windows:
    mask = (winter['date'] >= start) & (winter['date'] <= end)
    winter.loc[mask, 'in_ssw'] = True
    winter.loc[mask, 'ssw_event'] = onset.strftime('%Y-%m-%d')

# Classify regimes using ERA5 t2m and tp
winter_with_era5 = winter.dropna(subset=['t2m_K', 'tp_mm']).copy()
t2m_med = winter_with_era5['t2m_K'].median()
tp_med = winter_with_era5['tp_mm'].median()
conditions = [
    (winter_with_era5['t2m_K'] < t2m_med) & (winter_with_era5['tp_mm'] < tp_med),
    (winter_with_era5['t2m_K'] < t2m_med) & (winter_with_era5['tp_mm'] >= tp_med),
    (winter_with_era5['t2m_K'] >= t2m_med) & (winter_with_era5['tp_mm'] < tp_med),
    (winter_with_era5['t2m_K'] >= t2m_med) & (winter_with_era5['tp_mm'] >= tp_med),
]
labels = ['cold_dry', 'cold_wet', 'warm_dry', 'warm_wet']
winter_with_era5['regime'] = np.select(conditions, labels, default='unknown')
winter = winter_with_era5

# Avalanche column
aval_col = 'dry_natural_size_1234'
print(f'Using avalanche column: {aval_col}')
print(f'Winter days with ERA5: {len(winter)}')

# Event-level aggregation
events = winter[winter['ssw_event'].notna()].groupby('ssw_event').agg(
    cold_dry_frac=('regime', lambda x: (x == 'cold_dry').mean()),
    aval_rate=(aval_col, 'mean'),
    n_days=('date', 'count'),
    doy_mean=('doy', 'mean')
).reset_index()

# Baron-Kenny paths
ssw_cd_mean = events['cold_dry_frac'].mean()
ctl_cd_mean = (winter[~winter['in_ssw']]['regime'] == 'cold_dry').mean()

r_mediator, p_mediator = stats.pearsonr(events['cold_dry_frac'], events['aval_rate'])

ssw_aval = winter[winter['in_ssw']][aval_col].mean()
ctl_aval = winter[~winter['in_ssw']][aval_col].mean()

a = ssw_cd_mean - ctl_cd_mean
slope_b = np.polyfit(events['cold_dry_frac'], events['aval_rate'], 1)[0]
indirect = a * slope_b
total = ssw_aval - ctl_aval
prop_mediated = indirect / total if total != 0 else None

print('=== BARON-KENNY MEDIATION (Event-Level) ===')
print(f'Events analyzed: {len(events)}')
print(f'Path a (SSW -> cold_dry_frac): SSW={ssw_cd_mean:.3f}, Control={ctl_cd_mean:.3f}, effect=+{a:.3f}')
print(f'Path b (cold_dry_frac -> aval_rate): r={r_mediator:.3f}, P={p_mediator:.4f}')
print(f'Path c (Total): SSW={ssw_aval:.3f}, Control={ctl_aval:.3f}, diff={total:.3f}')
print(f'Indirect (a*b): {indirect:.4f}')
print(f'Proportion mediated: {prop_mediated:.1%}' if prop_mediated else 'N/A')

# Bootstrap mediation CI
n_boot = 5000
boot_prop = []
for _ in range(n_boot):
    idx = np.random.choice(len(events), len(events), replace=True)
    boot_events = events.iloc[idx]
    boot_cd = boot_events['cold_dry_frac'].mean()
    boot_a = boot_cd - ctl_cd_mean
    if boot_events['cold_dry_frac'].std() > 0:
        boot_slope = np.polyfit(boot_events['cold_dry_frac'], boot_events['aval_rate'], 1)[0]
    else:
        boot_slope = slope_b
    boot_indirect = boot_a * boot_slope
    boot_total = boot_events['aval_rate'].mean() - ctl_aval
    if boot_total != 0:
        boot_prop.append(boot_indirect / boot_total)

boot_prop = np.array(boot_prop)
med_ci_lo, med_ci_hi = np.percentile(boot_prop, [2.5, 97.5])
print(f'Bootstrap mediation CI: [{med_ci_lo:.1%}, {med_ci_hi:.1%}]')

# Effective sample size
icc_norway = 0.7
k_norway = 5
deff_norway = 1 + (k_norway - 1) * icc_norway
n_eff_norway = (4 * k_norway) / deff_norway
n_eff_swiss = 16
n_eff_utah = 4
n_eff_total = n_eff_swiss + n_eff_norway + n_eff_utah

print(f'\n=== EFFECTIVE SAMPLE SIZE ===')
print(f'Switzerland: {n_eff_swiss} independent events')
print(f'Norway: 4 events x 5 regions, ICC={icc_norway}, DEFF={deff_norway:.1f}, n_eff={n_eff_norway:.1f}')
print(f'Utah: {n_eff_utah} events')
print(f'Total effective n: {n_eff_total:.1f} (vs naive 36)')

n_eff_pairs = int(round(n_eff_total))
decrease_rate = 32 / 36
n_decrease_eff = int(round(decrease_rate * n_eff_pairs))
binom_p = stats.binomtest(n_decrease_eff, n_eff_pairs, 0.5, alternative='greater').pvalue
print(f'Adjusted sign test: {n_decrease_eff}/{n_eff_pairs} decrease, P={binom_p:.6f}')

results = {
    'mediation': {
        'n_events': len(events),
        'path_a_effect': float(a),
        'path_b_slope': float(slope_b),
        'path_b_r': float(r_mediator),
        'path_b_p': float(p_mediator),
        'path_c_total': float(total),
        'indirect_effect': float(indirect),
        'proportion_mediated': float(prop_mediated) if prop_mediated else None,
        'bootstrap_ci_lo': float(med_ci_lo),
        'bootstrap_ci_hi': float(med_ci_hi),
    },
    'effective_sample_size': {
        'swiss_n_eff': n_eff_swiss,
        'norway_n_eff': float(n_eff_norway),
        'utah_n_eff': n_eff_utah,
        'total_n_eff': float(n_eff_total),
        'naive_n': 36,
        'icc_norway': icc_norway,
        'deff_norway': float(deff_norway),
        'adjusted_sign_test_p': float(binom_p)
    }
}

with open('data/results/r28_mediation_analysis.json', 'w') as f:
    json.dump(results, f, indent=2)
print('\nResults saved.')
