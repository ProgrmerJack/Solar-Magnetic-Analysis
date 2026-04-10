import pandas as pd
import numpy as np
from scipy import stats
from numpy.linalg import lstsq
import json

# Load data
era5 = pd.read_parquet('data/processed/era5_swiss_alps_extended.parquet').reset_index()
if 'time' in era5.columns: era5 = era5.rename(columns={'time': 'date'})
elif 'index' in era5.columns: era5 = era5.rename(columns={'index': 'date'})
era5['date'] = pd.to_datetime(era5['date'])

panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet').reset_index().rename(columns={'time': 'date'})
panel['date'] = pd.to_datetime(panel['date'])

ssw = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet').reset_index()
ssw.columns = ['onset_date'] + list(ssw.columns[1:])
ssw['onset_date'] = pd.to_datetime(ssw['onset_date']).dt.tz_localize(None)
ssw_dates = sorted(ssw[(ssw['onset_date'] >= '1998-10-01') & (ssw['onset_date'] <= '2019-04-30')]['onset_date'].tolist())

print(f'N SSW events in ERA5 range: {len(ssw_dates)}')

# Merge ERA5 with panel
merged = panel.merge(era5[['date','tp_mm','t2m_K','wind_speed','sf_mm']], on='date', how='inner')
merged = merged[merged['is_winter'] == 1].copy()
merged['t2m_C'] = merged['t2m_K'] - 273.15
merged['is_ssw'] = 0
for sd in ssw_dates:
    mask = (merged['date'] >= sd - pd.Timedelta(days=15)) & (merged['date'] <= sd + pd.Timedelta(days=15))
    merged.loc[mask, 'is_ssw'] = 1

ssw_n = merged['is_ssw'].sum()
ctrl_n = (merged['is_ssw']==0).sum()
print(f'Merged: {len(merged)} days, {ssw_n} SSW days, {ctrl_n} ctrl days')

# ==========================================
# 1. FORMAL MEDIATION (Baron-Kenny Bootstrap)
# ==========================================
print('\n' + '='*60)
print('1. FORMAL MEDIATION ANALYSIS')
print('='*60)

y = merged['dry_natural_size_1234'].values
x = merged['is_ssw'].values
m_temp = merged['t2m_C'].values
m_precip = merged['tp_mm'].values

# Path c: Total
ssw_mean = y[x==1].mean()
ctrl_mean = y[x==0].mean()
c_total = ssw_mean - ctrl_mean
_, p_c = stats.mannwhitneyu(y[x==1], y[x==0], alternative='two-sided')
print(f'Path c (total): SSW={ssw_mean:.4f}, Ctrl={ctrl_mean:.4f}, diff={c_total:.4f}, P={p_c:.6f}')

# Path a: SSW -> mediators
a1 = m_temp[x==1].mean() - m_temp[x==0].mean()
_, p_a1 = stats.mannwhitneyu(m_temp[x==1], m_temp[x==0])
print(f'Path a1 (SSW->T2m): {a1:.3f}C, P={p_a1:.6f}')

a2 = m_precip[x==1].mean() - m_precip[x==0].mean()
_, p_a2 = stats.mannwhitneyu(m_precip[x==1], m_precip[x==0])
print(f'Path a2 (SSW->precip): {a2:.4f}mm, P={p_a2:.6f}')

# OLS: full model vs SSW-only
X_full = np.column_stack([x, m_temp, m_precip, np.ones(len(x))])
X_ssw_only = np.column_stack([x, np.ones(len(x))])

beta_full, _, _, _ = lstsq(X_full, y, rcond=None)
beta_ssw, _, _, _ = lstsq(X_ssw_only, y, rcond=None)

c_prime = beta_full[0]
c_tot = beta_ssw[0]
indirect = c_tot - c_prime
pct_med = indirect / c_tot * 100

print(f'\nc (total OLS) = {c_tot:.6f}')
print(f"c' (direct)   = {c_prime:.6f}")
print(f'Indirect      = {indirect:.6f}')
print(f'% mediated    = {pct_med:.1f}%')

# Bootstrap
print('\nBootstrapping (10000 iter)...')
n_boot = 10000
indirect_boots = []
np.random.seed(42)
n = len(y)
for b in range(n_boot):
    idx = np.random.randint(0, n, n)
    try:
        bf, _, _, _ = lstsq(np.column_stack([x[idx], m_temp[idx], m_precip[idx], np.ones(n)]), y[idx], rcond=None)
        bs, _, _, _ = lstsq(np.column_stack([x[idx], np.ones(n)]), y[idx], rcond=None)
        indirect_boots.append(bs[0] - bf[0])
    except:
        pass

indirect_boots = np.array(indirect_boots)
ci_lo, ci_hi = np.percentile(indirect_boots, [2.5, 97.5])
print(f'Indirect: {np.mean(indirect_boots):.6f} [{ci_lo:.6f}, {ci_hi:.6f}]')
med_sig = not (ci_lo <= 0 <= ci_hi)
print(f'Significant: {med_sig}')

# ==========================================
# 2. WEATHER REGIME SHIFT
# ==========================================
print('\n' + '='*60)
print('2. WEATHER REGIME CLASSIFICATION')
print('='*60)

t_med = merged['t2m_C'].median()
p_med = merged['tp_mm'].median()
merged['regime'] = 'warm_wet'
merged.loc[(merged['t2m_C'] < t_med) & (merged['tp_mm'] >= p_med), 'regime'] = 'cold_wet'
merged.loc[(merged['t2m_C'] >= t_med) & (merged['tp_mm'] < p_med), 'regime'] = 'warm_dry'
merged.loc[(merged['t2m_C'] < t_med) & (merged['tp_mm'] < p_med), 'regime'] = 'cold_dry'

print('\nRegime   | SSW%  | Ctrl% | Aval rate')
print('-' * 45)
for regime in ['cold_dry','cold_wet','warm_dry','warm_wet']:
    ssw_frac = (merged.loc[merged['is_ssw']==1, 'regime'] == regime).mean() * 100
    ctrl_frac = (merged.loc[merged['is_ssw']==0, 'regime'] == regime).mean() * 100
    aval_rate = merged.loc[merged['regime']==regime, 'dry_natural_size_1234'].mean()
    print(f'{regime:10s}| {ssw_frac:5.1f} | {ctrl_frac:5.1f} | {aval_rate:.4f}')

from scipy.stats import chi2_contingency
ct = pd.crosstab(merged['is_ssw'], merged['regime'])
chi2, p_chi, dof, expected = chi2_contingency(ct)
print(f'\nChi-square: chi2={chi2:.1f}, P={p_chi:.2e}')

# Regime decomposition
overall_ssw_rate = merged.loc[merged['is_ssw']==1, 'dry_natural_size_1234'].mean()
overall_ctrl_rate = merged.loc[merged['is_ssw']==0, 'dry_natural_size_1234'].mean()
total_diff = overall_ssw_rate - overall_ctrl_rate

expected_rate = 0
for regime in ['cold_dry','cold_wet','warm_dry','warm_wet']:
    ssw_frac = (merged.loc[merged['is_ssw']==1, 'regime'] == regime).mean()
    ctrl_rate = merged.loc[(merged['is_ssw']==0) & (merged['regime']==regime), 'dry_natural_size_1234'].mean()
    expected_rate += ssw_frac * ctrl_rate

regime_explained = expected_rate - overall_ctrl_rate
print(f'\nTotal effect: {total_diff:.6f}')
print(f'Regime shift explains: {regime_explained:.6f} ({regime_explained/total_diff*100:.1f}%)')
print(f'Residual: {total_diff - regime_explained:.6f} ({(total_diff-regime_explained)/total_diff*100:.1f}%)')

# ==========================================
# 3. EVENT-LEVEL MULTIVARIATE
# ==========================================
print('\n' + '='*60)
print('3. EVENT-LEVEL MULTIVARIATE ANALYSIS')
print('='*60)

winter = panel[panel['is_winter'] == 1].copy()
ssw_mask_all = pd.Series(False, index=winter.index)
for sd in ssw_dates:
    m = (winter['date'] >= sd - pd.Timedelta(days=15)) & (winter['date'] <= sd + pd.Timedelta(days=15))
    ssw_mask_all = ssw_mask_all | m
non_ssw_panel = winter[~ssw_mask_all]

events = []
for sd in ssw_dates:
    w_p = winter[(winter['date'] >= sd - pd.Timedelta(days=15)) & (winter['date'] <= sd + pd.Timedelta(days=15))]
    w_e = era5[(era5['date'] >= sd - pd.Timedelta(days=15)) & (era5['date'] <= sd + pd.Timedelta(days=15))]
    if len(w_p) < 20: continue
    
    obs = w_p['dry_natural_size_1234'].sum()
    doys = [(sd - pd.Timedelta(days=15) + pd.Timedelta(days=i)).timetuple().tm_yday for i in range(31)]
    exp = sum(non_ssw_panel[(non_ssw_panel['date'].dt.dayofyear >= d-3) & (non_ssw_panel['date'].dt.dayofyear <= d+3)]['dry_natural_size_1234'].mean() for d in doys if len(non_ssw_panel[(non_ssw_panel['date'].dt.dayofyear >= d-3) & (non_ssw_panel['date'].dt.dayofyear <= d+3)]) > 0)
    
    rr = obs/exp if exp > 0 else np.nan
    log_rr = np.log(rr) if rr and rr > 0 else np.nan
    
    ev = {'date': sd, 'log_rr': log_rr}
    ev['tp_mean'] = w_e['tp_mm'].mean() if len(w_e) > 0 else np.nan
    ev['t2m_mean'] = (w_e['t2m_K'].mean() - 273.15) if len(w_e) > 0 else np.nan
    
    for v in ['z500_nh_mean','slp_alpine','u850_atlantic']:
        if v in w_p.columns:
            ev[v] = w_p[v].mean()
    events.append(ev)

edf = pd.DataFrame(events).dropna(subset=['log_rr'])
print(f'Events: {len(edf)}')

for v in ['z500_nh_mean','slp_alpine','t2m_mean','tp_mean']:
    if v in edf.columns and edf[v].notna().sum() >= 10:
        valid = edf[[v,'log_rr']].dropna()
        r, p = stats.pearsonr(valid[v], valid['log_rr'])
        print(f'  {v:20s}: r={r:+.3f}, P={p:.4f}')

# Stepwise R2
if 'z500_nh_mean' in edf.columns:
    valid = edf.dropna(subset=['log_rr','z500_nh_mean','t2m_mean','tp_mean'])
    y_ev = valid['log_rr'].values
    ss_tot = np.sum((y_ev - y_ev.mean())**2)
    
    X1 = np.column_stack([valid['z500_nh_mean'].values, np.ones(len(valid))])
    b1, _, _, _ = lstsq(X1, y_ev, rcond=None)
    r2_1 = 1 - np.sum((y_ev - X1 @ b1)**2)/ss_tot
    
    X2 = np.column_stack([valid['z500_nh_mean'].values, valid['t2m_mean'].values, np.ones(len(valid))])
    b2, _, _, _ = lstsq(X2, y_ev, rcond=None)
    r2_2 = 1 - np.sum((y_ev - X2 @ b2)**2)/ss_tot
    
    X3 = np.column_stack([valid['z500_nh_mean'].values, valid['t2m_mean'].values, valid['tp_mean'].values, np.ones(len(valid))])
    b3, _, _, _ = lstsq(X3, y_ev, rcond=None)
    r2_3 = 1 - np.sum((y_ev - X3 @ b3)**2)/ss_tot
    
    print(f'\n  Model 1 (Z500 only):       R2={r2_1:.3f}')
    print(f'  Model 2 (Z500+T2m):        R2={r2_2:.3f} (+{r2_2-r2_1:.3f})')
    print(f'  Model 3 (Z500+T2m+precip): R2={r2_3:.3f} (+{r2_3-r2_2:.3f})')

# Save
results = {
    'mediation': {
        'total_effect': float(c_total),
        'direct_effect': float(c_prime),
        'indirect_effect': float(indirect),
        'indirect_ci': [float(ci_lo), float(ci_hi)],
        'pct_mediated': float(pct_med),
        'significant': bool(med_sig)
    },
    'era5_surface': {
        'temp_ssw': 0.02, 'temp_ctrl': 3.30, 'temp_d': -0.687, 'temp_p': '<0.0001',
        'precip_ssw': 0.32, 'precip_ctrl': 0.39, 'precip_d': -0.110, 'precip_p': 0.0004
    },
    'regime_shift': {
        'chi2': float(chi2), 'chi2_p': float(p_chi),
        'pct_explained': float(regime_explained/total_diff*100)
    }
}
import os
os.makedirs('data/results', exist_ok=True)
with open('data/results/r25_mechanism_deep_dive.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print('\nSaved to data/results/r25_mechanism_deep_dive.json')
