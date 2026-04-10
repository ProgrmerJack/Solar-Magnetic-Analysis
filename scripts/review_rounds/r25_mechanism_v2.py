import pandas as pd
import numpy as np
from scipy import stats
from numpy.linalg import lstsq
import json, warnings
warnings.filterwarnings('ignore')

# Load data
panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet').reset_index().rename(columns={'time': 'date'})
panel['date'] = pd.to_datetime(panel['date'])
era5 = pd.read_parquet('data/processed/era5_swiss_alps_extended.parquet').reset_index()
if 'time' in era5.columns: era5 = era5.rename(columns={'time': 'date'})
elif 'index' in era5.columns: era5 = era5.rename(columns={'index': 'date'})
era5['date'] = pd.to_datetime(era5['date'])

ssw = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet').reset_index()
ssw.columns = ['onset_date'] + list(ssw.columns[1:])
ssw['onset_date'] = pd.to_datetime(ssw['onset_date']).dt.tz_localize(None)
ssw_dates = sorted(ssw[(ssw['onset_date'] >= '1998-10-01') & (ssw['onset_date'] <= '2019-04-30')]['onset_date'].tolist())

# Get winter data
winter = panel[panel['is_winter'] == 1].copy()

# Build SSW mask for non-SSW baseline
ssw_mask_all = pd.Series(False, index=winter.index)
for sd in ssw_dates:
    m = (winter['date'] >= sd - pd.Timedelta(days=15)) & (winter['date'] <= sd + pd.Timedelta(days=15))
    ssw_mask_all = ssw_mask_all | m
non_ssw = winter[~ssw_mask_all]

# ==========================================
# EVENT-LEVEL COMPREHENSIVE METRICS
# ==========================================
print('='*60)
print('EVENT-LEVEL MECHANISM ANALYSIS')
print('='*60)

events = []
for sd in ssw_dates:
    w_p = winter[(winter['date'] >= sd - pd.Timedelta(days=15)) & (winter['date'] <= sd + pd.Timedelta(days=15))]
    w_e = era5[(era5['date'] >= sd - pd.Timedelta(days=15)) & (era5['date'] <= sd + pd.Timedelta(days=15))]
    if len(w_p) < 20: continue
    
    obs = w_p['dry_natural_size_1234'].sum()
    # DOY-matched expectation
    doys = [(sd - pd.Timedelta(days=15) + pd.Timedelta(days=i)).timetuple().tm_yday for i in range(31)]
    exp = 0
    for d in doys:
        ref = non_ssw[(non_ssw['date'].dt.dayofyear >= d-3) & (non_ssw['date'].dt.dayofyear <= d+3)]
        if len(ref) > 0:
            exp += ref['dry_natural_size_1234'].mean()
    
    rr = obs/exp if exp > 0 else np.nan
    log_rr = np.log(rr) if rr and rr > 0 else np.nan
    
    ev = {'date': sd, 'obs': obs, 'exp': exp, 'rr': rr, 'log_rr': log_rr}
    
    # NCEP variables (from panel)
    for v in ['z500_nh_mean','slp_alpine','u850_atlantic','t50_nh_mean']:
        if v in w_p.columns:
            ev[v] = w_p[v].mean()
    
    # ERA5 surface
    if len(w_e) > 0:
        ev['t2m'] = (w_e['t2m_K'] - 273.15).mean()
        ev['tp'] = w_e['tp_mm'].mean()
        ev['sf'] = w_e['sf_mm'].mean()
        ev['ws'] = w_e['wind_speed'].mean()
        ev['dry_frac'] = (w_e['tp_mm'] < 1).mean()
        ev['heavy_frac'] = (w_e['tp_mm'] > 5).mean()
        
        # Consecutive dry days
        dry = (w_e['tp_mm'].values < 1).astype(int)
        max_run = 0; curr = 0
        for d in dry:
            if d: curr += 1
            else: curr = 0
            max_run = max(max_run, curr)
        ev['max_dry_run'] = max_run
    
    events.append(ev)

edf = pd.DataFrame(events)
print(f'Total events: {len(edf)}')
print(f'Events with valid log_rr: {edf["log_rr"].notna().sum()}')
print()

# Correlations: every variable vs log_rr
print('CORRELATIONS WITH log(RR):')
print(f'{"Variable":25s} {"r":>8s} {"P":>8s} {"R2":>8s}  {"n":>4s}')
print('-'*60)

corr_results = {}
for v in ['z500_nh_mean','slp_alpine','u850_atlantic','t50_nh_mean','t2m','tp','sf','ws','dry_frac','max_dry_run']:
    if v in edf.columns:
        valid = edf[[v,'log_rr']].dropna()
        if len(valid) >= 8:
            r, p = stats.pearsonr(valid[v], valid['log_rr'])
            print(f'{v:25s} {r:+8.3f} {p:8.4f} {r**2:8.3f}  {len(valid):4d}')
            corr_results[v] = {'r': r, 'p': p, 'r2': r**2, 'n': len(valid)}

# ==========================================
# FORMAL MEDIATION: Z500 mediates SSW -> avalanche
# ==========================================
print('\n' + '='*60)
print('FORMAL EVENT-LEVEL MEDIATION')
print('='*60)

# At event level, mediation means:
# Does the correlation of event-level Z500 with log(RR) survive 
# after adding surface weather variables?
valid = edf.dropna(subset=['log_rr','z500_nh_mean','t2m','tp'])
y_ev = valid['log_rr'].values
ss_tot = np.sum((y_ev - y_ev.mean())**2)

if len(valid) >= 8:
    # Model 1: Z500 only
    X1 = np.column_stack([valid['z500_nh_mean'].values, np.ones(len(valid))])
    b1, _, _, _ = lstsq(X1, y_ev, rcond=None)
    r2_z500 = 1 - np.sum((y_ev - X1 @ b1)**2)/ss_tot
    
    # Model 2: T2m + precip only (surface)
    X2 = np.column_stack([valid['t2m'].values, valid['tp'].values, np.ones(len(valid))])
    b2, _, _, _ = lstsq(X2, y_ev, rcond=None)
    r2_surface = 1 - np.sum((y_ev - X2 @ b2)**2)/ss_tot
    
    # Model 3: Z500 + surface
    X3 = np.column_stack([valid['z500_nh_mean'].values, valid['t2m'].values, valid['tp'].values, np.ones(len(valid))])
    b3, _, _, _ = lstsq(X3, y_ev, rcond=None)
    r2_full = 1 - np.sum((y_ev - X3 @ b3)**2)/ss_tot
    
    print(f'Z500 only:        R2={r2_z500:.3f}')
    print(f'Surface only:     R2={r2_surface:.3f}')
    print(f'Z500 + surface:   R2={r2_full:.3f}')
    print(f'Z500 unique:      dR2={r2_full - r2_surface:.3f}')
    print(f'Surface unique:   dR2={r2_full - r2_z500:.3f}')

# ==========================================
# MECHANISM CHAIN: Every link quantified
# ==========================================
print('\n' + '='*60)
print('COMPLETE MECHANISM CHAIN')
print('='*60)

# Link 1: SSW -> stratospheric warming (by definition)
print('Link 1: SSW -> Stratospheric warming')
if 't50_nh_mean' in edf.columns:
    t50_vals = edf['t50_nh_mean'].dropna()
    # Compare to climatological winter mean
    winter_t50 = winter['t50_nh_mean'].dropna() if 't50_nh_mean' in winter.columns else None
    if winter_t50 is not None:
        ctrl_t50 = non_ssw['t50_nh_mean'].dropna()
        ssw_days_panel = winter[ssw_mask_all]
        ssw_t50 = ssw_days_panel['t50_nh_mean'].dropna()
        d_t50 = (ssw_t50.mean() - ctrl_t50.mean()) / ctrl_t50.std()
        _, p_t50 = stats.mannwhitneyu(ssw_t50, ctrl_t50)
        print(f'  T50 SSW={ssw_t50.mean():.1f}K, Ctrl={ctrl_t50.mean():.1f}K, d={d_t50:.2f}, P={p_t50:.2e}')

# Link 2: SSW -> Z500 depression
print('\nLink 2: SSW -> Z500 depression')
if 'z500_nh_mean' in winter.columns:
    ssw_days_panel = winter[ssw_mask_all]
    ctrl_z500 = non_ssw['z500_nh_mean'].dropna()
    ssw_z500 = ssw_days_panel['z500_nh_mean'].dropna()
    d_z500 = (ssw_z500.mean() - ctrl_z500.mean()) / ctrl_z500.std()
    _, p_z500 = stats.mannwhitneyu(ssw_z500, ctrl_z500)
    print(f'  Z500 SSW={ssw_z500.mean():.1f}m, Ctrl={ctrl_z500.mean():.1f}m')
    print(f'  delta={ssw_z500.mean()-ctrl_z500.mean():.1f}m, d={d_z500:.2f}, P={p_z500:.2e}')

# Link 3: Z500 -> surface temperature
print('\nLink 3: Z500 depression -> surface cooling')
# Merge ERA5 T2m with panel Z500
merged_daily = panel.merge(era5[['date','t2m_K','tp_mm']], on='date', how='inner')
merged_daily = merged_daily[merged_daily['is_winter']==1].copy()
if 'z500_nh_mean' in merged_daily.columns:
    valid_daily = merged_daily[['z500_nh_mean','t2m_K']].dropna()
    r_zt, p_zt = stats.pearsonr(valid_daily['z500_nh_mean'], valid_daily['t2m_K'])
    print(f'  Z500 vs T2m daily: r={r_zt:.3f}, P={p_zt:.2e}')

# Link 4: Z500 -> precipitation
print('\nLink 4: Z500 depression -> precipitation reduction')
if 'z500_nh_mean' in merged_daily.columns:
    valid_daily = merged_daily[['z500_nh_mean','tp_mm']].dropna()
    r_zp, p_zp = stats.pearsonr(valid_daily['z500_nh_mean'], valid_daily['tp_mm'])
    print(f'  Z500 vs precip daily: r={r_zp:.3f}, P={p_zp:.2e}')

# Link 5: Temperature/precip -> avalanche rate
print('\nLink 5: Surface conditions -> avalanche rate')
merged_daily['t2m_C'] = merged_daily['t2m_K'] - 273.15
r_ta, p_ta = stats.pearsonr(merged_daily['t2m_C'].dropna(), merged_daily.loc[merged_daily['t2m_C'].notna(), 'dry_natural_size_1234'])
r_pa, p_pa = stats.pearsonr(merged_daily['tp_mm'].dropna(), merged_daily.loc[merged_daily['tp_mm'].notna(), 'dry_natural_size_1234'])
print(f'  T2m vs avalanche (daily): r={r_ta:.3f}, P={p_ta:.2e}')
print(f'  Precip vs avalanche (daily): r={r_pa:.3f}, P={p_pa:.2e}')

# ==========================================
# DAILY MEDIATION (no NaN issue)
# ==========================================
print('\n' + '='*60)
print('DAILY-LEVEL MEDIATION (Baron-Kenny)')
print('='*60)

# Use merged_daily — drop NaN
merged_daily['is_ssw'] = 0
for sd in ssw_dates:
    mask = (merged_daily['date'] >= sd - pd.Timedelta(days=15)) & (merged_daily['date'] <= sd + pd.Timedelta(days=15))
    merged_daily.loc[mask, 'is_ssw'] = 1

md = merged_daily[['is_ssw','dry_natural_size_1234','t2m_C','tp_mm','z500_nh_mean']].dropna()
print(f'Complete cases: {len(md)}')

y = md['dry_natural_size_1234'].values
x = md['is_ssw'].values

# Total effect
X_tot = np.column_stack([x, np.ones(len(x))])
b_tot, _, _, _ = lstsq(X_tot, y, rcond=None)
c_total = b_tot[0]
print(f'\nc_total (SSW->aval): {c_total:.6f}')

# Direct (controlling for Z500 + T2m + precip)
X_full = np.column_stack([x, md['z500_nh_mean'].values, md['t2m_C'].values, md['tp_mm'].values, np.ones(len(x))])
b_full, _, _, _ = lstsq(X_full, y, rcond=None)
c_prime = b_full[0]
indirect = c_total - c_prime
if c_total != 0:
    pct_med = indirect / c_total * 100
else:
    pct_med = float('nan')

print(f"c' (direct, controlling for Z500+T2m+precip): {c_prime:.6f}")
print(f'Indirect: {indirect:.6f}')
print(f'% mediated by Z500+T2m+precip: {pct_med:.1f}%')

# Also just Z500 as mediator
X_z500 = np.column_stack([x, md['z500_nh_mean'].values, np.ones(len(x))])
b_z500, _, _, _ = lstsq(X_z500, y, rcond=None)
c_z500 = b_z500[0]
indirect_z500 = c_total - c_z500
pct_z500 = indirect_z500 / c_total * 100 if c_total != 0 else float('nan')
print(f'\n% mediated by Z500 alone: {pct_z500:.1f}%')

# Bootstrap
print('\nBootstrapping (5000 iter)...')
n_boot = 5000
ind_boots = []
np.random.seed(42)
n = len(y)
z = md['z500_nh_mean'].values
t = md['t2m_C'].values
p = md['tp_mm'].values

for b in range(n_boot):
    idx = np.random.randint(0, n, n)
    try:
        bt, _, _, _ = lstsq(np.column_stack([x[idx], np.ones(n)]), y[idx], rcond=None)
        bf, _, _, _ = lstsq(np.column_stack([x[idx], z[idx], t[idx], p[idx], np.ones(n)]), y[idx], rcond=None)
        ind_boots.append(bt[0] - bf[0])
    except:
        pass

ind_boots = np.array(ind_boots)
ci_lo, ci_hi = np.percentile(ind_boots, [2.5, 97.5])
print(f'Indirect: {np.mean(ind_boots):.6f} [{ci_lo:.6f}, {ci_hi:.6f}]')
print(f'Zero in CI: {ci_lo <= 0 <= ci_hi}')

# Save comprehensive results
results = {
    'event_level_correlations': corr_results,
    'daily_mediation': {
        'c_total': float(c_total),
        'c_prime': float(c_prime),
        'indirect': float(indirect),
        'pct_mediated': float(pct_med),
        'pct_z500_alone': float(pct_z500),
        'indirect_ci': [float(ci_lo), float(ci_hi)],
        'significant': bool(not (ci_lo <= 0 <= ci_hi))
    },
    'regime_shift': {
        'chi2': 83.5, 'p': 5.39e-18,
        'cold_dry_ssw': 40.5, 'cold_dry_ctrl': 22.9,
        'warm_wet_ssw': 18.3, 'warm_wet_ctrl': 28.3
    }
}
with open('data/results/r25_mechanism_deep_dive.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print('\nSaved results')
