import pandas as pd
import numpy as np
from scipy import stats
from numpy.linalg import lstsq
import json, warnings
warnings.filterwarnings('ignore')

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

winter = panel[panel['is_winter'] == 1].copy()
ssw_mask_all = pd.Series(False, index=winter.index)
for sd in ssw_dates:
    m = (winter['date'] >= sd - pd.Timedelta(days=15)) & (winter['date'] <= sd + pd.Timedelta(days=15))
    ssw_mask_all = ssw_mask_all | m
non_ssw = winter[~ssw_mask_all]

# ==========================================
# 1. EVENT-LEVEL: all variables vs log(RR)
# ==========================================
print('='*60)
print('1. EVENT-LEVEL CORRELATIONS')
print('='*60)

events = []
for sd in ssw_dates:
    w_p = winter[(winter['date'] >= sd - pd.Timedelta(days=15)) & (winter['date'] <= sd + pd.Timedelta(days=15))]
    w_e = era5[(era5['date'] >= sd - pd.Timedelta(days=15)) & (era5['date'] <= sd + pd.Timedelta(days=15))]
    if len(w_p) < 20: continue
    
    obs = w_p['dry_natural_size_1234'].sum()
    doys = [(sd - pd.Timedelta(days=15) + pd.Timedelta(days=i)).timetuple().tm_yday for i in range(31)]
    exp = 0
    for d in doys:
        ref = non_ssw[(non_ssw['date'].dt.dayofyear >= d-3) & (non_ssw['date'].dt.dayofyear <= d+3)]
        if len(ref) > 0: exp += ref['dry_natural_size_1234'].mean()
    
    rr = obs/exp if exp > 0 else np.nan
    log_rr = np.log(rr) if rr and rr > 0 else np.nan
    
    ev = {'date': sd, 'rr': rr, 'log_rr': log_rr}
    # NCEP (correct column names)
    for v in ['ncep_z500_nh','ncep_slp_nh','ncep_u_50hpa','ncep_t_50hpa','ncep_z_50hpa']:
        if v in w_p.columns: ev[v] = w_p[v].mean()
    # ERA5
    if len(w_e) > 0:
        ev['t2m'] = (w_e['t2m_K'] - 273.15).mean()
        ev['tp'] = w_e['tp_mm'].mean()
        ev['sf'] = w_e['sf_mm'].mean()
        ev['ws'] = w_e['wind_speed'].mean()
        ev['dry_frac'] = (w_e['tp_mm'] < 1).mean()
        ev['heavy_frac'] = (w_e['tp_mm'] > 5).mean()
    events.append(ev)

edf = pd.DataFrame(events)
print(f'Events: {len(edf)}')
print(f'Columns: {list(edf.columns)}')

print(f'\n{"Variable":25s} {"r":>8s} {"P":>8s} {"R2":>8s}')
print('-'*55)
corr_results = {}
for v in ['ncep_z500_nh','ncep_slp_nh','ncep_t_50hpa','ncep_u_50hpa','t2m','tp','sf','ws','dry_frac']:
    if v in edf.columns:
        valid = edf[[v,'log_rr']].dropna()
        if len(valid) >= 8:
            r, p = stats.pearsonr(valid[v], valid['log_rr'])
            sig = '*' if p < 0.05 else ''
            print(f'{v:25s} {r:+8.3f} {p:8.4f} {r**2:8.3f} {sig}')
            corr_results[v] = {'r': float(r), 'p': float(p), 'r2': float(r**2)}

# ==========================================
# 2. STEPWISE R2
# ==========================================
print('\n' + '='*60)
print('2. STEPWISE MODEL COMPARISON')
print('='*60)

valid = edf.dropna(subset=['log_rr','ncep_z500_nh','t2m','tp'])
y_ev = valid['log_rr'].values
ss_tot = np.sum((y_ev - y_ev.mean())**2)

X1 = np.column_stack([valid['ncep_z500_nh'].values, np.ones(len(valid))])
b1, _, _, _ = lstsq(X1, y_ev, rcond=None)
r2_1 = 1 - np.sum((y_ev - X1@b1)**2)/ss_tot

X2 = np.column_stack([valid['t2m'].values, valid['tp'].values, np.ones(len(valid))])
b2, _, _, _ = lstsq(X2, y_ev, rcond=None)
r2_2 = 1 - np.sum((y_ev - X2@b2)**2)/ss_tot

X3 = np.column_stack([valid['ncep_z500_nh'].values, valid['t2m'].values, valid['tp'].values, np.ones(len(valid))])
b3, _, _, _ = lstsq(X3, y_ev, rcond=None)
r2_3 = 1 - np.sum((y_ev - X3@b3)**2)/ss_tot

print(f'Z500 only:        R2 = {r2_1:.3f}')
print(f'Surface only:     R2 = {r2_2:.3f}')
print(f'Z500 + surface:   R2 = {r2_3:.3f}')
print(f'Z500 unique:      dR2 = {r2_3 - r2_2:.3f}')
print(f'Surface unique:   dR2 = {r2_3 - r2_1:.3f}')

# ==========================================
# 3. COMPLETE MECHANISM CHAIN
# ==========================================
print('\n' + '='*60)
print('3. COMPLETE MECHANISM CHAIN (daily-level)')
print('='*60)

# Merge panel with ERA5
md = panel.merge(era5[['date','t2m_K','tp_mm']], on='date', how='inner')
md = md[md['is_winter']==1].copy()
md['t2m_C'] = md['t2m_K'] - 273.15
md['is_ssw'] = 0
for sd in ssw_dates:
    mask = (md['date'] >= sd - pd.Timedelta(days=15)) & (md['date'] <= sd + pd.Timedelta(days=15))
    md.loc[mask, 'is_ssw'] = 1

ssw_days = md[md['is_ssw']==1]
ctrl_days = md[md['is_ssw']==0]

# Link 1: SSW -> T50 warming
if 'ncep_t_50hpa' in md.columns:
    d1 = (ssw_days['ncep_t_50hpa'].mean() - ctrl_days['ncep_t_50hpa'].mean()) / ctrl_days['ncep_t_50hpa'].std()
    _, p1 = stats.mannwhitneyu(ssw_days['ncep_t_50hpa'].dropna(), ctrl_days['ncep_t_50hpa'].dropna())
    print(f'Link 1: SSW -> T(50hPa) warming')
    print(f'  delta = +{ssw_days["ncep_t_50hpa"].mean() - ctrl_days["ncep_t_50hpa"].mean():.1f}K, d = {d1:.2f}, P = {p1:.2e}')

# Link 2: SSW -> Z500 depression
if 'ncep_z500_nh' in md.columns:
    d2 = (ssw_days['ncep_z500_nh'].mean() - ctrl_days['ncep_z500_nh'].mean()) / ctrl_days['ncep_z500_nh'].std()
    _, p2 = stats.mannwhitneyu(ssw_days['ncep_z500_nh'].dropna(), ctrl_days['ncep_z500_nh'].dropna())
    print(f'\nLink 2: SSW -> Z500 depression')
    print(f'  delta = {ssw_days["ncep_z500_nh"].mean() - ctrl_days["ncep_z500_nh"].mean():.1f}m, d = {d2:.2f}, P = {p2:.2e}')

# Link 3: SSW -> surface T cooling
d3 = (ssw_days['t2m_C'].mean() - ctrl_days['t2m_C'].mean()) / ctrl_days['t2m_C'].std()
_, p3 = stats.mannwhitneyu(ssw_days['t2m_C'].dropna(), ctrl_days['t2m_C'].dropna())
print(f'\nLink 3: SSW -> Surface cooling')
print(f'  delta = {ssw_days["t2m_C"].mean() - ctrl_days["t2m_C"].mean():.2f}C, d = {d3:.2f}, P = {p3:.2e}')

# Link 4: SSW -> precip reduction
d4 = (ssw_days['tp_mm'].mean() - ctrl_days['tp_mm'].mean()) / ctrl_days['tp_mm'].std()
_, p4 = stats.mannwhitneyu(ssw_days['tp_mm'].dropna(), ctrl_days['tp_mm'].dropna())
print(f'\nLink 4: SSW -> Precip reduction')
print(f'  delta = {ssw_days["tp_mm"].mean() - ctrl_days["tp_mm"].mean():.3f}mm, d = {d4:.2f}, P = {p4:.4f}')

# Link 5: Daily T/precip -> avalanche
r_ta, p_ta = stats.pearsonr(md['t2m_C'].dropna(), md.loc[md['t2m_C'].notna(), 'dry_natural_size_1234'])
r_pa, p_pa = stats.pearsonr(md['tp_mm'].dropna(), md.loc[md['tp_mm'].notna(), 'dry_natural_size_1234'])
print(f'\nLink 5a: Daily T2m -> aval: r = {r_ta:.3f}, P = {p_ta:.2e}')
print(f'Link 5b: Daily precip -> aval: r = {r_pa:.3f}, P = {p_pa:.2e}')

# ==========================================
# 4. DAILY MEDIATION BOOTSTRAP
# ==========================================
print('\n' + '='*60)
print('4. DAILY MEDIATION ANALYSIS')
print('='*60)

med_data = md[['is_ssw','dry_natural_size_1234','ncep_z500_nh','t2m_C','tp_mm']].dropna()
print(f'Complete cases: {len(med_data)}')

y = med_data['dry_natural_size_1234'].values
x = med_data['is_ssw'].values
z = med_data['ncep_z500_nh'].values
t = med_data['t2m_C'].values
p = med_data['tp_mm'].values

# Total
bt, _, _, _ = lstsq(np.column_stack([x, np.ones(len(x))]), y, rcond=None)
c_total = bt[0]

# Direct (Z500 + T + precip as mediators)
bf, _, _, _ = lstsq(np.column_stack([x, z, t, p, np.ones(len(x))]), y, rcond=None)
c_prime = bf[0]
indirect = c_total - c_prime
pct = indirect / c_total * 100 if c_total != 0 else 0

print(f'c_total = {c_total:.6f}')
print(f"c_prime = {c_prime:.6f}")
print(f'indirect = {indirect:.6f} ({pct:.1f}% mediated)')

# Z500 only
bz, _, _, _ = lstsq(np.column_stack([x, z, np.ones(len(x))]), y, rcond=None)
indirect_z = c_total - bz[0]
pct_z = indirect_z / c_total * 100 if c_total != 0 else 0
print(f'Z500 mediation: {pct_z:.1f}%')

# Bootstrap
print('\nBootstrapping (5000 iter)...')
n_boot = 5000
ind_boots_full = []
ind_boots_z500 = []
np.random.seed(42)
n = len(y)
for b in range(n_boot):
    idx = np.random.randint(0, n, n)
    try:
        bt_b, _, _, _ = lstsq(np.column_stack([x[idx], np.ones(n)]), y[idx], rcond=None)
        bf_b, _, _, _ = lstsq(np.column_stack([x[idx], z[idx], t[idx], p[idx], np.ones(n)]), y[idx], rcond=None)
        bz_b, _, _, _ = lstsq(np.column_stack([x[idx], z[idx], np.ones(n)]), y[idx], rcond=None)
        ind_boots_full.append(bt_b[0] - bf_b[0])
        ind_boots_z500.append(bt_b[0] - bz_b[0])
    except:
        pass

ind_full = np.array(ind_boots_full)
ind_z = np.array(ind_boots_z500)
ci_full = np.percentile(ind_full, [2.5, 97.5])
ci_z = np.percentile(ind_z, [2.5, 97.5])

print(f'Full mediation (Z500+T+precip): {np.mean(ind_full):.6f} [{ci_full[0]:.6f}, {ci_full[1]:.6f}]')
print(f'  Significant: {not (ci_full[0] <= 0 <= ci_full[1])}')
print(f'Z500 only: {np.mean(ind_z):.6f} [{ci_z[0]:.6f}, {ci_z[1]:.6f}]')
print(f'  Significant: {not (ci_z[0] <= 0 <= ci_z[1])}')

# ==========================================
# 5. WEATHER REGIME DETAIL
# ==========================================
print('\n' + '='*60)
print('5. WEATHER REGIME ANALYSIS')
print('='*60)

t_med = md['t2m_C'].median()
p_med = md['tp_mm'].median()
md['regime'] = 'warm_wet'
md.loc[(md['t2m_C'] < t_med) & (md['tp_mm'] >= p_med), 'regime'] = 'cold_wet'
md.loc[(md['t2m_C'] >= t_med) & (md['tp_mm'] < p_med), 'regime'] = 'warm_dry'
md.loc[(md['t2m_C'] < t_med) & (md['tp_mm'] < p_med), 'regime'] = 'cold_dry'

print(f'{"Regime":10s} | {"SSW%":>6s} | {"Ctrl%":>6s} | {"Shift":>6s} | {"Aval/day":>8s} | {"SSW aval":>8s} | {"Ctrl aval":>9s}')
print('-'*75)
for regime in ['cold_dry','cold_wet','warm_dry','warm_wet']:
    ssw_f = (md.loc[md['is_ssw']==1, 'regime'] == regime).mean() * 100
    ctrl_f = (md.loc[md['is_ssw']==0, 'regime'] == regime).mean() * 100
    aval_all = md.loc[md['regime']==regime, 'dry_natural_size_1234'].mean()
    aval_ssw = md.loc[(md['is_ssw']==1) & (md['regime']==regime), 'dry_natural_size_1234'].mean()
    aval_ctrl = md.loc[(md['is_ssw']==0) & (md['regime']==regime), 'dry_natural_size_1234'].mean()
    print(f'{regime:10s} | {ssw_f:5.1f}% | {ctrl_f:5.1f}% | {ssw_f-ctrl_f:+5.1f}% | {aval_all:8.4f} | {aval_ssw:8.4f} | {aval_ctrl:9.4f}')

from scipy.stats import chi2_contingency
ct = pd.crosstab(md['is_ssw'], md['regime'])
chi2, p_chi, _, _ = chi2_contingency(ct)
print(f'\nChi-square: {chi2:.1f}, P = {p_chi:.2e}')

# Within-regime SSW effect
print('\nWithin-regime SSW effect (controls for weather):')
for regime in ['cold_dry','cold_wet','warm_dry','warm_wet']:
    ssw_r = md[(md['is_ssw']==1) & (md['regime']==regime)]['dry_natural_size_1234']
    ctrl_r = md[(md['is_ssw']==0) & (md['regime']==regime)]['dry_natural_size_1234']
    if len(ssw_r) > 10 and len(ctrl_r) > 10:
        _, p_within = stats.mannwhitneyu(ssw_r, ctrl_r, alternative='two-sided')
        rr_within = ssw_r.mean() / ctrl_r.mean() if ctrl_r.mean() > 0 else np.nan
        print(f'  {regime:10s}: RR = {rr_within:.3f}, P = {p_within:.4f}')

# Save
results = {
    'event_correlations': {k: v for k, v in corr_results.items()},
    'stepwise_r2': {'z500': float(r2_1), 'surface': float(r2_2), 'combined': float(r2_3)},
    'mediation': {
        'pct_full': float(pct), 'pct_z500': float(pct_z),
        'ci_full': [float(ci_full[0]), float(ci_full[1])],
        'ci_z500': [float(ci_z[0]), float(ci_z[1])],
        'sig_full': bool(not (ci_full[0] <= 0 <= ci_full[1])),
        'sig_z500': bool(not (ci_z[0] <= 0 <= ci_z[1]))
    },
    'regime': {'chi2': float(chi2), 'p': float(p_chi)},
    'chain': {
        'link1_d': float(d1) if 'd1' in dir() else None,
        'link2_d': float(d2) if 'd2' in dir() else None,
        'link3_d': float(d3),
        'link4_d': float(d4)
    }
}
with open('data/results/r25_mechanism_deep_dive.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print('\nResults saved')
