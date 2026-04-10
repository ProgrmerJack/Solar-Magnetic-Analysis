"""
R11 Mechanism Upgrade Part 2: Daily Chain + Phase-Resolved Analysis
===================================================================
High-power daily-level tests and phase-resolved mechanism chain.
"""
import pandas as pd
import numpy as np
from scipy import stats
import json
import warnings
warnings.filterwarnings('ignore')

results = {}

# Load data
panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
era5 = pd.read_parquet('data/processed/era5_swiss_alps_extended.parquet')
ssw_cat = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')

panel.index = pd.to_datetime(panel.index)
era5.index = pd.to_datetime(era5.index)
ssw_cat.index = ssw_cat.index.tz_localize(None)

# Merge
era5_m = era5[['t2m_K', 'tp_mm', 'sf_mm', 'wind_speed']].copy()
era5_m.index.name = 'time'
merged = panel.join(era5_m, how='left')

# Climatologies
panel['doy'] = panel.index.dayofyear
era5['doy'] = era5.index.dayofyear
era5_clim = era5.groupby('doy')['t2m_K'].mean()

# Re-merge after adding doy to panel
era5_m = era5[['t2m_K', 'tp_mm', 'sf_mm', 'wind_speed']].copy()
era5_m.index.name = 'time'
merged = panel.join(era5_m, how='left')

for col in ['ncep_t_10hpa', 'ncep_t_50hpa', 'ncep_t_100hpa', 'ncep_u_10hpa',
            'ncep_z500_nh', 'ncep_u850_nh']:
    if col in panel.columns:
        clim = panel.groupby('doy')[col].mean()
        merged[f'{col}_anom'] = merged.apply(
            lambda r: r[col] - clim.get(r['doy'], r[col]) if not pd.isna(r[col]) else np.nan,
            axis=1
        )

# Compute surface T anomaly
merged['t2m_anom'] = merged.apply(
    lambda r: r['t2m_K'] - era5_clim.get(r['doy'], r['t2m_K']) if not pd.isna(r.get('t2m_K', np.nan)) else np.nan,
    axis=1
)

winter_mask = merged['doy'].apply(lambda x: x >= 305 or x <= 120)
winter = merged[winter_mask].copy()

ssw_dates = ssw_cat[(ssw_cat.index >= '1998-12-01') & (ssw_cat.index <= '2019-12-31')].index.tolist()

# ============================================================
# 1. DAILY-LEVEL CHAIN ANALYSIS (high power)
# ============================================================
print("=== 1. DAILY-LEVEL CHAIN ANALYSIS ===")

# Create lagged variables
for lag_d in [1, 3, 5, 7, 10, 14]:
    winter[f'strat_t10_anom_lag{lag_d}'] = winter['ncep_t_10hpa_anom'].shift(lag_d)
    winter[f't2m_anom_lag{lag_d}'] = winter['t2m_anom'].shift(lag_d)
    winter[f'ncep_u_10hpa_anom_lag{lag_d}'] = winter['ncep_u_10hpa_anom'].shift(lag_d)

# Daily chain correlations at various lags
print("\n--- Strat T (10 hPa) -> Surface T (2m) at various lags ---")
daily_chain = {}
for lag_d in [1, 3, 5, 7, 10, 14]:
    col_lag = f'strat_t10_anom_lag{lag_d}'
    valid = winter[['t2m_anom', col_lag]].dropna()
    if len(valid) > 50:
        r, p = stats.spearmanr(valid[col_lag], valid['t2m_anom'])
        daily_chain[f'strat_t10_to_t2m_lag{lag_d}'] = {'r': round(r, 4), 'p': round(p, 6), 'n': len(valid)}
        sig = '*' if p < 0.05 else ''
        print(f"  Lag {lag_d:2d}d: r={r:.4f}, P={p:.6f} (n={len(valid)}) {sig}")

# U10 deceleration -> Surface T
print("\n--- Strat U (10 hPa) deceleration -> Surface T at various lags ---")
winter['u10_decel'] = -winter['ncep_u_10hpa_anom'].diff()
for lag_d in [1, 3, 5, 7, 10, 14]:
    winter[f'u10_decel_lag{lag_d}'] = winter['u10_decel'].shift(lag_d)
    valid = winter[['t2m_anom', f'u10_decel_lag{lag_d}']].dropna()
    if len(valid) > 50:
        r, p = stats.spearmanr(valid[f'u10_decel_lag{lag_d}'], valid['t2m_anom'])
        sig = '*' if p < 0.05 else ''
        print(f"  Lag {lag_d:2d}d: r={r:.4f}, P={p:.6f} (n={len(valid)}) {sig}")
        daily_chain[f'u10_decel_to_t2m_lag{lag_d}'] = {'r': round(r, 4), 'p': round(p, 6)}

# Surface T -> Avalanche count at various lags
print("\n--- Surface T (2m) -> Dry Slab Count at various lags ---")
for lag_d in [0, 1, 3, 5, 7, 10, 14]:
    if lag_d == 0:
        t_col = 't2m_anom'
    else:
        t_col = f't2m_anom_lag{lag_d}'
    valid = winter[[t_col, 'aai_all_dry']].dropna()
    if len(valid) > 50:
        r, p = stats.spearmanr(valid[t_col], valid['aai_all_dry'])
        sig = '*' if p < 0.05 else ''
        print(f"  Lag {lag_d:2d}d: r={r:.4f}, P={p:.6f} (n={len(valid)}) {sig}")
        daily_chain[f't2m_to_aval_lag{lag_d}'] = {'r': round(r, 4), 'p': round(p, 6)}

# Strat T -> Avalanche (skip surface)
print("\n--- Strat T (10 hPa) -> Dry Slab Count at various lags ---")
for lag_d in [0, 3, 7, 10, 14, 21]:
    if lag_d == 0:
        col = 'ncep_t_10hpa_anom'
    else:
        if f'strat_t10_anom_lag{lag_d}' in winter.columns:
            col = f'strat_t10_anom_lag{lag_d}'
        else:
            winter[f'strat_t10_anom_lag{lag_d}'] = winter['ncep_t_10hpa_anom'].shift(lag_d)
            col = f'strat_t10_anom_lag{lag_d}'
    valid = winter[[col, 'aai_all_dry']].dropna()
    if len(valid) > 50:
        r, p = stats.spearmanr(valid[col], valid['aai_all_dry'])
        sig = '*' if p < 0.05 else ''
        print(f"  Lag {lag_d:2d}d: r={r:.4f}, P={p:.6f} (n={len(valid)}) {sig}")
        daily_chain[f'strat_t10_to_aval_lag{lag_d}'] = {'r': round(r, 4), 'p': round(p, 6)}

results['daily_chain'] = daily_chain

# ============================================================
# 2. PHASE-RESOLVED MECHANISM CHAIN
# ============================================================
print("\n=== 2. PHASE-RESOLVED MECHANISM CHAIN ===")

phases = [
    ('Early pre-SSW', -15, -8),
    ('Late pre-SSW', -7, -1),
    ('Onset', 0, 3),
    ('Early post-SSW', 4, 10),
    ('Late post-SSW', 11, 15),
    ('Cold reversal', 16, 30),
]

phase_data = []
for phase_name, lag_start, lag_end in phases:
    t2m_anoms = []
    strat_anoms = []
    aval_rates = []
    ctrl_aval_rates = []
    sf_anoms = []
    
    for ssw_date in ssw_dates:
        # Collect daily values in this phase
        phase_t2m = []
        phase_strat = []
        phase_aval = []
        phase_sf = []
        
        for offset in range(lag_start, lag_end + 1):
            day = ssw_date + pd.Timedelta(days=offset)
            doy = day.timetuple().tm_yday
            
            # Surface temperature anomaly
            if day in era5.index:
                t2m = era5.loc[day, 't2m_K']
                if isinstance(t2m, pd.Series):
                    t2m = t2m.iloc[0]
                clim_t = era5_clim.get(doy, t2m)
                phase_t2m.append(float(t2m) - float(clim_t))
                
                # Snowfall
                sf = era5.loc[day, 'sf_mm']
                if isinstance(sf, pd.Series):
                    sf = sf.iloc[0]
                sf_clim = era5.groupby('doy')['sf_mm'].mean().get(doy, sf)
                phase_sf.append(float(sf) - float(sf_clim))
            
            # Strat temperature anomaly
            if day in panel.index and 'ncep_t_10hpa_anom' in panel.columns:
                st = panel.loc[day, 'ncep_t_10hpa_anom']
                if isinstance(st, pd.Series):
                    st = st.iloc[0]
                if not pd.isna(st):
                    phase_strat.append(float(st))
            
            # Avalanche
            if day in panel.index:
                av = panel.loc[day, 'aai_all_dry']
                if isinstance(av, pd.Series):
                    av = av.iloc[0]
                if not pd.isna(av):
                    phase_aval.append(float(av))
        
        if phase_t2m:
            t2m_anoms.append(np.mean(phase_t2m))
        if phase_strat:
            strat_anoms.append(np.mean(phase_strat))
        if phase_aval:
            aval_rates.append(np.mean(phase_aval))
        if phase_sf:
            sf_anoms.append(np.mean(phase_sf))
    
    # Phase-level statistics
    mean_t2m = np.mean(t2m_anoms) if t2m_anoms else np.nan
    mean_strat = np.mean(strat_anoms) if strat_anoms else np.nan
    mean_aval = np.mean(aval_rates) if aval_rates else np.nan
    mean_sf = np.mean(sf_anoms) if sf_anoms else np.nan
    
    # T-test vs 0 for temperature anomaly
    t_t2m, p_t2m = stats.ttest_1samp(t2m_anoms, 0) if len(t2m_anoms) > 2 else (np.nan, np.nan)
    t_strat, p_strat = stats.ttest_1samp(strat_anoms, 0) if len(strat_anoms) > 2 else (np.nan, np.nan)
    
    phase_data.append({
        'phase': phase_name,
        'lag_range': f'{lag_start} to {lag_end}',
        't2m_anom_K': round(float(mean_t2m), 3),
        't2m_p': round(float(p_t2m), 4) if not np.isnan(p_t2m) else None,
        'strat_t10_anom_K': round(float(mean_strat), 2),
        'strat_p': round(float(p_strat), 4) if not np.isnan(p_strat) else None,
        'mean_aval_rate': round(float(mean_aval), 3),
        'sf_anom_mm': round(float(mean_sf), 3) if not np.isnan(mean_sf) else None,
        'n_events': len(t2m_anoms),
    })
    
    print(f"\n  {phase_name} (days {lag_start} to {lag_end}):")
    print(f"    Strat T10: {mean_strat:+.1f} K (P={p_strat:.4f})")
    print(f"    Surface T: {mean_t2m:+.3f} K (P={p_t2m:.4f})")
    print(f"    Snowfall:  {mean_sf:+.3f} mm/d")
    print(f"    Aval rate: {mean_aval:.3f}/day")

# Cross-phase correlations: does phase-level T predict phase-level avalanche rate?
print("\n--- Phase-level cross-correlations ---")
pdf = pd.DataFrame(phase_data)
# Exclude cold reversal for the mechanism test (it's the falsification)
pdf_main = pdf[pdf['phase'] != 'Cold reversal']

if len(pdf_main) >= 4:
    r_t2m_aval, p_t2m_aval = stats.spearmanr(pdf_main['t2m_anom_K'], pdf_main['mean_aval_rate'])
    print(f"  Phase T2m vs Aval rate (n={len(pdf_main)} phases): r={r_t2m_aval:.3f}, P={p_t2m_aval:.4f}")
    
    r_strat_aval, p_strat_aval = stats.spearmanr(pdf_main['strat_t10_anom_K'], pdf_main['mean_aval_rate'])
    print(f"  Phase Strat T vs Aval rate: r={r_strat_aval:.3f}, P={p_strat_aval:.4f}")
    
    r_strat_t2m, p_strat_t2m = stats.spearmanr(pdf_main['strat_t10_anom_K'], pdf_main['t2m_anom_K'])
    print(f"  Phase Strat T vs Surface T: r={r_strat_t2m:.3f}, P={p_strat_t2m:.4f}")

# Full 6-phase analysis
if len(pdf) >= 5:
    r_full, p_full = stats.spearmanr(pdf['t2m_anom_K'], pdf['mean_aval_rate'])
    print(f"  Full 6-phase T2m vs Aval (with cold reversal): r={r_full:.3f}, P={p_full:.4f}")

results['phase_resolved'] = {
    'phases': phase_data,
    'phase_t2m_vs_aval_r': round(float(r_t2m_aval), 3) if 'r_t2m_aval' in dir() else None,
    'phase_t2m_vs_aval_p': round(float(p_t2m_aval), 4) if 'p_t2m_aval' in dir() else None,
    'full_phase_r': round(float(r_full), 3) if 'r_full' in dir() else None,
    'full_phase_p': round(float(p_full), 4) if 'p_full' in dir() else None,
}

# ============================================================
# 3. SINTERING BY SSW TYPE (fixed)
# ============================================================
print("\n=== 3. SINTERING BY SSW TYPE ===")

published_types = {
    '1998-12-15': 'D', '1999-02-26': 'D', '2001-02-11': 'D',
    '2001-12-30': 'S', '2002-02-17': 'S', '2003-01-18': 'S',
    '2004-01-05': 'D', '2006-01-21': 'D', '2007-02-24': 'D',
    '2008-02-22': 'D', '2009-01-24': 'S', '2010-02-09': 'D',
    '2012-01-11': 'D', '2013-01-07': 'S', '2018-02-12': 'S',
    '2019-01-01': 'D',
}

with open('data/results/sintering_extended.json', 'r') as f:
    sintering = json.load(f)

for se in sintering['per_event']:
    se['type'] = published_types.get(se['ssw_date'], '?')
    se['warming'] = se['delta_T_K'] > 0

disp_s = [s['sintering_enhancement_pct'] for s in sintering['per_event'] if s['type'] == 'D']
split_s = [s['sintering_enhancement_pct'] for s in sintering['per_event'] if s['type'] == 'S']
warm_s = [s['sintering_enhancement_pct'] for s in sintering['per_event'] if s['warming']]
cool_s = [s['sintering_enhancement_pct'] for s in sintering['per_event'] if not s['warming']]

print(f"Displacement (n={len(disp_s)}): mean={np.mean(disp_s):.1f}%, "
      f"positive={sum(1 for x in disp_s if x > 0)}/{len(disp_s)}")
print(f"Split (n={len(split_s)}): mean={np.mean(split_s):.1f}%, "
      f"positive={sum(1 for x in split_s if x > 0)}/{len(split_s)}")
print(f"\nWarming events (n={len(warm_s)}): mean={np.mean(warm_s):.1f}%, "
      f"positive={sum(1 for x in warm_s if x > 0)}/{len(warm_s)}")
print(f"Cooling events (n={len(cool_s)}): mean={np.mean(cool_s):.1f}%, "
      f"positive={sum(1 for x in cool_s if x > 0)}/{len(cool_s)}")

# Warming-only sintering significance
if warm_s:
    t_ws, p_ws = stats.ttest_1samp(warm_s, 0)
    sign_ws = sum(1 for x in warm_s if x > 0)
    sign_p_ws = stats.binomtest(sign_ws, len(warm_s), 0.5).pvalue
    print(f"\nWarming-only sintering test:")
    print(f"  Mean: {np.mean(warm_s):.1f}%, Median: {np.median(warm_s):.1f}%")
    print(f"  t-test vs 0: P = {p_ws:.4f}")
    print(f"  Sign test: {sign_ws}/{len(warm_s)}, P = {sign_p_ws:.4f}")

# Displacement-only sintering
if disp_s:
    t_ds, p_ds = stats.ttest_1samp(disp_s, 0)
    sign_ds = sum(1 for x in disp_s if x > 0)
    sign_p_ds = stats.binomtest(sign_ds, len(disp_s), 0.5).pvalue
    print(f"\nDisplacement-only sintering test:")
    print(f"  Mean: {np.mean(disp_s):.1f}%, Median: {np.median(disp_s):.1f}%")
    print(f"  t-test vs 0: P = {p_ds:.4f}")
    print(f"  Sign test: {sign_ds}/{len(disp_s)}, P = {sign_p_ds:.4f}")

results['sintering_by_type'] = {
    'displacement': {'n': len(disp_s), 'mean_pct': round(np.mean(disp_s), 1), 
                     'positive': sum(1 for x in disp_s if x > 0),
                     'p_ttest': round(float(p_ds), 4) if 'p_ds' in dir() else None},
    'split': {'n': len(split_s), 'mean_pct': round(np.mean(split_s), 1),
              'positive': sum(1 for x in split_s if x > 0)},
    'warming': {'n': len(warm_s), 'mean_pct': round(np.mean(warm_s), 1),
                'positive': sum(1 for x in warm_s if x > 0),
                'p_ttest': round(float(p_ws), 4) if 'p_ws' in dir() else None,
                'p_sign': round(float(sign_p_ws), 4) if 'sign_p_ws' in dir() else None},
    'cooling': {'n': len(cool_s), 'mean_pct': round(np.mean(cool_s), 1),
                'positive': sum(1 for x in cool_s if x > 0)},
}

# ============================================================
# 4. EDDY HEAT FLUX PROXY (wave activity index)
# ============================================================
print("\n=== 4. EDDY HEAT FLUX PROXY ===")

# Newman et al. (2001) proxy: wave activity ~ -d[u_10]/dt
# Compute 5-day smoothed wave activity index
winter['wave_activity'] = -winter['ncep_u_10hpa_anom'].diff().rolling(5, center=True).mean()

# Does wave activity predict surface temperature with a lag?
print("\n--- Wave activity -> Surface T at various lags ---")
for lag_d in [0, 3, 5, 7, 10, 14, 21]:
    if lag_d > 0:
        winter[f'wave_lag{lag_d}'] = winter['wave_activity'].shift(lag_d)
        valid = winter[['t2m_anom', f'wave_lag{lag_d}']].dropna()
    else:
        valid = winter[['t2m_anom', 'wave_activity']].dropna()
    if len(valid) > 50:
        r, p = stats.spearmanr(valid.iloc[:, 1], valid.iloc[:, 0])
        sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else ''))
        print(f"  Lag {lag_d:2d}d: r={r:.4f}, P={p:.6f} (n={len(valid)}) {sig}")

# Does wave activity predict avalanche counts?
print("\n--- Wave activity -> Avalanche count at various lags ---")
for lag_d in [0, 3, 5, 7, 10, 14, 21]:
    if lag_d > 0:
        col = f'wave_lag{lag_d}'
        if col not in winter.columns:
            winter[col] = winter['wave_activity'].shift(lag_d)
    else:
        col = 'wave_activity'
    valid = winter[[col, 'aai_all_dry']].dropna()
    if len(valid) > 50:
        r, p = stats.spearmanr(valid[col], valid['aai_all_dry'])
        sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else ''))
        print(f"  Lag {lag_d:2d}d: r={r:.4f}, P={p:.6f} (n={len(valid)}) {sig}")

# Cumulative wave activity over pre-SSW window predicts event-level response?
print("\n--- Cumulative pre-SSW wave activity vs event response ---")
cum_wave = []
for ssw_date in ssw_dates:
    wave_vals = []
    for offset in range(-20, 0):
        day = ssw_date + pd.Timedelta(days=offset)
        if day in winter.index and 'wave_activity' in winter.columns:
            val = winter.loc[day, 'wave_activity']
            if isinstance(val, pd.Series):
                val = val.iloc[0]
            if not pd.isna(val):
                wave_vals.append(float(val))
    cum_wave.append(np.sum(wave_vals) if wave_vals else np.nan)

# Load event catalog from previous script
ec = pd.read_csv('data/results/ssw_event_catalog.csv')
ec['cum_wave'] = cum_wave

valid_cw = ec[['cum_wave', 'rr', 'surface_t_anom_K']].dropna()
if len(valid_cw) >= 5:
    r_cw_rr, p_cw_rr = stats.spearmanr(valid_cw['cum_wave'], valid_cw['rr'])
    r_cw_t, p_cw_t = stats.spearmanr(valid_cw['cum_wave'], valid_cw['surface_t_anom_K'])
    print(f"  Cumulative wave -> Aval RR: r={r_cw_rr:.3f}, P={p_cw_rr:.4f}")
    print(f"  Cumulative wave -> Surface T: r={r_cw_t:.3f}, P={p_cw_t:.4f}")

# ============================================================
# SAVE
# ============================================================
print("\n=== SAVING ===")
with open('data/results/r11_mechanism_upgrade_p2.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print("Saved to data/results/r11_mechanism_upgrade_p2.json")
print("\n=== COMPLETE ===")
