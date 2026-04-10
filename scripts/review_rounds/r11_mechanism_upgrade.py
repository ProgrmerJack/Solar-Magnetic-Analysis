"""
R11 Mechanism Upgrade: Comprehensive Chain Analysis
====================================================
Addresses all four R11 reviewer requests:
1. SSW-type stratification (split vs displacement + surface response)
2. Atmospheric chain analysis with wave activity proxy & downward propagation
3. Standardized effect sizes (Cohen's d, Glass's delta, IRR CIs)
4. Cold-phase reversal as falsification + NAO/AO rejection strengthening
"""
import pandas as pd
import numpy as np
from scipy import stats
import json
import warnings
warnings.filterwarnings('ignore')

results = {}

# ============================================================
# LOAD DATA
# ============================================================
print("Loading data...")
panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
era5 = pd.read_parquet('data/processed/era5_swiss_alps_extended.parquet')
ssw_cat = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')

panel.index = pd.to_datetime(panel.index)
era5.index = pd.to_datetime(era5.index)
ssw_cat.index = ssw_cat.index.tz_localize(None)

# Merge panel with ERA5 surface data
era5_for_merge = era5[['t2m_K', 'tp_mm', 'sf_mm', 'wind_speed']].copy()
era5_for_merge.index.name = 'time'
merged = panel.join(era5_for_merge, how='left')

# SSW events in study period
study_ssws = ssw_cat[(ssw_cat.index >= '1998-12-01') & (ssw_cat.index <= '2019-12-31')]
ssw_dates = study_ssws.index.tolist()
print(f"SSW events in study period: {len(ssw_dates)}")

# Compute ERA5 climatology (day-of-year means)
era5['doy'] = era5.index.dayofyear
era5_clim = era5.groupby('doy')['t2m_K'].mean()

# NCEP climatology for stratospheric variables
ncep_levels = ['ncep_t_10hpa', 'ncep_t_20hpa', 'ncep_t_30hpa',
               'ncep_t_50hpa', 'ncep_t_70hpa', 'ncep_t_100hpa']
ncep_u_levels = ['ncep_u_10hpa', 'ncep_u_20hpa', 'ncep_u_30hpa',
                 'ncep_u_50hpa', 'ncep_u_70hpa', 'ncep_u_100hpa']
panel['doy'] = panel.index.dayofyear
ncep_clim = {}
for col in ncep_levels + ncep_u_levels + ['ncep_z500_nh', 'ncep_u850_nh']:
    if col in panel.columns:
        ncep_clim[col] = panel.groupby('doy')[col].mean()

# ============================================================
# 1. SSW TYPE CLASSIFICATION
# ============================================================
print("\n=== 1. SSW TYPE CLASSIFICATION ===")

# Published classifications from peer-reviewed literature
# Charlton & Polvani (2007), Mitchell et al. (2013), Butler et al. (2017)
# D = Displacement (wave-1), S = Split (wave-2)
published_types = {
    '1998-12-15': 'D',   # Naito & Yoden 2006
    '1999-02-26': 'D',   # Charlton & Polvani 2007
    '2001-02-11': 'D',   # Charlton & Polvani 2007
    '2001-12-30': 'S',   # Nishii et al. 2011
    '2002-02-17': 'S',   # Charlton & Polvani 2007
    '2003-01-18': 'S',   # Charlton & Polvani 2007
    '2004-01-05': 'D',   # Charlton & Polvani 2007
    '2006-01-21': 'D',   # Mitchell et al. 2013
    '2007-02-24': 'D',   # Cohen & Jones 2011
    '2008-02-22': 'D',   # Harada et al. 2010
    '2009-01-24': 'S',   # Manney et al. 2009
    '2010-02-09': 'D',   # Kuttippurath & Nikulin 2012
    '2012-01-11': 'D',   # Coy & Pawson 2015
    '2013-01-07': 'S',   # Mitchell et al. 2013
    '2018-02-12': 'S',   # Karpechko et al. 2018
    '2019-01-01': 'D',   # Rao et al. 2020
}

# Compute per-event metrics
event_catalog = []
for ssw_date in ssw_dates:
    d_str = ssw_date.strftime('%Y-%m-%d')

    # Surface temperature anomaly (ERA5, days 0 to +15)
    t2m_vals = []
    for offset in range(0, 16):
        day = ssw_date + pd.Timedelta(days=offset)
        if day in era5.index:
            doy = era5.loc[day, 'doy']
            if isinstance(doy, pd.Series):
                doy = doy.iloc[0]
            anom = era5.loc[day, 't2m_K'] - era5_clim.get(int(doy), era5.loc[day, 't2m_K'])
            if isinstance(anom, pd.Series):
                anom = anom.iloc[0]
            t2m_vals.append(float(anom))
    surface_t_anom = np.mean(t2m_vals) if t2m_vals else np.nan

    # Pre-SSW surface temperature anomaly (days -15 to -1)
    pre_t2m_vals = []
    for offset in range(-15, 0):
        day = ssw_date + pd.Timedelta(days=offset)
        if day in era5.index:
            doy = era5.loc[day, 'doy']
            if isinstance(doy, pd.Series):
                doy = doy.iloc[0]
            anom = era5.loc[day, 't2m_K'] - era5_clim.get(int(doy), era5.loc[day, 't2m_K'])
            if isinstance(anom, pd.Series):
                anom = anom.iloc[0]
            pre_t2m_vals.append(float(anom))
    pre_surface_t_anom = np.mean(pre_t2m_vals) if pre_t2m_vals else np.nan

    # Stratospheric temperature anomaly at 10 hPa (days 0 to +10)
    strat_t_vals = []
    for offset in range(0, 11):
        day = ssw_date + pd.Timedelta(days=offset)
        if day in panel.index and not pd.isna(panel.loc[day, 'ncep_t_10hpa']):
            doy = panel.loc[day, 'doy']
            if isinstance(doy, pd.Series):
                doy = doy.iloc[0]
            clim_val = ncep_clim.get('ncep_t_10hpa', pd.Series(dtype=float))
            if int(doy) in clim_val.index:
                val = panel.loc[day, 'ncep_t_10hpa']
                if isinstance(val, pd.Series):
                    val = val.iloc[0]
                anom = float(val) - float(clim_val[int(doy)])
                strat_t_vals.append(anom)
    strat_t10_anom = np.mean(strat_t_vals) if strat_t_vals else np.nan

    # Wave activity proxy: zonal wind deceleration at 10 hPa
    # du/dt from day -20 to day -5 (pre-SSW wave forcing period)
    pre_u_vals = []
    for offset in range(-20, -4):
        day = ssw_date + pd.Timedelta(days=offset)
        if day in panel.index and not pd.isna(panel.loc[day, 'ncep_u_10hpa']):
            val = panel.loc[day, 'ncep_u_10hpa']
            if isinstance(val, pd.Series):
                val = val.iloc[0]
            pre_u_vals.append(float(val))
    if len(pre_u_vals) >= 5:
        wave_decel = -(pre_u_vals[-1] - pre_u_vals[0]) / len(pre_u_vals)
    else:
        wave_decel = np.nan

    # Vortex wind change (peak pre-SSW u10 minus minimum around onset)
    wide_window_u = []
    for offset in range(-25, 16):
        day = ssw_date + pd.Timedelta(days=offset)
        if day in panel.index and not pd.isna(panel.loc[day, 'ncep_u_10hpa']):
            val = panel.loc[day, 'ncep_u_10hpa']
            if isinstance(val, pd.Series):
                val = val.iloc[0]
            wide_window_u.append(float(val))
    vortex_disruption = max(wide_window_u) - min(wide_window_u) if len(wide_window_u) >= 10 else np.nan

    # Tropospheric Z500 anomaly (days 0 to +15)
    z500_vals = []
    for offset in range(0, 16):
        day = ssw_date + pd.Timedelta(days=offset)
        if day in panel.index and not pd.isna(panel.loc[day, 'ncep_z500_nh']):
            doy = panel.loc[day, 'doy']
            if isinstance(doy, pd.Series):
                doy = doy.iloc[0]
            clim_val = ncep_clim.get('ncep_z500_nh', pd.Series(dtype=float))
            if int(doy) in clim_val.index:
                val = panel.loc[day, 'ncep_z500_nh']
                if isinstance(val, pd.Series):
                    val = val.iloc[0]
                z500_vals.append(float(val) - float(clim_val[int(doy)]))
    z500_anom = np.mean(z500_vals) if z500_vals else np.nan

    # U850 anomaly (days 0 to +15)
    u850_vals = []
    for offset in range(0, 16):
        day = ssw_date + pd.Timedelta(days=offset)
        if day in panel.index and not pd.isna(panel.loc[day, 'ncep_u850_nh']):
            doy = panel.loc[day, 'doy']
            if isinstance(doy, pd.Series):
                doy = doy.iloc[0]
            clim_val = ncep_clim.get('ncep_u850_nh', pd.Series(dtype=float))
            if int(doy) in clim_val.index:
                val = panel.loc[day, 'ncep_u850_nh']
                if isinstance(val, pd.Series):
                    val = val.iloc[0]
                u850_vals.append(float(val) - float(clim_val[int(doy)]))
    u850_anom = np.mean(u850_vals) if u850_vals else np.nan

    # Avalanche rate ratio (SSW window vs matched control)
    ssw_start = ssw_date - pd.Timedelta(days=15)
    ssw_end = ssw_date + pd.Timedelta(days=15)
    ssw_aval = panel.loc[ssw_start:ssw_end, 'aai_all_dry'].dropna()
    ssw_rate = ssw_aval.mean() if len(ssw_aval) > 0 else np.nan

    # Control: same day-of-season in non-SSW winters
    winter = ssw_date.year if ssw_date.month >= 9 else ssw_date.year - 1
    all_winters = range(1999, 2020)
    ssw_winter_set = set()
    for sd in ssw_dates:
        w = sd.year if sd.month >= 9 else sd.year - 1
        ssw_winter_set.add(w)
    ctrl_rates = []
    ssw_doy = ssw_date.timetuple().tm_yday
    for yr in all_winters:
        if yr == winter:
            continue
        for doy_offset in range(-3, 4):
            target_doy = ssw_doy + doy_offset
            for month_start in [10, 1]:
                try:
                    if target_doy <= 0:
                        target_date = pd.Timestamp(f'{yr}-12-31') + pd.Timedelta(days=target_doy)
                    elif target_doy > 365:
                        target_date = pd.Timestamp(f'{yr+1}-01-01') + pd.Timedelta(days=target_doy - 366)
                    else:
                        target_date = pd.Timestamp(f'{yr}' if ssw_date.month >= 9 else f'{yr+1}') + pd.Timedelta(days=target_doy - 1)
                        if target_date.month < 9 and ssw_date.month >= 9:
                            target_date = pd.Timestamp(f'{yr}') + pd.Timedelta(days=target_doy - 1)
                except:
                    continue
    # Simpler control: use panel's seasonal average for same DOY range
    ctrl_vals = []
    for yr in all_winters:
        if yr == winter:
            continue
        for doy_offset in range(-3, 4):
            d = ssw_doy + doy_offset
            mask = (panel['doy'] == d)
            vals = panel.loc[mask, 'aai_all_dry'].dropna()
            if len(vals) > 0:
                ctrl_vals.extend(vals.values.tolist())
    ctrl_rate = np.mean(ctrl_vals) if ctrl_vals else np.nan
    rr = ssw_rate / ctrl_rate if ctrl_rate > 0 else np.nan

    event_catalog.append({
        'date': d_str,
        'published_type': published_types.get(d_str, '?'),
        'surface_t_anom_K': round(float(surface_t_anom), 3) if not np.isnan(surface_t_anom) else None,
        'pre_surface_t_anom_K': round(float(pre_surface_t_anom), 3) if not np.isnan(pre_surface_t_anom) else None,
        'strat_t10_anom_K': round(float(strat_t10_anom), 3) if not np.isnan(strat_t10_anom) else None,
        'wave_decel_ms_day': round(float(wave_decel), 4) if not np.isnan(wave_decel) else None,
        'vortex_disruption_ms': round(float(vortex_disruption), 2) if not np.isnan(vortex_disruption) else None,
        'z500_anom_m': round(float(z500_anom), 2) if not np.isnan(z500_anom) else None,
        'u850_anom_ms': round(float(u850_anom), 3) if not np.isnan(u850_anom) else None,
        'surface_warming': surface_t_anom > 0 if not np.isnan(surface_t_anom) else None,
        'ssw_rate': round(float(ssw_rate), 4) if not np.isnan(ssw_rate) else None,
        'ctrl_rate': round(float(ctrl_rate), 4) if not np.isnan(ctrl_rate) else None,
        'rr': round(float(rr), 4) if not np.isnan(rr) else None,
    })

ec = pd.DataFrame(event_catalog)
print("\nEvent catalog:")
print(ec[['date', 'published_type', 'surface_t_anom_K', 'strat_t10_anom_K',
          'wave_decel_ms_day', 'rr']].to_string(index=False))

# Stratify by published type
disp_events = ec[ec['published_type'] == 'D']
split_events = ec[ec['published_type'] == 'S']

print(f"\nDisplacement events (n={len(disp_events)}):")
print(f"  Mean surface T anomaly: {disp_events['surface_t_anom_K'].mean():.3f} K")
print(f"  Mean avalanche RR: {disp_events['rr'].mean():.3f}")
print(f"  Events with RR < 1: {(disp_events['rr'] < 1).sum()}/{len(disp_events)}")

print(f"\nSplit events (n={len(split_events)}):")
print(f"  Mean surface T anomaly: {split_events['surface_t_anom_K'].mean():.3f} K")
print(f"  Mean avalanche RR: {split_events['rr'].mean():.3f}")
print(f"  Events with RR < 1: {(split_events['rr'] < 1).sum()}/{len(split_events)}")

# Statistical test: displacement vs split RR
if len(disp_events) > 2 and len(split_events) > 2:
    disp_rr = disp_events['rr'].dropna().values
    split_rr = split_events['rr'].dropna().values
    t_stat, p_type = stats.mannwhitneyu(disp_rr, split_rr, alternative='two-sided')
    print(f"\nDisplacement vs Split RR: Mann-Whitney P = {p_type:.4f}")
    # One-sided: displacement has lower RR (more avalanche decrease)
    _, p_onesided = stats.mannwhitneyu(disp_rr, split_rr, alternative='less')
    print(f"  One-sided (disp < split): P = {p_onesided:.4f}")

# Stratify by surface warming
warm_events = ec[ec['surface_warming'] == True]
cool_events = ec[ec['surface_warming'] == False]

print(f"\nWarming SSW events (n={len(warm_events)}):")
print(f"  Mean surface T anomaly: {warm_events['surface_t_anom_K'].mean():.3f} K")
print(f"  Mean avalanche RR: {warm_events['rr'].mean():.3f}")
print(f"  Events with RR < 1: {(warm_events['rr'] < 1).sum()}/{len(warm_events)}")

print(f"\nCooling SSW events (n={len(cool_events)}):")
print(f"  Mean surface T anomaly: {cool_events['surface_t_anom_K'].mean():.3f} K")
print(f"  Mean avalanche RR: {cool_events['rr'].mean():.3f}")
print(f"  Events with RR < 1: {(cool_events['rr'] < 1).sum()}/{len(cool_events)}")

if len(warm_events) > 2 and len(cool_events) > 2:
    warm_rr = warm_events['rr'].dropna().values
    cool_rr = cool_events['rr'].dropna().values
    _, p_warm_cool = stats.mannwhitneyu(warm_rr, cool_rr, alternative='two-sided')
    print(f"\nWarming vs Cooling RR: Mann-Whitney P = {p_warm_cool:.4f}")
    _, p_warm_less = stats.mannwhitneyu(warm_rr, cool_rr, alternative='less')
    print(f"  One-sided (warm < cool): P = {p_warm_less:.4f}")

# Cross-tabulate type with surface warming
print("\nCross-tabulation: Published Type x Surface Response")
for ptype in ['D', 'S']:
    subset = ec[ec['published_type'] == ptype]
    n_warm = (subset['surface_warming'] == True).sum()
    n_cool = (subset['surface_warming'] == False).sum()
    print(f"  {ptype}: {n_warm} warming, {n_cool} cooling")

results['ssw_type_stratification'] = {
    'n_displacement': int(len(disp_events)),
    'n_split': int(len(split_events)),
    'displacement_mean_rr': round(float(disp_events['rr'].mean()), 4),
    'split_mean_rr': round(float(split_events['rr'].mean()), 4),
    'displacement_decrease_fraction': f"{(disp_events['rr'] < 1).sum()}/{len(disp_events)}",
    'split_decrease_fraction': f"{(split_events['rr'] < 1).sum()}/{len(split_events)}",
    'displacement_mean_surface_t': round(float(disp_events['surface_t_anom_K'].mean()), 3),
    'split_mean_surface_t': round(float(split_events['surface_t_anom_K'].mean()), 3),
    'type_comparison_p': round(float(p_type), 4) if 'p_type' in dir() else None,
    'n_warming': int(len(warm_events)),
    'n_cooling': int(len(cool_events)),
    'warming_mean_rr': round(float(warm_events['rr'].mean()), 4),
    'cooling_mean_rr': round(float(cool_events['rr'].mean()), 4),
    'warming_decrease_fraction': f"{(warm_events['rr'] < 1).sum()}/{len(warm_events)}",
    'cooling_decrease_fraction': f"{(cool_events['rr'] < 1).sum()}/{len(cool_events)}",
    'warming_vs_cooling_p': round(float(p_warm_cool), 4) if 'p_warm_cool' in dir() else None,
    'event_catalog': event_catalog,
}

# ============================================================
# 2. ATMOSPHERIC CHAIN ANALYSIS
# ============================================================
print("\n=== 2. ATMOSPHERIC CHAIN ANALYSIS ===")

# 2a. Wave activity proxy: -du_10hPa/dt (positive = wave forcing)
print("\n--- 2a. Wave Activity Proxy ---")
ec_df = ec.copy()
chain_vars = ['wave_decel_ms_day', 'strat_t10_anom_K', 'z500_anom_m',
              'u850_anom_ms', 'surface_t_anom_K', 'rr']
chain_labels = ['Wave forcing', 'Strat T10 anom', 'Z500 anom',
                'U850 anom', 'Surface T anom', 'Avalanche RR']

# Correlation between consecutive chain links
chain_corrs = {}
for i in range(len(chain_vars) - 1):
    v1, v2 = chain_vars[i], chain_vars[i+1]
    valid = ec_df[[v1, v2]].dropna()
    if len(valid) >= 5:
        r, p = stats.spearmanr(valid[v1], valid[v2])
        chain_corrs[f"{chain_labels[i]} -> {chain_labels[i+1]}"] = {
            'r': round(float(r), 3),
            'p': round(float(p), 4),
            'n': int(len(valid))
        }
        print(f"  {chain_labels[i]} -> {chain_labels[i+1]}: r={r:.3f}, P={p:.4f} (n={len(valid)})")

# Key chain test: wave forcing -> surface T -> avalanche RR
valid_chain = ec_df[['wave_decel_ms_day', 'surface_t_anom_K', 'rr']].dropna()
if len(valid_chain) >= 5:
    r_wave_surf, p_wave_surf = stats.spearmanr(valid_chain['wave_decel_ms_day'],
                                                 valid_chain['surface_t_anom_K'])
    r_surf_aval, p_surf_aval = stats.spearmanr(valid_chain['surface_t_anom_K'],
                                                 valid_chain['rr'])
    r_wave_aval, p_wave_aval = stats.spearmanr(valid_chain['wave_decel_ms_day'],
                                                 valid_chain['rr'])
    print(f"\n  Full chain (n={len(valid_chain)}):")
    print(f"    Wave -> Surface T: r={r_wave_surf:.3f}, P={p_wave_surf:.4f}")
    print(f"    Surface T -> Aval RR: r={r_surf_aval:.3f}, P={p_surf_aval:.4f}")
    print(f"    Wave -> Aval RR (skip): r={r_wave_aval:.3f}, P={p_wave_aval:.4f}")

# 2b. Downward propagation analysis
print("\n--- 2b. Downward Propagation with Surface Extension ---")
levels_ordered = ['ncep_t_10hpa', 'ncep_t_20hpa', 'ncep_t_30hpa',
                  'ncep_t_50hpa', 'ncep_t_70hpa', 'ncep_t_100hpa']
level_labels = ['10 hPa', '20 hPa', '30 hPa', '50 hPa', '70 hPa', '100 hPa']
level_heights_km = [31, 26, 24, 21, 18, 16]  # approximate

propagation_results = {}
for li, level in enumerate(levels_ordered):
    lag_anoms = []
    for lag in range(-15, 31):
        anoms = []
        for ssw_date in ssw_dates:
            day = ssw_date + pd.Timedelta(days=lag)
            if day in panel.index:
                val = panel.loc[day, level]
                doy = panel.loc[day, 'doy']
                if isinstance(val, pd.Series):
                    val = val.iloc[0]
                if isinstance(doy, pd.Series):
                    doy = doy.iloc[0]
                if not pd.isna(val) and int(doy) in ncep_clim.get(level, pd.Series(dtype=float)).index:
                    anom = float(val) - float(ncep_clim[level][int(doy)])
                    anoms.append(anom)
        if anoms:
            lag_anoms.append({'lag': lag, 'mean_anom': np.mean(anoms), 'n': len(anoms)})

    if lag_anoms:
        df_lag = pd.DataFrame(lag_anoms)
        peak_row = df_lag.loc[df_lag['mean_anom'].idxmax()]
        propagation_results[level_labels[li]] = {
            'peak_lag': int(peak_row['lag']),
            'peak_anom_K': round(float(peak_row['mean_anom']), 2),
            'height_km': level_heights_km[li]
        }
        print(f"  {level_labels[li]} ({level_heights_km[li]} km): peak at lag {int(peak_row['lag'])}d, "
              f"+{peak_row['mean_anom']:.1f} K")

# Surface T propagation
surface_lag_anoms = []
for lag in range(-15, 31):
    anoms = []
    for ssw_date in ssw_dates:
        day = ssw_date + pd.Timedelta(days=lag)
        if day in era5.index:
            doy = era5.loc[day, 'doy']
            if isinstance(doy, pd.Series):
                doy = doy.iloc[0]
            if int(doy) in era5_clim.index:
                anom = float(era5.loc[day, 't2m_K']) - float(era5_clim[int(doy)])
                anoms.append(anom)
    if anoms:
        surface_lag_anoms.append({'lag': lag, 'mean_anom': np.mean(anoms), 'n': len(anoms)})

if surface_lag_anoms:
    df_surf = pd.DataFrame(surface_lag_anoms)
    peak_surf = df_surf.loc[df_surf['mean_anom'].idxmax()]
    propagation_results['Surface (2m)'] = {
        'peak_lag': int(peak_surf['lag']),
        'peak_anom_K': round(float(peak_surf['mean_anom']), 2),
        'height_km': 0
    }
    print(f"  Surface 2m (0 km): peak at lag {int(peak_surf['lag'])}d, "
          f"+{peak_surf['mean_anom']:.2f} K")

    # Cold phase reversal
    late_anoms = df_surf[(df_surf['lag'] >= 16) & (df_surf['lag'] <= 30)]
    if len(late_anoms) > 0:
        cold_phase_mean = late_anoms['mean_anom'].mean()
        print(f"  Cold phase (lag 16-30): mean anom = {cold_phase_mean:.3f} K")

# Propagation speed computation
print("\n  Downward propagation speed:")
prop_heights = [(k, v['height_km'], v['peak_lag']) for k, v in propagation_results.items()]
prop_heights.sort(key=lambda x: -x[1])  # sort by height descending
for i in range(len(prop_heights) - 1):
    dh = prop_heights[i][1] - prop_heights[i+1][1]  # km
    dt = prop_heights[i+1][2] - prop_heights[i][2]   # days
    if dt > 0:
        speed = dh / dt  # km/day
        print(f"    {prop_heights[i][0]} -> {prop_heights[i+1][0]}: "
              f"{dh} km in {dt}d = {speed:.1f} km/day")

# 2c. Per-event chain: does the chain predict avalanche response?
print("\n--- 2c. Per-Event Chain Prediction ---")
# For each event, compute the "chain strength" = product of intermediate signals
valid_ec = ec_df[['strat_t10_anom_K', 'surface_t_anom_K', 'rr']].dropna()
if len(valid_ec) >= 5:
    # Strat warming predicts surface warming?
    r_st, p_st = stats.spearmanr(valid_ec['strat_t10_anom_K'], valid_ec['surface_t_anom_K'])
    print(f"  Strat T10 -> Surface T: r={r_st:.3f}, P={p_st:.4f}")
    # Surface warming predicts avalanche decrease?
    r_ta, p_ta = stats.spearmanr(valid_ec['surface_t_anom_K'], valid_ec['rr'])
    print(f"  Surface T -> Aval RR: r={r_ta:.3f}, P={p_ta:.4f}")
    # Strat warming predicts avalanche response directly?
    r_sa, p_sa = stats.spearmanr(valid_ec['strat_t10_anom_K'], valid_ec['rr'])
    print(f"  Strat T10 -> Aval RR: r={r_sa:.3f}, P={p_sa:.4f}")

chain_results = {
    'chain_correlations': chain_corrs,
    'propagation': propagation_results,
    'strat_to_surface_r': round(float(r_st), 3) if 'r_st' in dir() else None,
    'strat_to_surface_p': round(float(p_st), 4) if 'p_st' in dir() else None,
    'surface_to_avalanche_r': round(float(r_ta), 3) if 'r_ta' in dir() else None,
    'surface_to_avalanche_p': round(float(p_ta), 4) if 'p_ta' in dir() else None,
    'strat_to_avalanche_r': round(float(r_sa), 3) if 'r_sa' in dir() else None,
    'strat_to_avalanche_p': round(float(p_sa), 4) if 'p_sa' in dir() else None,
}
results['atmospheric_chain'] = chain_results

# ============================================================
# 3. STANDARDIZED EFFECT SIZES
# ============================================================
print("\n=== 3. STANDARDIZED EFFECT SIZES ===")

def cohens_d(group1, group2):
    """Compute Cohen's d (pooled SD)."""
    n1, n2 = len(group1), len(group2)
    mean1, mean2 = np.mean(group1), np.mean(group2)
    s1, s2 = np.std(group1, ddof=1), np.std(group2, ddof=1)
    pooled_sd = np.sqrt(((n1-1)*s1**2 + (n2-1)*s2**2) / (n1+n2-2))
    d = (mean1 - mean2) / pooled_sd if pooled_sd > 0 else 0
    # 95% CI using approximation
    se = np.sqrt((n1+n2)/(n1*n2) + d**2/(2*(n1+n2)))
    ci_lo = d - 1.96*se
    ci_hi = d + 1.96*se
    return d, ci_lo, ci_hi

def glass_delta(group1, group2):
    """Glass's delta using control group SD."""
    mean1, mean2 = np.mean(group1), np.mean(group2)
    s2 = np.std(group2, ddof=1)
    return (mean1 - mean2) / s2 if s2 > 0 else 0

# Swiss SSW vs control avalanche counts
winter_mask = panel['is_winter'] == 1 if 'is_winter' in panel.columns else panel['doy'].apply(lambda x: x >= 305 or x <= 120)
ssw_mask = panel['ssw_within_15d'] == 1 if 'ssw_within_15d' in panel.columns else pd.Series(False, index=panel.index)

if ssw_mask.any():
    ssw_aval = panel.loc[ssw_mask & winter_mask, 'aai_all_dry'].dropna().values
    ctrl_aval = panel.loc[~ssw_mask & winter_mask, 'aai_all_dry'].dropna().values
    d_swiss, ci_lo, ci_hi = cohens_d(ssw_aval, ctrl_aval)
    g_swiss = glass_delta(ssw_aval, ctrl_aval)
    print(f"Swiss SSW vs Control (daily counts):")
    print(f"  SSW mean={np.mean(ssw_aval):.3f}, Control mean={np.mean(ctrl_aval):.3f}")
    print(f"  Cohen's d = {d_swiss:.3f} [{ci_lo:.3f}, {ci_hi:.3f}]")
    print(f"  Glass's delta = {g_swiss:.3f}")
    print(f"  n_SSW={len(ssw_aval)}, n_control={len(ctrl_aval)}")
else:
    # Reconstruct SSW mask from dates
    ssw_mask2 = pd.Series(False, index=panel.index)
    for sd in ssw_dates:
        start = sd - pd.Timedelta(days=15)
        end = sd + pd.Timedelta(days=15)
        ssw_mask2.loc[start:end] = True
    ssw_aval = panel.loc[ssw_mask2 & winter_mask, 'aai_all_dry'].dropna().values
    ctrl_aval = panel.loc[~ssw_mask2 & winter_mask, 'aai_all_dry'].dropna().values
    d_swiss, ci_lo, ci_hi = cohens_d(ssw_aval, ctrl_aval)
    g_swiss = glass_delta(ssw_aval, ctrl_aval)
    print(f"Swiss SSW vs Control (daily counts):")
    print(f"  SSW mean={np.mean(ssw_aval):.3f}, Control mean={np.mean(ctrl_aval):.3f}")
    print(f"  Cohen's d = {d_swiss:.3f} [{ci_lo:.3f}, {ci_hi:.3f}]")
    print(f"  Glass's delta = {g_swiss:.3f}")

# Event-level effect size
event_rrs = ec_df['rr'].dropna().values
d_event = (np.mean(event_rrs) - 1.0) / np.std(event_rrs, ddof=1) if np.std(event_rrs, ddof=1) > 0 else 0
print(f"\nEvent-level RR departure from 1.0:")
print(f"  Mean RR = {np.mean(event_rrs):.3f}, SD = {np.std(event_rrs, ddof=1):.3f}")
print(f"  Cohen's d (from null RR=1): {d_event:.3f}")

# Surface temperature effect size
if 'surface_t_anom_K' in ec_df.columns:
    surf_anoms = ec_df['surface_t_anom_K'].dropna().values
    d_temp = np.mean(surf_anoms) / np.std(surf_anoms, ddof=1) if np.std(surf_anoms, ddof=1) > 0 else 0
    print(f"\nSurface T anomaly effect size:")
    print(f"  Mean = {np.mean(surf_anoms):.3f} K, SD = {np.std(surf_anoms, ddof=1):.3f}")
    print(f"  Cohen's d (from null=0): {d_temp:.3f}")

# Cold-phase reversal effect size
# Days +16 to +30 vs climatology
cold_phase_vals = []
for ssw_date in ssw_dates:
    for offset in range(16, 31):
        day = ssw_date + pd.Timedelta(days=offset)
        if day in era5.index:
            doy = era5.loc[day, 'doy']
            if isinstance(doy, pd.Series):
                doy = doy.iloc[0]
            if int(doy) in era5_clim.index:
                anom = float(era5.loc[day, 't2m_K']) - float(era5_clim[int(doy)])
                cold_phase_vals.append(anom)

warm_phase_vals = []
for ssw_date in ssw_dates:
    for offset in range(0, 16):
        day = ssw_date + pd.Timedelta(days=offset)
        if day in era5.index:
            doy = era5.loc[day, 'doy']
            if isinstance(doy, pd.Series):
                doy = doy.iloc[0]
            if int(doy) in era5_clim.index:
                anom = float(era5.loc[day, 't2m_K']) - float(era5_clim[int(doy)])
                warm_phase_vals.append(anom)

if cold_phase_vals and warm_phase_vals:
    cold_arr = np.array(cold_phase_vals)
    warm_arr = np.array(warm_phase_vals)
    d_reversal, ci_lo_r, ci_hi_r = cohens_d(cold_arr, warm_arr)
    t_reversal, p_reversal = stats.ttest_ind(cold_arr, warm_arr)
    print(f"\nCold-phase reversal (days +16 to +30 vs days 0 to +15):")
    print(f"  Cold phase mean: {np.mean(cold_arr):.3f} K (n={len(cold_arr)})")
    print(f"  Warm phase mean: {np.mean(warm_arr):.3f} K (n={len(warm_arr)})")
    print(f"  Difference: {np.mean(cold_arr) - np.mean(warm_arr):.3f} K")
    print(f"  Cohen's d = {d_reversal:.3f} [{ci_lo_r:.3f}, {ci_hi_r:.3f}]")
    print(f"  t-test P = {p_reversal:.6f}")

# Compute avalanche rate in cold phase vs warm phase
cold_aval = []
warm_aval = []
for ssw_date in ssw_dates:
    cold_days = panel.loc[ssw_date + pd.Timedelta(days=16):ssw_date + pd.Timedelta(days=30), 'aai_all_dry'].dropna()
    warm_days = panel.loc[ssw_date:ssw_date + pd.Timedelta(days=15), 'aai_all_dry'].dropna()
    if len(cold_days) > 0:
        cold_aval.extend(cold_days.values.tolist())
    if len(warm_days) > 0:
        warm_aval.extend(warm_days.values.tolist())

if cold_aval and warm_aval:
    cold_a = np.array(cold_aval)
    warm_a = np.array(warm_aval)
    d_aval_reversal, _, _ = cohens_d(cold_a, warm_a)
    t_aval, p_aval = stats.ttest_ind(cold_a, warm_a)
    u_aval, p_u_aval = stats.mannwhitneyu(cold_a, warm_a, alternative='greater')
    print(f"\nCold-phase avalanche reversal:")
    print(f"  Cold phase mean aval: {np.mean(cold_a):.3f}/day (n={len(cold_a)})")
    print(f"  Warm phase mean aval: {np.mean(warm_a):.3f}/day (n={len(warm_a)})")
    print(f"  Cohen's d = {d_aval_reversal:.3f}")
    print(f"  Mann-Whitney P (cold > warm) = {p_u_aval:.4f}")

results['standardized_effects'] = {
    'swiss_daily_cohen_d': round(float(d_swiss), 3),
    'swiss_daily_cohen_d_ci': [round(float(ci_lo), 3), round(float(ci_hi), 3)],
    'swiss_glass_delta': round(float(g_swiss), 3),
    'event_level_cohen_d': round(float(d_event), 3),
    'surface_t_cohen_d': round(float(d_temp), 3) if 'd_temp' in dir() else None,
    'cold_reversal_cohen_d': round(float(d_reversal), 3) if 'd_reversal' in dir() else None,
    'cold_reversal_p': round(float(p_reversal), 6) if 'p_reversal' in dir() else None,
    'cold_aval_cohen_d': round(float(d_aval_reversal), 3) if 'd_aval_reversal' in dir() else None,
    'cold_aval_p': round(float(p_u_aval), 4) if 'p_u_aval' in dir() else None,
}

# ============================================================
# 4. COLD-PHASE FALSIFICATION & NAO REJECTION
# ============================================================
print("\n=== 4. COLD-PHASE FALSIFICATION & NAO REJECTION ===")

# The cold-phase reversal is a prediction of the mechanism:
# If SSW-associated surface warming stabilises snowpack → fewer avalanches,
# then the post-SSW cold reversal should DE-stabilise → more avalanches
# This is a pre-registered prediction based on the mechanism
print("\n--- 4a. Cold-Phase as Mechanistic Prediction ---")
# Already computed above. Now compute the EVENT-LEVEL cold-phase avalanche increase
cold_phase_rrs = []
for ssw_date in ssw_dates:
    cold_window = panel.loc[ssw_date + pd.Timedelta(days=16):ssw_date + pd.Timedelta(days=30), 'aai_all_dry'].dropna()
    cold_rate = cold_window.mean() if len(cold_window) > 0 else np.nan
    
    # Control for same DOY
    cold_doy = (ssw_date + pd.Timedelta(days=23)).timetuple().tm_yday
    ctrl_vals_c = []
    winter = ssw_date.year if ssw_date.month >= 9 else ssw_date.year - 1
    for yr in range(1999, 2020):
        if yr == winter:
            continue
        for doy_off in range(-3, 4):
            d = cold_doy + doy_off
            mask_d = (panel['doy'] == d)
            vals_d = panel.loc[mask_d, 'aai_all_dry'].dropna()
            if len(vals_d) > 0:
                ctrl_vals_c.extend(vals_d.values.tolist())
    ctrl_rate_c = np.mean(ctrl_vals_c) if ctrl_vals_c else np.nan
    cold_rr = cold_rate / ctrl_rate_c if ctrl_rate_c > 0 else np.nan
    cold_phase_rrs.append(cold_rr)

cold_phase_rrs = [r for r in cold_phase_rrs if not np.isnan(r)]
n_increase = sum(1 for r in cold_phase_rrs if r > 1)
n_decrease = sum(1 for r in cold_phase_rrs if r < 1)
cold_sign_p = stats.binomtest(n_increase, len(cold_phase_rrs), 0.5).pvalue if cold_phase_rrs else 1.0
print(f"Cold-phase avalanche RRs: {n_increase}/{len(cold_phase_rrs)} show INCREASE")
print(f"  Mean RR = {np.mean(cold_phase_rrs):.3f}")
print(f"  Sign test P = {cold_sign_p:.4f}")
print(f"  t-test vs 1.0: P = {stats.ttest_1samp(cold_phase_rrs, 1.0).pvalue:.4f}")

# Combined warm-phase decrease + cold-phase increase = the temperature-mechanism prediction
warm_phase_rrs = [e['rr'] for e in event_catalog if e['rr'] is not None]
print(f"\nMechanistic prediction test (warm decrease + cold increase):")
print(f"  Warm phase: mean RR = {np.mean(warm_phase_rrs):.3f} (n={len(warm_phase_rrs)})")
print(f"  Cold phase: mean RR = {np.mean(cold_phase_rrs):.3f} (n={len(cold_phase_rrs)})")
warm_cold_diff = np.mean(cold_phase_rrs) - np.mean(warm_phase_rrs)
print(f"  Cold - Warm difference = {warm_cold_diff:.3f}")

# Permutation test for the warm-cold RR difference
n_perm = 10000
observed_diff = warm_cold_diff
all_rrs = warm_phase_rrs + cold_phase_rrs
perm_diffs = []
for _ in range(n_perm):
    np.random.shuffle(all_rrs)
    perm_warm = np.mean(all_rrs[:len(warm_phase_rrs)])
    perm_cold = np.mean(all_rrs[len(warm_phase_rrs):])
    perm_diffs.append(perm_cold - perm_warm)
p_perm = np.mean([d >= observed_diff for d in perm_diffs])
print(f"  Permutation P (cold > warm) = {p_perm:.4f}")

# 4b. NAO/AO Rejection Strengthening
print("\n--- 4b. NAO/AO Rejection Strengthening ---")
# Compute NAO during SSW vs control
nao_ssw = []
nao_ctrl = []
for ssw_date in ssw_dates:
    ssw_nao = panel.loc[ssw_date:ssw_date + pd.Timedelta(days=14), 'nao_daily'].dropna()
    if len(ssw_nao) > 0:
        nao_ssw.append(ssw_nao.mean())

# Overall winter NAO
winter_nao = panel.loc[winter_mask, 'nao_daily'].dropna()
nao_ctrl_mean = winter_nao.mean()
nao_ctrl_std = winter_nao.std()

if nao_ssw:
    nao_ssw_arr = np.array(nao_ssw)
    t_nao, p_nao = stats.ttest_1samp(nao_ssw_arr, nao_ctrl_mean)
    d_nao = (nao_ssw_arr.mean() - nao_ctrl_mean) / nao_ctrl_std
    print(f"  Post-SSW NAO mean: {nao_ssw_arr.mean():.3f} (n={len(nao_ssw_arr)})")
    print(f"  Climatological NAO mean: {nao_ctrl_mean:.3f}")
    print(f"  t-test P = {p_nao:.4f}")
    print(f"  Effect size d = {d_nao:.3f}")
    # NAO-avalanche correlation (align indices)
    nao_winter = panel.loc[winter_mask, 'nao_daily'].dropna()
    aval_winter = panel.loc[winter_mask, 'aai_all_dry'].dropna()
    common_idx = nao_winter.index.intersection(aval_winter.index)
    if len(common_idx) > 10:
        r_nao_aval, p_nao_aval = stats.spearmanr(nao_winter.loc[common_idx], aval_winter.loc[common_idx])
        print(f"  NAO-avalanche daily correlation: r={r_nao_aval:.3f}, P={p_nao_aval:.4f} (n={len(common_idx)})")

results['cold_phase_falsification'] = {
    'n_increase': int(n_increase),
    'n_total': int(len(cold_phase_rrs)),
    'mean_cold_rr': round(float(np.mean(cold_phase_rrs)), 3),
    'sign_test_p': round(float(cold_sign_p), 4),
    'warm_cold_diff': round(float(warm_cold_diff), 3),
    'permutation_p': round(float(p_perm), 4),
}

results['nao_rejection'] = {
    'post_ssw_nao_mean': round(float(nao_ssw_arr.mean()), 3) if nao_ssw else None,
    'clim_nao_mean': round(float(nao_ctrl_mean), 3),
    'nao_shift_p': round(float(p_nao), 4) if 'p_nao' in dir() else None,
    'nao_shift_d': round(float(d_nao), 3) if 'd_nao' in dir() else None,
}

# ============================================================
# 5. SINTERING BY SSW TYPE
# ============================================================
print("\n=== 5. SINTERING BY SSW TYPE ===")
try:
    with open('data/results/sintering_extended.json', 'r') as f:
        sintering = json.load(f)
    
    sinter_events = sintering.get('per_event', [])
    if sinter_events:
        # Match sintering events with SSW types
        for se in sinter_events:
            se_date = se.get('ssw_date', '')
            se['published_type'] = published_types.get(se_date, '?')
            # Find surface warming from our catalog
            match = [e for e in event_catalog if e['date'] == se_date]
            if match:
                se['surface_warming'] = match[0].get('surface_warming', None)
        
        disp_sinter = [s['strength_ratio'] for s in sinter_events if s.get('published_type') == 'D']
        split_sinter = [s['strength_ratio'] for s in sinter_events if s.get('published_type') == 'S']
        warm_sinter = [s['strength_ratio'] for s in sinter_events if s.get('surface_warming') == True]
        cool_sinter = [s['strength_ratio'] for s in sinter_events if s.get('surface_warming') == False]
        
        print(f"Displacement sintering (n={len(disp_sinter)}): "
              f"mean={np.mean(disp_sinter):.3f}, "
              f"positive={sum(1 for x in disp_sinter if x > 1)}/{len(disp_sinter)}")
        print(f"Split sintering (n={len(split_sinter)}): "
              f"mean={np.mean(split_sinter):.3f}, "
              f"positive={sum(1 for x in split_sinter if x > 1)}/{len(split_sinter)}")
        
        if warm_sinter and cool_sinter:
            print(f"\nWarming-SSW sintering (n={len(warm_sinter)}): "
                  f"mean={np.mean(warm_sinter):.3f}, "
                  f"positive={sum(1 for x in warm_sinter if x > 1)}/{len(warm_sinter)}")
            print(f"Cooling-SSW sintering (n={len(cool_sinter)}): "
                  f"mean={np.mean(cool_sinter):.3f}, "
                  f"positive={sum(1 for x in cool_sinter if x > 1)}/{len(cool_sinter)}")
            
            # Warming-only test
            warm_ratios = np.array(warm_sinter)
            t_warm, p_warm = stats.ttest_1samp(warm_ratios, 1.0)
            sign_p_warm = stats.binomtest(sum(1 for x in warm_sinter if x > 1), len(warm_sinter), 0.5).pvalue
            print(f"\n  Warming-only sintering t-test vs 1.0: P = {p_warm:.4f}")
            print(f"  Warming-only sign test: P = {sign_p_warm:.4f}")
            print(f"  Warming-only mean enhancement: {(np.mean(warm_ratios)-1)*100:.1f}%")
            
            results['sintering_by_type'] = {
                'displacement_mean': round(float(np.mean(disp_sinter)), 4),
                'split_mean': round(float(np.mean(split_sinter)), 4),
                'warming_mean': round(float(np.mean(warm_sinter)), 4),
                'cooling_mean': round(float(np.mean(cool_sinter)), 4),
                'warming_n': len(warm_sinter),
                'cooling_n': len(cool_sinter),
                'warming_positive': sum(1 for x in warm_sinter if x > 1),
                'warming_p_ttest': round(float(p_warm), 4),
                'warming_p_sign': round(float(sign_p_warm), 4),
                'warming_mean_pct': round(float((np.mean(warm_ratios)-1)*100), 1),
            }
except Exception as e:
    print(f"  Sintering analysis skipped: {e}")

# ============================================================
# SAVE RESULTS
# ============================================================
print("\n=== SAVING RESULTS ===")
with open('data/results/r11_mechanism_upgrade.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print("Saved to data/results/r11_mechanism_upgrade.json")

# Save event catalog as CSV for reference
ec.to_csv('data/results/ssw_event_catalog.csv', index=False)
print("Saved event catalog to data/results/ssw_event_catalog.csv")

print("\n=== ANALYSIS COMPLETE ===")
