"""
R14 Mechanism Breakthrough Analysis
====================================
Addresses ALL reviewer concerns about mechanism void:
1. ERA5 composite weather anomalies during SSW (T, precip, snowfall, wind)
2. Formal causal mediation analysis (SSW → Surface T → Avalanche)
3. Event-level dose-response (weather anomaly → avalanche reduction)
4. Stratospheric signal propagation (10hPa → 100hPa → surface)
5. Alternative mechanism quantification (wind vs T vs precip)
"""

import pandas as pd
import numpy as np
from scipy import stats
import json
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# =============================================================================
# DATA LOADING
# =============================================================================
print("=" * 70)
print("R14 MECHANISM BREAKTHROUGH ANALYSIS")
print("=" * 70)

panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
era5 = pd.read_parquet('data/processed/era5_swiss_alps_extended.parquet')
ssw_cat = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw_cat.index = ssw_cat.index.tz_localize(None)

# NCEP stratospheric data
ncep_strat = pd.read_parquet('data/processed/atmospheric/ncep_stratosphere.parquet')

# SSW events (major only, within ERA5 range)
ssw_dates = ssw_cat[ssw_cat['type'] == 'M'].index
ssw_dates = ssw_dates[(ssw_dates >= '1998-01-01') & (ssw_dates <= '2019-12-31')]
print(f"SSW events: {len(ssw_dates)}")

# SSW types
disp_dates_str = ['1998-12-15', '1999-02-26', '2001-02-11', '2004-01-05', '2006-01-21',
                  '2007-02-24', '2008-02-22', '2010-02-09', '2012-01-11', '2019-01-01']
split_dates_str = ['2001-12-30', '2002-02-17', '2003-01-18', '2009-01-24', '2013-01-07', '2018-02-12']

disp_dates = pd.to_datetime(disp_dates_str)
split_dates = pd.to_datetime(split_dates_str)

# Merge panel + ERA5
merged = panel.join(era5[['t2m_K', 'tp_mm', 'sf_mm', 'sd_m', 'u10', 'v10', 'wind_speed']], how='inner')
merged = merged[merged['is_winter'] == 1]  # winter only
print(f"Merged winter dataset: {len(merged)} days")

# Compute day-of-year climatology for anomalies
merged['doy'] = merged.index.dayofyear
for col in ['t2m_K', 'tp_mm', 'sf_mm', 'wind_speed', 'sd_m']:
    clim = merged.groupby('doy')[col].transform('mean')
    merged[f'{col}_anom'] = merged[col] - clim

# Avalanche anomaly
aval_clim = merged.groupby('doy')['aai_all_dry'].transform('mean')
merged['aval_anom'] = merged['aai_all_dry'] - aval_clim

results = {}

# =============================================================================
# PHASE 1: ERA5 COMPOSITE WEATHER ANOMALIES
# =============================================================================
print("\n" + "=" * 70)
print("PHASE 1: ERA5 COMPOSITE WEATHER ANOMALIES DURING SSW")
print("=" * 70)

def compute_composite(dates, merged_df, window_pre=20, window_post=30):
    """Compute composite weather anomalies around SSW events."""
    lags = range(-window_pre, window_post + 1)
    variables = ['t2m_K_anom', 'tp_mm_anom', 'sf_mm_anom', 'wind_speed_anom', 'sd_m_anom', 'aval_anom']
    
    composites = {var: {lag: [] for lag in lags} for var in variables}
    
    for ssw_date in dates:
        for lag in lags:
            target = ssw_date + pd.Timedelta(days=lag)
            if target in merged_df.index:
                for var in variables:
                    val = merged_df.loc[target, var]
                    if not np.isnan(val):
                        composites[var][lag].append(val)
    
    # Compute means and significance
    result = {}
    for var in variables:
        means = []
        p_values = []
        ci_lo = []
        ci_hi = []
        ns = []
        for lag in lags:
            vals = composites[var][lag]
            if len(vals) >= 3:
                m = np.mean(vals)
                se = np.std(vals, ddof=1) / np.sqrt(len(vals))
                t_stat, p_val = stats.ttest_1samp(vals, 0)
                means.append(m)
                p_values.append(p_val)
                ci_lo.append(m - 1.96 * se)
                ci_hi.append(m + 1.96 * se)
                ns.append(len(vals))
            else:
                means.append(np.nan)
                p_values.append(np.nan)
                ci_lo.append(np.nan)
                ci_hi.append(np.nan)
                ns.append(0)
        
        result[var] = {
            'lags': list(lags),
            'means': means,
            'p_values': p_values,
            'ci_lo': ci_lo,
            'ci_hi': ci_hi,
            'n': ns
        }
    return result

# All SSW events
comp_all = compute_composite(ssw_dates, merged)
comp_disp = compute_composite(disp_dates, merged)
comp_split = compute_composite(split_dates, merged)

# Print key results
print("\n--- ALL SSW EVENTS (n=%d) ---" % len(ssw_dates))
for var in ['t2m_K_anom', 'sf_mm_anom', 'tp_mm_anom', 'wind_speed_anom']:
    lags = comp_all[var]['lags']
    means = comp_all[var]['means']
    pvals = comp_all[var]['p_values']
    
    # Find mean over post-SSW window (0 to +15 days)
    post_idx = [i for i, l in enumerate(lags) if 0 <= l <= 15]
    post_mean = np.nanmean([means[i] for i in post_idx])
    post_sig = sum(1 for i in post_idx if pvals[i] < 0.05)
    
    pre_idx = [i for i, l in enumerate(lags) if -15 <= l < 0]
    pre_mean = np.nanmean([means[i] for i in pre_idx])
    
    print(f"  {var}: Pre[-15,0) mean={pre_mean:.4f}, Post[0,+15] mean={post_mean:.4f}, "
          f"sig days={post_sig}/{len(post_idx)}")

print("\n--- DISPLACEMENT SSW (n=%d) ---" % len(disp_dates))
for var in ['t2m_K_anom', 'sf_mm_anom', 'tp_mm_anom', 'wind_speed_anom']:
    lags = comp_disp[var]['lags']
    means = comp_disp[var]['means']
    pvals = comp_disp[var]['p_values']
    
    post_idx = [i for i, l in enumerate(lags) if 0 <= l <= 15]
    post_mean = np.nanmean([means[i] for i in post_idx])
    post_sig = sum(1 for i in post_idx if pvals[i] < 0.05)
    
    pre_idx = [i for i, l in enumerate(lags) if -15 <= l < 0]
    pre_mean = np.nanmean([means[i] for i in pre_idx])
    
    print(f"  {var}: Pre[-15,0) mean={pre_mean:.4f}, Post[0,+15] mean={post_mean:.4f}, "
          f"sig days={post_sig}/{len(post_idx)}")

# Mean post-SSW temperature anomaly for displacement
disp_t2m_post = []
for ssw_date in disp_dates:
    vals = []
    for lag in range(0, 16):
        target = ssw_date + pd.Timedelta(days=lag)
        if target in merged.index:
            vals.append(merged.loc[target, 't2m_K_anom'])
    if vals:
        disp_t2m_post.append(np.nanmean(vals))

t_val, p_val = stats.ttest_1samp(disp_t2m_post, 0)
print(f"\nDisplacement post-SSW T anomaly: mean={np.mean(disp_t2m_post):.3f} K, "
      f"t={t_val:.3f}, P={p_val:.4f}, n={len(disp_t2m_post)}")

results['phase1_composites'] = {
    'all_ssw': {
        'n_events': len(ssw_dates),
        't2m_post_mean': float(np.nanmean([comp_all['t2m_K_anom']['means'][i] 
                                            for i, l in enumerate(comp_all['t2m_K_anom']['lags']) if 0 <= l <= 15])),
        'sf_post_mean': float(np.nanmean([comp_all['sf_mm_anom']['means'][i] 
                                           for i, l in enumerate(comp_all['sf_mm_anom']['lags']) if 0 <= l <= 15])),
        'wind_post_mean': float(np.nanmean([comp_all['wind_speed_anom']['means'][i] 
                                             for i, l in enumerate(comp_all['wind_speed_anom']['lags']) if 0 <= l <= 15])),
    },
    'displacement': {
        'n_events': len(disp_dates),
        't2m_post_mean_K': float(np.mean(disp_t2m_post)),
        't2m_post_p': float(p_val),
    }
}

# =============================================================================
# PHASE 2: EVENT-LEVEL DOSE-RESPONSE (Weather → Avalanche)
# =============================================================================
print("\n" + "=" * 70)
print("PHASE 2: EVENT-LEVEL DOSE-RESPONSE")
print("=" * 70)

event_data = []
for ssw_date in ssw_dates:
    ssw_window = pd.date_range(ssw_date, ssw_date + pd.Timedelta(days=15))
    ctrl_start = ssw_date - pd.Timedelta(days=45)
    ctrl_end = ssw_date - pd.Timedelta(days=16)
    ctrl_window = pd.date_range(ctrl_start, ctrl_end)
    
    ssw_in = merged.index.isin(ssw_window)
    ctrl_in = merged.index.isin(ctrl_window)
    
    if ssw_in.sum() < 5 or ctrl_in.sum() < 10:
        continue
    
    ssw_data = merged.loc[ssw_in]
    ctrl_data = merged.loc[ctrl_in]
    
    event = {
        'ssw_date': str(ssw_date.date()),
        'type': 'displacement' if ssw_date in disp_dates else 'split',
        # Weather anomalies (SSW minus control)
        'delta_t2m': float(ssw_data['t2m_K'].mean() - ctrl_data['t2m_K'].mean()),
        'delta_sf': float(ssw_data['sf_mm'].mean() - ctrl_data['sf_mm'].mean()),
        'delta_tp': float(ssw_data['tp_mm'].mean() - ctrl_data['tp_mm'].mean()),
        'delta_wind': float(ssw_data['wind_speed'].mean() - ctrl_data['wind_speed'].mean()),
        'delta_sd': float(ssw_data['sd_m'].mean() - ctrl_data['sd_m'].mean()),
        # Avalanche anomaly
        'ssw_aval_rate': float(ssw_data['aai_all_dry'].mean()),
        'ctrl_aval_rate': float(ctrl_data['aai_all_dry'].mean()),
        'aval_rr': float(ssw_data['aai_all_dry'].mean() / max(ctrl_data['aai_all_dry'].mean(), 0.01)),
        'aval_diff': float(ssw_data['aai_all_dry'].mean() - ctrl_data['aai_all_dry'].mean()),
    }
    event_data.append(event)

edf = pd.DataFrame(event_data)
print(f"\nEvent-level dataset: {len(edf)} events")
print(edf[['ssw_date', 'type', 'delta_t2m', 'delta_sf', 'delta_wind', 'aval_diff']].to_string())

# Dose-response correlations
print("\n--- DOSE-RESPONSE CORRELATIONS (all events) ---")
dose_results = {}
for weather_var in ['delta_t2m', 'delta_sf', 'delta_tp', 'delta_wind']:
    r_s, p_s = stats.spearmanr(edf[weather_var], edf['aval_diff'])
    r_p, p_p = stats.pearsonr(edf[weather_var], edf['aval_diff'])
    print(f"  {weather_var} vs aval_diff: Spearman r={r_s:.3f} (P={p_s:.4f}), "
          f"Pearson r={r_p:.3f} (P={p_p:.4f})")
    dose_results[weather_var] = {
        'spearman_r': float(r_s), 'spearman_p': float(p_s),
        'pearson_r': float(r_p), 'pearson_p': float(p_p),
    }

# Displacement only
edf_disp = edf[edf['type'] == 'displacement']
print(f"\n--- DOSE-RESPONSE (displacement only, n={len(edf_disp)}) ---")
for weather_var in ['delta_t2m', 'delta_sf', 'delta_tp', 'delta_wind']:
    if len(edf_disp) >= 5:
        r_s, p_s = stats.spearmanr(edf_disp[weather_var], edf_disp['aval_diff'])
        print(f"  {weather_var} vs aval_diff: Spearman r={r_s:.3f} (P={p_s:.4f})")
        dose_results[f'disp_{weather_var}'] = {
            'spearman_r': float(r_s), 'spearman_p': float(p_s),
        }

results['phase2_dose_response'] = {
    'n_events': len(edf),
    'correlations': dose_results,
    'events': event_data,
}

# =============================================================================
# PHASE 3: FORMAL CAUSAL MEDIATION ANALYSIS
# =============================================================================
print("\n" + "=" * 70)
print("PHASE 3: FORMAL CAUSAL MEDIATION ANALYSIS")
print("=" * 70)

# Treatment: SSW window (1 = within 15 days of SSW, 0 = control winter day)
# Mediator: Surface temperature anomaly
# Outcome: Avalanche count (aai_all_dry)

# Create treatment variable
merged['ssw_treatment'] = 0
for ssw_date in ssw_dates:
    window = pd.date_range(ssw_date, ssw_date + pd.Timedelta(days=15))
    merged.loc[merged.index.isin(window), 'ssw_treatment'] = 1

# Displacement-only treatment
merged['disp_treatment'] = 0
for ssw_date in disp_dates:
    window = pd.date_range(ssw_date, ssw_date + pd.Timedelta(days=15))
    merged.loc[merged.index.isin(window), 'disp_treatment'] = 1

# Baron & Kenny mediation steps
# Step 1: Total effect (X → Y)
# Step 2: X → M (treatment → mediator)
# Step 3: X + M → Y (treatment + mediator → outcome)

from scipy.stats import mannwhitneyu

def mediation_analysis(treatment_col, mediator_col, outcome_col, df, label):
    """Baron & Kenny mediation with bootstrap Sobel test."""
    treat = df[df[treatment_col] == 1]
    ctrl = df[df[treatment_col] == 0]
    
    # Step 1: Total effect (X → Y)
    total_diff = treat[outcome_col].mean() - ctrl[outcome_col].mean()
    u_stat, p_total = mannwhitneyu(ctrl[outcome_col], treat[outcome_col], alternative='two-sided')
    
    # Step 2: X → M
    mediator_diff = treat[mediator_col].mean() - ctrl[mediator_col].mean()
    u_m, p_mediator = mannwhitneyu(ctrl[mediator_col], treat[mediator_col], alternative='two-sided')
    
    # Step 3: M → Y controlling for X (via regression)
    from numpy.linalg import lstsq
    X_mat = np.column_stack([
        df[treatment_col].values,
        df[mediator_col].values,
        np.ones(len(df))
    ])
    y_vec = df[outcome_col].values
    mask = ~np.isnan(X_mat).any(axis=1) & ~np.isnan(y_vec)
    X_clean = X_mat[mask]
    y_clean = y_vec[mask]
    
    beta, residuals, rank, sv = lstsq(X_clean, y_clean, rcond=None)
    
    # Indirect effect = a * b (product of coefficients)
    a = mediator_diff  # X → M path
    b = beta[1]  # M → Y path controlling for X
    indirect = a * b
    direct = beta[0]  # X → Y controlling for M
    
    # Residual computation for SE
    y_pred = X_clean @ beta
    resid = y_clean - y_pred
    n = len(y_clean)
    p_params = X_clean.shape[1]
    mse = np.sum(resid**2) / (n - p_params)
    XtX_inv = np.linalg.inv(X_clean.T @ X_clean)
    se_beta = np.sqrt(mse * np.diag(XtX_inv))
    
    se_a = np.std(treat[mediator_col].dropna()) / np.sqrt(len(treat))  # approximate
    se_b = se_beta[1]
    
    # Sobel test
    sobel_se = np.sqrt(a**2 * se_b**2 + b**2 * se_a**2)
    sobel_z = indirect / sobel_se if sobel_se > 0 else 0
    sobel_p = 2 * (1 - stats.norm.cdf(abs(sobel_z)))
    
    # Bootstrap mediation CI
    n_boot = 5000
    boot_indirect = []
    for _ in range(n_boot):
        idx = np.random.choice(len(df), len(df), replace=True)
        bdf = df.iloc[idx]
        bt = bdf[bdf[treatment_col] == 1]
        bc = bdf[bdf[treatment_col] == 0]
        if len(bt) < 5 or len(bc) < 5:
            continue
        ba = bt[mediator_col].mean() - bc[mediator_col].mean()
        
        bX = np.column_stack([bdf[treatment_col].values, bdf[mediator_col].values, np.ones(len(bdf))])
        by = bdf[outcome_col].values
        bmask = ~np.isnan(bX).any(axis=1) & ~np.isnan(by)
        if bmask.sum() < 10:
            continue
        try:
            bbeta, _, _, _ = lstsq(bX[bmask], by[bmask], rcond=None)
            boot_indirect.append(ba * bbeta[1])
        except:
            continue
    
    boot_ci = np.percentile(boot_indirect, [2.5, 97.5]) if len(boot_indirect) > 100 else [np.nan, np.nan]
    
    # Proportion mediated
    prop_mediated = indirect / total_diff if abs(total_diff) > 0.001 else np.nan
    
    print(f"\n--- {label} ---")
    print(f"  Total effect (X→Y): {total_diff:.4f}, P={p_total:.6f}")
    print(f"  Path a (X→M): {mediator_diff:.4f}, P={p_mediator:.6f}")
    print(f"  Path b (M→Y|X): {b:.4f}, SE={se_b:.4f}")
    print(f"  Direct effect (c'): {direct:.4f}, SE={se_beta[0]:.4f}")
    print(f"  Indirect effect (a*b): {indirect:.4f}")
    print(f"  Sobel test: z={sobel_z:.3f}, P={sobel_p:.6f}")
    print(f"  Bootstrap 95% CI: [{boot_ci[0]:.4f}, {boot_ci[1]:.4f}]")
    print(f"  Proportion mediated: {prop_mediated:.1%}" if not np.isnan(prop_mediated) else "  Proportion mediated: N/A")
    
    return {
        'total_effect': float(total_diff),
        'total_p': float(p_total),
        'path_a': float(mediator_diff),
        'path_a_p': float(p_mediator),
        'path_b': float(b),
        'path_b_se': float(se_b),
        'direct_effect': float(direct),
        'indirect_effect': float(indirect),
        'sobel_z': float(sobel_z),
        'sobel_p': float(sobel_p),
        'bootstrap_ci': [float(boot_ci[0]), float(boot_ci[1])],
        'proportion_mediated': float(prop_mediated) if not np.isnan(prop_mediated) else None,
        'n_treat': int(len(treat)),
        'n_ctrl': int(len(ctrl)),
    }

# Clean dataset for mediation
med_df = merged[['ssw_treatment', 'disp_treatment', 't2m_K_anom', 'aai_all_dry', 
                  'sf_mm_anom', 'wind_speed_anom']].dropna()

# Mediation 1: All SSW → Surface T → Avalanche
med1 = mediation_analysis('ssw_treatment', 't2m_K_anom', 'aai_all_dry', med_df, 
                           'All SSW → Surface T → Avalanche')

# Mediation 2: Displacement SSW → Surface T → Avalanche
med2 = mediation_analysis('disp_treatment', 't2m_K_anom', 'aai_all_dry', med_df,
                           'Displacement SSW → Surface T → Avalanche')

# Mediation 3: All SSW → Snowfall → Avalanche
med3 = mediation_analysis('ssw_treatment', 'sf_mm_anom', 'aai_all_dry', med_df,
                           'All SSW → Snowfall → Avalanche')

# Mediation 4: All SSW → Wind → Avalanche
med4 = mediation_analysis('ssw_treatment', 'wind_speed_anom', 'aai_all_dry', med_df,
                           'All SSW → Wind → Avalanche')

results['phase3_mediation'] = {
    'all_ssw_temp': med1,
    'disp_ssw_temp': med2,
    'all_ssw_snowfall': med3,
    'all_ssw_wind': med4,
}

# =============================================================================
# PHASE 4: STRATOSPHERIC SIGNAL PROPAGATION
# =============================================================================
print("\n" + "=" * 70)
print("PHASE 4: STRATOSPHERIC SIGNAL PROPAGATION")
print("=" * 70)

print("NCEP strat columns:", list(ncep_strat.columns))
print("Shape:", ncep_strat.shape)
print("Date range:", ncep_strat.index.min(), "-", ncep_strat.index.max())

# Compute anomalies for NCEP levels
ncep_levels = ['ncep_t_10hpa', 'ncep_t_20hpa', 'ncep_t_30hpa', 'ncep_t_50hpa', 
               'ncep_t_70hpa', 'ncep_t_100hpa']
available_levels = [c for c in ncep_levels if c in ncep_strat.columns]

if not available_levels:
    # Try panel data
    available_levels = [c for c in ncep_levels if c in panel.columns]
    strat_df = panel[available_levels].copy() if available_levels else None
    print(f"Using panel data for NCEP levels: {available_levels}")
else:
    strat_df = ncep_strat[available_levels].copy()
    print(f"Using NCEP strat data for levels: {available_levels}")

if strat_df is not None:
    strat_df['doy'] = strat_df.index.dayofyear
    for col in available_levels:
        clim = strat_df.groupby('doy')[col].transform('mean')
        strat_df[f'{col}_anom'] = strat_df[col] - clim
    
    # Add surface temperature
    strat_df = strat_df.join(era5[['t2m_K']], how='inner')
    strat_df['t2m_anom'] = strat_df['t2m_K'] - strat_df.groupby('doy')['t2m_K'].transform('mean')
    
    # Composite by lag
    levels_to_plot = [f'{c}_anom' for c in available_levels] + ['t2m_anom']
    lags = range(-20, 31)
    
    propagation = {}
    for level in levels_to_plot:
        lag_means = []
        lag_pvals = []
        for lag in lags:
            vals = []
            for ssw_date in disp_dates:  # displacement only for cleaner signal
                target = ssw_date + pd.Timedelta(days=lag)
                if target in strat_df.index:
                    v = strat_df.loc[target, level]
                    if not np.isnan(v):
                        vals.append(v)
            if len(vals) >= 3:
                lag_means.append(float(np.mean(vals)))
                _, p = stats.ttest_1samp(vals, 0)
                lag_pvals.append(float(p))
            else:
                lag_means.append(np.nan)
                lag_pvals.append(np.nan)
        propagation[level] = {'means': lag_means, 'p_values': lag_pvals}
    
    # Find peak timing at each level
    print("\n--- SIGNAL PROPAGATION TIMING (displacement events) ---")
    for level in levels_to_plot:
        means = propagation[level]['means']
        if any(not np.isnan(m) for m in means):
            valid = [(i, m) for i, m in enumerate(means) if not np.isnan(m)]
            peak_idx, peak_val = max(valid, key=lambda x: abs(x[1]))
            peak_lag = list(lags)[peak_idx]
            peak_p = propagation[level]['p_values'][peak_idx]
            print(f"  {level}: peak anomaly={peak_val:.2f} at lag {peak_lag:+d}d (P={peak_p:.4f})")
    
    results['phase4_propagation'] = {
        'lags': list(lags),
        'levels': {k: v for k, v in propagation.items()},
        'n_events': len(disp_dates),
    }

# =============================================================================
# PHASE 5: ALTERNATIVE MECHANISM QUANTIFICATION
# =============================================================================
print("\n" + "=" * 70)
print("PHASE 5: ALTERNATIVE MECHANISM QUANTIFICATION")
print("=" * 70)

# For each SSW event, quantify: Temperature, Snowfall, Wind, Snow Depth changes
print("\n--- WEATHER ANOMALIES BY SSW TYPE (0 to +15 days vs control) ---")

for label, dates in [('Displacement', disp_dates), ('Split', split_dates), ('All', ssw_dates)]:
    t2m_deltas = []
    sf_deltas = []
    wind_deltas = []
    sd_deltas = []
    
    for ssw_date in dates:
        ssw_win = pd.date_range(ssw_date, ssw_date + pd.Timedelta(days=15))
        ctrl_win = pd.date_range(ssw_date - pd.Timedelta(days=45), ssw_date - pd.Timedelta(days=16))
        
        ssw_in = merged.index.isin(ssw_win)
        ctrl_in = merged.index.isin(ctrl_win)
        
        if ssw_in.sum() < 5 or ctrl_in.sum() < 10:
            continue
        
        t2m_deltas.append(merged.loc[ssw_in, 't2m_K'].mean() - merged.loc[ctrl_in, 't2m_K'].mean())
        sf_deltas.append(merged.loc[ssw_in, 'sf_mm'].mean() - merged.loc[ctrl_in, 'sf_mm'].mean())
        wind_deltas.append(merged.loc[ssw_in, 'wind_speed'].mean() - merged.loc[ctrl_in, 'wind_speed'].mean())
        sd_deltas.append(merged.loc[ssw_in, 'sd_m'].mean() - merged.loc[ctrl_in, 'sd_m'].mean())
    
    print(f"\n  {label} (n={len(t2m_deltas)}):")
    
    for name, vals in [('Temperature (K)', t2m_deltas), ('Snowfall (mm/d)', sf_deltas),
                       ('Wind speed (m/s)', wind_deltas), ('Snow depth (m)', sd_deltas)]:
        arr = np.array(vals)
        if len(arr) >= 3:
            m = np.mean(arr)
            t, p = stats.ttest_1samp(arr, 0)
            pos = sum(1 for v in arr if v > 0)
            print(f"    {name}: mean={m:+.3f}, t={t:.2f}, P={p:.4f}, positive={pos}/{len(arr)}")

# Effect size comparison: which weather variable best predicts avalanche reduction?
print("\n--- MULTIPLE REGRESSION: Weather → Avalanche ---")
from numpy.linalg import lstsq

# Standardize predictors for comparable effect sizes
X_vars = ['t2m_K_anom', 'sf_mm_anom', 'wind_speed_anom']
y_var = 'aai_all_dry'

reg_df = merged[X_vars + [y_var, 'ssw_treatment']].dropna()
# Only SSW-affected days for this analysis
ssw_df = reg_df[reg_df['ssw_treatment'] == 1]

if len(ssw_df) > 20:
    X_std = (ssw_df[X_vars] - ssw_df[X_vars].mean()) / ssw_df[X_vars].std()
    X_mat = np.column_stack([X_std.values, np.ones(len(X_std))])
    y = ssw_df[y_var].values
    
    beta, _, _, _ = lstsq(X_mat, y, rcond=None)
    
    # Compute R-squared for each predictor
    y_pred = X_mat @ beta
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - y.mean())**2)
    r2 = 1 - ss_res / ss_tot
    
    print(f"\n  During SSW windows (n={len(ssw_df)}):")
    print(f"  Multiple R² = {r2:.4f}")
    for i, var in enumerate(X_vars):
        print(f"    β({var}) = {beta[i]:.4f} (standardized)")
    
    results['phase5_weather_regression'] = {
        'n': int(len(ssw_df)),
        'r_squared': float(r2),
        'standardized_betas': {var: float(beta[i]) for i, var in enumerate(X_vars)},
    }

# =============================================================================
# PHASE 6: ENHANCED DOSE-RESPONSE — Temperature thresholds
# =============================================================================
print("\n" + "=" * 70)
print("PHASE 6: TEMPERATURE THRESHOLD ANALYSIS")
print("=" * 70)

# Key question: is there a threshold of surface warming above which avalanche
# reduction becomes significant?

# Split events into "warmed" vs "cooled" based on surface T change
warmed_events = edf[edf['delta_t2m'] > 0]
cooled_events = edf[edf['delta_t2m'] <= 0]

print(f"\nWarmed events (n={len(warmed_events)}): mean aval_diff = {warmed_events['aval_diff'].mean():.3f}")
print(f"Cooled events (n={len(cooled_events)}): mean aval_diff = {cooled_events['aval_diff'].mean():.3f}")

if len(warmed_events) >= 3 and len(cooled_events) >= 3:
    u, p = mannwhitneyu(warmed_events['aval_diff'], cooled_events['aval_diff'], alternative='two-sided')
    print(f"Mann-Whitney warmed vs cooled: U={u:.1f}, P={p:.4f}")
    
    t, p_t = stats.ttest_ind(warmed_events['aval_diff'], cooled_events['aval_diff'])
    print(f"t-test warmed vs cooled: t={t:.3f}, P={p_t:.4f}")
    
    # Cohen's d
    pooled_std = np.sqrt(((len(warmed_events)-1)*warmed_events['aval_diff'].std()**2 + 
                           (len(cooled_events)-1)*cooled_events['aval_diff'].std()**2) / 
                          (len(warmed_events) + len(cooled_events) - 2))
    d = (warmed_events['aval_diff'].mean() - cooled_events['aval_diff'].mean()) / pooled_std
    print(f"Cohen's d = {d:.3f}")
    
    results['phase6_threshold'] = {
        'n_warmed': int(len(warmed_events)),
        'n_cooled': int(len(cooled_events)),
        'warmed_aval_diff': float(warmed_events['aval_diff'].mean()),
        'cooled_aval_diff': float(cooled_events['aval_diff'].mean()),
        'mw_p': float(p),
        'ttest_p': float(p_t),
        'cohens_d': float(d),
    }

# =============================================================================
# PHASE 7: CROSS-COUNTRY WEATHER ANALYSIS  
# =============================================================================
print("\n" + "=" * 70)
print("PHASE 7: NCEP HEMISPHERIC WEATHER DURING SSW")
print("=" * 70)

# NH-mean z500, slp, u850 during SSW — these are in the panel
nh_vars = ['ncep_z500_nh', 'ncep_slp_nh', 'ncep_u850_nh']
available_nh = [v for v in nh_vars if v in panel.columns]

for var in available_nh:
    panel[f'{var}_doy'] = panel.index.dayofyear
    clim = panel.groupby(f'{var}_doy')[var].transform('mean')
    panel[f'{var}_anom'] = panel[var] - clim

    # Composite for displacement SSW
    post_vals = []
    for ssw_date in disp_dates:
        for lag in range(0, 16):
            target = ssw_date + pd.Timedelta(days=lag)
            if target in panel.index:
                v = panel.loc[target, f'{var}_anom']
                if not np.isnan(v):
                    post_vals.append(v)
    
    if post_vals:
        m = np.mean(post_vals)
        t, p = stats.ttest_1samp(post_vals, 0)
        print(f"  {var} (displacement, post-SSW): mean_anom={m:.3f}, t={t:.2f}, P={p:.4f}")

# =============================================================================
# SAVE RESULTS
# =============================================================================
print("\n" + "=" * 70)
print("SAVING RESULTS")
print("=" * 70)

with open('data/results/r14_mechanism_breakthrough.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)

print("Results saved to data/results/r14_mechanism_breakthrough.json")
print("\nDONE.")
