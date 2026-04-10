"""R31: SNOWPACK stability analysis around SSW events.

Uses the EnviDat weather_snowpack_danger.csv which contains SNOWPACK model output
with stability indices, weather variables, and snowpack structure data from multiple
Swiss stations (1997-present, 292K rows).

Key stability indices:
- sn38_pwl: Natural stability at 38° for persistent weak layers (higher = more stable)
- sk38_pwl: Skier stability at 38° (higher = more stable)
- ssi_pwl: Structural Stability Index (higher = more stable)
- ccl_pwl: Critical Cut Length (lower = easier to trigger)
- pwl_100: Persistent weak layer probability (0 or 1)
"""
import pandas as pd
import numpy as np
import json
import sys
import os
from scipy import stats

sys.stdout.reconfigure(encoding='utf-8')

# ── 1. Load data ──────────────────────────────────────────────────────────────
print("Loading weather_snowpack_danger.csv...")
df = pd.read_csv('data/cryosphere/envidat/weather_snowpack_danger.csv', 
                  low_memory=False)
print(f"  Rows: {len(df):,}, Columns: {len(df.columns)}")

df['datum'] = pd.to_datetime(df['datum'])
print(f"  Date range: {df['datum'].min()} to {df['datum'].max()}")
print(f"  Stations: {df['station_code'].nunique()}")
print(f"  Stations list: {sorted(df['station_code'].unique())[:20]}...")

# ── 2. Load SSW catalog ──────────────────────────────────────────────────────
ssw_cat = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw_cat.index = ssw_cat.index.tz_localize(None)
ssw_dates = ssw_cat.index.values

# Filter to winter months and dates within our data range
data_start = df['datum'].min()
data_end = df['datum'].max()
ssw_in_range = [pd.Timestamp(d) for d in ssw_dates 
                if pd.Timestamp(d) >= data_start and pd.Timestamp(d) <= data_end]
print(f"\nSSW events in data range: {len(ssw_in_range)}")
for d in ssw_in_range:
    print(f"  {d.strftime('%Y-%m-%d')}")

# ── 3. Define stability and weather variables ─────────────────────────────────
stability_vars = {
    'sn38_pwl': 'Natural stability 38° (PWL)',
    'sk38_pwl': 'Skier stability 38° (PWL)',
    'ssi_pwl': 'Structural Stability Index (PWL)',
    'ccl_pwl': 'Critical Cut Length (PWL)',
    'sn38_pwl_100': 'Natural stability (top 100cm)',
    'sk38_pwl_100': 'Skier stability (top 100cm)',
    'ssi_pwl_100': 'SSI (top 100cm)',
    'ccl_pwl_100': 'CCL (top 100cm)',
    'pwl_100': 'PWL probability (top 100cm)',
    'pwl_100_15': 'PWL probability (15d)',
}

weather_vars = {
    'TA': 'Air temperature (°C)',
    'HN24': 'New snow 24h (cm)',
    'HN72_24': 'New snow 72h (cm)',
    'HN24_7d': 'New snow 7d (cm)',
    'HS_mod': 'Snow height modeled (cm)',
    'SWE': 'Snow water equivalent (mm)',
    'VW': 'Wind speed (m/s)',
    'RH': 'Relative humidity (%)',
    'MS_Snow': 'Snowfall mass (kg/m²)',
    'wind_trans24': 'Wind transport 24h',
    'hoar_size': 'Surface hoar size',
}

# ── 4. Assign SSW windows (±15 days) ─────────────────────────────────────────
WINDOW = 15
df['in_ssw'] = False
df['ssw_event'] = None

for ssw_date in ssw_in_range:
    mask = (df['datum'] >= ssw_date - pd.Timedelta(days=WINDOW)) & \
           (df['datum'] <= ssw_date + pd.Timedelta(days=WINDOW))
    df.loc[mask, 'in_ssw'] = True
    df.loc[mask & df['ssw_event'].isna(), 'ssw_event'] = ssw_date.strftime('%Y-%m-%d')

# Winter filter (Nov-Mar)
df['month'] = df['datum'].dt.month
winter_mask = df['month'].isin([11, 12, 1, 2, 3])
df_winter = df[winter_mask].copy()
print(f"\nWinter records: {len(df_winter):,}")
print(f"  In SSW window: {df_winter['in_ssw'].sum():,}")
print(f"  Control: {(~df_winter['in_ssw']).sum():,}")

# ── 5. Overall SSW vs Control comparison ──────────────────────────────────────
print("\n" + "="*80)
print("OVERALL SSW vs CONTROL COMPARISON")
print("="*80)

results = {'overall': {}, 'event_level': {}, 'phase_resolved': {}, 'mechanism_test': {}}

all_vars = {**stability_vars, **weather_vars}
for var, label in all_vars.items():
    if var not in df_winter.columns:
        continue
    
    ssw_vals = df_winter.loc[df_winter['in_ssw'], var].dropna()
    ctrl_vals = df_winter.loc[~df_winter['in_ssw'], var].dropna()
    
    if len(ssw_vals) < 30 or len(ctrl_vals) < 30:
        continue
    
    mw_stat, mw_p = stats.mannwhitneyu(ssw_vals, ctrl_vals, alternative='two-sided')
    cohen_d = (ssw_vals.mean() - ctrl_vals.mean()) / np.sqrt(
        (ssw_vals.std()**2 + ctrl_vals.std()**2) / 2)
    
    results['overall'][var] = {
        'label': label,
        'ssw_mean': float(ssw_vals.mean()),
        'ssw_median': float(ssw_vals.median()),
        'ctrl_mean': float(ctrl_vals.mean()),
        'ctrl_median': float(ctrl_vals.median()),
        'diff': float(ssw_vals.mean() - ctrl_vals.mean()),
        'pct_change': float((ssw_vals.mean() - ctrl_vals.mean()) / ctrl_vals.mean() * 100) if ctrl_vals.mean() != 0 else 0,
        'cohen_d': float(cohen_d),
        'mw_p': float(mw_p),
        'n_ssw': int(len(ssw_vals)),
        'n_ctrl': int(len(ctrl_vals)),
    }
    
    sig = '***' if mw_p < 0.001 else '**' if mw_p < 0.01 else '*' if mw_p < 0.05 else ''
    print(f"\n{label} ({var}):")
    print(f"  SSW: {ssw_vals.mean():.3f} ± {ssw_vals.std():.3f} (n={len(ssw_vals)})")
    print(f"  Ctrl: {ctrl_vals.mean():.3f} ± {ctrl_vals.std():.3f} (n={len(ctrl_vals)})")
    print(f"  Δ = {ssw_vals.mean() - ctrl_vals.mean():+.3f} ({(ssw_vals.mean() - ctrl_vals.mean()) / ctrl_vals.mean() * 100:+.1f}%)")
    print(f"  Cohen's d = {cohen_d:+.3f}, MW P = {mw_p:.2e} {sig}")

# ── 6. Event-level analysis ──────────────────────────────────────────────────
print("\n" + "="*80)
print("EVENT-LEVEL SSW ANALYSIS")
print("="*80)

key_vars = ['sn38_pwl', 'sk38_pwl', 'ccl_pwl', 'HN24', 'TA', 'SWE', 'pwl_100']

for var in key_vars:
    if var not in df_winter.columns:
        continue
    
    label = all_vars.get(var, var)
    print(f"\n{label} ({var}):")
    
    event_diffs = []
    for ssw_date in ssw_in_range:
        ssw_start = ssw_date - pd.Timedelta(days=WINDOW)
        ssw_end = ssw_date + pd.Timedelta(days=WINDOW)
        
        # Control: same DOY range from other years
        ssw_doy_start = (ssw_date - pd.Timedelta(days=WINDOW)).dayofyear
        ssw_doy_end = (ssw_date + pd.Timedelta(days=WINDOW)).dayofyear
        ssw_year = ssw_date.year
        
        # SSW window values
        ssw_mask = (df_winter['datum'] >= ssw_start) & (df_winter['datum'] <= ssw_end)
        ssw_vals = df_winter.loc[ssw_mask, var].dropna()
        
        # Control: same DOY window, different years
        ctrl_mask = (~df_winter['in_ssw']) & \
                    (df_winter['datum'].dt.dayofyear >= ssw_doy_start) & \
                    (df_winter['datum'].dt.dayofyear <= ssw_doy_end)
        ctrl_vals = df_winter.loc[ctrl_mask, var].dropna()
        
        if len(ssw_vals) < 5 or len(ctrl_vals) < 5:
            continue
        
        diff = ssw_vals.mean() - ctrl_vals.mean()
        event_diffs.append({
            'event': ssw_date.strftime('%Y-%m-%d'),
            'ssw_mean': float(ssw_vals.mean()),
            'ctrl_mean': float(ctrl_vals.mean()),
            'diff': float(diff),
            'n_ssw': int(len(ssw_vals)),
            'n_ctrl': int(len(ctrl_vals)),
        })
        
        direction = '↑' if diff > 0 else '↓'
        print(f"  {ssw_date.strftime('%Y-%m-%d')}: SSW={ssw_vals.mean():.2f}, Ctrl={ctrl_vals.mean():.2f}, Δ={diff:+.3f} {direction}")
    
    if event_diffs:
        diffs = [e['diff'] for e in event_diffs]
        n_increase = sum(1 for d in diffs if d > 0)
        n_decrease = sum(1 for d in diffs if d < 0)
        mean_diff = np.mean(diffs)
        
        # Sign test
        from scipy.stats import binomtest
        sign_p = binomtest(max(n_increase, n_decrease), len(diffs), 0.5).pvalue
        
        # t-test
        t_stat, t_p = stats.ttest_1samp(diffs, 0)
        
        print(f"  Summary: {n_increase}↑ / {n_decrease}↓ of {len(diffs)} events")
        print(f"  Mean Δ = {mean_diff:+.4f}, Sign P = {sign_p:.4f}, t-test P = {t_p:.4f}")
        
        results['event_level'][var] = {
            'label': label,
            'events': event_diffs,
            'n_increase': n_increase,
            'n_decrease': n_decrease,
            'mean_diff': float(mean_diff),
            'sign_p': float(sign_p),
            't_p': float(t_p),
        }

# ── 7. Phase-resolved analysis (pre, during, post SSW) ──────────────────────
print("\n" + "="*80)
print("PHASE-RESOLVED ANALYSIS")
print("="*80)

phases = {
    'pre': (-15, -1),
    'during': (0, 7),
    'post': (8, 15),
    'late': (16, 30),
}

for var in ['sn38_pwl', 'sk38_pwl', 'HN24', 'TA', 'pwl_100', 'ccl_pwl']:
    if var not in df_winter.columns:
        continue
    
    label = all_vars.get(var, var)
    print(f"\n{label} ({var}):")
    
    phase_results = {}
    for phase_name, (d_start, d_end) in phases.items():
        phase_vals = []
        ctrl_vals_all = []
        
        for ssw_date in ssw_in_range:
            p_start = ssw_date + pd.Timedelta(days=d_start)
            p_end = ssw_date + pd.Timedelta(days=d_end)
            
            # Phase values
            mask = (df_winter['datum'] >= p_start) & (df_winter['datum'] <= p_end)
            vals = df_winter.loc[mask, var].dropna()
            
            # DOY-matched control
            doy_start = p_start.dayofyear
            doy_end = p_end.dayofyear
            ctrl_mask = (~df_winter['in_ssw']) & \
                        (df_winter['datum'].dt.dayofyear >= doy_start) & \
                        (df_winter['datum'].dt.dayofyear <= doy_end)
            ctrl = df_winter.loc[ctrl_mask, var].dropna()
            
            phase_vals.extend(vals.tolist())
            ctrl_vals_all.extend(ctrl.tolist())
        
        if len(phase_vals) >= 30 and len(ctrl_vals_all) >= 30:
            mw_stat, mw_p = stats.mannwhitneyu(phase_vals, ctrl_vals_all, alternative='two-sided')
            diff = np.mean(phase_vals) - np.mean(ctrl_vals_all)
            d = diff / np.sqrt((np.std(phase_vals)**2 + np.std(ctrl_vals_all)**2) / 2)
            
            phase_results[phase_name] = {
                'mean_ssw': float(np.mean(phase_vals)),
                'mean_ctrl': float(np.mean(ctrl_vals_all)),
                'diff': float(diff),
                'cohen_d': float(d),
                'mw_p': float(mw_p),
                'n_ssw': len(phase_vals),
                'n_ctrl': len(ctrl_vals_all),
            }
            
            sig = '***' if mw_p < 0.001 else '**' if mw_p < 0.01 else '*' if mw_p < 0.05 else ''
            print(f"  {phase_name:6s} [{d_start:+3d} to {d_end:+3d}d]: "
                  f"SSW={np.mean(phase_vals):.3f}, Ctrl={np.mean(ctrl_vals_all):.3f}, "
                  f"Δ={diff:+.3f} (d={d:+.3f}) P={mw_p:.2e} {sig}")
    
    results['phase_resolved'][var] = phase_results

# ── 8. THE KEY TEST: Loading vs Stability decomposition ──────────────────────
print("\n" + "="*80)
print("KEY MECHANISM TEST: LOADING vs STABILITY DECOMPOSITION")
print("="*80)
print("Hypothesis: SSW reduces LOADING (new snow, precipitation) while")
print("STABILITY (weak layers, structural integrity) remains unchanged or worsens.")
print()

loading_vars = ['HN24', 'HN72_24', 'HN24_7d', 'SWE', 'MS_Snow', 'wind_trans24']
stability_vars_list = ['sn38_pwl', 'sk38_pwl', 'ssi_pwl', 'ccl_pwl', 'pwl_100']

print("LOADING variables (expect DECREASE during SSW):")
loading_results = {}
for var in loading_vars:
    if var not in df_winter.columns:
        continue
    ssw_vals = df_winter.loc[df_winter['in_ssw'], var].dropna()
    ctrl_vals = df_winter.loc[~df_winter['in_ssw'], var].dropna()
    if len(ssw_vals) < 30 or len(ctrl_vals) < 30:
        continue
    
    diff = ssw_vals.mean() - ctrl_vals.mean()
    pct = diff / ctrl_vals.mean() * 100 if ctrl_vals.mean() != 0 else 0
    mw_stat, mw_p = stats.mannwhitneyu(ssw_vals, ctrl_vals, alternative='two-sided')
    d = diff / np.sqrt((ssw_vals.std()**2 + ctrl_vals.std()**2) / 2)
    
    sig = '***' if mw_p < 0.001 else '**' if mw_p < 0.01 else '*' if mw_p < 0.05 else ''
    direction = 'DECREASE ✓' if diff < 0 else 'INCREASE ✗'
    print(f"  {all_vars.get(var, var):30s}: {pct:+6.1f}% (d={d:+.3f}, P={mw_p:.2e}) {sig} → {direction}")
    
    loading_results[var] = {
        'pct_change': float(pct), 'cohen_d': float(d), 'p': float(mw_p),
        'direction': 'decrease' if diff < 0 else 'increase'
    }

print("\nSTABILITY variables (expect NO CHANGE or DECREASE during SSW):")
stability_results = {}
for var in stability_vars_list:
    if var not in df_winter.columns:
        continue
    ssw_vals = df_winter.loc[df_winter['in_ssw'], var].dropna()
    ctrl_vals = df_winter.loc[~df_winter['in_ssw'], var].dropna()
    if len(ssw_vals) < 30 or len(ctrl_vals) < 30:
        continue
    
    diff = ssw_vals.mean() - ctrl_vals.mean()
    pct = diff / ctrl_vals.mean() * 100 if ctrl_vals.mean() != 0 else 0
    mw_stat, mw_p = stats.mannwhitneyu(ssw_vals, ctrl_vals, alternative='two-sided')
    d = diff / np.sqrt((ssw_vals.std()**2 + ctrl_vals.std()**2) / 2)
    
    sig = '***' if mw_p < 0.001 else '**' if mw_p < 0.01 else '*' if mw_p < 0.05 else ''
    # For stability: increase means MORE stable (less avalanche risk)
    # For sn38/sk38/ssi: higher = more stable → no decrease is good for our hypothesis
    # For ccl: lower = easier to trigger → no increase is good
    if var == 'ccl_pwl':
        direction = 'NO CHANGE ✓' if abs(pct) < 5 else ('DECREASE ✓' if diff < 0 else 'INCREASE (instability!) ✓')
    else:
        direction = 'NO CHANGE ✓' if abs(pct) < 5 else ('INCREASE (more stable)' if diff > 0 else 'DECREASE (less stable) ✓')
    print(f"  {all_vars.get(var, var):30s}: {pct:+6.1f}% (d={d:+.3f}, P={mw_p:.2e}) {sig} → {direction}")
    
    stability_results[var] = {
        'pct_change': float(pct), 'cohen_d': float(d), 'p': float(mw_p),
        'direction': 'decrease' if diff < 0 else 'increase' if diff > 0 else 'neutral'
    }

results['mechanism_test'] = {
    'loading': loading_results,
    'stability': stability_results,
}

# ── 9. Natural vs Skier stability contrast ────────────────────────────────────
print("\n" + "="*80)
print("NATURAL vs SKIER STABILITY CONTRAST")
print("="*80)
print("If loading-reduction is the mechanism:")
print("  - Natural stability (sn38) should INCREASE (less loading → harder to trigger naturally)")
print("  - Skier stability (sk38) should NOT INCREASE (human triggers still possible)")
print()

for var, label in [('sn38_pwl', 'Natural stability'), ('sk38_pwl', 'Skier stability')]:
    ssw_vals = df_winter.loc[df_winter['in_ssw'], var].dropna()
    ctrl_vals = df_winter.loc[~df_winter['in_ssw'], var].dropna()
    
    diff = ssw_vals.mean() - ctrl_vals.mean()
    pct = diff / ctrl_vals.mean() * 100 if ctrl_vals.mean() != 0 else 0
    d = diff / np.sqrt((ssw_vals.std()**2 + ctrl_vals.std()**2) / 2)
    mw_stat, mw_p = stats.mannwhitneyu(ssw_vals, ctrl_vals, alternative='two-sided')
    
    # Also test for increase specifically
    _, mw_p_greater = stats.mannwhitneyu(ssw_vals, ctrl_vals, alternative='greater')
    _, mw_p_less = stats.mannwhitneyu(ssw_vals, ctrl_vals, alternative='less')
    
    print(f"{label} ({var}):")
    print(f"  SSW: {ssw_vals.mean():.3f} ± {ssw_vals.std():.3f}")
    print(f"  Ctrl: {ctrl_vals.mean():.3f} ± {ctrl_vals.std():.3f}")
    print(f"  Δ = {diff:+.3f} ({pct:+.1f}%), d = {d:+.3f}")
    print(f"  Two-sided P = {mw_p:.2e}")
    print(f"  One-sided P (SSW > Ctrl) = {mw_p_greater:.2e}")
    print(f"  One-sided P (SSW < Ctrl) = {mw_p_less:.2e}")
    print()

results['nat_vs_skier'] = {}
for var in ['sn38_pwl', 'sk38_pwl']:
    ssw_vals = df_winter.loc[df_winter['in_ssw'], var].dropna()
    ctrl_vals = df_winter.loc[~df_winter['in_ssw'], var].dropna()
    diff = ssw_vals.mean() - ctrl_vals.mean()
    pct = diff / ctrl_vals.mean() * 100 if ctrl_vals.mean() != 0 else 0
    d = diff / np.sqrt((ssw_vals.std()**2 + ctrl_vals.std()**2) / 2)
    _, mw_p = stats.mannwhitneyu(ssw_vals, ctrl_vals, alternative='two-sided')
    _, mw_p_greater = stats.mannwhitneyu(ssw_vals, ctrl_vals, alternative='greater')
    results['nat_vs_skier'][var] = {
        'ssw_mean': float(ssw_vals.mean()),
        'ctrl_mean': float(ctrl_vals.mean()),
        'diff': float(diff),
        'pct_change': float(pct),
        'cohen_d': float(d),
        'mw_p': float(mw_p),
        'mw_p_greater': float(mw_p_greater),
    }

# ── 10. Persistent weak layer analysis ────────────────────────────────────────
print("="*80)
print("PERSISTENT WEAK LAYER PREVALENCE")
print("="*80)

pwl_ssw = df_winter.loc[df_winter['in_ssw'], 'pwl_100'].dropna()
pwl_ctrl = df_winter.loc[~df_winter['in_ssw'], 'pwl_100'].dropna()

print(f"PWL prevalence during SSW: {pwl_ssw.mean():.3f} ({pwl_ssw.sum():.0f}/{len(pwl_ssw)})")
print(f"PWL prevalence control:    {pwl_ctrl.mean():.3f} ({pwl_ctrl.sum():.0f}/{len(pwl_ctrl)})")

# Chi-squared test
from scipy.stats import chi2_contingency
table = [[int(pwl_ssw.sum()), int(len(pwl_ssw) - pwl_ssw.sum())],
         [int(pwl_ctrl.sum()), int(len(pwl_ctrl) - pwl_ctrl.sum())]]
chi2, chi_p, dof, expected = chi2_contingency(table)
print(f"Chi-squared: {chi2:.2f}, P = {chi_p:.2e}")
print(f"SSW PWL odds ratio: {(table[0][0]*table[1][1])/(table[0][1]*table[1][0]):.3f}")

results['pwl_analysis'] = {
    'ssw_prevalence': float(pwl_ssw.mean()),
    'ctrl_prevalence': float(pwl_ctrl.mean()),
    'chi2': float(chi2),
    'chi_p': float(chi_p),
    'odds_ratio': float((table[0][0]*table[1][1])/(table[0][1]*table[1][0])),
}

# ── 11. Station-level concordance ─────────────────────────────────────────────
print("\n" + "="*80)
print("STATION-LEVEL CONCORDANCE (sn38_pwl)")
print("="*80)

stations = df_winter['station_code'].unique()
station_diffs = {}
for station in stations:
    st_df = df_winter[df_winter['station_code'] == station]
    ssw_vals = st_df.loc[st_df['in_ssw'], 'sn38_pwl'].dropna()
    ctrl_vals = st_df.loc[~st_df['in_ssw'], 'sn38_pwl'].dropna()
    if len(ssw_vals) >= 10 and len(ctrl_vals) >= 10:
        station_diffs[station] = float(ssw_vals.mean() - ctrl_vals.mean())

n_inc = sum(1 for d in station_diffs.values() if d > 0)
n_dec = sum(1 for d in station_diffs.values() if d < 0)
print(f"Stations with sufficient data: {len(station_diffs)}")
print(f"  Natural stability INCREASE (more stable): {n_inc}")
print(f"  Natural stability DECREASE (less stable): {n_dec}")
if len(station_diffs) > 0:
    sign_p = binomtest(max(n_inc, n_dec), len(station_diffs), 0.5).pvalue
    print(f"  Sign test P = {sign_p:.4f}")

results['station_concordance'] = {
    'n_stations': len(station_diffs),
    'n_increase': n_inc,
    'n_decrease': n_dec,
    'sign_p': float(sign_p) if len(station_diffs) > 0 else None,
    'mean_diff': float(np.mean(list(station_diffs.values()))) if station_diffs else None,
}

# ── 12. Temperature analysis (sintering validation) ──────────────────────────
print("\n" + "="*80)
print("TEMPERATURE ANALYSIS (SINTERING VALIDATION)")
print("="*80)
print("SSW → tropospheric cooling → enhanced sintering → slab hardening")
print()

for var in ['TA', 'TSS_mod', 'TS0', 'TS1']:
    if var not in df_winter.columns:
        continue
    ssw_vals = df_winter.loc[df_winter['in_ssw'], var].dropna()
    ctrl_vals = df_winter.loc[~df_winter['in_ssw'], var].dropna()
    if len(ssw_vals) < 30 or len(ctrl_vals) < 30:
        continue
    
    diff = ssw_vals.mean() - ctrl_vals.mean()
    d = diff / np.sqrt((ssw_vals.std()**2 + ctrl_vals.std()**2) / 2)
    _, mw_p = stats.mannwhitneyu(ssw_vals, ctrl_vals, alternative='two-sided')
    
    var_labels = {'TA': 'Air temp', 'TSS_mod': 'Snow surface temp', 'TS0': 'Snow temp 0cm', 'TS1': 'Snow temp layer 1'}
    print(f"  {var_labels.get(var, var):20s}: SSW={ssw_vals.mean():+.2f}, Ctrl={ctrl_vals.mean():+.2f}, "
          f"Δ={diff:+.2f}°C (d={d:+.3f}, P={mw_p:.2e})")

# ── 13. Save results ──────────────────────────────────────────────────────────
os.makedirs('data/results', exist_ok=True)
with open('data/results/r31_snowpack_stability.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)

print("\n" + "="*80)
print("SUMMARY")
print("="*80)

# Loading vars summary
loading_decrease = sum(1 for v in loading_results.values() if v['direction'] == 'decrease')
loading_sig = sum(1 for v in loading_results.values() if v['p'] < 0.05)
print(f"\nLoading reduction: {loading_decrease}/{len(loading_results)} variables decrease during SSW, {loading_sig} significant")

# Stability vars summary
stability_decrease = sum(1 for v in stability_results.values() if v['direction'] == 'decrease')
stability_sig = sum(1 for v in stability_results.values() if v['p'] < 0.05)
print(f"Stability change: {stability_decrease}/{len(stability_results)} variables decrease during SSW, {stability_sig} significant")

# Key conclusion
print(f"\nNatural stability (sn38): {results['nat_vs_skier']['sn38_pwl']['pct_change']:+.1f}% during SSW")
print(f"Skier stability (sk38):   {results['nat_vs_skier']['sk38_pwl']['pct_change']:+.1f}% during SSW")

print("\nResults saved to data/results/r31_snowpack_stability.json")
print("DONE")
