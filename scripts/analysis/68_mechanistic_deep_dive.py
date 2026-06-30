"""
R68 Mechanistic Analysis: Wind, PWL, stability, temporal structure around SSW events.
Uses the 130-station SNOWPACK dataset to quantify:
1. Wind transport anomalies during SSW events
2. PWL (persistent weak layer) frequency
3. Hoar crystal formation
4. Temperature gradient (TS0-TS2)
5. Temporal lead-lag structure around SSW onset
6. Stability index changes
"""
import pandas as pd
import numpy as np
from scipy import stats
import json, os, warnings
warnings.filterwarnings('ignore')

base_data = r'C:\Users\Jack0\Solar-Magnetic-Analysis\data'
base_results = os.path.join(base_data, 'results')

# Load SSW events
ssw_events = pd.read_csv(os.path.join(base_results, 'ssw_event_catalog.csv'))
ssw_dates = pd.to_datetime(ssw_events['date'])
print(f"SSW events: {len(ssw_dates)}")

# Load weather/snowpack data (huge file, read efficiently)
print("Loading weather/snowpack data...")
cols_needed = [
    'datum', 'station_code', 'elevation_station',
    'TA', 'VW', 'VW_drift', 'DW',
    'wind_trans24', 'wind_trans24_3d', 'wind_trans24_7d',
    'HN24', 'SWE', 'MS_Snow', 'MS_Wind', 'MS_Rain',
    'hoar_size', 'pwl_100', 'ssi_pwl', 'sk38_pwl',
    'TS0', 'TS1', 'TS2',
    'Sclass2', 'zSd_mean', 'Sd', 'dangerLevel',
    'HS_mod', 'min_ccl_pen'
]
df = pd.read_csv(
    os.path.join(base_data, 'cryosphere', 'envidat', 'weather_snowpack_danger.csv'),
    usecols=[c for c in cols_needed if c in pd.read_csv(
        os.path.join(base_data, 'cryosphere', 'envidat', 'weather_snowpack_danger.csv'),
        nrows=0
    ).columns],
    parse_dates=['datum']
)
print(f"Loaded {len(df)} rows, {len(df.columns)} columns")
print(f"Date range: {df['datum'].min()} to {df['datum'].max()}")

# Create SSW window flags
# For each SSW event, flag days 0-30 after onset as "SSW period"
# Also flag days -30 to -1 as "pre-SSW" for comparison
df['is_ssw'] = False
df['is_pre_ssw'] = False
df['ssw_day'] = np.nan  # day relative to SSW onset

for onset in ssw_dates:
    mask_ssw = (df['datum'] >= onset) & (df['datum'] <= onset + pd.Timedelta(days=30))
    mask_pre = (df['datum'] >= onset - pd.Timedelta(days=30)) & (df['datum'] < onset)
    df.loc[mask_ssw, 'is_ssw'] = True
    df.loc[mask_pre, 'is_pre_ssw'] = True
    # Day relative to onset
    for day_offset in range(-30, 31):
        day_mask = df['datum'] == onset + pd.Timedelta(days=day_offset)
        df.loc[day_mask, 'ssw_day'] = day_offset

# Filter to winter months only (Nov-Apr)
df['month'] = df['datum'].dt.month
df_winter = df[df['month'].isin([11, 12, 1, 2, 3, 4])].copy()
df_ctrl = df_winter[~df_winter['is_ssw'] & ~df_winter['is_pre_ssw']].copy()
df_ssw = df_winter[df_winter['is_ssw']].copy()
df_pre = df_winter[df_winter['is_pre_ssw']].copy()

print(f"\nWinter data: {len(df_winter)} rows")
print(f"SSW period: {len(df_ssw)} rows")
print(f"Pre-SSW: {len(df_pre)} rows")
print(f"Control: {len(df_ctrl)} rows")

results = {}

# === 1. Wind Transport Analysis ===
print("\n=== WIND TRANSPORT ===")
for col in ['wind_trans24', 'VW', 'VW_drift']:
    ssw_mean = df_ssw[col].mean()
    ctrl_mean = df_ctrl[col].mean()
    pre_mean = df_pre[col].mean()
    
    # Effect size
    pooled_std = df_winter[col].std()
    d_ssw = (ssw_mean - ctrl_mean) / pooled_std if pooled_std > 0 else 0
    
    # Mann-Whitney U
    ssw_vals = df_ssw[col].dropna()
    ctrl_vals = df_ctrl[col].dropna()
    if len(ssw_vals) > 10 and len(ctrl_vals) > 10:
        u, p = stats.mannwhitneyu(ssw_vals, ctrl_vals, alternative='two-sided')
        print(f"  {col}: SSW={ssw_mean:.3f}, Pre={pre_mean:.3f}, Ctrl={ctrl_mean:.3f}, d={d_ssw:.3f}, P={p:.4f}")
    
    results[f'wind_{col}'] = {
        'ssw_mean': round(ssw_mean, 3),
        'pre_mean': round(pre_mean, 3),
        'ctrl_mean': round(ctrl_mean, 3),
        'cohen_d': round(d_ssw, 3),
        'pct_change': round((ssw_mean - ctrl_mean) / ctrl_mean * 100, 1) if ctrl_mean != 0 else None
    }

# === 2. PWL Analysis ===
print("\n=== PERSISTENT WEAK LAYERS ===")
for col in ['pwl_100', 'ssi_pwl', 'sk38_pwl', 'hoar_size']:
    ssw_mean = df_ssw[col].mean()
    ctrl_mean = df_ctrl[col].mean()
    pre_mean = df_pre[col].mean()
    pooled_std = df_winter[col].std()
    d_ssw = (ssw_mean - ctrl_mean) / pooled_std if pooled_std > 0 else 0
    
    ssw_vals = df_ssw[col].dropna()
    ctrl_vals = df_ctrl[col].dropna()
    if len(ssw_vals) > 10 and len(ctrl_vals) > 10:
        u, p = stats.mannwhitneyu(ssw_vals, ctrl_vals, alternative='two-sided')
        pct = (ssw_mean - ctrl_mean) / ctrl_mean * 100 if ctrl_mean != 0 else 0
        print(f"  {col}: SSW={ssw_mean:.4f}, Ctrl={ctrl_mean:.4f}, Δ={pct:+.1f}%, d={d_ssw:.3f}, P={p:.4f}")
    
    results[f'pwl_{col}'] = {
        'ssw_mean': round(ssw_mean, 4),
        'ctrl_mean': round(ctrl_mean, 4),
        'cohen_d': round(d_ssw, 3),
        'pct_change': round((ssw_mean - ctrl_mean) / ctrl_mean * 100, 1) if ctrl_mean != 0 else None
    }

# === 3. Temperature and Gradient ===
print("\n=== TEMPERATURE ===")
for col in ['TA', 'TS0', 'TS1', 'TS2']:
    ssw_mean = df_ssw[col].mean()
    ctrl_mean = df_ctrl[col].mean()
    pooled_std = df_winter[col].std()
    d = (ssw_mean - ctrl_mean) / pooled_std if pooled_std > 0 else 0
    
    ssw_vals = df_ssw[col].dropna()
    ctrl_vals = df_ctrl[col].dropna()
    if len(ssw_vals) > 10 and len(ctrl_vals) > 10:
        u, p = stats.mannwhitneyu(ssw_vals, ctrl_vals, alternative='two-sided')
        print(f"  {col}: SSW={ssw_mean:.2f}, Ctrl={ctrl_mean:.2f}, Δ={ssw_mean-ctrl_mean:+.2f}, d={d:.3f}, P={p:.4f}")
    
    results[f'temp_{col}'] = {
        'ssw_mean': round(ssw_mean, 2),
        'ctrl_mean': round(ctrl_mean, 2),
        'delta': round(ssw_mean - ctrl_mean, 2),
        'cohen_d': round(d, 3)
    }

# Temperature gradient proxy: TS0 - TS2 (surface to depth)
df_winter['TG_proxy'] = (df_winter['TS0'] - df_winter['TS2']).abs()
tg_ssw = df_winter.loc[df_winter['is_ssw'], 'TG_proxy'].mean()
tg_ctrl = df_winter.loc[~df_winter['is_ssw'] & ~df_winter['is_pre_ssw'], 'TG_proxy'].mean()
tg_d = (tg_ssw - tg_ctrl) / df_winter['TG_proxy'].std()
print(f"  TG proxy (|TS0-TS2|): SSW={tg_ssw:.2f}, Ctrl={tg_ctrl:.2f}, d={tg_d:.3f}")
results['temp_gradient'] = {'ssw': round(tg_ssw, 2), 'ctrl': round(tg_ctrl, 2), 'cohen_d': round(tg_d, 3)}

# === 4. SWE and New Snow ===
print("\n=== SNOW ===")
for col in ['SWE', 'HN24', 'MS_Snow', 'MS_Rain']:
    ssw_mean = df_ssw[col].mean()
    ctrl_mean = df_ctrl[col].mean()
    pooled_std = df_winter[col].std()
    d = (ssw_mean - ctrl_mean) / pooled_std if pooled_std > 0 else 0
    pct = (ssw_mean - ctrl_mean) / ctrl_mean * 100 if ctrl_mean != 0 else 0
    print(f"  {col}: SSW={ssw_mean:.3f}, Ctrl={ctrl_mean:.3f}, Δ={pct:+.1f}%, d={d:.3f}")
    results[f'snow_{col}'] = {
        'ssw_mean': round(ssw_mean, 3),
        'ctrl_mean': round(ctrl_mean, 3),
        'pct_change': round(pct, 1),
        'cohen_d': round(d, 3)
    }

# === 5. Stability Metrics ===
print("\n=== STABILITY ===")
for col in ['Sclass2', 'Sd', 'zSd_mean', 'dangerLevel', 'min_ccl_pen']:
    ssw_mean = df_ssw[col].mean()
    ctrl_mean = df_ctrl[col].mean()
    pooled_std = df_winter[col].std()
    d = (ssw_mean - ctrl_mean) / pooled_std if pooled_std > 0 else 0
    ssw_vals = df_ssw[col].dropna()
    ctrl_vals = df_ctrl[col].dropna()
    if len(ssw_vals) > 10 and len(ctrl_vals) > 10:
        u, p = stats.mannwhitneyu(ssw_vals, ctrl_vals, alternative='two-sided')
        print(f"  {col}: SSW={ssw_mean:.3f}, Ctrl={ctrl_mean:.3f}, d={d:.3f}, P={p:.4f}")
    results[f'stability_{col}'] = {
        'ssw_mean': round(ssw_mean, 3),
        'ctrl_mean': round(ctrl_mean, 3),
        'cohen_d': round(d, 3)
    }

# === 6. Temporal Lead-Lag Structure ===
print("\n=== TEMPORAL LEAD-LAG ===")
# Compute daily means across all stations for key variables relative to SSW onset
lag_data = df[df['ssw_day'].notna()].copy()
temporal = lag_data.groupby('ssw_day').agg({
    'wind_trans24': 'mean',
    'TA': 'mean',
    'HN24': 'mean',
    'SWE': 'mean',
    'pwl_100': 'mean',
    'hoar_size': 'mean',
    'MS_Rain': 'mean',
    'dangerLevel': 'mean'
}).reset_index()

# Compute pre-SSW baseline (days -30 to -1)
baseline = temporal[temporal['ssw_day'] < 0]
onset_period = temporal[(temporal['ssw_day'] >= 0) & (temporal['ssw_day'] <= 15)]
late_period = temporal[(temporal['ssw_day'] > 15) & (temporal['ssw_day'] <= 30)]

print("Variable: Pre-SSW → Onset(0-15d) → Late(15-30d)")
for col in ['wind_trans24', 'TA', 'HN24', 'pwl_100', 'MS_Rain', 'dangerLevel']:
    pre = baseline[col].mean()
    onset = onset_period[col].mean()
    late = late_period[col].mean()
    print(f"  {col}: {pre:.3f} → {onset:.3f} → {late:.3f}")
    results[f'temporal_{col}'] = {
        'pre_ssw': round(pre, 3), 'onset_0_15': round(onset, 3), 'late_15_30': round(late, 3)
    }

# Peak day for PWL
if not temporal.empty:
    peak_pwl_day = temporal.loc[temporal['pwl_100'].idxmax(), 'ssw_day']
    peak_rain_min_day = temporal.loc[temporal['MS_Rain'].idxmin(), 'ssw_day']
    print(f"  Peak PWL at day: {peak_pwl_day}")
    print(f"  Min rain at day: {peak_rain_min_day}")
    results['temporal_peaks'] = {
        'peak_pwl_day': int(peak_pwl_day),
        'min_rain_day': int(peak_rain_min_day)
    }

# === 7. Elevation-stratified analysis ===
print("\n=== ELEVATION STRATIFICATION ===")
df_winter['elev_band'] = pd.cut(df_winter['elevation_station'], bins=[0, 1500, 2000, 2500, 4000], labels=['<1500m', '1500-2000m', '2000-2500m', '>2500m'])
for band in ['<1500m', '1500-2000m', '2000-2500m', '>2500m']:
    band_df = df_winter[df_winter['elev_band'] == band]
    ssw_swe = band_df.loc[band_df['is_ssw'], 'SWE'].mean()
    ctrl_swe = band_df.loc[~band_df['is_ssw'] & ~band_df['is_pre_ssw'], 'SWE'].mean()
    pct = (ssw_swe - ctrl_swe) / ctrl_swe * 100 if ctrl_swe != 0 else 0
    
    ssw_pwl = band_df.loc[band_df['is_ssw'], 'pwl_100'].mean()
    ctrl_pwl = band_df.loc[~band_df['is_ssw'] & ~band_df['is_pre_ssw'], 'pwl_100'].mean()
    pwl_pct = (ssw_pwl - ctrl_pwl) / ctrl_pwl * 100 if ctrl_pwl != 0 else 0
    
    n_stations = band_df['station_code'].nunique()
    print(f"  {band} (n={n_stations}): SWE Δ={pct:+.1f}%, PWL Δ={pwl_pct:+.1f}%")
    results[f'elev_{band}'] = {
        'n_stations': int(n_stations),
        'swe_pct_change': round(pct, 1),
        'pwl_pct_change': round(pwl_pct, 1)
    }

# Save results
with open(os.path.join(base_results, 'r68_mechanistic_analysis.json'), 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nResults saved to r68_mechanistic_analysis.json")
