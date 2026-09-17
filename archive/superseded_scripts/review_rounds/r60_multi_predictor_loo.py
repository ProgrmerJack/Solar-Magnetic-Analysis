"""
Multi-predictor leave-one-out Brier Skill Score analysis for SSW-avalanche prediction.

Analysis at EVENT LEVEL (n=16 SSW events):
1. SSW-only: base rate predictor (14/16 = 87.5% suppression)
2. Z500-only: logistic regression on mean Z500 anomaly during SSW window, LOO
3. Z500+SSW-lag: Z500 anomaly + lag between SSW onset and peak blocking
4. Z500+Precip: Z500 + mean precipitation during SSW window

Targets:
- Binary: is this SSW event a suppression event? (RR < 1)
- Continuous: event-level RR (observed/expected avalanche count)
"""

import pandas as pd
import numpy as np
import json
from datetime import timedelta
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# LOAD DATA
# ============================================================================
print("Loading data...")
analysis_df = pd.read_parquet(
    r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\analysis_panel_v2.parquet'
)
ssw_catalog = pd.read_parquet(
    r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\atmospheric\ssw_catalog.parquet'
)

# Filter analysis data
df = analysis_df[
    (analysis_df['winter_id'].notna()) & 
    (analysis_df['dry_natural_size_1234'].notna())
].copy()

print(f"Analysis data shape after filtering: {df.shape}")
print(f"SSW catalog shape: {ssw_catalog.shape}")

# Ensure datetime index is tz-aware
if df.index.tz is None:
    df.index = df.index.tz_localize('UTC')
if ssw_catalog.index.tz is None:
    ssw_catalog.index = ssw_catalog.index.tz_localize('UTC')

# ============================================================================
# BUILD EVENT-LEVEL DATASET
# ============================================================================
print("\nBuilding event-level dataset...")

def get_event_features_and_target(ssw_onset_date, df, ssw_catalog, window_days=15):
    """
    For a single SSW event, compute:
    - observed: avalanche count in ±window_days around SSW onset
    - expected: DOY-matched non-SSW control (same DOY ±window_days in other years)
    - RR = observed / expected
    - Z500 anomaly during SSW window
    - Lag between SSW onset and peak blocking (proxy: std of Z500 during window)
    - Precipitation during SSW window
    """
    # SSW window: ±15 days around onset
    ssw_start = ssw_onset_date - timedelta(days=window_days)
    ssw_end = ssw_onset_date + timedelta(days=window_days)
    
    # Get observations in SSW window
    ssw_window_mask = (df.index >= ssw_start) & (df.index <= ssw_end)
    ssw_window_data = df[ssw_window_mask]
    observed = ssw_window_data['dry_natural_size_1234'].sum()
    
    # Get DOY (day of year) for SSW onset
    ssw_doy = ssw_onset_date.timetuple().tm_yday
    
    # DOY-matched control: same DOY ±window_days in OTHER years (non-SSW years)
    ssw_year = ssw_onset_date.year
    expected_sum = 0
    control_count = 0
    
    for year in df.index.year.unique():
        if year == ssw_year:
            continue  # Skip the SSW year itself
        
        # Check if this year is a non-SSW year (no SSW events that overlap with our window)
        is_ssw_year = False
        for other_onset in ssw_catalog.index:
            if other_onset.year == year:
                other_ssw_start = other_onset - timedelta(days=window_days)
                other_ssw_end = other_onset + timedelta(days=window_days)
                if (ssw_start <= other_ssw_end) and (other_ssw_start <= ssw_end):
                    is_ssw_year = True
                    break
        
        if is_ssw_year:
            continue
        
        # Get data for this DOY ±window_days in this year
        try:
            target_date = pd.Timestamp(year=year, month=ssw_onset_date.month, 
                                      day=ssw_onset_date.day, tz='UTC')
            control_start = target_date - timedelta(days=window_days)
            control_end = target_date + timedelta(days=window_days)
            
            control_mask = (df.index >= control_start) & (df.index <= control_end)
            control_data = df[control_mask]
            expected_sum += control_data['dry_natural_size_1234'].sum()
            control_count += 1
        except:
            continue
    
    # Expected: average from control years
    expected = expected_sum / control_count if control_count > 0 else observed
    
    # Avoid division by zero
    if expected == 0:
        expected = 1
    
    rr = observed / expected
    
    # Z500 anomaly during SSW window
    z500_window = ssw_window_data['ncep_z500_nh'].dropna()
    z500_mean = z500_window.mean() if len(z500_window) > 0 else np.nan
    z500_overall_mean = df['ncep_z500_nh'].mean()
    z500_anomaly = z500_mean - z500_overall_mean if not np.isnan(z500_mean) else np.nan
    
    # SSW-lag proxy: std of Z500 during window (as proxy for lag between onset and peak blocking)
    z500_std = z500_window.std() if len(z500_window) > 1 else np.nan
    
    # Precipitation during SSW window
    precip_window = ssw_window_data['snotel_prec_mean'].dropna()
    precip_mean = precip_window.mean() if len(precip_window) > 0 else np.nan
    
    return {
        'ssw_onset_date': ssw_onset_date,
        'observed': observed,
        'expected': expected,
        'rr': rr,
        'is_suppression': 1 if rr < 1 else 0,
        'z500_anomaly': z500_anomaly,
        'z500_std': z500_std,  # proxy for SSW-lag
        'precip_mean': precip_mean,
        'n_control_years': control_count
    }

# Build event-level data
event_list = []
for onset_date in ssw_catalog.index:
    features = get_event_features_and_target(onset_date, df, ssw_catalog)
    event_list.append(features)

event_df = pd.DataFrame(event_list)
print(f"\nEvent-level dataset shape: {event_df.shape}")
print(f"Suppression events (RR < 1): {event_df['is_suppression'].sum()} / {len(event_df)}")
print(f"Base rate of suppression: {event_df['is_suppression'].mean():.3f}")

# ============================================================================
# REMOVE ROWS WITH NaN IN KEY FEATURES
# ============================================================================
print("\nChecking for NaN values...")
print(event_df.isnull().sum())

event_df = event_df.dropna(subset=['z500_anomaly', 'precip_mean'])
print(f"\nEvent-level dataset after removing NaN: {event_df.shape}")

# ============================================================================
# DEFINE MODELS AND LOO BSS COMPUTATION
# ============================================================================

def brier_skill_score_binary(y_true, y_pred_proba, base_rate=None):
    """
    Brier Skill Score for binary prediction.
    BSS = 1 - BS(model) / BS(climatology)
    
    y_true: binary targets (0/1)
    y_pred_proba: predicted probabilities [0, 1]
    base_rate: climatological base rate (if None, use mean of y_true)
    """
    if base_rate is None:
        base_rate = y_true.mean()
    
    # Brier score
    bs_model = np.mean((y_pred_proba - y_true) ** 2)
    bs_climatology = base_rate * (1 - base_rate)
    
    # Brier Skill Score
    bss = 1 - (bs_model / bs_climatology) if bs_climatology > 0 else np.nan
    return bss, bs_model, bs_climatology

def brier_skill_score_continuous(y_true, y_pred, base_rate=None):
    """
    Continuous BSS = 1 - MSE(model) / MSE(climatology)
    
    y_true: continuous targets
    y_pred: predicted values
    base_rate: climatological value (if None, use mean of y_true)
    """
    if base_rate is None:
        base_rate = y_true.mean()
    
    mse_model = np.mean((y_pred - y_true) ** 2)
    mse_climatology = np.mean((base_rate - y_true) ** 2)
    
    bss = 1 - (mse_model / mse_climatology) if mse_climatology > 0 else np.nan
    return bss, mse_model, mse_climatology

# ============================================================================
# 1. SSW-ONLY MODEL (base rate constant predictor)
# ============================================================================
print("\n" + "="*80)
print("MODEL 1: SSW-ONLY (base rate constant predictor)")
print("="*80)

y_binary = event_df['is_suppression'].values
base_rate = y_binary.mean()
y_pred_ssw_only = np.full_like(y_binary, base_rate, dtype=float)

bss_ssw_only_binary, bs_ssw, bs_clim_ssw = brier_skill_score_binary(
    y_binary, y_pred_ssw_only, base_rate=base_rate
)
print(f"Base rate (suppression): {base_rate:.3f}")
print(f"LOO BSS (binary): {bss_ssw_only_binary:.4f}")
print(f"Brier score: {bs_ssw:.4f}, Climatology: {bs_clim_ssw:.4f}")

# Continuous: RR
y_continuous = event_df['rr'].values
y_pred_rr_ssw_only = np.full_like(y_continuous, y_continuous.mean())
bss_ssw_only_continuous, mse_ssw, mse_clim_ssw = brier_skill_score_continuous(
    y_continuous, y_pred_rr_ssw_only
)
print(f"LOO BSS (continuous RR): {bss_ssw_only_continuous:.4f}")

# ============================================================================
# 2. Z500-ONLY MODEL (logistic regression, LOO)
# ============================================================================
print("\n" + "="*80)
print("MODEL 2: Z500-ONLY (logistic regression, LOO)")
print("="*80)

X_z500_only = event_df[['z500_anomaly']].values
scaler_z500 = StandardScaler()
X_z500_only_scaled = scaler_z500.fit_transform(X_z500_only)

y_pred_z500_loo = np.zeros(len(event_df))

for i in range(len(event_df)):
    # Leave-one-out: use all but i
    X_train = np.delete(X_z500_only_scaled, i, axis=0)
    y_train = np.delete(y_binary, i)
    
    # Fit logistic regression
    clf = LogisticRegression(random_state=42, max_iter=1000)
    clf.fit(X_train, y_train)
    
    # Predict on held-out sample
    X_test = X_z500_only_scaled[i:i+1]
    y_pred_z500_loo[i] = clf.predict_proba(X_test)[0, 1]

bss_z500_only, bs_z500, bs_clim_z500 = brier_skill_score_binary(
    y_binary, y_pred_z500_loo, base_rate=base_rate
)
print(f"LOO BSS (binary): {bss_z500_only:.4f}")
print(f"Brier score: {bs_z500:.4f}, Climatology: {bs_clim_z500:.4f}")

# Continuous: fit linear regression
from sklearn.linear_model import LinearRegression

y_pred_z500_loo_continuous = np.zeros(len(event_df))
for i in range(len(event_df)):
    X_train = np.delete(X_z500_only_scaled, i, axis=0)
    y_train = np.delete(y_continuous, i)
    
    reg = LinearRegression()
    reg.fit(X_train, y_train)
    
    X_test = X_z500_only_scaled[i:i+1]
    y_pred_z500_loo_continuous[i] = reg.predict(X_test)[0]

bss_z500_only_cont, mse_z500, mse_clim_z500 = brier_skill_score_continuous(
    y_continuous, y_pred_z500_loo_continuous
)
print(f"LOO BSS (continuous RR): {bss_z500_only_cont:.4f}")

# ============================================================================
# 3. Z500+SSW-LAG MODEL (logistic regression, LOO)
# ============================================================================
print("\n" + "="*80)
print("MODEL 3: Z500+SSW-LAG (logistic regression, LOO)")
print("="*80)

X_z500_lag = event_df[['z500_anomaly', 'z500_std']].values
scaler_z500_lag = StandardScaler()
X_z500_lag_scaled = scaler_z500_lag.fit_transform(X_z500_lag)

y_pred_z500_lag_loo = np.zeros(len(event_df))

for i in range(len(event_df)):
    X_train = np.delete(X_z500_lag_scaled, i, axis=0)
    y_train = np.delete(y_binary, i)
    
    clf = LogisticRegression(random_state=42, max_iter=1000)
    clf.fit(X_train, y_train)
    
    X_test = X_z500_lag_scaled[i:i+1]
    y_pred_z500_lag_loo[i] = clf.predict_proba(X_test)[0, 1]

bss_z500_lag, bs_z500_lag, bs_clim_z500_lag = brier_skill_score_binary(
    y_binary, y_pred_z500_lag_loo, base_rate=base_rate
)
print(f"LOO BSS (binary): {bss_z500_lag:.4f}")
print(f"Brier score: {bs_z500_lag:.4f}, Climatology: {bs_clim_z500_lag:.4f}")

# Continuous
y_pred_z500_lag_loo_continuous = np.zeros(len(event_df))
for i in range(len(event_df)):
    X_train = np.delete(X_z500_lag_scaled, i, axis=0)
    y_train = np.delete(y_continuous, i)
    
    reg = LinearRegression()
    reg.fit(X_train, y_train)
    
    X_test = X_z500_lag_scaled[i:i+1]
    y_pred_z500_lag_loo_continuous[i] = reg.predict(X_test)[0]

bss_z500_lag_cont, mse_z500_lag, mse_clim_z500_lag = brier_skill_score_continuous(
    y_continuous, y_pred_z500_lag_loo_continuous
)
print(f"LOO BSS (continuous RR): {bss_z500_lag_cont:.4f}")

# ============================================================================
# 4. Z500+PRECIP MODEL (logistic regression, LOO)
# ============================================================================
print("\n" + "="*80)
print("MODEL 4: Z500+PRECIP (logistic regression, LOO)")
print("="*80)

X_z500_precip = event_df[['z500_anomaly', 'precip_mean']].values
scaler_z500_precip = StandardScaler()
X_z500_precip_scaled = scaler_z500_precip.fit_transform(X_z500_precip)

y_pred_z500_precip_loo = np.zeros(len(event_df))

for i in range(len(event_df)):
    X_train = np.delete(X_z500_precip_scaled, i, axis=0)
    y_train = np.delete(y_binary, i)
    
    clf = LogisticRegression(random_state=42, max_iter=1000)
    clf.fit(X_train, y_train)
    
    X_test = X_z500_precip_scaled[i:i+1]
    y_pred_z500_precip_loo[i] = clf.predict_proba(X_test)[0, 1]

bss_z500_precip, bs_z500_precip, bs_clim_z500_precip = brier_skill_score_binary(
    y_binary, y_pred_z500_precip_loo, base_rate=base_rate
)
print(f"LOO BSS (binary): {bss_z500_precip:.4f}")
print(f"Brier score: {bs_z500_precip:.4f}, Climatology: {bs_clim_z500_precip:.4f}")

# Continuous
y_pred_z500_precip_loo_continuous = np.zeros(len(event_df))
for i in range(len(event_df)):
    X_train = np.delete(X_z500_precip_scaled, i, axis=0)
    y_train = np.delete(y_continuous, i)
    
    reg = LinearRegression()
    reg.fit(X_train, y_train)
    
    X_test = X_z500_precip_scaled[i:i+1]
    y_pred_z500_precip_loo_continuous[i] = reg.predict(X_test)[0]

bss_z500_precip_cont, mse_z500_precip, mse_clim_z500_precip = brier_skill_score_continuous(
    y_continuous, y_pred_z500_precip_loo_continuous
)
print(f"LOO BSS (continuous RR): {bss_z500_precip_cont:.4f}")

# ============================================================================
# SUMMARY AND COMPARISON
# ============================================================================
print("\n" + "="*80)
print("SUMMARY: MODEL COMPARISON")
print("="*80)

results_summary = {
    'n_events': len(event_df),
    'base_rate_suppression': float(base_rate),
    'models': {
        'ssw_only': {
            'description': 'Base rate constant predictor',
            'loo_bss_binary': float(bss_ssw_only_binary),
            'brier_score': float(bs_ssw),
            'loo_bss_continuous': float(bss_ssw_only_continuous),
            'mse': float(mse_ssw)
        },
        'z500_only': {
            'description': 'Z500 anomaly only (logistic regression)',
            'loo_bss_binary': float(bss_z500_only),
            'brier_score': float(bs_z500),
            'loo_bss_continuous': float(bss_z500_only_cont),
            'mse': float(mse_z500)
        },
        'z500_ssw_lag': {
            'description': 'Z500 anomaly + Z500 std (SSW-lag proxy)',
            'loo_bss_binary': float(bss_z500_lag),
            'brier_score': float(bs_z500_lag),
            'loo_bss_continuous': float(bss_z500_lag_cont),
            'mse': float(mse_z500_lag)
        },
        'z500_precip': {
            'description': 'Z500 anomaly + precipitation',
            'loo_bss_binary': float(bss_z500_precip),
            'brier_score': float(bs_z500_precip),
            'loo_bss_continuous': float(bss_z500_precip_cont),
            'mse': float(mse_z500_precip)
        }
    }
}

print("\nBINARY PREDICTION (suppression vs. non-suppression):")
print(f"{'Model':<25} {'LOO BSS':<12} {'Brier Score':<12}")
print("-" * 50)
for model_name, model_results in results_summary['models'].items():
    print(f"{model_name:<25} {model_results['loo_bss_binary']:<12.4f} {model_results['brier_score']:<12.4f}")

print("\nCONTINUOUS PREDICTION (event-level RR):")
print(f"{'Model':<25} {'LOO BSS':<12} {'MSE':<12}")
print("-" * 50)
for model_name, model_results in results_summary['models'].items():
    print(f"{model_name:<25} {model_results['loo_bss_continuous']:<12.4f} {model_results['mse']:<12.4f}")

print("\nKEY FINDING: Does SSW add value to Z500 predictor?")
z500_bss = results_summary['models']['z500_only']['loo_bss_binary']
z500_ssw_lag_bss = results_summary['models']['z500_ssw_lag']['loo_bss_binary']
improvement = z500_ssw_lag_bss - z500_bss

print(f"Z500-only BSS: {z500_bss:.4f}")
print(f"Z500+SSW-lag BSS: {z500_ssw_lag_bss:.4f}")
print(f"Improvement: {improvement:+.4f}")

if improvement > 0:
    print(f"✓ SSW-lag adds VALUE to Z500 predictor")
else:
    print(f"✗ SSW-lag does NOT add value to Z500 predictor")

# ============================================================================
# SAVE RESULTS
# ============================================================================
output_dir = r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\results'
import os
os.makedirs(output_dir, exist_ok=True)

output_path = os.path.join(output_dir, 'r60_multi_predictor_loo.json')
with open(output_path, 'w') as f:
    json.dump(results_summary, f, indent=2)

print(f"\n✓ Results saved to: {output_path}")

print("\nDone!")
