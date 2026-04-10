"""
R32: Mediation analysis — Do SNOWPACK trigger proxies predict avalanche counts?
Closes the quantitative gap: regime shift → trigger suppression → avalanche reduction.
"""
import pandas as pd
import numpy as np
from scipy import stats
import json, warnings
warnings.filterwarnings('ignore')

print("=" * 70)
print("R32: SNOWPACK TRIGGER → AVALANCHE MEDIATION ANALYSIS")
print("=" * 70)

# Load SNOWPACK data
snow = pd.read_csv('data/cryosphere/envidat/weather_snowpack_danger.csv', low_memory=False)
snow['date'] = pd.to_datetime(snow['datum'])
print(f"\nSNOWPACK: {len(snow):,} station-days, {snow['station_code'].nunique()} stations")

# Load avalanche panel
panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
panel = panel.reset_index()
panel['date'] = pd.to_datetime(panel['time'])
print(f"Panel: {len(panel):,} days")

# Load SSW catalog
ssw = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw = ssw.reset_index()
ssw['onset'] = pd.to_datetime(ssw['onset_date']).dt.tz_localize(None)

# Winter filter for SNOWPACK
snow_winter = snow[(snow['date'].dt.month >= 11) | (snow['date'].dt.month <= 4)].copy()

# ============================================================
# SECTION 1: Daily national trigger index from SNOWPACK
# ============================================================
print("\n" + "=" * 70)
print("SECTION 1: NATIONAL TRIGGER INDEX")
print("=" * 70)

# For each day, compute fraction of stations with each trigger type
daily_triggers = snow_winter.groupby('date').agg(
    n_stations=('station_code', 'nunique'),
    frac_surface_warm=('TSS_mod', lambda x: (x > -2).mean()),
    frac_rain=('MS_Rain', lambda x: (x > 0).mean()),
    mean_iswr=('ISWR', 'mean'),
    frac_high_iswr=('ISWR', lambda x: (x > 120).mean()),
    mean_sn38=('sn38_pwl', 'mean'),
    mean_hn24=('HN24', 'mean'),
    mean_tss=('TSS_mod', 'mean'),
).reset_index()

# Composite trigger index: fraction of stations with ANY surface trigger
def any_trigger(group):
    warm = group['TSS_mod'] > -2
    rain = group['MS_Rain'] > 0
    return (warm | rain).mean()

daily_any_trigger = snow_winter.groupby('date').apply(any_trigger, include_groups=False).reset_index()
daily_any_trigger.columns = ['date', 'frac_any_trigger']
daily_triggers = daily_triggers.merge(daily_any_trigger, on='date')

# Merge with avalanche panel
merged = panel[['date', 'dry_natural_size_1234']].merge(daily_triggers, on='date', how='inner')
merged = merged.dropna(subset=['dry_natural_size_1234', 'frac_surface_warm'])
print(f"Merged dataset: {len(merged)} days")

# ============================================================
# SECTION 2: Trigger proxies predict avalanche counts
# ============================================================
print("\n" + "=" * 70)
print("SECTION 2: DO TRIGGERS PREDICT AVALANCHE COUNTS?")
print("=" * 70)

aval = merged['dry_natural_size_1234'].values
predictors = {
    'frac_surface_warm': merged['frac_surface_warm'].values,
    'frac_rain': merged['frac_rain'].values,
    'mean_iswr': merged['mean_iswr'].values,
    'frac_high_iswr': merged['frac_high_iswr'].values,
    'frac_any_trigger': merged['frac_any_trigger'].values,
    'mean_sn38': merged['mean_sn38'].values,
    'mean_hn24': merged['mean_hn24'].values,
}

print(f"\nCorrelation of daily trigger indices with avalanche counts:")
print(f"{'Trigger proxy':<25} {'r':>8} {'P':>12} {'R²':>8}")
print("-" * 55)

for name, vals in predictors.items():
    mask = np.isfinite(vals) & np.isfinite(aval)
    r, p = stats.pearsonr(vals[mask], aval[mask])
    print(f"{name:<25} {r:>8.3f} {p:>12.2e} {r**2:>8.3f}")

# ============================================================
# SECTION 3: SSW effect on national trigger indices
# ============================================================
print("\n" + "=" * 70)
print("SECTION 3: SSW EFFECT ON NATIONAL TRIGGER INDICES")
print("=" * 70)

# Mark SSW windows
merged['in_ssw'] = False
ssw_study = ssw[(ssw['onset'] >= '1997-01-01') & (ssw['onset'] <= '2020-12-31')]
for _, row in ssw_study.iterrows():
    onset = row['onset']
    start = onset - pd.Timedelta(days=15)
    end = onset + pd.Timedelta(days=15)
    merged.loc[(merged['date'] >= start) & (merged['date'] <= end), 'in_ssw'] = True

n_ssw = merged['in_ssw'].sum()
n_ctrl = (~merged['in_ssw']).sum()
print(f"\nSSW days: {n_ssw}, Control days: {n_ctrl}")

print(f"\n{'Variable':<25} {'SSW mean':>10} {'Ctrl mean':>10} {'Ratio':>8} {'Cohen d':>8} {'P':>12}")
print("-" * 75)

for name in ['frac_surface_warm', 'frac_rain', 'mean_iswr', 'frac_high_iswr', 
             'frac_any_trigger', 'mean_sn38', 'mean_hn24']:
    ssw_vals = merged.loc[merged['in_ssw'], name].dropna()
    ctrl_vals = merged.loc[~merged['in_ssw'], name].dropna()
    ssw_m = ssw_vals.mean()
    ctrl_m = ctrl_vals.mean()
    ratio = ssw_m / ctrl_m if ctrl_m != 0 else np.nan
    pooled_std = np.sqrt((ssw_vals.std()**2 + ctrl_vals.std()**2) / 2)
    d = (ssw_m - ctrl_m) / pooled_std if pooled_std > 0 else 0
    _, p = stats.mannwhitneyu(ssw_vals, ctrl_vals, alternative='two-sided')
    print(f"{name:<25} {ssw_m:>10.4f} {ctrl_m:>10.4f} {ratio:>8.3f} {d:>8.3f} {p:>12.2e}")

# Avalanche counts
ssw_aval = merged.loc[merged['in_ssw'], 'dry_natural_size_1234'].dropna()
ctrl_aval = merged.loc[~merged['in_ssw'], 'dry_natural_size_1234'].dropna()
print(f"\n{'Avalanche counts':<25} {ssw_aval.mean():>10.4f} {ctrl_aval.mean():>10.4f} "
      f"{ssw_aval.mean()/ctrl_aval.mean():>8.3f}")

# ============================================================
# SECTION 4: EXPECTED AVALANCHE REDUCTION FROM TRIGGER CHANGES
# ============================================================
print("\n" + "=" * 70)
print("SECTION 4: EXPECTED VS OBSERVED AVALANCHE REDUCTION")
print("=" * 70)

# Linear model: avalanche = f(triggers)
from numpy.polynomial import polynomial as P
from sklearn.linear_model import LinearRegression

# Use multiple triggers
X_cols = ['frac_surface_warm', 'frac_rain', 'mean_iswr', 'mean_hn24', 'mean_sn38']
X = merged[X_cols].values
y = merged['dry_natural_size_1234'].values
mask = np.all(np.isfinite(X), axis=1) & np.isfinite(y)
X_clean = X[mask]
y_clean = y[mask]

lr = LinearRegression()
lr.fit(X_clean, y_clean)
r2 = lr.score(X_clean, y_clean)
print(f"\nMultiple regression: avalanche ~ triggers")
print(f"R² = {r2:.4f}")
for i, col in enumerate(X_cols):
    print(f"  {col:<25} coef = {lr.coef_[i]:>10.4f}")

# Predict avalanche counts under SSW vs control trigger distributions
X_ssw = merged.loc[merged['in_ssw'], X_cols].dropna()
X_ctrl = merged.loc[~merged['in_ssw'], X_cols].dropna()

pred_ssw = lr.predict(X_ssw.values).mean()
pred_ctrl = lr.predict(X_ctrl.values).mean()
expected_rr = pred_ssw / pred_ctrl if pred_ctrl != 0 else np.nan

print(f"\nPredicted avalanche rate (SSW triggers): {pred_ssw:.3f}")
print(f"Predicted avalanche rate (Ctrl triggers): {pred_ctrl:.3f}")
print(f"Expected RR from trigger changes: {expected_rr:.3f}")
print(f"Observed RR: {ssw_aval.mean() / ctrl_aval.mean():.3f}")
print(f"Trigger-mediated fraction: {(1 - expected_rr) / (1 - ssw_aval.mean()/ctrl_aval.mean()) * 100:.1f}%")

# ============================================================
# SECTION 5: REGIME-TRIGGER DECOMPOSITION
# ============================================================
print("\n" + "=" * 70)
print("SECTION 5: REGIME × TRIGGER INTERACTION")
print("=" * 70)

# Within cold-dry days only: do trigger proxies still differ SSW vs control?
# Use ERA5 data from panel for regime classification
panel_merged = panel[['date', 'dry_natural_size_1234']].copy()
panel_merged['date'] = pd.to_datetime(panel_merged['date'])

# Check if ERA5 columns exist
era5_cols = [c for c in panel.columns if 'era5' in c.lower() or 't2m' in c.lower() or 'tp' in c.lower()]
print(f"ERA5-related columns: {era5_cols[:10]}")

# Use temperature and precipitation from SNOWPACK data for regime classification
daily_weather = snow_winter.groupby('date').agg(
    mean_ta=('TA', 'mean'),
    mean_precip=('HN24', 'mean'),
).reset_index()

merged2 = merged.merge(daily_weather[['date', 'mean_ta', 'mean_precip']], on='date', how='left',
                        suffixes=('', '_weather'))

# Regime classification
ta_med = merged2['mean_ta'].median()
precip_med = merged2['mean_precip'].median()
merged2['regime'] = 'unknown'
merged2.loc[(merged2['mean_ta'] < ta_med) & (merged2['mean_precip'] < precip_med), 'regime'] = 'cold-dry'
merged2.loc[(merged2['mean_ta'] < ta_med) & (merged2['mean_precip'] >= precip_med), 'regime'] = 'cold-wet'
merged2.loc[(merged2['mean_ta'] >= ta_med) & (merged2['mean_precip'] < precip_med), 'regime'] = 'warm-dry'
merged2.loc[(merged2['mean_ta'] >= ta_med) & (merged2['mean_precip'] >= precip_med), 'regime'] = 'warm-wet'

print("\nRegime distribution (SSW vs Control):")
for regime in ['cold-dry', 'cold-wet', 'warm-dry', 'warm-wet']:
    ssw_frac = (merged2.loc[merged2['in_ssw'], 'regime'] == regime).mean()
    ctrl_frac = (merged2.loc[~merged2['in_ssw'], 'regime'] == regime).mean()
    print(f"  {regime:<12} SSW: {ssw_frac:.1%}  Ctrl: {ctrl_frac:.1%}  Shift: {(ssw_frac-ctrl_frac)*100:+.1f}pp")

# Within cold-dry: avalanche rates SSW vs control
for regime in ['cold-dry', 'cold-wet', 'warm-dry', 'warm-wet']:
    subset = merged2[merged2['regime'] == regime]
    ssw_rate = subset.loc[subset['in_ssw'], 'dry_natural_size_1234'].mean()
    ctrl_rate = subset.loc[~subset['in_ssw'], 'dry_natural_size_1234'].mean()
    rr = ssw_rate / ctrl_rate if ctrl_rate > 0 else np.nan
    print(f"  {regime:<12} Within-regime RR: {rr:.3f} (SSW: {ssw_rate:.2f}, Ctrl: {ctrl_rate:.2f})")

# ============================================================
# SECTION 6: TRIGGER-DAY VS NO-TRIGGER-DAY AVALANCHE RATES
# ============================================================
print("\n" + "=" * 70)
print("SECTION 6: AVALANCHE RATES ON TRIGGER VS NO-TRIGGER DAYS")
print("=" * 70)

# Define "trigger day" = surface warming OR rain-on-snow at >10% of stations
merged['trigger_day'] = (merged['frac_surface_warm'] > 0.1) | (merged['frac_rain'] > 0.1)

trig_rate = merged.loc[merged['trigger_day'], 'dry_natural_size_1234'].mean()
no_trig_rate = merged.loc[~merged['trigger_day'], 'dry_natural_size_1234'].mean()
print(f"Trigger-day avalanche rate: {trig_rate:.3f}")
print(f"No-trigger-day avalanche rate: {no_trig_rate:.3f}")
print(f"Rate ratio (trigger/no-trigger): {trig_rate/no_trig_rate:.2f}")

# SSW effect on trigger-day frequency
ssw_trig_frac = merged.loc[merged['in_ssw'], 'trigger_day'].mean()
ctrl_trig_frac = merged.loc[~merged['in_ssw'], 'trigger_day'].mean()
print(f"\nTrigger-day frequency: SSW={ssw_trig_frac:.1%}, Control={ctrl_trig_frac:.1%}")
print(f"Reduction: {(1 - ssw_trig_frac/ctrl_trig_frac)*100:.1f}%")

# Expected avalanche change from trigger-day frequency change alone
expected_rate_ssw = ssw_trig_frac * trig_rate + (1 - ssw_trig_frac) * no_trig_rate
expected_rate_ctrl = ctrl_trig_frac * trig_rate + (1 - ctrl_trig_frac) * no_trig_rate
expected_rr2 = expected_rate_ssw / expected_rate_ctrl
print(f"\nExpected RR from trigger-day frequency change: {expected_rr2:.3f}")
print(f"Observed RR: {ssw_aval.mean()/ctrl_aval.mean():.3f}")

# ============================================================
# SECTION 7: HIGH-ISWR DAYS ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("SECTION 7: SOLAR RADIATION AND AVALANCHE RATES")
print("=" * 70)

# Tercile split on ISWR
iswr_terciles = merged['mean_iswr'].quantile([0.33, 0.67]).values
merged['iswr_cat'] = 'low'
merged.loc[merged['mean_iswr'] >= iswr_terciles[0], 'iswr_cat'] = 'mid'
merged.loc[merged['mean_iswr'] >= iswr_terciles[1], 'iswr_cat'] = 'high'

for cat in ['low', 'mid', 'high']:
    subset = merged[merged['iswr_cat'] == cat]
    rate = subset['dry_natural_size_1234'].mean()
    ssw_frac = subset['in_ssw'].mean()
    print(f"  ISWR {cat:<5}: avalanche rate = {rate:.3f}, SSW fraction = {ssw_frac:.1%}")

# SSW shifts days toward low-ISWR
print(f"\nISWR category distribution:")
for cat in ['low', 'mid', 'high']:
    ssw_frac = (merged.loc[merged['in_ssw'], 'iswr_cat'] == cat).mean()
    ctrl_frac = (merged.loc[~merged['in_ssw'], 'iswr_cat'] == cat).mean()
    print(f"  {cat:<5}: SSW={ssw_frac:.1%}, Ctrl={ctrl_frac:.1%}, Shift={ssw_frac-ctrl_frac:+.1%}")

# ============================================================
# SAVE RESULTS
# ============================================================
results = {
    'trigger_avalanche_correlations': {},
    'expected_vs_observed_rr': {
        'expected_rr_multivariate': float(expected_rr),
        'observed_rr': float(ssw_aval.mean() / ctrl_aval.mean()),
    },
    'trigger_day_analysis': {
        'trigger_day_rate': float(trig_rate),
        'no_trigger_day_rate': float(no_trig_rate),
        'trigger_nontrigger_ratio': float(trig_rate / no_trig_rate),
        'ssw_trigger_frac': float(ssw_trig_frac),
        'ctrl_trigger_frac': float(ctrl_trig_frac),
        'expected_rr_trigger_frequency': float(expected_rr2),
    }
}

for name, vals in predictors.items():
    mask = np.isfinite(vals) & np.isfinite(aval)
    r, p = stats.pearsonr(vals[mask], aval[mask])
    results['trigger_avalanche_correlations'][name] = {'r': float(r), 'p': float(p), 'r2': float(r**2)}

with open('data/results/r32_mediation.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\n" + "=" * 70)
print("RESULTS SAVED to data/results/r32_mediation.json")
print("=" * 70)
