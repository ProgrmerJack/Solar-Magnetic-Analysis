"""
Lagged temperature-avalanche analysis within SSW windows.
Tests: does warmth on day N predict fewer avalanches on day N+2 to N+5?
This directly tests the cumulative sintering hypothesis.
"""
import pandas as pd
import numpy as np
from scipy import stats

# Load ERA5 and panel
era5 = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\era5_swiss_alps_daily.parquet')
panel = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\analysis_panel_v2.parquet')
ssw = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\atmospheric\ssw_catalog.parquet')

# Align dates
era5.index = pd.to_datetime(era5.index)
panel_dates = panel.index
ssw_dates = ssw.index.tz_localize(None)

# Get SSW windows in ERA5 period
era5_ssw_mask = (ssw_dates >= era5.index.min()) & (ssw_dates <= era5.index.max())
era5_ssw_dates = ssw_dates[era5_ssw_mask]

# Merge ERA5 temperature with panel avalanche data
common_dates = era5.index.intersection(panel_dates)
merged = pd.DataFrame({
    't2m_anom': era5.loc[common_dates, 't2m_K_anom'] if 't2m_K_anom' in era5.columns else era5.loc[common_dates, 't2m_K'] - era5.loc[common_dates, 't2m_K'].mean(),
    'aval': panel.loc[common_dates, 'dry_natural_size_1234']
})

# Mark SSW window days
in_ssw = np.zeros(len(merged), dtype=bool)
for d in era5_ssw_dates:
    mask = (merged.index >= d - pd.Timedelta(days=15)) & (merged.index <= d + pd.Timedelta(days=15))
    in_ssw |= mask

ssw_data = merged[in_ssw].copy()
ctrl_data = merged[~in_ssw].copy()

print(f"SSW window days: {len(ssw_data)}, Control days: {len(ctrl_data)}")

# Lagged correlation within SSW windows
print("\n=== LAGGED T2m_anom → Avalanche correlation (SSW windows only) ===")
print(f"{'Lag(d)':<10} {'r':>8} {'P':>10} {'n':>6} {'Interpretation'}")
print("-" * 60)
for lag in range(0, 11):
    if lag == 0:
        r, p = stats.spearmanr(ssw_data['t2m_anom'], ssw_data['aval'])
        n = len(ssw_data)
    else:
        # Temperature on day N vs avalanches on day N+lag
        temp_shifted = ssw_data['t2m_anom'].shift(lag).dropna()
        aval_aligned = ssw_data['aval'].loc[temp_shifted.index]
        valid = temp_shifted.notna() & aval_aligned.notna()
        r, p = stats.spearmanr(temp_shifted[valid], aval_aligned[valid])
        n = valid.sum()
    
    interpretation = ""
    if p < 0.05 and r < 0:
        interpretation = "*** SIGNIFICANT: warmth → fewer avalanches"
    elif p < 0.05 and r > 0:
        interpretation = "warmth → more avalanches"
    elif p < 0.10:
        interpretation = "marginal"
    
    print(f"  {lag:<8} {r:>+8.3f} {p:>10.4f} {n:>6} {interpretation}")

# Cumulative temperature: 3-day and 5-day rolling mean
print("\n=== CUMULATIVE T2m_anom (rolling mean) → Avalanche ===")
for window in [3, 5, 7]:
    ssw_data_copy = ssw_data.copy()
    ssw_data_copy[f't2m_roll{window}'] = ssw_data_copy['t2m_anom'].rolling(window, min_periods=window).mean()
    valid = ssw_data_copy[f't2m_roll{window}'].notna() & ssw_data_copy['aval'].notna()
    r, p = stats.spearmanr(ssw_data_copy.loc[valid, f't2m_roll{window}'], ssw_data_copy.loc[valid, 'aval'])
    print(f"  {window}-day rolling mean: r={r:+.3f}, P={p:.4f}, n={valid.sum()}")

# Also test: cumulative T leading avalanches by 2 days
print("\n=== CUMULATIVE T2m_anom (5-day rolling) leading avalanche by 2 days ===")
ssw_data_copy = ssw_data.copy()
ssw_data_copy['t2m_roll5'] = ssw_data_copy['t2m_anom'].rolling(5, min_periods=5).mean()
ssw_data_copy['t2m_roll5_lag2'] = ssw_data_copy['t2m_roll5'].shift(2)
valid = ssw_data_copy['t2m_roll5_lag2'].notna() & ssw_data_copy['aval'].notna()
if valid.sum() > 10:
    r, p = stats.spearmanr(ssw_data_copy.loc[valid, 't2m_roll5_lag2'], ssw_data_copy.loc[valid, 'aval'])
    print(f"  r={r:+.3f}, P={p:.4f}, n={valid.sum()}")
else:
    print(f"  Insufficient data: n={valid.sum()}")

# Compare: same analysis on control days
print("\n=== CONTROL: Lagged T2m_anom → Avalanche (non-SSW days) ===")
for lag in [0, 2, 5]:
    if lag == 0:
        r, p = stats.spearmanr(ctrl_data['t2m_anom'], ctrl_data['aval'])
        n = len(ctrl_data)
    else:
        temp_shifted = ctrl_data['t2m_anom'].shift(lag).dropna()
        aval_aligned = ctrl_data['aval'].loc[temp_shifted.index]
        valid = temp_shifted.notna() & aval_aligned.notna()
        r, p = stats.spearmanr(temp_shifted[valid], aval_aligned[valid])
        n = valid.sum()
    print(f"  Lag {lag}d: r={r:+.3f}, P={p:.4f}, n={n}")
