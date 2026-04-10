"""Test if ERA5 temperature anomaly exists BEFORE SSW onset - confirming common-cause."""
import pandas as pd, numpy as np
from scipy import stats

era5 = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\era5_swiss_alps_daily.parquet')
era5.index = pd.to_datetime(era5.index)
ssw = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\atmospheric\ssw_catalog.parquet')
ssw.index = pd.to_datetime(ssw.index).tz_localize(None)

era5_w = era5[(era5.index.month >= 11) | (era5.index.month <= 4)].copy()
era5_ssw = ssw[(ssw.index >= '2004-01-01') & (ssw.index <= '2013-12-31')]
print(f'SSW events: {len(era5_ssw)}')

# Non-SSW control
all_ssw_mask = pd.Series(False, index=era5_w.index)
for d in era5_ssw.index:
    win = (era5_w.index >= d - pd.Timedelta(days=30)) & (era5_w.index <= d + pd.Timedelta(days=30))
    all_ssw_mask = all_ssw_mask | win
ctrl = era5_w[~all_ssw_mask]

print('\n=== TEMPERATURE ANOMALY BY WINDOW RELATIVE TO SSW ONSET ===')
windows = [
    ('Pre-30 to Pre-16', -30, -16),
    ('Pre-15 to Pre-1', -15, -1),
    ('Post-0 to Post-15', 0, 15),
    ('Post-16 to Post-30', 16, 30),
    ('Full +/-15d', -15, 15),
]

for label, start, end in windows:
    mask = pd.Series(False, index=era5_w.index)
    for d in era5_ssw.index:
        win = (era5_w.index >= d + pd.Timedelta(days=start)) & (era5_w.index <= d + pd.Timedelta(days=end))
        mask = mask | win
    ssw_subset = era5_w[mask]
    
    s = ssw_subset['t2m_K_anom'].dropna()
    c = ctrl['t2m_K_anom'].dropna()
    diff = s.mean() - c.mean()
    _, p = stats.ttest_ind(s, c)
    _, p_mw = stats.mannwhitneyu(s, c, alternative='two-sided')
    print(f'{label:25s}: mean_anom={s.mean():+.3f} K  ctrl={c.mean():+.3f} K  diff={diff:+.3f} K  t-P={p:.4f}  MW-P={p_mw:.4f}  n={len(s)}')

# Also check wind direction by window
print('\n=== U10 (ZONAL WIND) BY WINDOW ===')
for label, start, end in windows:
    mask = pd.Series(False, index=era5_w.index)
    for d in era5_ssw.index:
        win = (era5_w.index >= d + pd.Timedelta(days=start)) & (era5_w.index <= d + pd.Timedelta(days=end))
        mask = mask | win
    ssw_subset = era5_w[mask]
    
    s = ssw_subset['u10'].dropna()
    c = ctrl['u10'].dropna()
    diff = s.mean() - c.mean()
    _, p = stats.ttest_ind(s, c)
    print(f'{label:25s}: u10_ssw={s.mean():+.3f}  u10_ctrl={c.mean():+.3f}  diff={diff:+.3f}  P={p:.4f}  n={len(s)}')

# Composite time series at daily lag
print('\n=== COMPOSITE T2m ANOMALY AT DAILY LAGS ===')
print('Lag(d)  mean_T2m_anom(K)  n_days  P_value')
for lag in range(-20, 21, 5):
    temps = []
    for d in era5_ssw.index:
        target = d + pd.Timedelta(days=lag)
        if target in era5.index:
            val = era5.loc[target, 't2m_K_anom']
            if pd.notna(val):
                temps.append(val)
    if len(temps) >= 3:
        mean_t = np.mean(temps)
        _, p = stats.ttest_1samp(temps, 0)
        print(f'{lag:+4d}    {mean_t:+.3f}              {len(temps)}       {p:.4f}')
