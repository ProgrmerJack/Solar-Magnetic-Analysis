"""ERA5 wind mechanism analysis for SSW-avalanche pathway."""
import pandas as pd, numpy as np
from scipy import stats

era5 = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\era5_swiss_alps_daily.parquet')
ssw = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\atmospheric\ssw_catalog.parquet')
ssw.index = pd.to_datetime(ssw.index).tz_localize(None)

era5_w = era5[(era5.index.month >= 11) | (era5.index.month <= 4)].copy()
era5_w.index = pd.to_datetime(era5_w.index)

era5_ssw = ssw[(ssw.index >= '2004-01-01') & (ssw.index <= '2013-12-31')]
print(f'SSW events in ERA5 period: {len(era5_ssw)}')
for d in era5_ssw.index:
    print(f'  {d.date()}')

ssw_mask = pd.Series(False, index=era5_w.index)
for d in era5_ssw.index:
    win = (era5_w.index >= d - pd.Timedelta(days=15)) & (era5_w.index <= d + pd.Timedelta(days=15))
    ssw_mask = ssw_mask | win

ssw_days = era5_w[ssw_mask]
ctrl_days = era5_w[~ssw_mask]
print(f'SSW days: {len(ssw_days)}, Control days: {len(ctrl_days)}')

print('\n=== ERA5 MECHANISM ANALYSIS (raw values) ===')
for var in ['wind_speed', 'u10', 'v10', 't2m_K', 'tp_mm', 'sf_mm', 'sd_m']:
    s = ssw_days[var].dropna()
    c = ctrl_days[var].dropna()
    diff = s.mean() - c.mean()
    d_cohen = diff / c.std() if c.std() > 0 else 0
    _, p_t = stats.ttest_ind(s, c)
    _, p_mw = stats.mannwhitneyu(s, c, alternative='two-sided')
    print(f'{var:15s}: SSW={s.mean():.4f}  Ctrl={c.mean():.4f}  Diff={diff:+.4f}  d={d_cohen:+.3f}  t-P={p_t:.4f}  MW-P={p_mw:.4f}')

print('\n=== DESEASONALIZED ANOMALIES ===')
for var in ['wind_speed_anom', 't2m_K_anom', 'sf_mm_anom', 'tp_mm_anom']:
    s = ssw_days[var].dropna()
    c = ctrl_days[var].dropna()
    diff = s.mean() - c.mean()
    _, p_t = stats.ttest_ind(s, c)
    print(f'{var:20s}: SSW={s.mean():+.4f}  Ctrl={c.mean():+.4f}  Diff={diff:+.4f}  P={p_t:.4f}')

# Wind direction
ssw_dir = np.degrees(np.arctan2(ssw_days['v10'].mean(), ssw_days['u10'].mean())) % 360
ctrl_dir = np.degrees(np.arctan2(ctrl_days['v10'].mean(), ctrl_days['u10'].mean())) % 360
print(f'\nWind direction: SSW={ssw_dir:.1f}deg  Ctrl={ctrl_dir:.1f}deg  Shift={ssw_dir-ctrl_dir:+.1f}deg')

# Matched-control analysis (same DOY, other years)
print('\n=== MATCHED-CONTROL WIND ANALYSIS (event-level) ===')
event_diffs = []
for d in era5_ssw.index:
    win_s = d - pd.Timedelta(days=15)
    win_e = d + pd.Timedelta(days=15)
    ssw_data = era5_w[(era5_w.index >= win_s) & (era5_w.index <= win_e)]
    
    # Control: same DOY range, other years
    ctrl_data = []
    for yr_offset in range(-5, 6):
        if yr_offset == 0:
            continue
        cs = win_s + pd.DateOffset(years=yr_offset)
        ce = win_e + pd.DateOffset(years=yr_offset)
        c = era5_w[(era5_w.index >= cs) & (era5_w.index <= ce)]
        if len(c) > 0:
            ctrl_data.append(c)
    if not ctrl_data:
        continue
    ctrl_all = pd.concat(ctrl_data)
    
    ws_diff = ssw_data['wind_speed'].mean() - ctrl_all['wind_speed'].mean()
    t2m_diff = ssw_data['t2m_K'].mean() - ctrl_all['t2m_K'].mean()
    sf_diff = ssw_data['sf_mm'].mean() - ctrl_all['sf_mm'].mean()
    event_diffs.append({
        'date': d.date(), 'wind_diff': ws_diff, 't2m_diff': t2m_diff, 'sf_diff': sf_diff,
        'ws_ssw': ssw_data['wind_speed'].mean(), 'ws_ctrl': ctrl_all['wind_speed'].mean()
    })
    print(f'  {d.date()}: wind={ws_diff:+.3f} m/s  T2m={t2m_diff:+.2f} K  SF={sf_diff:+.4f} mm')

edf = pd.DataFrame(event_diffs)
print(f'\nMatched event-level wind speed diff:')
print(f'  Mean: {edf["wind_diff"].mean():+.4f} m/s')
print(f'  Positive (SSW windier): {(edf["wind_diff"] > 0).sum()}/{len(edf)}')
t, p = stats.ttest_1samp(edf['wind_diff'], 0)
print(f'  Paired t-test: P={p:.4f}')
_, p_w = stats.wilcoxon(edf['wind_diff'])
print(f'  Wilcoxon: P={p_w:.4f}')

print(f'\nMatched event-level T2m diff:')
print(f'  Mean: {edf["t2m_diff"].mean():+.4f} K')
print(f'  Positive (SSW warmer): {(edf["t2m_diff"] > 0).sum()}/{len(edf)}')
t, p = stats.ttest_1samp(edf['t2m_diff'], 0)
print(f'  Paired t-test: P={p:.4f}')

print(f'\nMatched event-level snowfall diff:')
print(f'  Mean: {edf["sf_diff"].mean():+.4f} mm')
print(f'  Positive (SSW more snow): {(edf["sf_diff"] > 0).sum()}/{len(edf)}')
t, p = stats.ttest_1samp(edf['sf_diff'], 0)
print(f'  Paired t-test: P={p:.4f}')

# Additional: wind speed variability (std)
print('\n=== WIND VARIABILITY ===')
ws_ssw_std = ssw_days['wind_speed'].std()
ws_ctrl_std = ctrl_days['wind_speed'].std()
print(f'Wind speed std: SSW={ws_ssw_std:.4f}  Ctrl={ws_ctrl_std:.4f}')
_, p = stats.levene(ssw_days['wind_speed'].dropna(), ctrl_days['wind_speed'].dropna())
print(f'Levene test for equal variance: P={p:.4f}')

# High-wind day frequency
high_wind_thresh = ctrl_days['wind_speed'].quantile(0.90)
ssw_high = (ssw_days['wind_speed'] > high_wind_thresh).mean()
ctrl_high = (ctrl_days['wind_speed'] > high_wind_thresh).mean()
print(f'Fraction of high-wind days (>P90={high_wind_thresh:.2f}): SSW={ssw_high:.3f}  Ctrl={ctrl_high:.3f}')
print(f'Ratio: {ssw_high/ctrl_high:.3f}')
