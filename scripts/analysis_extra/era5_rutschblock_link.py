"""Test if ERA5 temperature anomaly predicts Rutschblock stability at event level."""
import pandas as pd, numpy as np
from scipy import stats

era5 = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\era5_swiss_alps_daily.parquet')
ssw = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\atmospheric\ssw_catalog.parquet')
ssw.index = pd.to_datetime(ssw.index).tz_localize(None)
rb = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\cryosphere\slf_stability-tests-avalanche_Rutschblock_data_(Switzerland).parquet')
rb.index = pd.to_datetime(rb.index).tz_localize(None)
rb = rb[np.isfinite(rb['stabclass'])]

# ERA5 period SSW events
era5_ssw = ssw[(ssw.index >= '2004-01-01') & (ssw.index <= '2013-12-31')]
# Rutschblock period starts 2001 but ERA5 starts 2004, so overlap is 2004-2013
print(f'SSW events in ERA5+Rutschblock overlap: {len(era5_ssw)}')

results = []
for d in era5_ssw.index:
    win_s = d - pd.Timedelta(days=15)
    win_e = d + pd.Timedelta(days=15)
    
    # ERA5 temperature anomaly in window
    era5_win = era5[(era5.index >= win_s) & (era5.index <= win_e)]
    if len(era5_win) == 0:
        continue
    t2m_anom = era5_win['t2m_K_anom'].mean()
    ws_anom = era5_win['wind_speed_anom'].mean()
    
    # Rutschblock stability in window
    rb_ssw = rb[(rb.index >= win_s) & (rb.index <= win_e)]
    if len(rb_ssw) < 5:
        continue
    
    # Control Rutschblock: same month, other years
    month = d.month
    ctrl_mask = ((rb.index.month == month) | (rb.index.month == (month % 12 + 1)))
    ctrl_mask = ctrl_mask & ~((rb.index >= win_s) & (rb.index <= win_e))
    rb_ctrl = rb[ctrl_mask]
    if len(rb_ctrl) < 10:
        continue
    
    stab_shift = rb_ssw['stabclass'].mean() - rb_ctrl['stabclass'].mean()
    results.append({
        'date': d.date(), 't2m_anom': t2m_anom, 'ws_anom': ws_anom,
        'stab_shift': stab_shift, 'n_rb': len(rb_ssw)
    })

df = pd.DataFrame(results)
print(f'Events with both ERA5 + Rutschblock: {len(df)}')
print(df.to_string(index=False))

# Temperature anomaly -> Rutschblock shift
r, p = stats.spearmanr(df['t2m_anom'], df['stab_shift'])
print(f'\nT2m anomaly -> Rutschblock stability shift:')
print(f'  Spearman r = {r:.3f}, P = {p:.4f}')

r2, p2 = stats.pearsonr(df['t2m_anom'], df['stab_shift'])
print(f'  Pearson r = {r2:.3f}, P = {p2:.4f}')

# Wind anomaly -> Rutschblock shift
r3, p3 = stats.spearmanr(df['ws_anom'], df['stab_shift'])
print(f'\nWind speed anomaly -> Rutschblock stability shift:')
print(f'  Spearman r = {r3:.3f}, P = {p3:.4f}')

# Interpretation
print(f'\nPositive t2m_anom + positive stab_shift = warmer -> more stable')
print(f'Events with both warmer + more stable: '
      f'{((df["t2m_anom"] > 0) & (df["stab_shift"] > 0)).sum()}/{len(df)}')
