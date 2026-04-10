import pandas as pd, numpy as np
from scipy import stats

era5 = pd.read_parquet('data/processed/era5_swiss_alps_daily.parquet')
try:
    era5_ext = pd.read_parquet('data/processed/era5_swiss_alps_extended.parquet')
    era5 = pd.concat([era5, era5_ext]).sort_index()
    era5 = era5[~era5.index.duplicated(keep='last')]
except: pass

ssw_cat = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw_dates = ssw_cat.index.tz_localize(None)
ssw_in = ssw_dates[(ssw_dates >= era5.index.min()) & (ssw_dates <= era5.index.max())]

winter = era5[era5['doy'].isin(list(range(1,121)) + list(range(305,366)))].copy()

def ssw_mask(idx, dates, w=15):
    m = pd.Series(False, index=idx)
    for d in dates:
        m |= (idx >= d - pd.Timedelta(days=w)) & (idx <= d + pd.Timedelta(days=w))
    return m

sm = ssw_mask(winter.index, ssw_in)

winter['rain_mm'] = (winter['tp_mm'] - winter['sf_mm']).clip(lower=0)
winter['snow_frac'] = winter['sf_mm'] / winter['tp_mm'].clip(lower=0.01)

print('=== Precipitation Partitioning During SSW ===')
for col, name in [('tp_mm', 'Total precip'), ('sf_mm', 'Snowfall'), ('rain_mm', 'Rainfall'), ('snow_frac', 'Snow fraction')]:
    ssw_v = winter.loc[sm, col].dropna()
    ctrl_v = winter.loc[~sm, col].dropna()
    d = (ssw_v.mean() - ctrl_v.mean()) / ctrl_v.std()
    _, p = stats.mannwhitneyu(ssw_v, ctrl_v, alternative='two-sided')
    pct = (ssw_v.mean() - ctrl_v.mean()) / ctrl_v.mean() * 100
    print(f'{name:20s}: SSW={ssw_v.mean():.4f} Ctrl={ctrl_v.mean():.4f} chg={pct:+.1f}% d={d:+.3f} P={p:.2e}')

# Wind direction sectors
winter['wind_dir'] = np.degrees(np.arctan2(-winter['u10'], -winter['v10'])) % 360
bins = [0, 45, 90, 135, 180, 225, 270, 315, 360]
labels = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
winter['sector'] = pd.cut(winter['wind_dir'], bins=bins, labels=labels, include_lowest=True)

ssw_sect = winter.loc[sm, 'sector'].value_counts(normalize=True).sort_index()
ctrl_sect = winter.loc[~sm, 'sector'].value_counts(normalize=True).sort_index()

print('\n=== Wind Direction Shift ===')
print(f'{"Sector":>8} {"SSW%":>8} {"Ctrl%":>8} {"Delta":>8}')
for s in labels:
    sv = ssw_sect.get(s, 0) * 100
    cv = ctrl_sect.get(s, 0) * 100
    print(f'{s:>8} {sv:>8.1f} {cv:>8.1f} {sv-cv:>+8.1f}')

# Easterly component enhancement (continental vs maritime)
winter['u_east'] = -winter['u10']  # positive = from east (continental)
ssw_e = winter.loc[sm, 'u_east'].mean()
ctrl_e = winter.loc[~sm, 'u_east'].mean()
d_e = (ssw_e - ctrl_e) / winter.loc[~sm, 'u_east'].std()
_, p_e = stats.mannwhitneyu(winter.loc[sm, 'u_east'], winter.loc[~sm, 'u_east'])
print(f'\nEasterly component: SSW={ssw_e:.4f} Ctrl={ctrl_e:.4f} d={d_e:.3f} P={p_e:.2e}')
print('Positive = more continental (easterly), negative = more maritime (westerly)')
