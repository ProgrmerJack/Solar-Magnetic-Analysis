import pandas as pd
import numpy as np
from scipy import stats

era5 = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\era5_swiss_alps_daily.parquet')
panel = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\analysis_panel_v2.parquet')
ssw_cat = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\atmospheric\ssw_catalog.parquet')

ssw_cat.index = ssw_cat.index.tz_localize(None)
aval = panel['dry_natural_size_1234']
aval.index = pd.to_datetime(aval.index)

# SSW events in ERA5 period (2004-2013)
era5_events = ssw_cat[(ssw_cat.index >= '2004-01-01') & (ssw_cat.index <= '2013-12-31')]
print(f'SSW events in ERA5 period: {len(era5_events)}')
for d in era5_events.index:
    print(f'  {d.strftime("%Y-%m-%d")}')

print('\n=== ERA5 REGIONAL RESPONSE DURING SSW EVENTS ===')
variables = ['tp_mm', 'sf_mm', 't2m_K', 'wind_speed', 'sd_m']
var_labels = ['Total precip', 'Snowfall', 'Temperature', 'Wind speed', 'Snow depth']

for var, label in zip(variables, var_labels):
    pre_vals = []
    post_vals = []
    for onset in era5_events.index:
        onset_ts = pd.Timestamp(onset)
        pre = era5.loc[onset_ts - pd.Timedelta(days=30):onset_ts - pd.Timedelta(days=1), var]
        post = era5.loc[onset_ts - pd.Timedelta(days=15):onset_ts + pd.Timedelta(days=15), var]
        if len(pre) > 10 and len(post) > 10:
            pre_vals.append(pre.mean())
            post_vals.append(post.mean())
    
    pre_arr = np.array(pre_vals)
    post_arr = np.array(post_vals)
    n = min(len(pre_arr), len(post_arr))
    if n > 3:
        t, p = stats.ttest_rel(pre_arr[:n], post_arr[:n])
        change_pct = 100 * (post_arr[:n].mean() - pre_arr[:n].mean()) / pre_arr[:n].mean() if pre_arr[:n].mean() != 0 else 0
        print(f'{label:15s}: pre={pre_arr[:n].mean():.3f}, SSW_window={post_arr[:n].mean():.3f}, change={change_pct:+.1f}%, P={p:.4f}')

# Use anomaly columns for deseasonalized comparison
print('\n=== DESEASONALIZED ANOMALIES DURING SSW ===')
anom_vars = ['tp_mm_anom', 'sf_mm_anom', 't2m_K_anom', 'wind_speed_anom']
anom_labels = ['Precip anomaly', 'Snowfall anomaly', 'Temp anomaly', 'Wind anomaly']

for var, label in zip(anom_vars, anom_labels):
    ssw_anoms = []
    for onset in era5_events.index:
        onset_ts = pd.Timestamp(onset)
        window = era5.loc[onset_ts - pd.Timedelta(days=15):onset_ts + pd.Timedelta(days=15), var]
        if len(window) > 10:
            ssw_anoms.append(window.mean())
    
    anoms = np.array(ssw_anoms)
    anoms = anoms[np.isfinite(anoms)]
    if len(anoms) > 3:
        t, p = stats.ttest_1samp(anoms, 0)
        direction = '+' if anoms.mean() > 0 else '-'
        sign_count = (anoms > 0).sum() if anoms.mean() > 0 else (anoms < 0).sum()
        print(f'{label:20s}: mean={anoms.mean():+.4f}, {direction} in {sign_count}/{len(anoms)} events, t={t:.3f}, P={p:.4f}')

# Key test: per-event snowfall anomaly vs avalanche ratio
print('\n=== PER-EVENT: Snowfall anomaly vs Avalanche ratio ===')
sf_anoms = []
aval_ratios = []
for onset in era5_events.index:
    onset_ts = pd.Timestamp(onset)
    window = pd.date_range(onset_ts - pd.Timedelta(days=15), onset_ts + pd.Timedelta(days=15))
    
    sf = era5.loc[era5.index.isin(window), 'sf_mm_anom']
    a = aval.reindex(window)
    
    doy = onset_ts.timetuple().tm_yday
    doys = set()
    for d in range(-15, 16):
        doys.add((doy + d) % 366)
    clim = aval[aval.index.dayofyear.isin(doys)].mean()
    
    if len(sf) > 10 and a.notna().sum() > 10 and clim > 0:
        sf_anoms.append(sf.mean())
        aval_ratios.append(a.mean() / clim)
        print(f'{onset_ts.strftime("%Y-%m-%d")}: sf_anom={sf.mean():+.4f}, aval_ratio={a.mean()/clim:.3f}')

r, p = stats.spearmanr(sf_anoms, aval_ratios)
print(f'\nSnowfall anomaly vs avalanche ratio: r={r:.3f}, P={p:.4f} (n={len(sf_anoms)})')
