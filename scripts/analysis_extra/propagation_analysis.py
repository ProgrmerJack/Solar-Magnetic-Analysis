import pandas as pd
import numpy as np
from scipy import stats

strat = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\atmospheric\ncep_stratosphere.parquet')
trop = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\atmospheric\ncep_troposphere.parquet')
panel = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\analysis_panel_v2.parquet')
ssw_cat = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\atmospheric\ssw_catalog.parquet')

strat.index = strat.index.tz_localize(None)
trop.index = trop.index.tz_localize(None)
ssw_cat.index = ssw_cat.index.tz_localize(None)

# Compute NAM-like index at each level: standardized zonal wind anomaly
# This shows the downward propagation of stratospheric signals

levels = [10, 20, 30, 50, 70, 100]
u_cols = [f'uwnd_ms_{l}hPa' for l in levels]

# Show SSW composite at each level: zonal wind anomaly from day -30 to +60
swiss_events = ssw_cat[(ssw_cat.index >= '1998-11-01') & (ssw_cat.index <= '2019-04-30')].index

print('=== MULTI-LEVEL PROPAGATION COMPOSITE (n=16 events) ===')
print(f'{"Day":>5s}', end='')
for l in levels:
    print(f' {l:>7d}hPa', end='')
print()

for day_offset in range(-15, 61, 5):
    means = []
    for l in levels:
        col = f'uwnd_ms_{l}hPa'
        vals = []
        for onset in swiss_events:
            target = onset + pd.Timedelta(days=day_offset)
            if target in strat.index:
                # Compute anomaly: subtract climatological value for this DOY
                doy = target.timetuple().tm_yday
                clim = strat[col][strat.index.dayofyear == doy].mean()
                vals.append(strat.loc[target, col] - clim)
        means.append(np.nanmean(vals) if vals else np.nan)
    print(f'{day_offset:+5d}', end='')
    for m in means:
        print(f' {m:+8.1f}', end='')
    print()

# Now test: does the tropospheric signal predict avalanches?
# Use SLP and 500hPa height as tropospheric proxies
print()
print('=== TROPOSPHERIC RESPONSE AROUND SSW ===')
print('Columns in troposphere:', list(trop.columns))

# Test tropospheric variables around SSW events
aval = panel['dry_natural_size_1234']
aval.index = pd.to_datetime(aval.index)

for col in trop.columns:
    pre_vals = []
    post_vals = []
    for onset in swiss_events:
        onset_ts = pd.Timestamp(onset)
        pre_window = pd.date_range(onset_ts - pd.Timedelta(days=30), onset_ts - pd.Timedelta(days=1))
        post_window = pd.date_range(onset_ts, onset_ts + pd.Timedelta(days=30))
        pre_vals.append(trop[col].reindex(pre_window).mean())
        post_vals.append(trop[col].reindex(post_window).mean())
    
    pre_arr = np.array([v for v in pre_vals if not np.isnan(v)])
    post_arr = np.array([v for v in post_vals if not np.isnan(v)])
    if len(pre_arr) > 3 and len(post_arr) > 3:
        t, p = stats.ttest_rel(pre_arr[:min(len(pre_arr),len(post_arr))], 
                               post_arr[:min(len(pre_arr),len(post_arr))])
        print(f'{col}: pre={np.mean(pre_arr):.2f}, post={np.mean(post_arr):.2f}, change={np.mean(post_arr)-np.mean(pre_arr):+.2f}, paired-t P={p:.4f}')

# Also check: AO/NAO indices if available
ao_cols = [c for c in panel.columns if 'ao' in c.lower() or 'nao' in c.lower() or 'annular' in c.lower()]
print(f'\nAO/NAO columns in panel: {ao_cols}')
