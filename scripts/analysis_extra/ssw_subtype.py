import pandas as pd
import numpy as np
from scipy import stats

ssw = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\atmospheric\ssw_catalog.parquet')
strat = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\atmospheric\ncep_stratosphere.parquet')
panel = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\analysis_panel_v2.parquet')

strat.index = strat.index.tz_localize(None)
ssw.index = ssw.index.tz_localize(None)

swiss_ssw = ssw[(ssw.index >= '1998-11-01') & (ssw.index <= '2019-04-30')]
u10 = strat['uwnd_ms_10hPa']
u100 = strat['uwnd_ms_100hPa']
aval = panel['dry_natural_size_1234']
aval.index = pd.to_datetime(aval.index)

events = []
for onset in swiss_ssw.index:
    onset_ts = pd.Timestamp(onset)
    pre_window = pd.date_range(onset_ts - pd.Timedelta(days=30), onset_ts - pd.Timedelta(days=10))
    pre_u10 = u10.reindex(pre_window).mean()
    post_window = pd.date_range(onset_ts, onset_ts + pd.Timedelta(days=30))
    min_u10 = u10.reindex(post_window).min()
    rev_days = int((u10.reindex(post_window) < 0).sum())
    decel = pre_u10 - min_u10
    
    ssw_window = pd.date_range(onset_ts - pd.Timedelta(days=15), onset_ts + pd.Timedelta(days=15))
    aval_ssw = aval.reindex(ssw_window).mean()
    doy_center = onset_ts.timetuple().tm_yday
    doy_range = set()
    for d in range(-15, 16):
        doy_range.add((doy_center + d) % 366)
    clim = aval[aval.index.dayofyear.isin(doy_range)].mean()
    ratio = aval_ssw / clim if clim > 0 else np.nan
    
    pre_u100 = u100.reindex(pre_window).mean()
    post_u100_window = pd.date_range(onset_ts + pd.Timedelta(days=10), onset_ts + pd.Timedelta(days=60))
    min_u100 = u100.reindex(post_u100_window).min()
    surface_prop = pre_u100 - min_u100 if not pd.isna(min_u100) else np.nan
    
    events.append({
        'onset': onset_ts.strftime('%Y-%m-%d'),
        'type': swiss_ssw.loc[onset, 'type'],
        'decel': decel, 'min_u10': min_u10,
        'rev_days': rev_days, 'surface_prop': surface_prop,
        'aval_ratio': ratio
    })

df = pd.DataFrame(events)
print('=== SSW Events with Avalanche Response ===')
for _, r in df.iterrows():
    marker = 'v' if r['aval_ratio'] < 1 else '^'
    print(f"{r['onset']} type={r['type']} decel={r['decel']:5.1f} min_u={r['min_u10']:6.1f} rev={r['rev_days']:2d}d surf={r['surface_prop']:5.1f} aval_ratio={r['aval_ratio']:5.3f} {marker}")

print()
vals = df.dropna(subset=['decel','aval_ratio'])
r1, p1 = stats.spearmanr(vals['decel'], vals['aval_ratio'])
print(f'Deceleration vs aval_ratio: r={r1:.3f}, P={p1:.4f}')
r2, p2 = stats.spearmanr(vals['rev_days'], vals['aval_ratio'])
print(f'Reversal duration vs aval_ratio: r={r2:.3f}, P={p2:.4f}')
r3, p3 = stats.spearmanr(vals['surface_prop'], vals['aval_ratio'])
print(f'Surface propagation vs aval_ratio: r={r3:.3f}, P={p3:.4f}')

for t in df['type'].unique():
    sub = df[df['type'] == t]
    below = (sub['aval_ratio'] < 1).sum()
    print(f"\nType {t}: n={len(sub)}, mean_ratio={sub['aval_ratio'].mean():.3f}, median={sub['aval_ratio'].median():.3f}, <1: {below}/{len(sub)}")
