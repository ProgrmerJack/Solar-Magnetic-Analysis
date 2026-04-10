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

swiss_events = ssw_cat[(ssw_cat.index >= '1998-11-01') & (ssw_cat.index <= '2019-04-30')].index

aval = panel['dry_natural_size_1234']
aval.index = pd.to_datetime(aval.index)

nao = panel['nao_daily']
nao.index = pd.to_datetime(nao.index)

u850 = trop['uwnd_850hPa_ms']

# 1. NAO response during SSW events
print('=== NAO RESPONSE DURING SSW EVENTS ===')
nao_pre = []
nao_post = []
for onset in swiss_events:
    onset_ts = pd.Timestamp(onset)
    pre_window = pd.date_range(onset_ts - pd.Timedelta(days=30), onset_ts - pd.Timedelta(days=1))
    post_window = pd.date_range(onset_ts, onset_ts + pd.Timedelta(days=30))
    nao_pre.append(nao.reindex(pre_window).mean())
    nao_post.append(nao.reindex(post_window).mean())

nao_pre = np.array(nao_pre)
nao_post = np.array(nao_post)
mask = np.isfinite(nao_pre) & np.isfinite(nao_post)
print(f'NAO pre-SSW: {nao_pre[mask].mean():.3f} +/- {nao_pre[mask].std():.3f}')
print(f'NAO post-SSW: {nao_post[mask].mean():.3f} +/- {nao_post[mask].std():.3f}')
t, p = stats.ttest_rel(nao_pre[mask], nao_post[mask])
print(f'Paired t-test: t={t:.3f}, P={p:.4f}')
print(f'NAO decrease: {(nao_post[mask] < nao_pre[mask]).sum()}/{mask.sum()} events')

# 2. Per-event: NAO change vs avalanche change
print('\n=== PER-EVENT: NAO change vs Avalanche anomaly ===')
nao_changes = []
aval_ratios = []
u850_changes = []
for onset in swiss_events:
    onset_ts = pd.Timestamp(onset)
    pre = pd.date_range(onset_ts - pd.Timedelta(days=30), onset_ts - pd.Timedelta(days=1))
    post = pd.date_range(onset_ts - pd.Timedelta(days=15), onset_ts + pd.Timedelta(days=15))
    
    n_pre = nao.reindex(pre).mean()
    n_post = nao.reindex(post).mean()
    
    a_ssw = aval.reindex(post).mean()
    doy = onset_ts.timetuple().tm_yday
    doys = set()
    for d in range(-15, 16):
        doys.add((doy + d) % 366)
    a_clim = aval[aval.index.dayofyear.isin(doys)].mean()
    a_ratio = a_ssw / a_clim if a_clim > 0 else np.nan
    
    u_pre = u850.reindex(pre).mean()
    u_post = u850.reindex(post).mean()
    
    nao_changes.append(n_post - n_pre)
    aval_ratios.append(a_ratio)
    u850_changes.append(u_post - u_pre)
    
    print(f'{onset_ts.strftime("%Y-%m-%d")}: NAO_chg={n_post-n_pre:+.2f}, U850_chg={u_post-u_pre:+.2f}, aval_ratio={a_ratio:.3f}')

# Correlations
nao_ch = np.array(nao_changes)
aval_r = np.array(aval_ratios)
u850_ch = np.array(u850_changes)
mask = np.isfinite(nao_ch) & np.isfinite(aval_r)

print(f'\n=== MEDIATION TEST ===')
r1, p1 = stats.spearmanr(nao_ch[mask], aval_r[mask])
print(f'NAO change vs aval ratio: r={r1:.3f}, P={p1:.4f}')

r2, p2 = stats.spearmanr(u850_ch[mask], aval_r[mask])
print(f'U850 change vs aval ratio: r={r2:.3f}, P={p2:.4f}')

# 3. Direct test: daily NAO vs daily avalanche (deseasonalized)
print('\n=== DAILY NAO vs DAILY AVALANCHES (winter only) ===')
common = nao.index.intersection(aval.index)
n_c = nao.loc[common].values.astype(float)
a_c = aval.loc[common].values.astype(float)
winter = np.array([d.month in [11,12,1,2,3,4] for d in common])
n_w = n_c[winter]
a_w = a_c[winter]
valid = np.isfinite(n_w) & np.isfinite(a_w)
n_w = n_w[valid]
a_w = a_w[valid]

# Deseasonalize
doy_arr = np.array([d.timetuple().tm_yday for d in common[winter][valid]])
n_anom = n_w.copy()
a_anom = a_w.copy()
for d in set(doy_arr):
    dm = doy_arr == d
    n_anom[dm] -= n_w[dm].mean()
    a_anom[dm] -= a_w[dm].mean()

for lag in [0, 5, 10, 14]:
    if lag == 0:
        nn, aa = n_anom, a_anom
    else:
        nn = n_anom[:-lag]
        aa = a_anom[lag:]
    r, p = stats.spearmanr(nn, aa)
    print(f'Lag {lag:2d}d: r={r:.4f}, P={p:.2e}')

# NAO negative vs positive: avalanche rates
nao_neg = a_w[n_w < -0.5]
nao_pos = a_w[n_w > 0.5]
print(f'\nNAO negative (<-0.5): mean aval={nao_neg.mean():.3f}/day (n={len(nao_neg)})')
print(f'NAO positive (>+0.5): mean aval={nao_pos.mean():.3f}/day (n={len(nao_pos)})')
u, p = stats.mannwhitneyu(nao_neg, nao_pos, alternative='two-sided')
print(f'MWU P={p:.4f}')
