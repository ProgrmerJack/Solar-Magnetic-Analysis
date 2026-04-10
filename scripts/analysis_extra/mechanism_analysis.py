"""
Event-level dose-response and phase decomposition analysis.
Tests whether temperature anomaly predicts avalanche rate ratio,
and whether pre-SSW vs post-SSW phases differ.
"""
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import binom
import json

panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
era5 = pd.read_parquet('data/processed/era5_swiss_alps_daily.parquet')
ssw = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw.index = ssw.index.tz_localize(None)
aval = 'dry_natural_size_1234'

panel.index.name = 'date'
era5.index.name = 'date'
merged = panel.join(era5[['t2m_K_anom', 'sf_mm_anom']], how='inner').dropna(subset=[aval])

ssw_in_range = ssw[(ssw.index >= merged.index.min()) & (ssw.index <= merged.index.max())]
print(f'SSW events in ERA5 range: {len(ssw_in_range)}')

results = []
for sd in ssw_in_range.index:
    w = merged[(merged.index >= sd - pd.Timedelta(days=15)) & (merged.index <= sd + pd.Timedelta(days=15))]
    ctrl_parts = []
    for yr in range(-5, 6):
        if yr == 0:
            continue
        c = merged[(merged.index >= sd - pd.Timedelta(days=15) + pd.DateOffset(years=yr)) &
                    (merged.index <= sd + pd.Timedelta(days=15) + pd.DateOffset(years=yr))]
        ctrl_parts.append(c)
    ctrl = pd.concat(ctrl_parts)
    if len(w) < 5 or len(ctrl) < 10:
        continue

    rr = w[aval].mean() / ctrl[aval].mean() if ctrl[aval].mean() > 0 else np.nan
    pre_w = merged[(merged.index >= sd - pd.Timedelta(days=15)) & (merged.index < sd)]
    post_w = merged[(merged.index >= sd) & (merged.index <= sd + pd.Timedelta(days=15))]

    results.append({
        'date': str(sd.date()),
        'rr': float(rr) if not np.isnan(rr) else None,
        'log_rr': float(np.log(rr)) if (not np.isnan(rr) and rr > 0) else None,
        't_anom': float(w['t2m_K_anom'].mean()),
        'sf_anom': float(w['sf_mm_anom'].mean()),
        'pre_aval': float(pre_w[aval].mean()) if len(pre_w) > 0 else None,
        'post_aval': float(post_w[aval].mean()) if len(post_w) > 0 else None,
        'pre_t': float(pre_w['t2m_K_anom'].mean()) if len(pre_w) > 0 else None,
        'post_t': float(post_w['t2m_K_anom'].mean()) if len(post_w) > 0 else None,
    })

df = pd.DataFrame(results)
print()
print('=== EVENT-LEVEL RESULTS ===')
for _, r in df.iterrows():
    d = r['date']
    rr_v = r['rr'] if r['rr'] is not None else float('nan')
    t_v = r['t_anom'] if r['t_anom'] is not None else float('nan')
    sf_v = r['sf_anom'] if r['sf_anom'] is not None else float('nan')
    pre_v = r['pre_aval'] if r['pre_aval'] is not None else float('nan')
    post_v = r['post_aval'] if r['post_aval'] is not None else float('nan')
    print(f"  {d}: RR={rr_v:.2f} T={t_v:+.2f}K SF={sf_v:+.1f}mm PRE={pre_v:.2f} POST={post_v:.2f}")

# Dose-response: T anomaly -> avalanche RR
v = df.dropna(subset=['log_rr', 't_anom'])
x_arr = v['t_anom'].values.astype(float)
y_arr = v['log_rr'].values.astype(float)
r_s, p_s = stats.spearmanr(x_arr, y_arr)
r_p, p_p = stats.pearsonr(x_arr, y_arr)
slope, intercept, rv, pv, se = stats.linregress(x_arr, y_arr)
print()
print('=== DOSE-RESPONSE: T anomaly -> avalanche RR ===')
print(f'  Spearman r={r_s:.3f} P={p_s:.4f}')
print(f'  Pearson  r={r_p:.3f} P={p_p:.4f}')
print(f'  Linear: log(RR) = {intercept:.3f} + {slope:.3f}*T  R2={rv**2:.3f}')
print(f'  n events: {len(v)}')

# Phase decomposition
print()
print('=== PHASE DECOMPOSITION ===')
both_valid = df.dropna(subset=['pre_aval', 'post_aval'])
print(f'  n events with both phases: {len(both_valid)}')
print(f'  Mean PRE aval rate:  {both_valid["pre_aval"].mean():.3f}')
print(f'  Mean POST aval rate: {both_valid["post_aval"].mean():.3f}')
print(f'  PRE T mean:  {both_valid["pre_t"].mean():+.3f}K')
print(f'  POST T mean: {both_valid["post_t"].mean():+.3f}K')

t_stat, p_phase = stats.ttest_rel(both_valid['pre_aval'].astype(float), both_valid['post_aval'].astype(float))
print(f'  Paired t-test PRE vs POST aval: t={t_stat:.3f}, P={p_phase:.4f}')

t_stat2, p_phase2 = stats.ttest_rel(both_valid['pre_t'].astype(float), both_valid['post_t'].astype(float))
print(f'  Paired t-test PRE vs POST temp: t={t_stat2:.3f}, P={p_phase2:.4f}')

# Snowfall dose-response
sf_arr = v['sf_anom'].values.astype(float)
r_sf, p_sf = stats.spearmanr(sf_arr, y_arr)
print()
print(f'=== SNOWFALL DOSE-RESPONSE ===')
print(f'  Spearman r={r_sf:.3f} P={p_sf:.4f}')

# Bidirectional response analysis
print()
print('=== BIDIRECTIONAL RESPONSE ===')
bidi_results = []
for sd in ssw_in_range.index:
    pre_w = merged[(merged.index >= sd - pd.Timedelta(days=15)) & (merged.index < sd)]
    post_w = merged[(merged.index >= sd) & (merged.index <= sd + pd.Timedelta(days=15))]
    
    pre_ctrls = []
    post_ctrls = []
    for yr in range(-5, 6):
        if yr == 0:
            continue
        pc = merged[(merged.index >= sd - pd.Timedelta(days=15) + pd.DateOffset(years=yr)) &
                     (merged.index < sd + pd.DateOffset(years=yr))]
        oc = merged[(merged.index >= sd + pd.DateOffset(years=yr)) &
                     (merged.index <= sd + pd.Timedelta(days=15) + pd.DateOffset(years=yr))]
        pre_ctrls.append(pc)
        post_ctrls.append(oc)
    
    pre_ctrl = pd.concat(pre_ctrls)
    post_ctrl = pd.concat(post_ctrls)
    
    if len(pre_w) < 3 or len(post_w) < 3:
        continue
    if len(pre_ctrl) < 5 or len(post_ctrl) < 5:
        continue
    if pre_ctrl[aval].mean() == 0 or post_ctrl[aval].mean() == 0:
        continue
    
    pre_rr = pre_w[aval].mean() / pre_ctrl[aval].mean()
    post_rr = post_w[aval].mean() / post_ctrl[aval].mean()
    
    bidi_results.append({
        'date': str(sd.date()),
        'pre_rr': float(pre_rr),
        'post_rr': float(post_rr),
        'pre_t': float(pre_w['t2m_K_anom'].mean()),
        'post_t': float(post_w['t2m_K_anom'].mean()),
    })

bdf = pd.DataFrame(bidi_results)
print(f'Events with bidirectional data: {len(bdf)}')
print(f'  Mean PRE rate ratio:  {bdf["pre_rr"].mean():.3f}')
print(f'  Mean POST rate ratio: {bdf["post_rr"].mean():.3f}')
print(f'  PRE < 1 in {(bdf["pre_rr"] < 1).sum()}/{len(bdf)} events')
print(f'  POST < 1 in {(bdf["post_rr"] < 1).sum()}/{len(bdf)} events')

n_pre_decrease = int((bdf['pre_rr'] < 1).sum())
n_post_decrease = int((bdf['post_rr'] < 1).sum())
n_total = len(bdf)
p_pre_sign = float(binom.sf(n_pre_decrease - 1, n_total, 0.5))
p_post_sign = float(binom.sf(n_post_decrease - 1, n_total, 0.5))
print(f'  PRE sign test: {n_pre_decrease}/{n_total} decrease, P={p_pre_sign:.4f}')
print(f'  POST sign test: {n_post_decrease}/{n_total} decrease, P={p_post_sign:.4f}')

# Temperature correlation with each phase
r_pre, p_pre = stats.spearmanr(bdf['pre_t'].values, bdf['pre_rr'].values)
r_post, p_post = stats.spearmanr(bdf['post_t'].values, bdf['post_rr'].values)
print(f'  PRE: T-RR correlation r={r_pre:.3f} P={p_pre:.4f}')
print(f'  POST: T-RR correlation r={r_post:.3f} P={p_post:.4f}')

# Save results
output = {
    'dose_response': {
        'spearman_r': float(r_s), 'spearman_p': float(p_s),
        'pearson_r': float(r_p), 'pearson_p': float(p_p),
        'slope': float(slope), 'intercept': float(intercept), 'r2': float(rv**2),
        'n_events': len(v)
    },
    'phase_decomposition': {
        'n_events': len(both_valid),
        'pre_mean_rate': float(both_valid['pre_aval'].mean()),
        'post_mean_rate': float(both_valid['post_aval'].mean()),
        'pre_mean_t': float(both_valid['pre_t'].mean()),
        'post_mean_t': float(both_valid['post_t'].mean()),
        'paired_t_aval': float(t_stat), 'paired_p_aval': float(p_phase),
    },
    'bidirectional': {
        'n_events': len(bdf),
        'pre_mean_rr': float(bdf['pre_rr'].mean()),
        'post_mean_rr': float(bdf['post_rr'].mean()),
        'pre_n_decrease': n_pre_decrease,
        'post_n_decrease': n_post_decrease,
        'pre_sign_p': p_pre_sign,
        'post_sign_p': p_post_sign,
    },
    'events': results
}

with open('data/results/mechanism_analysis.json', 'w') as f:
    json.dump(output, f, indent=2, default=str)
print('\nSaved to data/results/mechanism_analysis.json')
