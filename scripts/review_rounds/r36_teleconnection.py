"""
R36 Upgrade Part 2: NAO/AO-mediated mechanism analysis
Uses freshly downloaded daily teleconnection indices to quantify
the stratosphere-surface pathway through NAO/AO.
"""
import pandas as pd, numpy as np, json
from scipy import stats

# Load NAO/AO daily
def load_cpc_index(path):
    """Parse CPC daily teleconnection index format"""
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    y, m, d = int(parts[0]), int(parts[1]), int(parts[2])
                    val = float(parts[3])
                    rows.append({'date': pd.Timestamp(y, m, d), 'value': val})
                except:
                    continue
    return pd.DataFrame(rows).set_index('date')

nao = load_cpc_index('data/atmospheric/nao_daily.txt')
nao.columns = ['nao']
ao = load_cpc_index('data/atmospheric/ao_daily.txt')
ao.columns = ['ao']
pna = load_cpc_index('data/atmospheric/pna_daily.txt')
pna.columns = ['pna']

print(f"NAO: {nao.index.min().date()} to {nao.index.max().date()} ({len(nao)} days)")
print(f"AO:  {ao.index.min().date()} to {ao.index.max().date()} ({len(ao)} days)")
print(f"PNA: {pna.index.min().date()} to {pna.index.max().date()} ({len(pna)} days)")

# Merge
idx = nao.join(ao, how='outer').join(pna, how='outer')
print(f"Combined: {idx.index.min().date()} to {idx.index.max().date()}")

# Load SSW catalog
ssw_cat = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw_dates = ssw_cat.index.tz_localize(None)

# Load avalanche panel
panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')

# ============================================================
# PART 1: NAO/AO during SSW events
# ============================================================
print("\n" + "=" * 60)
print("PART 1: Teleconnection Indices During SSW Events")
print("=" * 60)

def get_window(series, date, w=15):
    mask = (series.index >= date - pd.Timedelta(days=w)) & (series.index <= date + pd.Timedelta(days=w))
    return series[mask]

# Winter background
winter_mask = idx.index.month.isin([1, 2, 3, 11, 12])
winter_idx = idx[winter_mask]

def ssw_mask_fn(index, dates, w=15):
    m = pd.Series(False, index=index)
    for d in dates:
        m |= (index >= d - pd.Timedelta(days=w)) & (index <= d + pd.Timedelta(days=w))
    return m

sm = ssw_mask_fn(winter_idx.index, ssw_dates)

for col in ['nao', 'ao', 'pna']:
    ssw_vals = winter_idx.loc[sm, col].dropna()
    ctrl_vals = winter_idx.loc[~sm, col].dropna()
    d = (ssw_vals.mean() - ctrl_vals.mean()) / ctrl_vals.std()
    _, p = stats.mannwhitneyu(ssw_vals, ctrl_vals, alternative='two-sided')
    print(f"{col.upper():>4}: SSW={ssw_vals.mean():+.3f} Ctrl={ctrl_vals.mean():+.3f} d={d:+.3f} P={p:.2e}")

# Event-level NAO/AO
print("\nEvent-level NAO/AO (mean during +/-15d window):")
ssw_in_range = ssw_dates[(ssw_dates >= idx.index.min()) & (ssw_dates <= idx.index.max())]
event_data = []
for d in ssw_in_range:
    w = get_window(idx, d, 15)
    if len(w) > 10:
        nao_val = w['nao'].mean()
        ao_val = w['ao'].mean()
        # Get avalanche RR from panel if available
        av_w = (panel.index >= d - pd.Timedelta(days=15)) & (panel.index <= d + pd.Timedelta(days=15))
        if av_w.sum() > 0 and 'dry_natural_size_1234' in panel.columns:
            av_ssw = panel.loc[av_w, 'dry_natural_size_1234'].mean()
            # Control
            ctrl_doy = panel.index.dayofyear
            target_doy = d.dayofyear
            ctrl_mask = (ctrl_doy >= target_doy - 10) & (ctrl_doy <= target_doy + 10) & ~av_w
            av_ctrl = panel.loc[ctrl_mask, 'dry_natural_size_1234'].mean()
            log_rr = np.log(max(av_ssw, 0.01) / max(av_ctrl, 0.01))
        else:
            log_rr = np.nan
        
        event_data.append({
            'onset': d,
            'nao': nao_val,
            'ao': ao_val,
            'log_rr': log_rr,
        })
        sign_nao = "neg" if nao_val < 0 else "POS"
        sign_ao = "neg" if ao_val < 0 else "POS"
        print(f"  {d.date()}: NAO={nao_val:+.2f}({sign_nao}) AO={ao_val:+.2f}({sign_ao}) log(RR)={log_rr:+.2f}")

edf = pd.DataFrame(event_data)

# ============================================================
# PART 2: NAO/AO correlation with avalanche RR
# ============================================================
print("\n" + "=" * 60)
print("PART 2: NAO/AO vs Avalanche Response")
print("=" * 60)

valid = edf.dropna(subset=['log_rr', 'nao', 'ao'])
if len(valid) >= 5:
    r_nao, p_nao = stats.pearsonr(valid['nao'], valid['log_rr'])
    r_ao, p_ao = stats.pearsonr(valid['ao'], valid['log_rr'])
    print(f"NAO vs log(RR): r={r_nao:.3f}, P={p_nao:.3f}")
    print(f"AO  vs log(RR): r={r_ao:.3f}, P={p_ao:.3f}")
    
    # NAO negative = stronger avalanche suppression?
    neg_nao = valid[valid['nao'] < 0]
    pos_nao = valid[valid['nao'] >= 0]
    if len(neg_nao) >= 3 and len(pos_nao) >= 3:
        print(f"\nNAO-stratified avalanche response:")
        print(f"  NAO<0 events (n={len(neg_nao)}): mean log(RR)={neg_nao['log_rr'].mean():.3f}")
        print(f"  NAO>=0 events (n={len(pos_nao)}): mean log(RR)={pos_nao['log_rr'].mean():.3f}")
        _, p_diff = stats.mannwhitneyu(neg_nao['log_rr'], pos_nao['log_rr'], alternative='less')
        print(f"  Mann-Whitney P (NAO<0 more negative)={p_diff:.3f}")

# ============================================================
# PART 3: NAO/AO Phase Cascade Timing
# For each SSW, track NAO/AO from -30 to +30 days
# ============================================================
print("\n" + "=" * 60)
print("PART 3: NAO/AO Phase Evolution Around SSW")
print("=" * 60)

lags = range(-30, 31)
nao_composite = []
ao_composite = []

for d in ssw_in_range:
    for lag in lags:
        target = d + pd.Timedelta(days=lag)
        if target in idx.index:
            nao_val = idx.loc[target, 'nao']
            ao_val = idx.loc[target, 'ao']
            if not np.isnan(nao_val):
                nao_composite.append({'lag': lag, 'value': nao_val})
            if not np.isnan(ao_val):
                ao_composite.append({'lag': lag, 'value': ao_val})

nao_comp = pd.DataFrame(nao_composite).groupby('lag')['value'].agg(['mean', 'std', 'count'])
ao_comp = pd.DataFrame(ao_composite).groupby('lag')['value'].agg(['mean', 'std', 'count'])

print("Lag   NAO_mean   AO_mean")
for lag in [-30, -20, -15, -10, -5, 0, 5, 10, 15, 20, 30]:
    nao_m = nao_comp.loc[lag, 'mean'] if lag in nao_comp.index else np.nan
    ao_m = ao_comp.loc[lag, 'mean'] if lag in ao_comp.index else np.nan
    print(f"{lag:+4d}   {nao_m:+.3f}    {ao_m:+.3f}")

# When does AO reach minimum?
ao_min_lag = ao_comp['mean'].idxmin()
nao_min_lag = nao_comp['mean'].idxmin()
print(f"\nAO minimum at lag={ao_min_lag}d (value={ao_comp.loc[ao_min_lag, 'mean']:+.3f})")
print(f"NAO minimum at lag={nao_min_lag}d (value={nao_comp.loc[nao_min_lag, 'mean']:+.3f})")

# ============================================================
# PART 4: Proportion of variance explained by NAO/AO
# ============================================================
print("\n" + "=" * 60)
print("PART 4: Mediation through NAO/AO")
print("=" * 60)

if len(valid) >= 5:
    from sklearn.linear_model import LinearRegression
    X = valid[['nao', 'ao']].values
    y = valid['log_rr'].values
    reg = LinearRegression().fit(X, y)
    r2 = reg.score(X, y)
    print(f"Multiple regression R2 (NAO+AO -> log(RR)): {r2:.3f}")
    print(f"  NAO beta: {reg.coef_[0]:.4f}")
    print(f"  AO beta:  {reg.coef_[1]:.4f}")

# ============================================================
# Save
# ============================================================
results = {
    'nao_ssw_d': float((winter_idx.loc[sm, 'nao'].dropna().mean() - winter_idx.loc[~sm, 'nao'].dropna().mean()) / winter_idx.loc[~sm, 'nao'].dropna().std()),
    'ao_ssw_d': float((winter_idx.loc[sm, 'ao'].dropna().mean() - winter_idx.loc[~sm, 'ao'].dropna().mean()) / winter_idx.loc[~sm, 'ao'].dropna().std()),
    'nao_vs_rr_r': float(r_nao) if 'r_nao' in dir() else None,
    'nao_vs_rr_p': float(p_nao) if 'p_nao' in dir() else None,
    'ao_vs_rr_r': float(r_ao) if 'r_ao' in dir() else None,
    'ao_min_lag': int(ao_min_lag),
    'nao_min_lag': int(nao_min_lag),
    'nao_ao_r2': float(r2) if 'r2' in dir() else None,
}

with open('data/results/r36_teleconnection.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nSaved to data/results/r36_teleconnection.json")
