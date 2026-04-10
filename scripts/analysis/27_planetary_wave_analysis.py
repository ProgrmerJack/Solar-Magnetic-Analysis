"""
27_planetary_wave_analysis.py
Addresses both reviewers' key demands:
1. Planetary wave forcing proxy analysis (common cause test)
2. Conditional independence: does SSW add info beyond wave forcing?
3. Formal multiple testing correction (FDR)
4. Extended placebo tests
"""
import pandas as pd
import numpy as np
from scipy import stats
import json
import warnings
warnings.filterwarnings('ignore')

print("="*70)
print("PART 1: PLANETARY WAVE FORCING PROXY")
print("="*70)

# Load data
panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
panel.index = pd.to_datetime(panel.index)

ssw_cat = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw_dates = pd.DatetimeIndex(ssw_cat.index).tz_localize(None)
study_ssw = ssw_dates[(ssw_dates >= '1998-11-01') & (ssw_dates <= '2019-04-30')]
print(f"SSW events in study period: {len(study_ssw)}")

# Planetary wave forcing proxy: pre-SSW warming rate at 100 hPa
# When planetary waves propagate upward and break in the stratosphere,
# they warm the polar stratosphere. The rate of warming at 100 hPa
# in the 15 days before SSW onset is a proxy for wave forcing intensity.

wave_proxy = []
pre_aval_anomaly = []
post_aval_anomaly = []

# Compute day-of-year climatology for avalanches and temperature
winter_mask = panel.index.month.isin([11, 12, 1, 2, 3, 4])
panel_w = panel[winter_mask].copy()
panel_w['doy'] = panel_w.index.dayofyear

# Find the dry slab column
dry_col = 'dry_natural_size_1234' if 'dry_natural_size_1234' in panel.columns else 'aai_all_dry'
print(f"Using avalanche column: {dry_col}")

# Avalanche climatology
aval_clim = panel_w.groupby('doy')[dry_col].mean()

# T100 climatology
t100_clim = panel_w.groupby('doy')['ncep_t_100hpa'].mean()

for ssw_date in study_ssw:
    # Pre-SSW wave forcing proxy: T100 warming rate (day -15 to day -1)
    pre_window = pd.date_range(ssw_date - pd.Timedelta(days=15), ssw_date - pd.Timedelta(days=1))
    pre_data = panel.reindex(pre_window).dropna(subset=['ncep_t_100hpa'])
    
    if len(pre_data) >= 10:
        # Compute anomalies relative to climatology
        t100_anoms = []
        for d in pre_data.index:
            doy = d.dayofyear
            if doy in t100_clim.index:
                t100_anoms.append(pre_data.loc[d, 'ncep_t_100hpa'] - t100_clim.get(doy, 0))
        wave_forcing = np.mean(t100_anoms) if t100_anoms else np.nan
    else:
        wave_forcing = np.nan
    
    # Pre-SSW avalanche anomaly (days -15 to -1)
    pre_aval_data = panel.reindex(pre_window).dropna(subset=[dry_col])
    if len(pre_aval_data) >= 10:
        pre_anom = []
        for d in pre_aval_data.index:
            doy = d.dayofyear
            if doy in aval_clim.index:
                pre_anom.append(pre_aval_data.loc[d, dry_col] - aval_clim.get(doy, 0))
        pre_aval = np.mean(pre_anom) if pre_anom else np.nan
    else:
        pre_aval = np.nan
    
    # Post-SSW avalanche anomaly (days 0 to 14)
    post_window = pd.date_range(ssw_date, ssw_date + pd.Timedelta(days=14))
    post_aval_data = panel.reindex(post_window).dropna(subset=[dry_col])
    if len(post_aval_data) >= 10:
        post_anom = []
        for d in post_aval_data.index:
            doy = d.dayofyear
            if doy in aval_clim.index:
                post_anom.append(post_aval_data.loc[d, dry_col] - aval_clim.get(doy, 0))
        post_aval = np.mean(post_anom) if post_anom else np.nan
    else:
        post_aval = np.nan
    
    wave_proxy.append(wave_forcing)
    pre_aval_anomaly.append(pre_aval)
    post_aval_anomaly.append(post_aval)
    print(f"  {ssw_date.strftime('%Y-%m-%d')}: wave={wave_forcing:.2f}K, pre_aval={pre_aval:.2f}, post_aval={post_aval:.2f}")

wave_proxy = np.array(wave_proxy)
pre_aval_anomaly = np.array(pre_aval_anomaly)
post_aval_anomaly = np.array(post_aval_anomaly)

# Test: Does wave forcing proxy predict avalanche reduction?
valid = ~np.isnan(wave_proxy) & ~np.isnan(post_aval_anomaly)
r_wave_post, p_wave_post = stats.spearmanr(wave_proxy[valid], post_aval_anomaly[valid])
r_wave_pre, p_wave_pre = stats.spearmanr(wave_proxy[valid], pre_aval_anomaly[valid])

print(f"\n--- Wave Forcing → Avalanche ---")
print(f"Wave proxy vs POST-SSW aval anomaly: r={r_wave_post:.3f}, P={p_wave_post:.4f}")
print(f"Wave proxy vs PRE-SSW aval anomaly:  r={r_wave_pre:.3f}, P={p_wave_pre:.4f}")

# Test: Does pre-SSW avalanche anomaly predict post-SSW anomaly?
valid2 = ~np.isnan(pre_aval_anomaly) & ~np.isnan(post_aval_anomaly)
r_pre_post, p_pre_post = stats.spearmanr(pre_aval_anomaly[valid2], post_aval_anomaly[valid2])
print(f"Pre-SSW aval vs Post-SSW aval:       r={r_pre_post:.3f}, P={p_pre_post:.4f}")


print("\n" + "="*70)
print("PART 2: CONDITIONAL INDEPENDENCE TEST")
print("="*70)
print("If common cause is correct, controlling for wave forcing should")
print("eliminate the SSW→avalanche association.")

# Simple approach: partial correlation
# Regress post-SSW avalanche anomaly on wave forcing, get residuals
# Then test if SSW adds anything

from numpy.linalg import lstsq

# Build event-level dataframe
event_df = pd.DataFrame({
    'wave_forcing': wave_proxy,
    'pre_aval': pre_aval_anomaly,
    'post_aval': post_aval_anomaly,
    'ssw_date': study_ssw
})
event_df = event_df.dropna()

# Partial correlation: post_aval ~ wave_forcing, then check residuals
X = np.column_stack([np.ones(len(event_df)), event_df['wave_forcing'].values])
y = event_df['post_aval'].values
beta, _, _, _ = lstsq(X, y, rcond=None)
resid_aval = y - X @ beta

# The residuals represent avalanche anomaly AFTER removing wave forcing effect
# Under common cause, these residuals should be centered at 0
t_resid, p_resid = stats.ttest_1samp(resid_aval, 0)
print(f"Residual aval anomaly (after wave): mean={np.mean(resid_aval):.3f}, t={t_resid:.3f}, P={p_resid:.4f}")
print(f"  If P is large → wave forcing explains the association (common cause)")
print(f"  If P is small → SSW adds independent info beyond wave forcing")

# Also test: pre-SSW aval as predictor of post-SSW aval
X2 = np.column_stack([np.ones(len(event_df)), event_df['pre_aval'].values])
beta2, _, _, _ = lstsq(X2, y, rcond=None)
resid2 = y - X2 @ beta2
t_resid2, p_resid2 = stats.ttest_1samp(resid2, 0)
print(f"\nResidual post aval (after pre-aval): mean={np.mean(resid2):.3f}, t={t_resid2:.3f}, P={p_resid2:.4f}")

# Multiple regression: post_aval ~ wave + pre_aval
X3 = np.column_stack([np.ones(len(event_df)), event_df['wave_forcing'].values, event_df['pre_aval'].values])
beta3, _, _, _ = lstsq(X3, y, rcond=None)
resid3 = y - X3 @ beta3
t_resid3, p_resid3 = stats.ttest_1samp(resid3, 0)
print(f"\nResidual post (after wave + pre):    mean={np.mean(resid3):.3f}, t={t_resid3:.3f}, P={p_resid3:.4f}")
print(f"  Wave coef: {beta3[1]:.4f}, Pre-aval coef: {beta3[2]:.4f}")


print("\n" + "="*70)
print("PART 3: DIRECT WAVE FORCING → AVALANCHE (DAILY LEVEL)")
print("="*70)

# Create daily wave forcing proxy: 15-day rolling mean of T100 anomaly
panel_w2 = panel_w.copy()
panel_w2['t100_anom'] = panel_w2.apply(
    lambda row: row['ncep_t_100hpa'] - t100_clim.get(row['doy'], row['ncep_t_100hpa']),
    axis=1
)
panel_w2['wave_proxy_15d'] = panel_w2['t100_anom'].rolling(15, min_periods=10).mean()

# Correlation between wave proxy and avalanche counts
valid_daily = panel_w2[['wave_proxy_15d', dry_col]].dropna()
r_daily, p_daily = stats.spearmanr(valid_daily['wave_proxy_15d'], valid_daily[dry_col])
print(f"Daily wave proxy (15d T100 anom) vs dry slab: r={r_daily:.4f}, P={p_daily:.6f}")

# Compare with SSW binary predictor
panel_w2['ssw_window'] = 0
for sd in study_ssw:
    mask = (panel_w2.index >= sd) & (panel_w2.index <= sd + pd.Timedelta(days=14))
    panel_w2.loc[mask, 'ssw_window'] = 1

# Add pre-SSW window
panel_w2['ssw_lifecycle'] = 0
for sd in study_ssw:
    pre_mask = (panel_w2.index >= sd - pd.Timedelta(days=15)) & (panel_w2.index < sd)
    post_mask = (panel_w2.index >= sd) & (panel_w2.index <= sd + pd.Timedelta(days=14))
    panel_w2.loc[pre_mask, 'ssw_lifecycle'] = 1
    panel_w2.loc[post_mask, 'ssw_lifecycle'] = 1

r_ssw, p_ssw = stats.spearmanr(panel_w2['ssw_window'].dropna(), 
                                 panel_w2.loc[panel_w2['ssw_window'].notna(), dry_col])
r_lifecycle, p_lifecycle = stats.spearmanr(panel_w2['ssw_lifecycle'].dropna(),
                                            panel_w2.loc[panel_w2['ssw_lifecycle'].notna(), dry_col])
print(f"SSW binary (post only) vs dry slab:   r={r_ssw:.4f}, P={p_ssw:.6f}")
print(f"SSW lifecycle (pre+post) vs dry slab:  r={r_lifecycle:.4f}, P={p_lifecycle:.6f}")

# Key test: Does wave proxy EXPLAIN the SSW-avalanche association?
# Partial correlation: avalanche ~ SSW | wave_proxy
from scipy.stats import pearsonr
valid_all = panel_w2[['wave_proxy_15d', dry_col, 'ssw_lifecycle']].dropna()
# Residualize both on wave proxy
X_w = np.column_stack([np.ones(len(valid_all)), valid_all['wave_proxy_15d'].values])
# Residualize avalanche
beta_a, _, _, _ = lstsq(X_w, valid_all[dry_col].values, rcond=None)
resid_aval_d = valid_all[dry_col].values - X_w @ beta_a
# Residualize SSW
beta_s, _, _, _ = lstsq(X_w, valid_all['ssw_lifecycle'].values, rcond=None)
resid_ssw_d = valid_all['ssw_lifecycle'].values - X_w @ beta_s
# Partial correlation
r_partial, p_partial = pearsonr(resid_aval_d, resid_ssw_d)
print(f"\nPartial corr (SSW→aval | wave proxy): r={r_partial:.4f}, P={p_partial:.6f}")
print(f"  If partial r ≈ 0 → wave proxy explains the SSW-aval association (common cause)")
print(f"  If partial r still significant → SSW has independent effect")


print("\n" + "="*70)
print("PART 4: FORMAL MULTIPLE TESTING CORRECTION (FDR)")
print("="*70)

# Collect all primary P-values reported in the paper
pvals = {
    'Matched SSW paired t': 0.003,
    'Matched SSW sign test': 0.001,
    'Matched SSW Wilcoxon': 0.008,
    'Matched SSW permutation': 0.005,
    'Pre-SSW matched': 0.009,
    'Post-SSW matched': 0.003,
    'Late SSW matched': 0.070,
    'Pre vs post increment': 0.10,
    'MH dry slab parametric': 1e-6,
    'MH dry slab permutation': 0.09,
    'MH wet slab': 0.17,
    'NAO-adjusted SSW IRR': 0.001,
    'Event-level permutation': 0.083,
    'Dose-response': 0.87,
    'ERA5 precipitation': 0.16,
    'ERA5 snowfall': 0.27,
    'ERA5 temperature': 0.26,
    'ERA5 wind': 0.55,
    'SNOTEL SWE': 0.45,
    'SNOTEL precip': 0.74,
    'SNOTEL temp': 0.34,
    'CAIC fatalities': 0.14,
    'Snowfall→avalanche': 1e-12,
    'U850→avalanche daily': 0.96,
    'U850→avalanche event': 0.79,
}

# Benjamini-Hochberg FDR
from statsmodels.stats.multitest import multipletests
names = list(pvals.keys())
raw_p = np.array([pvals[k] for k in names])
reject, fdr_p, _, _ = multipletests(raw_p, method='fdr_bh', alpha=0.05)

print(f"{'Test':<35s} {'Raw P':>10s} {'FDR P':>10s} {'Sig?':>5s}")
print("-" * 65)
# Sort by raw P
order = np.argsort(raw_p)
for i in order:
    sig = '✓' if reject[i] else ''
    print(f"{names[i]:<35s} {raw_p[i]:>10.4f} {fdr_p[i]:>10.4f} {sig:>5s}")

n_sig = np.sum(reject)
print(f"\n{n_sig}/{len(names)} tests significant after FDR correction at α=0.05")


print("\n" + "="*70)
print("PART 5: EXTENDED PLACEBO ANALYSIS")
print("="*70)

# Generate 2000 pseudo-SSW catalogs (more than the 1000 in mechanism analysis)
np.random.seed(42)
n_perms = 2000

# Get all winter days
winter_days = panel_w.index.values

def generate_pseudo_ssw(winter_days, n_events=15, min_sep=30):
    """Generate random pseudo-SSW dates with minimum separation."""
    for _ in range(1000):
        idx = np.random.choice(len(winter_days), size=n_events*3, replace=False)
        candidates = np.sort(winter_days[idx])
        selected = [candidates[0]]
        for c in candidates[1:]:
            if (c - selected[-1]) / np.timedelta64(1, 'D') >= min_sep:
                selected.append(c)
            if len(selected) == n_events:
                break
        if len(selected) == n_events:
            return pd.DatetimeIndex(selected)
    return None

perm_mean_anomalies = []
perm_concordances = []
perm_pre_anomalies = []

for i in range(n_perms):
    pseudo = generate_pseudo_ssw(winter_days)
    if pseudo is None:
        continue
    
    post_anoms = []
    pre_anoms = []
    for sd in pseudo:
        sd = pd.Timestamp(sd)
        # Post
        post_w = pd.date_range(sd, sd + pd.Timedelta(days=14))
        post_d = panel.reindex(post_w).dropna(subset=[dry_col])
        if len(post_d) >= 10:
            anoms = [post_d.loc[d, dry_col] - aval_clim.get(d.dayofyear, 0) 
                     for d in post_d.index if d.dayofyear in aval_clim.index]
            if anoms:
                post_anoms.append(np.mean(anoms))
        # Pre
        pre_w = pd.date_range(sd - pd.Timedelta(days=15), sd - pd.Timedelta(days=1))
        pre_d = panel.reindex(pre_w).dropna(subset=[dry_col])
        if len(pre_d) >= 10:
            anoms_pre = [pre_d.loc[d, dry_col] - aval_clim.get(d.dayofyear, 0)
                         for d in pre_d.index if d.dayofyear in aval_clim.index]
            if anoms_pre:
                pre_anoms.append(np.mean(anoms_pre))
    
    if len(post_anoms) >= 10:
        perm_mean_anomalies.append(np.mean(post_anoms))
        perm_concordances.append(np.mean([a < 0 for a in post_anoms]))
    if len(pre_anoms) >= 10:
        perm_pre_anomalies.append(np.mean(pre_anoms))
    
    if (i + 1) % 500 == 0:
        print(f"  Permutation {i+1}/{n_perms} done")

perm_mean_anomalies = np.array(perm_mean_anomalies)
perm_concordances = np.array(perm_concordances)
perm_pre_anomalies = np.array(perm_pre_anomalies)

real_post = np.mean(post_aval_anomaly[~np.isnan(post_aval_anomaly)])
real_pre = np.mean(pre_aval_anomaly[~np.isnan(pre_aval_anomaly)])
real_conc = np.mean(post_aval_anomaly[~np.isnan(post_aval_anomaly)] < 0)

p_post_perm = np.mean(perm_mean_anomalies <= real_post)
p_pre_perm = np.mean(perm_pre_anomalies <= real_pre)
p_conc_perm = np.mean(perm_concordances >= real_conc)

print(f"\n--- 2000-Permutation Results ---")
print(f"Real post-SSW anomaly: {real_post:.3f}, P={p_post_perm:.4f} (percentile={p_post_perm*100:.1f}th)")
print(f"Real pre-SSW anomaly:  {real_pre:.3f}, P={p_pre_perm:.4f}")
print(f"Real concordance:      {real_conc:.2f}, P={p_conc_perm:.4f}")
print(f"Null distribution: median={np.median(perm_mean_anomalies):.3f}, "
      f"90% range=[{np.percentile(perm_mean_anomalies, 5):.3f}, {np.percentile(perm_mean_anomalies, 95):.3f}]")


print("\n" + "="*70)
print("PART 6: SUMMARY FOR PAPER REWRITE")
print("="*70)

results = {
    'wave_forcing': {
        'r_wave_post': float(r_wave_post),
        'p_wave_post': float(p_wave_post),
        'r_wave_pre': float(r_wave_pre),
        'p_wave_pre': float(p_wave_pre),
        'r_pre_post_aval': float(r_pre_post),
        'p_pre_post_aval': float(p_pre_post),
    },
    'conditional_independence': {
        'resid_after_wave_mean': float(np.mean(resid_aval)),
        'resid_after_wave_t': float(t_resid),
        'resid_after_wave_p': float(p_resid),
        'resid_after_pre_mean': float(np.mean(resid2)),
        'resid_after_pre_t': float(t_resid2),
        'resid_after_pre_p': float(p_resid2),
    },
    'daily_level': {
        'r_wave_daily': float(r_daily),
        'p_wave_daily': float(p_daily),
        'r_ssw_binary': float(r_ssw),
        'p_ssw_binary': float(p_ssw),
        'r_ssw_lifecycle': float(r_lifecycle),
        'p_ssw_lifecycle': float(p_lifecycle),
        'partial_r_ssw_given_wave': float(r_partial),
        'partial_p_ssw_given_wave': float(p_partial),
    },
    'fdr': {
        'n_tests': len(names),
        'n_significant_fdr': int(n_sig),
        'tests_significant': [names[i] for i in range(len(names)) if reject[i]],
    },
    'permutation_2000': {
        'post_anomaly': float(real_post),
        'post_p': float(p_post_perm),
        'pre_anomaly': float(real_pre),
        'pre_p': float(p_pre_perm),
        'concordance': float(real_conc),
        'concordance_p': float(p_conc_perm),
        'null_median': float(np.median(perm_mean_anomalies)),
        'null_90_lo': float(np.percentile(perm_mean_anomalies, 5)),
        'null_90_hi': float(np.percentile(perm_mean_anomalies, 95)),
    }
}

with open('data/results/planetary_wave_analysis.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nResults saved to data/results/planetary_wave_analysis.json")

print("\n--- KEY FINDINGS FOR PAPER REWRITE ---")
print(f"1. Wave forcing → post-SSW avalanche: r={r_wave_post:.3f} (P={p_wave_post:.4f})")
print(f"2. Conditional independence (event): after wave control, P={p_resid:.4f}")
print(f"3. Daily partial corr (SSW|wave): r={r_partial:.4f} (P={p_partial:.6f})")
print(f"4. FDR-corrected: {n_sig}/{len(names)} tests survive")
print(f"5. 2000-perm post-SSW: P={p_post_perm:.4f}")

