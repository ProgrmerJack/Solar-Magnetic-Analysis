"""
Comprehensive causal mediation analysis: quantify how much of the SSW→avalanche
effect flows through each identified pathway.

Pathways tested:
1. SSW → Z500 depression → weather regime shift → avalanche decrease
2. SSW → precipitation phase shift → avalanche decrease
3. SSW → temperature cooling → sintering change → avalanche decrease
4. SSW → trigger-day reduction → avalanche decrease

Uses Baron-Kenny 4-step approach + bootstrapped indirect effects.
"""
import pandas as pd, numpy as np, json, warnings
from scipy import stats
from scipy.stats import binomtest
warnings.filterwarnings('ignore')

panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
era5 = pd.read_parquet('data/processed/era5_swiss_alps_extended.parquet')
ssw_cat = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw_dates = ssw_cat.index.tz_localize(None)

# Align dates
common_dates = panel.index.intersection(era5.index)
panel_a = panel.loc[common_dates]
era5_a = era5.loc[common_dates]

# SSW binary indicator
ssw_in = ssw_dates[(ssw_dates >= panel.index.min()) & (ssw_dates <= panel.index.max())]
ssw_window = pd.Series(False, index=common_dates)
for d in ssw_in:
    mask = (common_dates >= d - pd.Timedelta(days=15)) & (common_dates <= d + pd.Timedelta(days=15))
    ssw_window[mask] = True

# Winter only
doy = common_dates.dayofyear
winter = (doy <= 120) | (doy >= 305)
idx = winter & common_dates.isin(panel_a.index) & common_dates.isin(era5_a.index)

df = pd.DataFrame({
    'ssw': ssw_window[idx].astype(float).values,
    'avcount': panel_a.loc[idx, 'dry_natural_size_1234'].values,
    't2m': era5_a.loc[idx, 't2m_K'].values if 't2m_K' in era5_a.columns else np.nan,
    'tp': era5_a.loc[idx, 'tp_mm'].values if 'tp_mm' in era5_a.columns else np.nan,
    'sf': era5_a.loc[idx, 'sf_mm'].values if 'sf_mm' in era5_a.columns else np.nan,
}, index=common_dates[idx])

# Compute derived variables
df['rainfall'] = df['tp'] - df['sf']
df['snow_frac'] = np.where(df['tp'] > 0.01, df['sf'] / df['tp'], 0)
df['doy'] = df.index.dayofyear
df['doy_c'] = df['doy'] - df['doy'].mean()
df['doy_c2'] = df['doy_c'] ** 2

# Z500 if available
if 'ncep_z500_nh' in panel_a.columns:
    df['z500'] = panel_a.loc[idx, 'ncep_z500_nh'].values
    
# Temperature anomaly
if 'ncep_t_10hpa' in panel_a.columns:
    df['t10'] = panel_a.loc[idx, 'ncep_t_10hpa'].values

df = df.dropna(subset=['avcount', 'ssw', 'tp', 'sf', 't2m'])

print(f"Sample: {len(df)} winter days, {df['ssw'].sum():.0f} SSW-window days")
print(f"{'='*70}")

# ===== Step 1: Total effect (c path) =====
from numpy.linalg import lstsq

def ols(X, y):
    """Simple OLS with intercept"""
    X = np.column_stack([np.ones(len(X)), X])
    beta, _, _, _ = lstsq(X, y, rcond=None)
    yhat = X @ beta
    resid = y - yhat
    n, k = X.shape
    se = np.sqrt(np.sum(resid**2) / (n - k) * np.diag(np.linalg.inv(X.T @ X)))
    t_stats = beta / se
    p_vals = 2 * stats.t.sf(np.abs(t_stats), n - k)
    return beta, se, t_stats, p_vals, np.corrcoef(y, yhat)[0,1]**2

# Control for DOY
controls = df[['doy_c', 'doy_c2']].values

print("\n=== TOTAL EFFECT: SSW → Avalanche count ===")
X_total = np.column_stack([df['ssw'].values, controls])
y = df['avcount'].values
beta, se, t, p, r2 = ols(X_total, y)
total_effect = beta[1]
print(f"  SSW coefficient: {total_effect:.4f} (SE={se[1]:.4f}, t={t[1]:.2f}, P={p[1]:.4e})")
print(f"  Interpretation: SSW window reduces daily count by {-total_effect:.2f} avalanches")
print(f"  R² = {r2:.4f}")

# ===== Step 2: SSW → Mediators (a paths) =====
mediators = {}
print(f"\n=== A PATHS: SSW → Mediators ===")

for name, col in [('Temperature (K)', 't2m'), ('Total precip (mm)', 'tp'),
                   ('Snowfall (mm)', 'sf'), ('Rainfall (mm)', 'rainfall'),
                   ('Snow fraction', 'snow_frac')]:
    if col in df.columns:
        X_a = np.column_stack([df['ssw'].values, controls])
        y_m = df[col].values
        beta_a, se_a, t_a, p_a, r2_a = ols(X_a, y_m)
        mediators[col] = {'name': name, 'a': beta_a[1], 'se_a': se_a[1], 'p_a': p_a[1]}
        sig = "***" if p_a[1] < 0.001 else "**" if p_a[1] < 0.01 else "*" if p_a[1] < 0.05 else "ns"
        print(f"  {name}: a={beta_a[1]:.4f} (P={p_a[1]:.4e}) {sig}")

if 'z500' in df.columns:
    X_a = np.column_stack([df['ssw'].values, controls])
    y_m = df['z500'].values[~np.isnan(df['z500'].values)]
    if len(y_m) == len(X_a):
        beta_a, se_a, t_a, p_a, r2_a = ols(X_a, y_m)
        mediators['z500'] = {'name': 'Z500 (gpm)', 'a': beta_a[1], 'se_a': se_a[1], 'p_a': p_a[1]}
        sig = "***" if p_a[1] < 0.001 else "ns"
        print(f"  Z500 (gpm): a={beta_a[1]:.4f} (P={p_a[1]:.4e}) {sig}")

# ===== Step 3: Mediator → Avalanche controlling for SSW (b paths) =====
print(f"\n=== B PATHS: Mediator → Avalanche (controlling SSW) ===")

for col, info in mediators.items():
    X_b = np.column_stack([df['ssw'].values, df[col].values, controls])
    y = df['avcount'].values
    
    # Handle NaN
    valid = ~(np.isnan(X_b).any(axis=1) | np.isnan(y))
    if valid.sum() < 100:
        continue
    
    beta_b, se_b, t_b, p_b, r2_b = ols(X_b[valid], y[valid])
    info['b'] = beta_b[2]  # mediator coefficient
    info['se_b'] = se_b[2]
    info['p_b'] = p_b[2]
    info['c_prime'] = beta_b[1]  # direct effect after controlling mediator
    info['p_c_prime'] = p_b[1]
    info['indirect'] = info['a'] * info['b']
    info['pct_mediated'] = (info['indirect'] / total_effect * 100) if total_effect != 0 else 0
    
    sig = "***" if p_b[2] < 0.001 else "**" if p_b[2] < 0.01 else "*" if p_b[2] < 0.05 else "ns"
    print(f"  {info['name']}: b={info['b']:.4f} (P={p_b[2]:.4e}) {sig}")
    print(f"    Indirect a×b = {info['indirect']:.4f} ({info['pct_mediated']:.1f}% of total)")

# ===== Step 4: Bootstrap confidence intervals for indirect effects =====
print(f"\n=== BOOTSTRAP INDIRECT EFFECTS (1000 iterations) ===")
np.random.seed(42)
n_boot = 1000

for col, info in mediators.items():
    if 'indirect' not in info:
        continue
    
    boot_indirect = []
    valid = ~(np.isnan(df[col].values) | np.isnan(df['avcount'].values))
    ssw_v = df['ssw'].values[valid]
    med_v = df[col].values[valid]
    av_v = df['avcount'].values[valid]
    ctrl_v = controls[valid]
    n_valid = valid.sum()
    
    for _ in range(n_boot):
        idx = np.random.choice(n_valid, n_valid, replace=True)
        try:
            # a path
            X_a = np.column_stack([np.ones(n_valid), ssw_v[idx], ctrl_v[idx]])
            beta_a, _, _, _ = lstsq(X_a, med_v[idx], rcond=None)
            
            # b path  
            X_b = np.column_stack([np.ones(n_valid), ssw_v[idx], med_v[idx], ctrl_v[idx]])
            beta_b, _, _, _ = lstsq(X_b, av_v[idx], rcond=None)
            
            boot_indirect.append(beta_a[1] * beta_b[2])
        except:
            continue
    
    if len(boot_indirect) > 100:
        ci_lo, ci_hi = np.percentile(boot_indirect, [2.5, 97.5])
        sig = "YES" if (ci_lo > 0 or ci_hi < 0) else "NO"
        info['boot_ci_lo'] = ci_lo
        info['boot_ci_hi'] = ci_hi
        info['boot_sig'] = sig
        print(f"  {info['name']}: indirect = {info['indirect']:.4f} "
              f"[95% CI: {ci_lo:.4f}, {ci_hi:.4f}] "
              f"Significant: {sig}")

# ===== Step 5: Combined multi-mediator model =====
print(f"\n=== MULTI-MEDIATOR MODEL ===")

# Use all significant mediators
sig_mediators = [col for col, info in mediators.items() 
                 if 'p_a' in info and info['p_a'] < 0.05]
print(f"Significant a-path mediators: {sig_mediators}")

if len(sig_mediators) > 0:
    med_cols = np.column_stack([df[col].values for col in sig_mediators])
    X_multi = np.column_stack([df['ssw'].values, med_cols, controls])
    valid = ~(np.isnan(X_multi).any(axis=1) | np.isnan(df['avcount'].values))
    
    beta_m, se_m, t_m, p_m, r2_m = ols(X_multi[valid], df['avcount'].values[valid])
    
    c_prime_multi = beta_m[1]
    pct_remaining = (c_prime_multi / total_effect * 100) if total_effect != 0 else 0
    pct_mediated_total = 100 - pct_remaining
    
    print(f"  Direct effect (c'): {c_prime_multi:.4f} (P={p_m[1]:.4e})")
    print(f"  Total effect (c): {total_effect:.4f}")
    print(f"  % mediated by all pathways: {pct_mediated_total:.1f}%")
    print(f"  % unexplained (direct): {pct_remaining:.1f}%")

# ===== Step 6: Formal trigger-day analysis =====
print(f"\n=== TRIGGER-DAY MEDIATION ===")

# A "trigger day" proxy: day where conditions favor triggering
# Rapid warming or rain-on-snow
df['warming_trigger'] = (df['t2m'] - df['t2m'].rolling(3, center=True).mean()) > 1.0
df['rain_trigger'] = df['rainfall'] > 0.5  # mm
df['any_trigger'] = df['warming_trigger'] | df['rain_trigger']

trigger_rate_ssw = df.loc[df['ssw']==1, 'any_trigger'].mean()
trigger_rate_ctrl = df.loc[df['ssw']==0, 'any_trigger'].mean()
trigger_rr = trigger_rate_ssw / max(trigger_rate_ctrl, 0.001)

# Count on trigger vs non-trigger days
count_trigger = df.loc[df['any_trigger'], 'avcount'].mean()
count_no_trigger = df.loc[~df['any_trigger'], 'avcount'].mean()

# Effect through trigger pathway
trigger_indirect = (trigger_rate_ssw - trigger_rate_ctrl) * (count_trigger - count_no_trigger)
trigger_pct = (trigger_indirect / total_effect * 100) if total_effect != 0 else 0

print(f"  Trigger day rate: SSW={trigger_rate_ssw:.3f}, Control={trigger_rate_ctrl:.3f}")
print(f"  Trigger RR = {trigger_rr:.3f}")
print(f"  Avalanche count: trigger days={count_trigger:.2f}, non-trigger={count_no_trigger:.2f}")
print(f"  Trigger pathway: {trigger_pct:.1f}% of total effect")

# Chi-square test for trigger frequency
from scipy.stats import chi2_contingency
ct = pd.crosstab(df['ssw'].astype(int), df['any_trigger'].astype(int))
if ct.shape == (2,2):
    chi2, p_chi, _, _ = chi2_contingency(ct)
    print(f"  Chi-square for trigger frequency: χ²={chi2:.2f}, P={p_chi:.4e}")

# ===== Summary =====
print(f"\n{'='*70}")
print(f"MEDIATION SUMMARY")
print(f"{'='*70}")
print(f"Total SSW→Avalanche effect: {total_effect:.4f} avalanches/day")
print(f"\nPathway decomposition:")

pathway_results = {}
for col, info in mediators.items():
    if 'indirect' in info:
        sig_marker = "✓" if info.get('boot_sig') == 'YES' else "✗"
        print(f"  {info['name']}: {info['pct_mediated']:.1f}% {sig_marker}")
        pathway_results[col] = {
            'name': info['name'],
            'a': float(info['a']),
            'b': float(info['b']),
            'indirect': float(info['indirect']),
            'pct_mediated': float(info['pct_mediated']),
            'boot_ci': [float(info.get('boot_ci_lo', 0)), float(info.get('boot_ci_hi', 0))],
            'significant': info.get('boot_sig', 'NA')
        }

if 'pct_mediated_total' in dir():
    print(f"\n  TOTAL mediated: {pct_mediated_total:.1f}%")
    print(f"  Unexplained (direct SSW effect): {pct_remaining:.1f}%")

print(f"\n  Trigger-day pathway: {trigger_pct:.1f}%")

# Save
results = {
    'total_effect': float(total_effect),
    'pathways': pathway_results,
    'trigger_pathway_pct': float(trigger_pct),
    'n_days': int(len(df)),
    'n_ssw_days': int(df['ssw'].sum()),
}

with open('data/results/r37_mediation.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nSaved to data/results/r37_mediation.json")
