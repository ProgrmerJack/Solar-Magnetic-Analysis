"""
Threshold vs Continuous model comparison for SSW-avalanche association.
Compares: (1) continuous vortex metric, (2) binary SSW indicator, (3) threshold model.
Uses AIC/BIC to assess which fits best.
Also tests precursor timing by finding vortex weakening onset before canonical SSW date.
"""
import pandas as pd
import numpy as np
from scipy import stats
import json, warnings
warnings.filterwarnings('ignore')

print("Loading data...")
panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
strat = pd.read_parquet('data/processed/atmospheric/ncep_stratosphere.parquet')
ssw_cat = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')

# Align dates
panel_dates = panel.index
strat['date'] = pd.to_datetime(strat.index)
if strat['date'].dt.tz is not None:
    strat['date'] = strat['date'].dt.tz_localize(None)
strat = strat.set_index('date')

# Get Swiss-period SSW events
ssw_dates = pd.to_datetime(ssw_cat.index).tz_localize(None)
swiss_mask = (ssw_dates >= panel_dates.min()) & (ssw_dates <= panel_dates.max())
swiss_ssw = ssw_dates[swiss_mask].sort_values()
swiss_ssw = pd.Series(swiss_ssw)
print("Swiss SSW events:", len(swiss_ssw))

# Winter days only (Nov-Apr)
winter_mask = panel.index.month.isin([11, 12, 1, 2, 3, 4])
winter_panel = panel[winter_mask].copy()

# Merge vortex data
vortex_col = 'uwnd_ms_10hPa'
common_dates = winter_panel.index.intersection(strat.index)
winter_panel = winter_panel.loc[common_dates]
winter_panel['vortex'] = strat.loc[common_dates, vortex_col].values

# Day-of-year climatology for deseasonalizing
doy = winter_panel.index.dayofyear
clim = winter_panel.groupby(doy)['vortex'].transform('mean')
winter_panel['vortex_anom'] = winter_panel['vortex'] - clim

# Avalanche column
aval_col = 'dry_natural_size_1234'

# Create SSW window indicators
winter_panel['ssw_window'] = 0  # binary: within 15d of SSW
winter_panel['days_since_ssw'] = np.nan
for sd in swiss_ssw:
    delta = (winter_panel.index - sd).days
    mask_window = (delta >= -15) & (delta <= 15)
    winter_panel.loc[mask_window, 'ssw_window'] = 1
    for idx in winter_panel.index[mask_window]:
        d = (idx - sd).days
        if pd.isna(winter_panel.loc[idx, 'days_since_ssw']) or abs(d) < abs(winter_panel.loc[idx, 'days_since_ssw']):
            winter_panel.at[idx, 'days_since_ssw'] = d

print("\n=== MODEL COMPARISON: Threshold vs Continuous ===")
print("(Predicting daily dry slab avalanche counts)")

y = winter_panel[aval_col].values.astype(float)
valid = ~np.isnan(y) & ~np.isnan(winter_panel['vortex_anom'].values)
y_v = y[valid]
n = len(y_v)

# Model 1: Intercept only (null)
mu0 = np.mean(y_v)
ll_null = np.sum(stats.poisson.logpmf(y_v.astype(int), mu0))
aic_null = -2 * ll_null + 2 * 1
bic_null = -2 * ll_null + np.log(n) * 1

# Model 2: Continuous vortex anomaly (linear Poisson)
from numpy.polynomial import polynomial as P
# Simple approach: exp(a + b*x) Poisson
# Use iterative reweighted least squares approximation
x_cont = winter_panel['vortex_anom'].values[valid]
# Fit via GLM-like: log(mu) = a + b*x
# Quick Newton-Raphson for Poisson GLM
def poisson_glm_fit(X, y, max_iter=50):
    n_feat = X.shape[1]
    beta = np.zeros(n_feat)
    for _ in range(max_iter):
        eta = X @ beta
        mu = np.exp(np.clip(eta, -20, 20))
        W = np.diag(mu)
        z = eta + (y - mu) / mu
        try:
            beta_new = np.linalg.solve(X.T @ W @ X, X.T @ W @ z)
        except np.linalg.LinAlgError:
            break
        if np.max(np.abs(beta_new - beta)) < 1e-8:
            beta = beta_new
            break
        beta = beta_new
    eta = X @ beta
    mu = np.exp(np.clip(eta, -20, 20))
    ll = np.sum(stats.poisson.logpmf(y.astype(int), mu))
    return beta, ll

# Too memory-intensive with full diagonal W for large n. Use simpler approach.
# Statsmodels-free Poisson GLM via scipy minimize
from scipy.optimize import minimize

def neg_poisson_ll(beta, X, y):
    eta = X @ beta
    mu = np.exp(np.clip(eta, -20, 20))
    return -np.sum(y * np.log(mu + 1e-10) - mu)

X_cont = np.column_stack([np.ones(n), x_cont])
res_cont = minimize(neg_poisson_ll, [np.log(mu0), 0], args=(X_cont, y_v), method='Nelder-Mead')
ll_cont = -res_cont.fun
aic_cont = -2 * ll_cont + 2 * 2
bic_cont = -2 * ll_cont + np.log(n) * 2

# Model 3: Binary SSW window indicator
x_ssw = winter_panel['ssw_window'].values[valid].astype(float)
X_ssw = np.column_stack([np.ones(n), x_ssw])
res_ssw = minimize(neg_poisson_ll, [np.log(mu0), 0], args=(X_ssw, y_v), method='Nelder-Mead')
ll_ssw = -res_ssw.fun
aic_ssw = -2 * ll_ssw + 2 * 2
bic_ssw = -2 * ll_ssw + np.log(n) * 2

# Model 4: Threshold (SSW binary) + continuous interaction
X_both = np.column_stack([np.ones(n), x_ssw, x_cont, x_ssw * x_cont])
res_both = minimize(neg_poisson_ll, [np.log(mu0), 0, 0, 0], args=(X_both, y_v), method='Nelder-Mead')
ll_both = -res_both.fun
aic_both = -2 * ll_both + 2 * 4
bic_both = -2 * ll_both + np.log(n) * 4

# Model 5: Vortex quintile (stepped threshold)
quintiles = pd.qcut(x_cont, 5, labels=False)
X_quint = np.column_stack([np.ones(n)] + [((quintiles == q).astype(float)) for q in range(1, 5)])
res_quint = minimize(neg_poisson_ll, np.zeros(5), args=(X_quint, y_v), method='Nelder-Mead')
ll_quint = -res_quint.fun
aic_quint = -2 * ll_quint + 2 * 5
bic_quint = -2 * ll_quint + np.log(n) * 5

print(f"\nN = {n} winter days")
print(f"{'Model':<35} {'LL':>10} {'AIC':>10} {'BIC':>10} {'dAIC':>8} {'dBIC':>8}")
print("-" * 85)
models = [
    ("Null (intercept only)", ll_null, aic_null, bic_null, 1),
    ("Continuous vortex anomaly", ll_cont, aic_cont, bic_cont, 2),
    ("Binary SSW indicator", ll_ssw, aic_ssw, bic_ssw, 2),
    ("SSW + vortex + interaction", ll_both, aic_both, bic_both, 4),
    ("Vortex quintiles", ll_quint, aic_quint, bic_quint, 5),
]

best_aic = min(m[2] for m in models)
best_bic = min(m[3] for m in models)

for name, ll, aic, bic, k in models:
    print(f"{name:<35} {ll:>10.1f} {aic:>10.1f} {bic:>10.1f} {aic-best_aic:>8.1f} {bic-best_bic:>8.1f}")

# LRT: SSW indicator vs null
lr_stat = 2 * (ll_ssw - ll_null)
lr_p = 1 - stats.chi2.cdf(lr_stat, 1)
print(f"\nLRT SSW vs Null: chi2 = {lr_stat:.2f}, P = {lr_p:.4f}")

# LRT: Continuous vs null
lr_stat2 = 2 * (ll_cont - ll_null)
lr_p2 = 1 - stats.chi2.cdf(lr_stat2, 1)
print(f"LRT Continuous vs Null: chi2 = {lr_stat2:.2f}, P = {lr_p2:.4f}")

# LRT: SSW vs Continuous (non-nested, use AIC/BIC comparison)
print(f"\nSSW vs Continuous (same df):")
print(f"  AIC difference: {aic_cont - aic_ssw:+.1f} (positive = SSW better)")
print(f"  BIC difference: {bic_cont - bic_ssw:+.1f} (positive = SSW better)")

# IRR for SSW indicator
beta_ssw = res_ssw.x
irr_ssw = np.exp(beta_ssw[1])
print(f"\nSSW indicator IRR: {irr_ssw:.3f} (rate ratio during SSW window)")

print("\n=== PRECURSOR TIMING ANALYSIS ===")
print("Finding vortex weakening onset before canonical SSW date...")

precursor_lags = []
for sd in swiss_ssw:
    # Look at vortex data -60 to 0 days relative to SSW
    pre_dates = pd.date_range(sd - pd.Timedelta(days=60), sd, freq='D')
    pre_dates = pre_dates.intersection(strat.index)
    if len(pre_dates) < 30:
        continue
    
    vortex_series = strat.loc[pre_dates, vortex_col]
    # Find first day where vortex drops below 1 std below climatology
    doy_vals = pre_dates.dayofyear
    clim_vals = strat.groupby(strat.index.dayofyear)[vortex_col].mean()
    std_vals = strat.groupby(strat.index.dayofyear)[vortex_col].std()
    
    onset_day = None
    for i, d in enumerate(pre_dates):
        dy = d.dayofyear
        if dy in clim_vals.index:
            threshold = clim_vals[dy] - 1.0 * std_vals[dy]
            if vortex_series.iloc[i] < threshold:
                onset_day = d
                break
    
    if onset_day is not None:
        lag = (sd - onset_day).days
        precursor_lags.append(lag)
        print(f"  {sd.strftime('%Y-%m-%d')}: vortex weakening onset {lag}d before SSW")
    else:
        precursor_lags.append(0)
        print(f"  {sd.strftime('%Y-%m-%d')}: no clear precursor found")

precursor_lags = np.array(precursor_lags)
print(f"\nPrecursor lag statistics:")
print(f"  Mean: {np.mean(precursor_lags):.1f} days before SSW")
print(f"  Median: {np.median(precursor_lags):.1f} days")
print(f"  Range: {np.min(precursor_lags)}-{np.max(precursor_lags)} days")

# Test if anchoring to precursor onset changes the avalanche anomaly timing
print("\n=== AVALANCHE ANOMALY RELATIVE TO PRECURSOR ONSET ===")
# Compute mean avalanche anomaly in windows relative to precursor onset
doy_aval_clim = panel.groupby(panel.index.dayofyear)[aval_col].mean()

pre_onset_anom = []
post_onset_anom = []
for i, sd in enumerate(swiss_ssw):
    if i >= len(precursor_lags):
        break
    onset = sd - pd.Timedelta(days=int(precursor_lags[i]))
    
    # Pre-onset: -15 to -1 days
    pre_dates = pd.date_range(onset - pd.Timedelta(days=15), onset - pd.Timedelta(days=1))
    pre_dates = pre_dates.intersection(panel.index)
    if len(pre_dates) > 0:
        vals = panel.loc[pre_dates, aval_col]
        clim_v = [doy_aval_clim.get(d.dayofyear, np.nan) for d in pre_dates]
        anom = np.nanmean(vals) - np.nanmean(clim_v)
        pre_onset_anom.append(anom)
    
    # Post-onset: 0 to 14 days  
    post_dates = pd.date_range(onset, onset + pd.Timedelta(days=14))
    post_dates = post_dates.intersection(panel.index)
    if len(post_dates) > 0:
        vals = panel.loc[post_dates, aval_col]
        clim_v = [doy_aval_clim.get(d.dayofyear, np.nan) for d in post_dates]
        anom = np.nanmean(vals) - np.nanmean(clim_v)
        post_onset_anom.append(anom)

pre_onset_anom = np.array(pre_onset_anom)
post_onset_anom = np.array(post_onset_anom)

print(f"Relative to PRECURSOR onset:")
print(f"  Pre-onset anomaly: {np.mean(pre_onset_anom):.3f} aval/day")
print(f"  Post-onset anomaly: {np.mean(post_onset_anom):.3f} aval/day")
t_pre, p_pre = stats.ttest_1samp(pre_onset_anom, 0)
t_post, p_post = stats.ttest_1samp(post_onset_anom, 0)
print(f"  Pre-onset P: {p_pre:.4f}")
print(f"  Post-onset P: {p_post:.4f}")

# Compare pre vs post (is there now asymmetry?)
t_diff, p_diff = stats.ttest_rel(pre_onset_anom[:len(post_onset_anom)], post_onset_anom[:len(pre_onset_anom)])
print(f"  Pre vs Post difference P: {p_diff:.4f}")

# Save results
results = {
    "model_comparison": {
        "n_days": int(n),
        "null": {"LL": float(ll_null), "AIC": float(aic_null), "BIC": float(bic_null)},
        "continuous_vortex": {"LL": float(ll_cont), "AIC": float(aic_cont), "BIC": float(bic_cont)},
        "binary_ssw": {"LL": float(ll_ssw), "AIC": float(aic_ssw), "BIC": float(bic_ssw), "IRR": float(irr_ssw)},
        "ssw_plus_interaction": {"LL": float(ll_both), "AIC": float(aic_both), "BIC": float(bic_both)},
        "vortex_quintiles": {"LL": float(ll_quint), "AIC": float(aic_quint), "BIC": float(bic_quint)},
        "lrt_ssw_vs_null": {"chi2": float(lr_stat), "P": float(lr_p)},
        "lrt_cont_vs_null": {"chi2": float(lr_stat2), "P": float(lr_p2)},
        "aic_ssw_minus_cont": float(aic_ssw - aic_cont),
        "bic_ssw_minus_cont": float(bic_ssw - bic_cont),
        "winner": "binary_ssw" if aic_ssw < aic_cont else "continuous"
    },
    "precursor_timing": {
        "mean_lag_days": float(np.mean(precursor_lags)),
        "median_lag_days": float(np.median(precursor_lags)),
        "range": [int(np.min(precursor_lags)), int(np.max(precursor_lags))],
        "individual_lags": [int(x) for x in precursor_lags],
        "pre_onset_anomaly": float(np.mean(pre_onset_anom)),
        "post_onset_anomaly": float(np.mean(post_onset_anom)),
        "pre_onset_P": float(p_pre),
        "post_onset_P": float(p_post),
        "pre_vs_post_P": float(p_diff)
    }
}

with open('data/results/threshold_vs_continuous.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nResults saved to data/results/threshold_vs_continuous.json")
print("DONE")
