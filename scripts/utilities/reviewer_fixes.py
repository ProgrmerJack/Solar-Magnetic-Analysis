"""
Reviewer-requested analyses:
1. Zero-inflated Poisson (ZIP) vs standard Poisson
2. Window sensitivity (10, 15, 20, 30 day windows)
3. Seasonal distribution of SSW events vs controls
"""
import pandas as pd
import numpy as np
from scipy import stats
from scipy.optimize import minimize
from scipy.special import gammaln
import json, warnings
warnings.filterwarnings('ignore')

panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
ssw_cat = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw_dates = ssw_cat.index.tz_localize(None)

# Load NCEP vortex data separately
ncep = pd.read_parquet('data/processed/atmospheric/ncep_stratosphere.parquet')
ncep.index = ncep.index.tz_localize(None)

panel['month'] = panel.index.month
winter = panel[panel['month'].isin([11,12,1,2,3,4])].copy()

# Merge vortex data
winter = winter.join(ncep[['uwnd_ms_10hPa']], how='left')

y_raw = winter['dry_natural_size_1234'].values
y = np.nan_to_num(y_raw, nan=0).astype(int)
vortex = winter['uwnd_ms_10hPa'].values

print(f"Winter days: {len(y)}, zeros: {(y==0).sum()} ({(y==0).mean()*100:.1f}%)")
print(f"Count mean: {y.mean():.3f}, var: {y.var():.3f}, dispersion: {y.var()/max(y.mean(),1e-10):.2f}")

# ================================================================
# 1. ZERO-INFLATED MODELS
# ================================================================
print("\n=== ZERO-INFLATED MODEL COMPARISON ===")

# Deseasonalize vortex
doy = winter.index.dayofyear.values
vortex_series = pd.Series(vortex, index=winter.index)
doy_mean = np.array([vortex_series[doy == d].mean() for d in range(1, 367)])
vortex_anom = np.array([v - doy_mean[min(d-1, 365)] if not np.isnan(v) else np.nan 
                        for v, d in zip(vortex, doy)])

# SSW indicator
ssw_indicator = np.zeros(len(winter))
for sd in ssw_dates:
    mask = (winter.index >= sd - pd.Timedelta(days=15)) & (winter.index <= sd + pd.Timedelta(days=15))
    ssw_indicator[mask] = 1

valid = ~np.isnan(vortex_anom) & (y >= 0)
y_v = y[valid]
x_v = vortex_anom[valid]
ssw_v = ssw_indicator[valid]
n = len(y_v)
print(f"Valid obs: {n}, SSW days: {int(ssw_v.sum())}")
print(f"Count zeros: {(y_v==0).sum()} ({(y_v==0).mean()*100:.1f}%)")
print(f"Count mean: {y_v.mean():.3f}, var: {y_v.var():.3f}")

# --- Standard Poisson with SSW ---
X_ssw = np.column_stack([np.ones(n), ssw_v])

def poisson_nll(params, y, X):
    eta = X @ params
    mu = np.exp(np.clip(eta, -20, 20))
    return -np.sum(y * np.log(mu + 1e-300) - mu - gammaln(y + 1))

res_pois_null = minimize(poisson_nll, [np.log(max(y_v.mean(), 0.01))], 
                         args=(y_v, np.ones((n,1))), method='Nelder-Mead')
ll_pois_null = -res_pois_null.fun
aic_pois_null = 2*1 - 2*ll_pois_null

res_pois_ssw = minimize(poisson_nll, [np.log(max(y_v.mean(), 0.01)), -0.2], 
                        args=(y_v, X_ssw), method='Nelder-Mead')
ll_pois_ssw = -res_pois_ssw.fun
aic_pois_ssw = 2*2 - 2*ll_pois_ssw
pois_irr = np.exp(res_pois_ssw.x[1])

print(f"\nPoisson null: AIC={aic_pois_null:.1f}, LL={ll_pois_null:.1f}")
print(f"Poisson SSW:  AIC={aic_pois_ssw:.1f}, LL={ll_pois_ssw:.1f}, IRR={pois_irr:.3f}")

# --- Zero-Inflated Poisson (ZIP) with SSW ---
def zip_nll(params, y, X):
    k = X.shape[1]
    beta = params[:k]
    gamma = params[k:]
    eta = X @ beta
    mu = np.exp(np.clip(eta, -20, 20))
    eta_z = X @ gamma
    pi = 1.0 / (1.0 + np.exp(-np.clip(eta_z, -20, 20)))
    
    ll = np.zeros(len(y))
    zero_mask = y == 0
    ll[zero_mask] = np.log(pi[zero_mask] + (1 - pi[zero_mask]) * np.exp(-mu[zero_mask]) + 1e-300)
    nonzero = ~zero_mask
    ll[nonzero] = (np.log(1 - pi[nonzero] + 1e-300) + y[nonzero] * np.log(mu[nonzero] + 1e-300) 
                   - mu[nonzero] - gammaln(y[nonzero] + 1))
    return -np.sum(ll)

# ZIP null (intercept only in both parts)
def zip_nll_null(params, y):
    beta0, gamma0 = params
    mu = np.exp(np.clip(beta0, -20, 20))
    pi = 1.0 / (1.0 + np.exp(-np.clip(gamma0, -20, 20)))
    ll = np.zeros(len(y))
    zero_mask = y == 0
    ll[zero_mask] = np.log(pi + (1 - pi) * np.exp(-mu) + 1e-300)
    nonzero = ~zero_mask
    ll[nonzero] = (np.log(1 - pi + 1e-300) + y[nonzero] * np.log(mu + 1e-300) 
                   - mu - gammaln(y[nonzero] + 1))
    return -np.sum(ll)

res_zip_null = minimize(zip_nll_null, [np.log(max(y_v.mean(), 0.01)), -1.0], 
                        args=(y_v,), method='Nelder-Mead', options={'maxiter': 50000})
ll_zip_null = -res_zip_null.fun
aic_zip_null = 2*2 - 2*ll_zip_null

# ZIP with SSW in both count and zero parts
init_zip = [np.log(max(y_v.mean(), 0.01)), -0.2, -1.0, 0.0]
res_zip_ssw = minimize(zip_nll, init_zip, args=(y_v, X_ssw), method='Nelder-Mead',
                       options={'maxiter': 100000, 'xatol': 1e-10, 'fatol': 1e-10})
ll_zip_ssw = -res_zip_ssw.fun
aic_zip_ssw = 2*4 - 2*ll_zip_ssw

beta_zip = res_zip_ssw.x[:2]
gamma_zip = res_zip_ssw.x[2:]
zip_irr = np.exp(beta_zip[1])
zip_pi_base = 1/(1+np.exp(-gamma_zip[0]))
zip_pi_ssw = 1/(1+np.exp(-(gamma_zip[0]+gamma_zip[1])))

print(f"\nZIP null:     AIC={aic_zip_null:.1f}, LL={ll_zip_null:.1f}")
print(f"ZIP SSW:      AIC={aic_zip_ssw:.1f}, LL={ll_zip_ssw:.1f}")
print(f"  ZIP count IRR = {zip_irr:.3f}")
print(f"  ZIP zero-inflation: base={zip_pi_base:.3f}, SSW={zip_pi_ssw:.3f}")

# Vuong test
eta_p = X_ssw @ res_pois_ssw.x
mu_p = np.exp(np.clip(eta_p, -20, 20))
ll_p_i = y_v * np.log(mu_p + 1e-300) - mu_p - gammaln(y_v + 1)

eta_c = X_ssw @ beta_zip
mu_z = np.exp(np.clip(eta_c, -20, 20))
eta_zi = X_ssw @ gamma_zip
pi_z = 1.0 / (1.0 + np.exp(-np.clip(eta_zi, -20, 20)))
ll_z_i = np.zeros(n)
zero_mask = y_v == 0
ll_z_i[zero_mask] = np.log(pi_z[zero_mask] + (1 - pi_z[zero_mask]) * np.exp(-mu_z[zero_mask]) + 1e-300)
nonzero = ~zero_mask
ll_z_i[nonzero] = (np.log(1 - pi_z[nonzero] + 1e-300) + y_v[nonzero] * np.log(mu_z[nonzero] + 1e-300) 
                   - mu_z[nonzero] - gammaln(y_v[nonzero] + 1))

m = ll_z_i - ll_p_i
if m.std() > 0:
    vuong_stat = np.sqrt(n) * m.mean() / m.std()
    vuong_p = 2 * (1 - stats.norm.cdf(abs(vuong_stat)))
else:
    vuong_stat, vuong_p = 0.0, 1.0
print(f"\nVuong test (ZIP vs Poisson): z={vuong_stat:.3f}, P={vuong_p:.4f}")

# Overdispersion
mu_hat = np.exp(np.clip(X_ssw @ res_pois_ssw.x, -20, 20))
pearson_chi2 = np.sum((y_v - mu_hat)**2 / (mu_hat + 1e-300))
disp_ratio = pearson_chi2 / (n - 2)
print(f"Overdispersion ratio: {disp_ratio:.2f}")

# LRT: ZIP SSW vs ZIP null
lrt_zip = 2 * (ll_zip_ssw - ll_zip_null)
lrt_zip_p = 1 - stats.chi2.cdf(max(lrt_zip, 0), df=2)
print(f"LRT ZIP SSW vs ZIP null: chi2={lrt_zip:.1f}, P={lrt_zip_p:.6f}")

# Key question: does the SSW effect survive in ZIP?
print(f"\n*** KEY RESULT: SSW effect in ZIP ***")
print(f"  Poisson IRR = {pois_irr:.3f}")
print(f"  ZIP count IRR = {zip_irr:.3f}")
print(f"  → SSW effect {'survives' if abs(zip_irr - 1) > 0.05 else 'disappears in'} zero-inflated model")

# ================================================================
# 2. WINDOW SENSITIVITY
# ================================================================
print("\n\n=== WINDOW SENSITIVITY ===")

swiss_ssw = sorted([d for d in ssw_dates if 1998 <= d.year <= 2019])

for window in [10, 15, 20, 30]:
    decreases = 0
    total = 0
    deltas = []
    for sd in swiss_ssw:
        post_mask = (winter.index >= sd) & (winter.index < sd + pd.Timedelta(days=window))
        if post_mask.sum() == 0:
            continue
        ssw_mean = winter.loc[post_mask, 'dry_natural_size_1234'].fillna(0).mean()
        
        # Matched control
        control_means = []
        for yr in range(1998, 2020):
            if yr == sd.year:
                continue
            for d_off in range(-3, 4):
                try:
                    ctrl_start = pd.Timestamp(year=yr, month=sd.month, day=sd.day) + pd.Timedelta(days=d_off)
                    ctrl_mask = (winter.index >= ctrl_start) & (winter.index < ctrl_start + pd.Timedelta(days=window))
                    if ctrl_mask.sum() > 0:
                        control_means.append(winter.loc[ctrl_mask, 'dry_natural_size_1234'].fillna(0).mean())
                except:
                    pass
        
        if control_means:
            ctrl_mean = np.mean(control_means)
            delta = ssw_mean - ctrl_mean
            deltas.append(delta)
            if delta < 0:
                decreases += 1
            total += 1
    
    if total > 0:
        sign_res = stats.binomtest(decreases, total, 0.5)
        print(f"Window {window:2d}d: {decreases}/{total} decrease, sign P={sign_res.pvalue:.4f}, "
              f"median Δ={np.median(deltas):.2f}, mean Δ={np.mean(deltas):.2f}")

# ================================================================
# 3. SEASONAL DISTRIBUTION
# ================================================================
print("\n\n=== SEASONAL DISTRIBUTION ===")

ssw_months = [d.month for d in swiss_ssw]
ssw_doys = [d.dayofyear for d in swiss_ssw]
month_names = ['','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
print("SSW events by month:")
for m in [11, 12, 1, 2, 3]:
    count = ssw_months.count(m)
    if count > 0:
        print(f"  {month_names[m]}: {count}")

print(f"\nBaseline dry slab rate by month (winter):")
for m in [11, 12, 1, 2, 3, 4]:
    mdata = winter[winter['month'] == m]['dry_natural_size_1234'].fillna(0)
    print(f"  {month_names[m]}: mean={mdata.mean():.3f}, n={len(mdata)}")

# Save results
results = {
    'zip_comparison': {
        'n_obs': n,
        'n_zeros': int((y_v==0).sum()),
        'pct_zeros': round((y_v==0).mean()*100, 1),
        'poisson_null_aic': round(aic_pois_null, 1),
        'poisson_ssw_aic': round(aic_pois_ssw, 1),
        'poisson_irr': round(float(pois_irr), 3),
        'zip_null_aic': round(aic_zip_null, 1),
        'zip_ssw_aic': round(aic_zip_ssw, 1),
        'zip_count_irr': round(float(zip_irr), 3),
        'zip_zero_infl_base': round(float(zip_pi_base), 3),
        'zip_zero_infl_ssw': round(float(zip_pi_ssw), 3),
        'vuong_z': round(float(vuong_stat), 3),
        'vuong_p': round(float(vuong_p), 4),
        'overdispersion': round(float(disp_ratio), 2),
        'lrt_zip_ssw_chi2': round(float(lrt_zip), 1),
        'lrt_zip_ssw_p': round(float(lrt_zip_p), 6),
    }
}

with open('data/results/reviewer_fixes.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved.")
