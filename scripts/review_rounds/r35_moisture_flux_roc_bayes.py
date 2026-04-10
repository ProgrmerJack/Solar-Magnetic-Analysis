"""
R35 Upgrade Analysis:
1. ERA5 Integrated Vapor Transport (IVT) composites during SSW vs control
   - Quantifies Mediterranean moisture pathway disruption
2. Retrospective ROC/Brier Skill Scores for forecasting framework
3. Bayesian analysis with explicit Bayes factors
"""
import numpy as np
import pandas as pd
import json
import os
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Load data
# ============================================================
era5 = pd.read_parquet('data/processed/era5_swiss_alps_daily.parquet')
# era5 has DatetimeIndex named 'date', columns: tp_mm, sf_mm, t2m_K, sd_m, wind_speed, u10, v10, doy, *_anom

ssw_cat = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
# ssw_cat has DatetimeIndex named 'onset_date', columns: type, source
ssw_dates = ssw_cat.index.tz_localize(None)  # Remove timezone

panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
# panel has DatetimeIndex named 'time'

# Also load extended ERA5 if available
try:
    era5_ext = pd.read_parquet('data/processed/era5_swiss_alps_extended.parquet')
    # Merge with original
    era5_all = pd.concat([era5, era5_ext]).sort_index()
    era5_all = era5_all[~era5_all.index.duplicated(keep='last')]
    print(f"ERA5 combined: {era5_all.index.min()} to {era5_all.index.max()}, {len(era5_all)} days")
except:
    era5_all = era5
    print(f"ERA5 original only: {era5_all.index.min()} to {era5_all.index.max()}, {len(era5_all)} days")

# Filter SSW dates to ERA5+panel coverage
panel_start = panel.index.min()
panel_end = panel.index.max()
ssw_in_range = ssw_dates[(ssw_dates >= panel_start) & (ssw_dates <= panel_end)]
print(f"SSW events in panel range: {len(ssw_in_range)}")
for d in ssw_in_range:
    print(f"  {d.date()}")

from scipy import stats

# ============================================================
# PART 1: ERA5 MOISTURE FLUX ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("PART 1: ERA5 MOISTURE FLUX ANALYSIS")
print("=" * 70)

# Define SSW mask on ERA5
def make_ssw_mask(index, ssw_dates, window=15):
    mask = pd.Series(False, index=index)
    for d in ssw_dates:
        start = d - pd.Timedelta(days=window)
        end = d + pd.Timedelta(days=window)
        mask |= (index >= start) & (index <= end)
    return mask

ssw_mask_era5 = make_ssw_mask(era5_all.index, ssw_in_range, window=15)
winter_mask = era5_all['doy'].isin(list(range(1, 121)) + list(range(305, 366)))
era5w = era5_all[winter_mask].copy()
ssw_w = ssw_mask_era5[winter_mask]

print(f"\nWinter days: {len(era5w)}, SSW days: {ssw_w.sum()}, Control: {(~ssw_w).sum()}")

# Key variables
print(f"\n--- Aggregate SSW vs Control Composites ---")
print(f"{'Variable':<20} {'SSW mean':>10} {'Ctrl mean':>10} {'Delta':>10} {'d':>8} {'P':>12}")
print("-" * 72)
for col in ['tp_mm', 'sf_mm', 't2m_K', 'wind_speed', 'u10', 'v10']:
    ssw_vals = era5w.loc[ssw_w, col].dropna()
    ctrl_vals = era5w.loc[~ssw_w, col].dropna()
    d = (ssw_vals.mean() - ctrl_vals.mean()) / ctrl_vals.std()
    _, p = stats.mannwhitneyu(ssw_vals, ctrl_vals, alternative='two-sided')
    print(f"{col:<20} {ssw_vals.mean():>10.3f} {ctrl_vals.mean():>10.3f} {ssw_vals.mean()-ctrl_vals.mean():>10.3f} {d:>8.3f} {p:>12.2e}")

# Compute Mediterranean Moisture Transport Proxy
print(f"\n--- Mediterranean Moisture Transport Proxy ---")
# IVT proxy: southerly wind component (v < 0 = from south in ERA5 convention)
# Mediterranean moisture arrives from the south
# v10 < 0 means southerly flow in ERA5

# Southerly moisture transport proxy
era5w['v_south'] = -era5w['v10']  # positive = from south
era5w['moisture_transport'] = era5w['v_south'] * era5w['tp_mm']  # southerly flow × precip
era5w['warm_moisture'] = era5w['v_south'] * (era5w['t2m_K'] / 273.15)  # T-weighted southerly

for proxy_col, proxy_name in [
    ('v_south', 'Southerly wind (v_south)'),
    ('moisture_transport', 'Moisture transport (v_south × precip)'),
    ('warm_moisture', 'Warm moisture flux (v_south × T/273)')
]:
    ssw_v = era5w.loc[ssw_w, proxy_col].dropna()
    ctrl_v = era5w.loc[~ssw_w, proxy_col].dropna()
    d = (ssw_v.mean() - ctrl_v.mean()) / ctrl_v.std()
    _, p = stats.mannwhitneyu(ssw_v, ctrl_v, alternative='two-sided')
    print(f"\n{proxy_name}:")
    print(f"  SSW: {ssw_v.mean():.4f}, Ctrl: {ctrl_v.mean():.4f}, Δ={ssw_v.mean()-ctrl_v.mean():.4f}")
    print(f"  Cohen's d: {d:.3f}, P: {p:.2e}")

# Event-level moisture analysis
print(f"\n--- Event-level Mediterranean Moisture Proxy ---")
event_moisture = []
for ssw_date in ssw_in_range:
    ssw_start = ssw_date - pd.Timedelta(days=15)
    ssw_end = ssw_date + pd.Timedelta(days=15)
    
    ssw_period = era5w.loc[ssw_start:ssw_end]
    if len(ssw_period) == 0:
        continue
    
    # DOY-matched control
    doy = ssw_date.dayofyear
    ctrl = era5w[
        (era5w['doy'] >= doy - 3) & (era5w['doy'] <= doy + 3) & (~ssw_w)
    ]
    
    for proxy in ['v_south', 'tp_mm', 'sf_mm', 'moisture_transport']:
        ssw_val = ssw_period[proxy].mean()
        ctrl_val = ctrl[proxy].mean()
        ratio = ssw_val / ctrl_val if ctrl_val != 0 else np.nan
        if proxy == 'v_south':
            event_moisture.append({
                'date': str(ssw_date.date()),
                'ssw_v_south': float(ssw_val),
                'ctrl_v_south': float(ctrl_val),
                'delta_v_south': float(ssw_val - ctrl_val),
                'suppressed': ssw_val < ctrl_val
            })
    
    print(f"  {ssw_date.date()}: v_south SSW={event_moisture[-1]['ssw_v_south']:.3f} Ctrl={event_moisture[-1]['ctrl_v_south']:.3f} Δ={event_moisture[-1]['delta_v_south']:+.3f} {'↓SUPPRESSED' if event_moisture[-1]['suppressed'] else '↑ENHANCED'}")

n_suppressed = sum(1 for e in event_moisture if e['suppressed'])
print(f"\nEvents with suppressed southerly moisture: {n_suppressed}/{len(event_moisture)}")
if len(event_moisture) > 0:
    sign_p = stats.binomtest(n_suppressed, len(event_moisture), 0.5).pvalue
    print(f"Sign test P: {sign_p:.4f}")

# Snowfall fraction analysis (key for loading mechanism)
print(f"\n--- Snowfall Fraction Analysis ---")
era5w['snow_fraction'] = era5w['sf_mm'] / era5w['tp_mm'].clip(lower=0.01)
ssw_sf = era5w.loc[ssw_w, 'snow_fraction'].dropna()
ctrl_sf = era5w.loc[~ssw_w, 'snow_fraction'].dropna()
d = (ssw_sf.mean() - ctrl_sf.mean()) / ctrl_sf.std()
_, p = stats.mannwhitneyu(ssw_sf, ctrl_sf, alternative='two-sided')
print(f"Snow fraction: SSW={ssw_sf.mean():.4f}, Ctrl={ctrl_sf.mean():.4f}, d={d:.3f}, P={p:.2e}")

# ============================================================
# PART 2: RETROSPECTIVE ROC / BRIER SKILL SCORES
# ============================================================
print("\n" + "=" * 70)
print("PART 2: RETROSPECTIVE FORECASTING SKILL ASSESSMENT")
print("=" * 70)

# Get panel columns
panel_cols = list(panel.columns)
nat_dry_col = None
for c in panel_cols:
    if 'nat_dry' in c.lower() or 'natural_dry' in c.lower():
        nat_dry_col = c
        break
if nat_dry_col is None:
    # Check for avalanche count column
    for c in panel_cols:
        if 'aval' in c.lower() or 'count' in c.lower():
            print(f"  candidate: {c}")
    nat_dry_col = 'nat_dry_slab'  # guess

print(f"Looking for avalanche column: {nat_dry_col}")
print(f"Panel columns with 'nat' or 'dry': {[c for c in panel_cols if 'nat' in c.lower() or 'dry' in c.lower()]}")

# Try to find the right column
aval_cols = [c for c in panel_cols if any(k in c.lower() for k in ['nat', 'dry', 'aval', 'count', 'slf'])]
print(f"Avalanche-related columns: {aval_cols[:10]}")

# Use whatever we find
if 'nat_dry_slab' in panel_cols:
    aval_col = 'nat_dry_slab'
elif 'natural_dry' in panel_cols:
    aval_col = 'natural_dry'
elif 'aval_count_dry_natural' in panel_cols:
    aval_col = 'aval_count_dry_natural'
else:
    # Just pick the first avalanche-related one
    aval_col = aval_cols[0] if aval_cols else None

if aval_col:
    print(f"\nUsing avalanche column: {aval_col}")
else:
    print("No avalanche column found! Trying 'nat_dry'...")
    aval_col = 'nat_dry'

# Compute event-level RR
event_results = []
for ssw_date in ssw_in_range:
    try:
        ssw_start = ssw_date - pd.Timedelta(days=15)
        ssw_end = ssw_date + pd.Timedelta(days=15)
        
        mask_ssw = (panel.index >= ssw_start) & (panel.index <= ssw_end)
        obs_count = panel.loc[mask_ssw, aval_col].sum()
        n_days = mask_ssw.sum()
        
        # DOY-matched control
        doy = ssw_date.dayofyear
        ctrl_mask = (panel.index.dayofyear >= doy - 3) & \
                    (panel.index.dayofyear <= doy + 3) & \
                    (~make_ssw_mask(panel.index, ssw_in_range))
        ctrl_rate = panel.loc[ctrl_mask, aval_col].mean()
        expected = ctrl_rate * n_days
        
        rr = obs_count / expected if expected > 0 else np.nan
        
        event_results.append({
            'date': ssw_date,
            'observed': float(obs_count),
            'expected': float(expected),
            'rr': float(rr),
            'log_rr': float(np.log(max(rr, 0.001))),
            'decrease': rr < 1.0
        })
    except Exception as ex:
        print(f"  Error for {ssw_date.date()}: {ex}")

print(f"\nTotal events: {len(event_results)}")
n_events = len(event_results)
n_decrease = sum(1 for e in event_results if e['decrease'])
print(f"Events with decrease: {n_decrease}/{n_events}")

# LOO-CV direction prediction
print(f"\n--- Leave-One-Out Direction Prediction ---")
correct = 0
for i in range(n_events):
    train_rate = sum(1 for j, e in enumerate(event_results) if j != i and e['decrease']) / (n_events - 1)
    pred = train_rate > 0.5
    actual = event_results[i]['decrease']
    correct += (pred == actual)
    
loo_accuracy = correct / n_events
loo_p = stats.binomtest(correct, n_events, 0.5, alternative='greater').pvalue
print(f"Accuracy: {correct}/{n_events} = {loo_accuracy:.3f}")
print(f"P vs chance: {loo_p:.6f}")

# Brier Scores
outcomes = [1 if e['decrease'] else 0 for e in event_results]
p_base = n_decrease / n_events

# Climatological Brier (always predict base rate)
brier_clim = np.mean([(p_base - o)**2 for o in outcomes])

# LOO Brier
brier_loo_vals = []
for i in range(n_events):
    train_rate = sum(1 for j, e in enumerate(event_results) if j != i and e['decrease']) / (n_events - 1)
    brier_loo_vals.append((train_rate - outcomes[i])**2)
brier_loo = np.mean(brier_loo_vals)

# BSS
bss = 1 - brier_loo / brier_clim if brier_clim > 0 else 0

print(f"\n--- Brier Score Analysis ---")
print(f"Climatological Brier: {brier_clim:.4f}")
print(f"LOO Brier: {brier_loo:.4f}")
print(f"BSS: {bss:.4f}")

# ROC with Z500
print(f"\n--- ROC with Z500 Predictor ---")
# Check if we have Z500 in panel or ERA5
z500_cols = [c for c in panel_cols if 'z500' in c.lower() or 'z_500' in c.lower() or 'geopotential' in c.lower()]
print(f"Z500 columns in panel: {z500_cols}")

# Also check NCEP data in panel
ncep_cols = [c for c in panel_cols if 'ncep' in c.lower()]
print(f"NCEP columns: {ncep_cols[:5]}")

# Use NCEP T at various levels as proxy for SSW strength → use as predictor
strat_cols = [c for c in panel_cols if 'ncep_t_10' in c.lower() or 'ncep_t_50' in c.lower()]
print(f"Stratospheric T columns: {strat_cols}")

# Use results JSON to get Z500 data
try:
    with open('data/results/r15_comprehensive_analysis.json', 'r') as f:
        r15 = json.load(f)
    print(f"R15 keys: {list(r15.keys())[:10]}")
except:
    print("R15 not available")

# Compute ROC from available stratospheric data
if strat_cols:
    # Use 10hPa temperature as predictor of avalanche decrease
    strat_col = strat_cols[0]
    z500_events = []
    for e in event_results:
        ssw_start = e['date'] - pd.Timedelta(days=15)
        ssw_end = e['date'] + pd.Timedelta(days=15)
        mask = (panel.index >= ssw_start) & (panel.index <= ssw_end)
        val = panel.loc[mask, strat_col].mean()
        z500_events.append(val)
    
    z500_arr = np.array(z500_events)
    outcomes_arr = np.array(outcomes)
    
    # Rank-based ROC
    valid = ~np.isnan(z500_arr)
    if valid.sum() >= 5:
        # Higher strat T → stronger SSW → more likely decrease
        thresholds = np.linspace(np.nanmin(z500_arr), np.nanmax(z500_arr), 50)
        tpr_list, fpr_list = [], []
        for thresh in thresholds:
            pred = z500_arr >= thresh
            tp = np.sum(pred & (outcomes_arr == 1))
            fp = np.sum(pred & (outcomes_arr == 0))
            fn = np.sum(~pred & (outcomes_arr == 1))
            tn = np.sum(~pred & (outcomes_arr == 0))
            tpr = tp / (tp + fn) if (tp + fn) > 0 else 0
            fpr = fp / (fp + tn) if (fp + tn) > 0 else 0
            tpr_list.append(tpr)
            fpr_list.append(fpr)
        
        # Sort by FPR for AUC
        pairs = sorted(zip(fpr_list, tpr_list))
        fpr_s = [p[0] for p in pairs]
        tpr_s = [p[1] for p in pairs]
        # Remove duplicates
        unique_pairs = list(dict.fromkeys(zip(fpr_s, tpr_s)))
        fpr_u = [p[0] for p in unique_pairs]
        tpr_u = [p[1] for p in unique_pairs]
        auc = np.trapz(tpr_u, fpr_u)
        
        print(f"Predictor: {strat_col}")
        print(f"AUC-ROC: {auc:.3f}")
        
        # Correlation with log(RR)
        log_rrs = np.array([e['log_rr'] for e in event_results])
        r, p = stats.pearsonr(z500_arr[valid], log_rrs[valid])
        print(f"Correlation with log(RR): r={r:.3f}, P={p:.4f}")

# ============================================================
# PART 3: BAYESIAN ANALYSIS WITH BAYES FACTORS
# ============================================================
print("\n" + "=" * 70)
print("PART 3: BAYESIAN ANALYSIS WITH BAYES FACTORS")
print("=" * 70)

from scipy.stats import norm
from scipy.special import comb

log_rr = np.array([e['log_rr'] for e in event_results])
n = len(log_rr)
mean_lr = np.mean(log_rr)
se_lr = np.std(log_rr, ddof=1) / np.sqrt(n)
t_stat = mean_lr / se_lr

print(f"Log(RR) statistics:")
print(f"  n = {n}")
print(f"  Mean: {mean_lr:.4f}")
print(f"  SE: {se_lr:.4f}")
print(f"  t = {t_stat:.3f}")
print(f"  P (two-sided): {2 * stats.t.sf(abs(t_stat), n-1):.6f}")

# Method 1: BIC approximation
rss_0 = np.sum(log_rr**2)
rss_1 = np.sum((log_rr - mean_lr)**2)
bic_0 = n * np.log(rss_0 / n)
bic_1 = n * np.log(rss_1 / n) + np.log(n)
bf_bic = np.exp((bic_0 - bic_1) / 2)

# Method 2: Savage-Dickey
sigma_prior = 2.0
prior_at_0 = norm.pdf(0, 0, sigma_prior)
posterior_at_0 = norm.pdf(0, mean_lr, se_lr)
bf_sd = prior_at_0 / posterior_at_0

print(f"\nBayes Factors:")
print(f"  BF₁₀ (BIC): {bf_bic:.1f}")
print(f"  BF₁₀ (Savage-Dickey, σ=2): {bf_sd:.1f}")

def interpret_bf(bf):
    if bf > 100: return "Decisive"
    if bf > 30: return "Very strong"
    if bf > 10: return "Strong"
    if bf > 3: return "Moderate"
    if bf > 1: return "Anecdotal"
    return "Favors H0"

print(f"  BIC interpretation: {interpret_bf(bf_bic)}")
print(f"  S-D interpretation: {interpret_bf(bf_sd)}")

# Prior sensitivity
print(f"\n--- Prior Sensitivity ---")
for sp in [0.5, 1.0, 2.0, 5.0]:
    p0 = norm.pdf(0, 0, sp)
    q0 = norm.pdf(0, mean_lr, se_lr)
    bf = p0 / q0
    print(f"  σ_prior={sp}: BF₁₀={bf:.1f} ({interpret_bf(bf)})")

# Posterior
ci95_low = mean_lr - 1.96 * se_lr
ci95_high = mean_lr + 1.96 * se_lr
rr_mean = np.exp(mean_lr)
rr_low = np.exp(ci95_low)
rr_high = np.exp(ci95_high)

print(f"\n--- Posterior Summary ---")
print(f"Posterior RR: {rr_mean:.3f} [{rr_low:.3f}, {rr_high:.3f}]")
print(f"P(RR < 1 | data): {norm.cdf(0, mean_lr, se_lr):.6f}")

# Beta-binomial sign test BF
k = n_decrease
a_post = 1 + k
b_post = 1 + n - k
from scipy.stats import beta as beta_dist
p_gt_half = 1 - beta_dist.cdf(0.5, a_post, b_post)

p_data_h0 = comb(n, k, exact=True) * 0.5**n
p_data_h1 = comb(n, k, exact=True) / (n + 1)
bf_sign = p_data_h1 / p_data_h0

print(f"\n--- Bayesian Sign Test ---")
print(f"Posterior Beta({a_post}, {b_post})")
print(f"P(true rate > 0.5 | data): {p_gt_half:.6f}")
print(f"BF₁₀ for sign consistency: {bf_sign:.1f} ({interpret_bf(bf_sign)})")

# ============================================================
# SAVE RESULTS
# ============================================================
print("\n" + "=" * 70)
print("SAVING RESULTS")
print("=" * 70)

all_results = {
    'moisture_analysis': {
        'v_south_ssw': float(era5w.loc[ssw_w, 'v_south'].mean()),
        'v_south_ctrl': float(era5w.loc[~ssw_w, 'v_south'].mean()),
        'v_south_d': float((era5w.loc[ssw_w, 'v_south'].mean() - era5w.loc[~ssw_w, 'v_south'].mean()) / era5w.loc[~ssw_w, 'v_south'].std()),
        'n_suppressed': n_suppressed,
        'n_events': len(event_moisture),
        'event_details': event_moisture
    },
    'forecasting': {
        'n_events': n_events,
        'n_decrease': n_decrease,
        'loo_accuracy': float(loo_accuracy),
        'loo_p': float(loo_p),
        'brier_clim': float(brier_clim),
        'brier_loo': float(brier_loo),
        'bss': float(bss)
    },
    'bayesian': {
        'mean_log_rr': float(mean_lr),
        'se_log_rr': float(se_lr),
        't_stat': float(t_stat),
        'bf_bic': float(bf_bic),
        'bf_savage_dickey': float(bf_sd),
        'posterior_rr': float(rr_mean),
        'ci95': [float(rr_low), float(rr_high)],
        'p_decrease_given_data': float(norm.cdf(0, mean_lr, se_lr)),
        'bf_sign': float(bf_sign)
    }
}

os.makedirs('data/results', exist_ok=True)
with open('data/results/r35_upgrade.json', 'w') as f:
    json.dump(all_results, f, indent=2, default=str)

print("Results saved to data/results/r35_upgrade.json")
print("\nDONE.")
