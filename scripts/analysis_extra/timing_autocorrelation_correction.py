"""
Correct the cross-country temporal correlation for autocorrelation.
R10-A flagged that 5-day rolling windows create ~6-8 effective df, not 29.
Implements block bootstrap and effective-n corrections.

Regenerates profiles from source data, then applies corrections.
"""
import json
import numpy as np
import pandas as pd
from scipy import stats
import os

BASE = os.path.join(os.path.dirname(__file__), '..')

# === Regenerate Swiss profile ===
panel = pd.read_parquet(os.path.join(BASE, 'data', 'processed', 'analysis_panel_v2.parquet'))
ssw_cat = pd.read_parquet(os.path.join(BASE, 'data', 'processed', 'atmospheric', 'ssw_catalog.parquet'))

panel.index = pd.to_datetime(panel.index).tz_localize(None)
ssw_dates = pd.to_datetime(ssw_cat.index).tz_localize(None)

# Swiss dry slab rate
aval_col = [c for c in panel.columns if 'dry' in c.lower() and 'nat' in c.lower()][0]
panel_clean = panel[[aval_col]].dropna()
baseline = panel_clean[aval_col].mean()

lag_range = range(-15, 16)
swiss_rr = {}
for lag in lag_range:
    vals = []
    for ssw_d in ssw_dates:
        target = ssw_d + pd.Timedelta(days=lag)
        window = pd.date_range(target - pd.Timedelta(days=2), target + pd.Timedelta(days=2))
        hits = panel_clean.reindex(window).dropna()
        if len(hits) > 0:
            vals.append(hits[aval_col].mean())
    if vals:
        swiss_rr[lag] = np.mean(vals) / baseline

# === Regenerate Norway profile ===
nve = pd.read_csv(os.path.join(BASE, 'data', 'cryosphere', 'norway_nve', 'nve_ssw_analysis.csv'))
nve_ssw = nve[nve['ssw_date'] != 'control'].copy()
nve_ssw['ssw_date'] = pd.to_datetime(nve_ssw['ssw_date'])
nve_ssw['date'] = pd.to_datetime(nve_ssw['date'])
nve_control = nve[nve['ssw_date'] == 'control']
norway_baseline = nve_control['danger_level'].mean()

norway_rr = {}
for lag in lag_range:
    vals = []
    for ssw_d in nve_ssw['ssw_date'].unique():
        target = ssw_d + pd.Timedelta(days=lag)
        window_start = target - pd.Timedelta(days=2)
        window_end = target + pd.Timedelta(days=2)
        hits = nve_ssw[(nve_ssw['ssw_date'] == ssw_d) & 
                       (nve_ssw['date'] >= window_start) & 
                       (nve_ssw['date'] <= window_end)]
        if len(hits) > 0:
            vals.append(hits['danger_level'].mean())
    if vals:
        norway_rr[lag] = np.mean(vals) / norway_baseline

# Build aligned arrays
common_lags = sorted(set(swiss_rr.keys()) & set(norway_rr.keys()))
swiss_vals = np.array([swiss_rr[l] for l in common_lags])
norway_vals = np.array([norway_rr[l] for l in common_lags])
lags_arr = np.array(common_lags)

print(f"Common lags: {len(common_lags)} (range {min(lags_arr)} to {max(lags_arr)})")

# Raw Spearman correlation
r_raw, p_raw = stats.spearmanr(swiss_vals, norway_vals)
print(f"\nRaw Spearman: r={r_raw:.3f}, P={p_raw:.6f}")

# === Method 1: Effective-N correction (Bretherton et al. 1999) ===
# Compute lag-1 autocorrelation for both series
def lag1_autocorr(x):
    n = len(x)
    xm = x - np.mean(x)
    return np.sum(xm[:-1] * xm[1:]) / np.sum(xm**2)

r1_swiss = lag1_autocorr(swiss_vals)
r1_norway = lag1_autocorr(norway_vals)
print(f"\nLag-1 autocorrelation: Swiss={r1_swiss:.3f}, Norway={r1_norway:.3f}")

# Effective degrees of freedom (Bretherton et al. 1999)
n = len(swiss_vals)
n_eff = n * (1 - r1_swiss * r1_norway) / (1 + r1_swiss * r1_norway)
n_eff = max(n_eff, 3)  # floor
print(f"Effective N: {n_eff:.1f} (raw N={n})")

# Recompute P-value with effective df
t_stat = r_raw * np.sqrt((n_eff - 2) / (1 - r_raw**2))
p_corrected = 2 * stats.t.sf(abs(t_stat), df=n_eff - 2)
print(f"Corrected P (Bretherton): t={t_stat:.3f}, df={n_eff-2:.1f}, P={p_corrected:.6f}")

# === Method 2: Block bootstrap ===
np.random.seed(42)
n_boot = 10000
block_sizes = [3, 5, 7]

for block_size in block_sizes:
    n_blocks = n // block_size
    boot_correlations = []
    
    for _ in range(n_boot):
        # Sample block starting positions with replacement
        starts = np.random.randint(0, n - block_size + 1, size=n_blocks)
        indices = np.concatenate([np.arange(s, s + block_size) for s in starts])[:n]
        
        if len(indices) >= 5:
            s_boot = swiss_vals[indices]
            n_boot_vals = norway_vals[indices]
            r_boot, _ = stats.spearmanr(s_boot, n_boot_vals)
            if not np.isnan(r_boot):
                boot_correlations.append(r_boot)
    
    boot_arr = np.array(boot_correlations)
    # P-value: proportion of bootstrap samples with r <= 0
    p_boot = np.mean(boot_arr <= 0) * 2  # two-sided
    ci_lo = np.percentile(boot_arr, 2.5)
    ci_hi = np.percentile(boot_arr, 97.5)
    
    print(f"\nBlock bootstrap (block={block_size}, n_boot={len(boot_correlations)}):")
    print(f"  Mean r={np.mean(boot_arr):.3f}, median={np.median(boot_arr):.3f}")
    print(f"  95% CI: [{ci_lo:.3f}, {ci_hi:.3f}]")
    print(f"  P(r<=0, two-sided): {p_boot:.4f}")

# === Method 3: Non-overlapping windows only ===
# Use every 5th lag to eliminate overlap completely
step = 5
independent_swiss = swiss_vals[::step]
independent_norway = norway_vals[::step]
independent_lags = lags_arr[::step]
print(f"\nNon-overlapping (every 5th lag): n={len(independent_swiss)}")
print(f"  Lags used: {independent_lags}")
r_indep, p_indep = stats.spearmanr(independent_swiss, independent_norway)
print(f"  Spearman r={r_indep:.3f}, P={p_indep:.4f}")

# === Method 4: Phase-level correlation (fully independent) ===
# 5 phases: early pre, late pre, onset, early post, late post
swiss_phases = [0.69, 0.77, 0.85, 0.93, 1.10]  # from analysis
norway_phases = [0.87, 0.88, 0.88, 0.92, 0.98]
r_phase, p_phase = stats.spearmanr(swiss_phases, norway_phases)
print(f"\nPhase-level correlation (n=5 independent phases):")
print(f"  Swiss phases: {swiss_phases}")
print(f"  Norway phases: {norway_phases}")
print(f"  Spearman r={r_phase:.3f}, P={p_phase:.4f}")

# === Summary ===
print("\n" + "="*60)
print("SUMMARY: Cross-country temporal correlation corrections")
print("="*60)
print(f"Raw:                    r={r_raw:.3f}, P={p_raw:.6f} (N={n})")
print(f"Bretherton corrected:   r={r_raw:.3f}, P={p_corrected:.6f} (N_eff={n_eff:.1f})")
print(f"Block bootstrap (b=5):  95% CI reported above")
print(f"Non-overlapping (Δ=5):  r={r_indep:.3f}, P={p_indep:.4f} (N={len(independent_swiss)})")
print(f"Phase-level:            r={r_phase:.3f}, P={p_phase:.4f} (N=5)")

# Save corrected results
results = {
    'raw': {'r': float(r_raw), 'p': float(p_raw), 'n': int(n)},
    'bretherton_corrected': {
        'r': float(r_raw), 'p': float(p_corrected),
        'n_eff': float(n_eff),
        'r1_swiss': float(r1_swiss), 'r1_norway': float(r1_norway)
    },
    'non_overlapping': {'r': float(r_indep), 'p': float(p_indep), 'n': int(len(independent_swiss))},
    'phase_level': {'r': float(r_phase), 'p': float(p_phase), 'n': 5}
}

out_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'results', 'timing_autocorrelation_correction.json')
with open(out_path, 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nResults saved to {out_path}")
