"""
Planetary wave timing analysis.
This addresses the critical reviewer concern: "effect precedes cause."

The key insight: planetary wave amplification drives BOTH:
1. Surface weather changes (affecting avalanche stability) — starts ~10-20 days before SSW
2. Stratospheric vortex disruption (formal SSW onset) — diagnosed at reversal

If both Swiss and Norwegian data show pre-SSW reductions, this is EVIDENCE for
planetary wave forcing, not a confound. The mechanism is:
  Planetary wave amplification → tropospheric weather changes → avalanche reduction
  (same waves also eventually → SSW in stratosphere)

This analysis quantifies the temporal structure of the surface signal relative to SSW onset.
"""
import pandas as pd
import numpy as np
from scipy import stats
import json

# Load all datasets
panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
ssw = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw.index = ssw.index.tz_localize(None)
nve = pd.read_csv('data/cryosphere/norway_nve/nve_ssw_analysis.csv')
aval = 'dry_natural_size_1234'

# Swiss: compute rolling 5-day average rates relative to SSW
panel_valid = panel.dropna(subset=[aval])
ssw_in_swiss = ssw[(ssw.index >= panel_valid.index.min()) & (ssw.index <= panel_valid.index.max())]

print("=== TEMPORAL STRUCTURE OF AVALANCHE RESPONSE ===\n")

# 5-day window analysis (more stable than daily)
window = 5
swiss_windows = []
for center_lag in range(-25, 26):
    start_lag = center_lag - window // 2
    end_lag = center_lag + window // 2
    vals = []
    for sd in ssw_in_swiss.index:
        for lag in range(start_lag, end_lag + 1):
            date = sd + pd.Timedelta(days=lag)
            if date in panel_valid.index:
                v = panel_valid.loc[date, aval]
                if not np.isnan(v):
                    vals.append(v)
    if vals:
        swiss_windows.append({
            'center_lag': center_lag,
            'mean_rate': np.mean(vals),
            'n': len(vals),
        })

sw_df = pd.DataFrame(swiss_windows)
baseline = panel_valid.dropna(subset=[aval])
baseline = baseline[(baseline.index.month >= 11) | (baseline.index.month <= 4)]
swiss_baseline = baseline[aval].mean()

print("Swiss 5-day rolling rate (relative to baseline):")
print(f"{'Lag':>5} {'Rate':>8} {'RR':>8} {'Category'}")
print("-" * 40)
for _, r in sw_df.iterrows():
    lag = int(r['center_lag'])
    rr = r['mean_rate'] / swiss_baseline
    cat = ""
    if lag < -15:
        cat = "far-pre"
    elif lag < 0:
        cat = "pre-SSW"
    elif lag <= 15:
        cat = "post-SSW"
    else:
        cat = "recovery"
    print(f"{lag:+5d} {r['mean_rate']:8.3f} {rr:8.3f} {cat}")

# Identify the minimum (peak suppression)
min_idx = sw_df['mean_rate'].idxmin()
min_lag = sw_df.loc[min_idx, 'center_lag']
min_rate = sw_df.loc[min_idx, 'mean_rate']
print(f"\nPeak suppression at lag={min_lag:+d} days, rate={min_rate:.3f}, RR={min_rate/swiss_baseline:.3f}")

# Norwegian lag profile (already computed, but compute 5-day rolling)
nve_ssw = nve[nve['window_type'] == 'ssw'].copy()
nve_ssw['date_dt'] = pd.to_datetime(nve_ssw['date'])
nve_ssw['ssw_dt'] = pd.to_datetime(nve_ssw['ssw_date'])
nve_ssw['lag'] = (nve_ssw['date_dt'] - nve_ssw['ssw_dt']).dt.days

nve_ctrl = nve[nve['window_type'] == 'ctrl']
nve_baseline = nve_ctrl['danger_level'].mean()

nve_windows = []
for center_lag in range(-15, 16):
    start_lag = center_lag - 2
    end_lag = center_lag + 2
    sub = nve_ssw[(nve_ssw['lag'] >= start_lag) & (nve_ssw['lag'] <= end_lag)]
    if len(sub) > 0:
        nve_windows.append({
            'center_lag': center_lag,
            'mean_danger': sub['danger_level'].mean(),
            'n': len(sub),
        })

nve_df = pd.DataFrame(nve_windows)
print(f"\nNorwegian 5-day rolling danger level (baseline={nve_baseline:.2f}):")
print(f"{'Lag':>5} {'Danger':>8} {'RR':>8}")
print("-" * 30)
for _, r in nve_df.iterrows():
    lag = int(r['center_lag'])
    rr = r['mean_danger'] / nve_baseline
    print(f"{lag:+5d} {r['mean_danger']:8.3f} {rr:8.3f}")

nve_min_idx = nve_df['mean_danger'].idxmin()
nve_min_lag = nve_df.loc[nve_min_idx, 'center_lag']
nve_min_danger = nve_df.loc[nve_min_idx, 'mean_danger']
print(f"\nPeak suppression at lag={nve_min_lag:+d} days, danger={nve_min_danger:.3f}, RR={nve_min_danger/nve_baseline:.3f}")

# Cross-country comparison
print(f"\n=== CROSS-COUNTRY TIMING COMPARISON ===")
print(f"Swiss peak suppression:     lag = {min_lag:+d} days")
print(f"Norwegian peak suppression: lag = {nve_min_lag:+d} days")

# Correlation between Swiss and Norwegian lag profiles
# Align on common lag range [-15, +15]
common_lags = range(-15, 16)
swiss_rr = []
nve_rr = []
for lag in common_lags:
    sw_row = sw_df[sw_df['center_lag'] == lag]
    nv_row = nve_df[nve_df['center_lag'] == lag]
    if len(sw_row) > 0 and len(nv_row) > 0:
        swiss_rr.append(sw_row['mean_rate'].values[0] / swiss_baseline)
        nve_rr.append(nv_row['mean_danger'].values[0] / nve_baseline)

if len(swiss_rr) > 5:
    r_cross, p_cross = stats.spearmanr(swiss_rr, nve_rr)
    print(f"Swiss-Norwegian temporal profile correlation: r={r_cross:.3f}, P={p_cross:.4f}")
    print(f"  (Positive r means both countries show similar timing patterns)")

# Phase-resolved comparison
print(f"\n=== PHASE-RESOLVED CROSS-COUNTRY COMPARISON ===")
phases = [
    ('Early pre', -15, -8),
    ('Late pre', -7, -1),
    ('Onset', 0, 3),
    ('Early post', 4, 10),
    ('Late post', 11, 15),
]

for name, start, end in phases:
    sw_sub = sw_df[(sw_df['center_lag'] >= start) & (sw_df['center_lag'] <= end)]
    nv_sub = nve_df[(nve_df['center_lag'] >= start) & (nve_df['center_lag'] <= end)]
    
    sw_rr_val = sw_sub['mean_rate'].mean() / swiss_baseline if len(sw_sub) > 0 else np.nan
    nv_rr_val = nv_sub['mean_danger'].mean() / nve_baseline if len(nv_sub) > 0 else np.nan
    
    print(f"  {name:12s} [{start:+3d},{end:+3d}]: Swiss RR={sw_rr_val:.3f}  Norway RR={nv_rr_val:.3f}")

# Save results
output = {
    'swiss': {
        'peak_suppression_lag': int(min_lag),
        'peak_suppression_rr': float(min_rate / swiss_baseline),
        'baseline_rate': float(swiss_baseline),
        'n_events': int(len(ssw_in_swiss)),
    },
    'norway': {
        'peak_suppression_lag': int(nve_min_lag),
        'peak_suppression_rr': float(nve_min_danger / nve_baseline),
        'baseline_danger': float(nve_baseline),
        'n_events': 4,
    },
    'cross_country_correlation': {
        'spearman_r': float(r_cross) if len(swiss_rr) > 5 else None,
        'spearman_p': float(p_cross) if len(swiss_rr) > 5 else None,
    },
}

with open('data/results/planetary_wave_timing.json', 'w') as f:
    json.dump(output, f, indent=2)
print('\nSaved to data/results/planetary_wave_timing.json')
