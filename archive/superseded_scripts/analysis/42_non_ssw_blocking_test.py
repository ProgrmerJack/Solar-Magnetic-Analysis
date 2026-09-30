"""
Script 42: Non-SSW Blocking Equivalence Test

The decisive mechanistic test (identified by stratospheric dynamics reviewer):
Does equivalent Z500 blocking magnitude OUTSIDE of SSW windows produce 
equivalent avalanche suppression?

If yes → blocking is the mechanism, SSW is just a marker
If no → SSW adds something beyond blocking alone
"""
import pandas as pd
import numpy as np
from scipy import stats
import json
from pathlib import Path

DATA = Path("data/processed")
OUT = Path("data/results")

# Load data
activity = pd.read_parquet(DATA / "cryosphere/slf_activity.parquet")
ssw_cat = pd.read_parquet(DATA / "atmospheric/ssw_catalog.parquet")
ncep_trop = pd.read_parquet(DATA / "atmospheric/ncep_troposphere.parquet")

# Get Z500 (Alpine sector mean from NCEP)
z500 = ncep_trop['hgt_500hPa_m'].copy()
z500.index = z500.index.tz_localize(None)

# Get activity
act = activity['aai_dry_natural'].copy()
if hasattr(act.index, 'tz'):
    act.index = act.index.tz_localize(None)

# Align dates
common_idx = z500.index.intersection(act.index)
z500 = z500.loc[common_idx]
act = act.loc[common_idx]

# Compute DOY-based Z500 anomalies (deviation from DOY mean)
z500_df = pd.DataFrame({'z500': z500, 'act': act, 'doy': z500.index.dayofyear})
# Only winter days (Nov-Apr)
winter_mask = z500_df.index.month.isin([11, 12, 1, 2, 3, 4])
z500_winter = z500_df[winter_mask].copy()

doy_mean = z500_winter.groupby('doy')['z500'].mean()
z500_winter['z500_anom'] = z500_winter.apply(lambda r: r['z500'] - doy_mean.get(r['doy'], r['z500']), axis=1)

# Mark SSW windows (±15 days from onset)
ssw_dates = set()
for onset in ssw_cat.index:
    onset_dt = pd.Timestamp(onset).tz_localize(None)
    for d in range(-15, 16):
        ssw_dates.add(onset_dt + pd.Timedelta(days=d))

z500_winter['is_ssw'] = z500_winter.index.isin(ssw_dates)

# Split into SSW and non-SSW days
ssw_days = z500_winter[z500_winter['is_ssw']]
non_ssw_days = z500_winter[~z500_winter['is_ssw']]

print(f"Total winter days: {len(z500_winter)}")
print(f"SSW window days: {len(ssw_days)}")
print(f"Non-SSW days: {len(non_ssw_days)}")
print(f"\nZ500 anomaly during SSW: {ssw_days['z500_anom'].mean():.1f} m (median: {ssw_days['z500_anom'].median():.1f})")
print(f"Z500 anomaly non-SSW: {non_ssw_days['z500_anom'].mean():.1f} m (median: {non_ssw_days['z500_anom'].median():.1f})")

# Key test: Find non-SSW days with EQUIVALENT Z500 blocking
# SSW days tend to have HIGH Z500 anomaly (blocking)
ssw_z500_thresh = ssw_days['z500_anom'].quantile(0.25)  # 25th percentile of SSW Z500
print(f"\nSSW Z500 anomaly 25th percentile: {ssw_z500_thresh:.1f} m")

# Non-SSW blocking days: Z500 anomaly >= SSW 25th percentile
non_ssw_blocking = non_ssw_days[non_ssw_days['z500_anom'] >= ssw_z500_thresh]
non_ssw_normal = non_ssw_days[non_ssw_days['z500_anom'] < ssw_z500_thresh]

print(f"Non-SSW blocking days (Z500 >= {ssw_z500_thresh:.0f}m): {len(non_ssw_blocking)}")
print(f"Non-SSW normal days: {len(non_ssw_normal)}")

# Compare avalanche rates
ssw_rate = ssw_days['act'].mean()
non_ssw_block_rate = non_ssw_blocking['act'].mean()
non_ssw_normal_rate = non_ssw_normal['act'].mean()
all_non_ssw_rate = non_ssw_days['act'].mean()

print(f"\n=== AVALANCHE RATES ===")
print(f"SSW window days: {ssw_rate:.2f} aval/day")
print(f"Non-SSW blocking days (equivalent Z500): {non_ssw_block_rate:.2f} aval/day")  
print(f"Non-SSW normal days: {non_ssw_normal_rate:.2f} aval/day")
print(f"All non-SSW days: {all_non_ssw_rate:.2f} aval/day")

# RR: SSW vs all non-SSW
rr_ssw = ssw_rate / all_non_ssw_rate if all_non_ssw_rate > 0 else np.nan
# RR: non-SSW blocking vs non-SSW normal
rr_blocking = non_ssw_block_rate / non_ssw_normal_rate if non_ssw_normal_rate > 0 else np.nan

print(f"\nRR (SSW vs non-SSW): {rr_ssw:.3f}")
print(f"RR (non-SSW blocking vs non-SSW normal): {rr_blocking:.3f}")

# Statistical test
u_stat, u_p = stats.mannwhitneyu(ssw_days['act'], non_ssw_blocking['act'], alternative='two-sided')
print(f"\nMann-Whitney SSW vs non-SSW-blocking: U={u_stat:.0f}, P={u_p:.4f}")

# THE DECISIVE COMPARISON:
# If SSW days have LOWER avalanche rate than non-SSW-blocking days with equivalent Z500,
# then SSW adds suppression beyond Z500 alone.
# If rates are similar, blocking fully explains the effect.
print(f"\n=== DECISIVE TEST ===")
print(f"SSW avalanche rate: {ssw_rate:.2f}")
print(f"Non-SSW equivalent-blocking rate: {non_ssw_block_rate:.2f}")
if ssw_rate < non_ssw_block_rate:
    excess = (1 - ssw_rate/non_ssw_block_rate) * 100
    print(f"SSW days have {excess:.1f}% LOWER rate than equivalent-blocking non-SSW days")
    print(f"→ SSW adds suppression BEYOND blocking alone")
else:
    excess = (ssw_rate/non_ssw_block_rate - 1) * 100
    print(f"SSW days have {excess:.1f}% HIGHER rate than equivalent-blocking non-SSW days")
    print(f"→ Blocking alone explains the effect; SSW adds nothing")

# Quartile analysis: bin Z500 anomaly and compare SSW vs non-SSW within bins
print(f"\n=== QUARTILE ANALYSIS ===")
z500_quartiles = z500_winter['z500_anom'].quantile([0.25, 0.5, 0.75])
bins = [-np.inf, z500_quartiles[0.25], z500_quartiles[0.5], z500_quartiles[0.75], np.inf]
labels = ['Q1 (low Z500)', 'Q2', 'Q3', 'Q4 (high Z500/blocking)']
z500_winter['z500_q'] = pd.cut(z500_winter['z500_anom'], bins=bins, labels=labels)

results = {}
for q in labels:
    q_data = z500_winter[z500_winter['z500_q'] == q]
    q_ssw = q_data[q_data['is_ssw']]
    q_non = q_data[~q_data['is_ssw']]
    
    if len(q_ssw) >= 5 and len(q_non) >= 5:
        rr = q_ssw['act'].mean() / q_non['act'].mean() if q_non['act'].mean() > 0 else np.nan
        u, p = stats.mannwhitneyu(q_ssw['act'], q_non['act'], alternative='two-sided')
        print(f"{q}: SSW={q_ssw['act'].mean():.2f} (n={len(q_ssw)}), non-SSW={q_non['act'].mean():.2f} (n={len(q_non)}), RR={rr:.3f}, P={p:.4f}")
        results[q] = {'ssw_rate': float(q_ssw['act'].mean()), 'non_ssw_rate': float(q_non['act'].mean()),
                       'rr': float(rr), 'p': float(p), 'n_ssw': int(len(q_ssw)), 'n_non': int(len(q_non))}
    else:
        print(f"{q}: SSW n={len(q_ssw)}, non-SSW n={len(q_non)} — insufficient")

# Save results
output = {
    'ssw_rate': float(ssw_rate),
    'non_ssw_blocking_rate': float(non_ssw_block_rate),
    'non_ssw_normal_rate': float(non_ssw_normal_rate),
    'rr_ssw_vs_all': float(rr_ssw),
    'rr_blocking_vs_normal': float(rr_blocking),
    'mann_whitney_ssw_vs_blocking': {'U': float(u_stat), 'P': float(u_p)},
    'z500_threshold': float(ssw_z500_thresh),
    'n_ssw_days': int(len(ssw_days)),
    'n_blocking_days': int(len(non_ssw_blocking)),
    'quartile_analysis': results
}

OUT.mkdir(parents=True, exist_ok=True)
with open(OUT / '42_non_ssw_blocking_test.json', 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nSaved to {OUT / '42_non_ssw_blocking_test.json'}")

