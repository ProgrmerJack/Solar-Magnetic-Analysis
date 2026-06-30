"""
Script 56: Blocking-without-SSW comparison
==========================================
THE single most important analysis demanded by ALL 5 reviewers:
Does Alpine blocking alone suppress avalanches, or does the SSW add unique value?

If blocking WITHOUT SSW doesn't suppress avalanches → SSW pathway is vindicated.
If blocking WITHOUT SSW does suppress avalanches → SSW is redundant, just a blocking proxy.

Method:
1. Compute daily Z500 anomaly from NCEP (1998-2019)
2. Define "blocking episodes": Z500 anomaly > 1 SD for ≥5 consecutive days
3. Cross-reference with SSW catalog (±30 day window)
4. Classify: SSW-blocking, non-SSW-blocking, no-blocking
5. Compare avalanche activity (natural dry slab) across classes
"""

import pandas as pd
import numpy as np
from scipy import stats
import json, os

OUT = "data/results/56_blocking_without_ssw.json"
os.makedirs("data/results", exist_ok=True)

# Load data
panel = pd.read_parquet("data/processed/analysis_panel.parquet")
ssw_cat = pd.read_parquet("data/processed/atmospheric/ssw_catalog.parquet")
slf = pd.read_parquet("data/processed/cryosphere/slf_activity.parquet")

# Ensure timezone consistency
panel.index = panel.index.tz_localize(None) if panel.index.tz is None else panel.index.tz_convert("UTC").tz_localize(None)
ssw_cat.index = ssw_cat.index.tz_localize(None) if ssw_cat.index.tz is None else ssw_cat.index.tz_convert("UTC").tz_localize(None)
slf.index = slf.index.tz_localize(None) if slf.index.tz is None else slf.index.tz_convert("UTC").tz_localize(None)

# Get SSW onset dates
ssw_dates = ssw_cat.index.to_list()

# Z500 from panel
z500 = panel["ncep_z500_nh"].dropna()

# Restrict to winter season (Nov-Apr) for climatological relevance
z500_winter = z500[z500.index.month.isin([11, 12, 1, 2, 3, 4])]

# Compute anomalies (subtract monthly climatology)
z500_clim = z500_winter.groupby(z500_winter.index.month).transform("mean")
z500_std = z500_winter.groupby(z500_winter.index.month).transform("std")
z500_anom = (z500_winter - z500_clim) / z500_std  # standardized anomaly

# Define blocking: Z500 anomaly > 1 SD for ≥5 consecutive days
blocking_threshold = 1.0
min_duration = 5

is_blocking_day = (z500_anom > blocking_threshold).astype(int)

# Find blocking episodes (runs of consecutive blocking days)
blocking_episodes = []
in_episode = False
ep_start = None

for i, (date, val) in enumerate(is_blocking_day.items()):
    if val == 1 and not in_episode:
        in_episode = True
        ep_start = date
    elif val == 0 and in_episode:
        in_episode = False
        ep_end = is_blocking_day.index[i-1]
        duration = (ep_end - ep_start).days + 1
        if duration >= min_duration:
            blocking_episodes.append({
                "start": ep_start,
                "end": ep_end,
                "duration": duration,
                "mean_z500_anom": float(z500_anom.loc[ep_start:ep_end].mean())
            })

if in_episode:
    ep_end = is_blocking_day.index[-1]
    duration = (ep_end - ep_start).days + 1
    if duration >= min_duration:
        blocking_episodes.append({
            "start": ep_start,
            "end": ep_end,
            "duration": duration,
            "mean_z500_anom": float(z500_anom.loc[ep_start:ep_end].mean())
        })

print(f"Found {len(blocking_episodes)} blocking episodes (Z500 > 1σ, ≥5 days)")

# Classify: SSW-associated vs non-SSW blocking
ssw_window = 30  # days before/after SSW onset

for ep in blocking_episodes:
    ep["ssw_associated"] = False
    ep["ssw_onset"] = None
    for ssw_date in ssw_dates:
        if abs((ep["start"] - ssw_date).days) <= ssw_window or abs((ep["end"] - ssw_date).days) <= ssw_window:
            ep["ssw_associated"] = True
            ep["ssw_onset"] = str(ssw_date.date())
            break

ssw_blocking = [ep for ep in blocking_episodes if ep["ssw_associated"]]
non_ssw_blocking = [ep for ep in blocking_episodes if not ep["ssw_associated"]]

print(f"SSW-associated blocking: {len(ssw_blocking)}")
print(f"Non-SSW blocking: {len(non_ssw_blocking)}")

# Get avalanche activity for each episode type
# Use aai_dry_natural from SLF
aai_col = "aai_dry_natural"
if aai_col not in slf.columns:
    candidates = [c for c in slf.columns if "dry" in c.lower() and "nat" in c.lower()]
    if candidates:
        aai_col = candidates[0]
    else:
        aai_col = slf.columns[0]
        print(f"WARNING: Using fallback column {aai_col}")

print(f"Using avalanche column: {aai_col}")

def get_aai_for_episodes(episodes, slf_data, col):
    """Get mean daily AAI during episodes and matched climatology."""
    episode_aai = []
    clim_aai = []
    
    for ep in episodes:
        start, end = ep["start"], ep["end"]
        # Episode AAI
        mask = (slf_data.index >= start) & (slf_data.index <= end)
        ep_data = slf_data.loc[mask, col].dropna()
        if len(ep_data) > 0:
            episode_aai.append(ep_data.mean())
            
            # Climatological comparison: same DOY range across all years
            doy_start = start.timetuple().tm_yday
            doy_end = end.timetuple().tm_yday
            if doy_end >= doy_start:
                clim_mask = (slf_data.index.dayofyear >= doy_start) & (slf_data.index.dayofyear <= doy_end)
            else:  # wraps around year
                clim_mask = (slf_data.index.dayofyear >= doy_start) | (slf_data.index.dayofyear <= doy_end)
            clim_data = slf_data.loc[clim_mask, col].dropna()
            clim_aai.append(clim_data.mean())
    
    return np.array(episode_aai), np.array(clim_aai)

# Compute AAI for SSW-blocking episodes
ssw_ep_aai, ssw_clim_aai = get_aai_for_episodes(ssw_blocking, slf, aai_col)
# Compute AAI for non-SSW blocking episodes
nonsw_ep_aai, nonsw_clim_aai = get_aai_for_episodes(non_ssw_blocking, slf, aai_col)

# Also get all-days climatology for reference
winter_mask = slf.index.month.isin([11, 12, 1, 2, 3, 4])
all_winter_aai = slf.loc[winter_mask, aai_col].dropna().mean()

print(f"\n=== RESULTS ===")
print(f"All-winter mean AAI: {all_winter_aai:.3f}")

# SSW-blocking
if len(ssw_ep_aai) > 0:
    ssw_rr = ssw_ep_aai.mean() / ssw_clim_aai.mean() if ssw_clim_aai.mean() > 0 else np.nan
    print(f"\nSSW-blocking episodes (n={len(ssw_ep_aai)}):")
    print(f"  Mean AAI during: {ssw_ep_aai.mean():.3f}")
    print(f"  Climatological AAI: {ssw_clim_aai.mean():.3f}")
    print(f"  Ratio (RR): {ssw_rr:.3f}")
    if len(ssw_ep_aai) >= 3:
        t_ssw, p_ssw = stats.ttest_1samp(ssw_ep_aai - ssw_clim_aai, 0)
        print(f"  Paired t-test vs climatology: t={t_ssw:.3f}, P={p_ssw:.4f}")
    else:
        t_ssw, p_ssw = np.nan, np.nan

# Non-SSW blocking
if len(nonsw_ep_aai) > 0:
    nonsw_rr = nonsw_ep_aai.mean() / nonsw_clim_aai.mean() if nonsw_clim_aai.mean() > 0 else np.nan
    print(f"\nNon-SSW blocking episodes (n={len(nonsw_ep_aai)}):")
    print(f"  Mean AAI during: {nonsw_ep_aai.mean():.3f}")
    print(f"  Climatological AAI: {nonsw_clim_aai.mean():.3f}")
    print(f"  Ratio (RR): {nonsw_rr:.3f}")
    if len(nonsw_ep_aai) >= 3:
        t_non, p_non = stats.ttest_1samp(nonsw_ep_aai - nonsw_clim_aai, 0)
        print(f"  Paired t-test vs climatology: t={t_non:.3f}, P={p_non:.4f}")
    else:
        t_non, p_non = np.nan, np.nan

# Direct comparison: SSW-blocking vs non-SSW blocking
if len(ssw_ep_aai) >= 3 and len(nonsw_ep_aai) >= 3:
    # Compare raw AAI
    t_comp, p_comp = stats.ttest_ind(ssw_ep_aai, nonsw_ep_aai)
    # Compare RRs (episode/climatology ratios)
    ssw_ratios = ssw_ep_aai / ssw_clim_aai
    nonsw_ratios = nonsw_ep_aai / nonsw_clim_aai
    t_rr, p_rr = stats.ttest_ind(ssw_ratios, nonsw_ratios)
    # Mann-Whitney (non-parametric)
    u_stat, p_mw = stats.mannwhitneyu(ssw_ep_aai, nonsw_ep_aai, alternative='less')
    
    print(f"\n=== DIRECT COMPARISON: SSW vs Non-SSW blocking ===")
    print(f"  SSW-blocking mean AAI: {ssw_ep_aai.mean():.3f} (n={len(ssw_ep_aai)})")
    print(f"  Non-SSW blocking mean AAI: {nonsw_ep_aai.mean():.3f} (n={len(nonsw_ep_aai)})")
    print(f"  t-test (raw AAI): t={t_comp:.3f}, P={p_comp:.4f}")
    print(f"  Mann-Whitney (SSW < non-SSW): U={u_stat:.0f}, P={p_mw:.4f}")
    print(f"  Ratio comparison: SSW RR={np.median(ssw_ratios):.3f}, non-SSW RR={np.median(nonsw_ratios):.3f}")
    print(f"  t-test (RR): t={t_rr:.3f}, P={p_rr:.4f}")
    
    # Effect size (Cohen's d)
    pooled_std = np.sqrt(((len(ssw_ep_aai)-1)*ssw_ep_aai.std()**2 + 
                          (len(nonsw_ep_aai)-1)*nonsw_ep_aai.std()**2) / 
                         (len(ssw_ep_aai)+len(nonsw_ep_aai)-2))
    cohens_d = (ssw_ep_aai.mean() - nonsw_ep_aai.mean()) / pooled_std if pooled_std > 0 else np.nan
    print(f"  Cohen's d (SSW vs non-SSW): {cohens_d:.3f}")

# Sign test: what fraction of SSW-blocking episodes show suppression?
if len(ssw_ep_aai) > 0:
    ssw_suppressed = np.sum(ssw_ep_aai < ssw_clim_aai)
    ssw_sign_p = stats.binomtest(ssw_suppressed, len(ssw_ep_aai), 0.5).pvalue
    print(f"\nSSW-blocking sign test: {ssw_suppressed}/{len(ssw_ep_aai)} suppressed, P={ssw_sign_p:.4f}")

if len(nonsw_ep_aai) > 0:
    nonsw_suppressed = np.sum(nonsw_ep_aai < nonsw_clim_aai)
    nonsw_sign_p = stats.binomtest(nonsw_suppressed, len(nonsw_ep_aai), 0.5).pvalue
    print(f"Non-SSW blocking sign test: {nonsw_suppressed}/{len(nonsw_ep_aai)} suppressed, P={nonsw_sign_p:.4f}")

# Summary interpretation
print(f"\n=== INTERPRETATION ===")
if len(ssw_ep_aai) > 0 and len(nonsw_ep_aai) > 0:
    if ssw_rr < nonsw_rr and (p_comp < 0.1 or p_mw < 0.1):
        print("SSW-blocking produces STRONGER suppression than non-SSW blocking.")
        print("→ SSW adds unique value beyond mere blocking.")
        interpretation = "SSW_adds_value"
    elif abs(ssw_rr - nonsw_rr) < 0.15:
        print("SSW-blocking and non-SSW blocking produce SIMILAR suppression.")
        print("→ Blocking is the proximate mechanism; SSW framing may be redundant.")
        interpretation = "blocking_sufficient"
    else:
        print(f"SSW RR={ssw_rr:.3f} vs Non-SSW RR={nonsw_rr:.3f}")
        if nonsw_rr > 0.9:
            print("Non-SSW blocking does NOT suppress avalanches.")
            print("→ SSW-specific atmospheric state is required for suppression.")
            interpretation = "SSW_required"
        else:
            interpretation = "ambiguous"

# Save results
results = {
    "n_blocking_episodes": len(blocking_episodes),
    "n_ssw_blocking": len(ssw_blocking),
    "n_non_ssw_blocking": len(non_ssw_blocking),
    "ssw_blocking_mean_aai": float(ssw_ep_aai.mean()) if len(ssw_ep_aai) > 0 else None,
    "non_ssw_blocking_mean_aai": float(nonsw_ep_aai.mean()) if len(nonsw_ep_aai) > 0 else None,
    "all_winter_mean_aai": float(all_winter_aai),
    "ssw_blocking_RR": float(ssw_rr) if len(ssw_ep_aai) > 0 else None,
    "non_ssw_blocking_RR": float(nonsw_rr) if len(nonsw_ep_aai) > 0 else None,
    "direct_comparison_t": float(t_comp) if len(ssw_ep_aai) >= 3 and len(nonsw_ep_aai) >= 3 else None,
    "direct_comparison_p": float(p_comp) if len(ssw_ep_aai) >= 3 and len(nonsw_ep_aai) >= 3 else None,
    "mann_whitney_p": float(p_mw) if len(ssw_ep_aai) >= 3 and len(nonsw_ep_aai) >= 3 else None,
    "cohens_d": float(cohens_d) if len(ssw_ep_aai) >= 3 and len(nonsw_ep_aai) >= 3 else None,
    "ssw_sign_test": f"{ssw_suppressed}/{len(ssw_ep_aai)}" if len(ssw_ep_aai) > 0 else None,
    "ssw_sign_p": float(ssw_sign_p) if len(ssw_ep_aai) > 0 else None,
    "non_ssw_sign_test": f"{nonsw_suppressed}/{len(nonsw_ep_aai)}" if len(nonsw_ep_aai) > 0 else None,
    "non_ssw_sign_p": float(nonsw_sign_p) if len(nonsw_ep_aai) > 0 else None,
    "interpretation": interpretation,
    "blocking_threshold_sigma": blocking_threshold,
    "min_duration_days": min_duration,
    "ssw_window_days": ssw_window,
    "aai_column": aai_col
}

with open(OUT, "w") as f:
    json.dump(results, f, indent=2, default=str)

print(f"\nResults saved to {OUT}")
