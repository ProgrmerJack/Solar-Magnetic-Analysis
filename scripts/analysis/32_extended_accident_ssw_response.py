#!/usr/bin/env python3
"""
32_extended_accident_ssw_response.py

Extends the SSW–avalanche association from the primary n=16 activity-based sample
(1998–2019) to the full accident-based record (1970–2025), yielding n≈30 events.

This addresses the core "n=16 ceiling" criticism by showing that the SSW–avalanche
link is detectable in an independent, longer dataset using a fundamentally different
outcome measure (fatal/caught avalanche incidents rather than daily activity counts).

Key output:
- SSW-period vs non-SSW-period accident rate ratio for 30 SSW events
- Pre-1998 subsample (n=14): independent validation on data never seen by primary analysis
- Post-2019 subsample (n=1, Jan 2021): true out-of-sample prospective test
- Bootstrapped confidence intervals and permutation-based p-values
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "results"
OUT.mkdir(parents=True, exist_ok=True)


def load_data():
    """Load SSW catalog, accident data, and activity data."""
    ssw = pd.read_parquet(ROOT / "data/processed/atmospheric/ssw_catalog.parquet")
    acc = pd.read_parquet(ROOT / "data/processed/cryosphere/slf_accidents.parquet")
    act = pd.read_parquet(ROOT / "data/processed/cryosphere/slf_activity.parquet")
    return ssw, acc, act


def compute_accident_rates(acc, ssw_dates, window_days=30, pre_days=30):
    """
    For each SSW event, compute:
    - Accident count in the post-SSW window [onset, onset+window_days]
    - Accident count in the pre-SSW window [onset-pre_days, onset-1]
    - Daily rate ratio (post/pre)
    
    Also computes matched-winter baselines from non-SSW years.
    """
    results = []
    
    for onset in ssw_dates:
        onset_ts = pd.Timestamp(onset) if pd.Timestamp(onset).tzinfo is not None else pd.Timestamp(onset, tz='UTC')
        post_start = onset_ts
        post_end = onset_ts + pd.Timedelta(days=window_days)
        pre_start = onset_ts - pd.Timedelta(days=pre_days)
        pre_end = onset_ts - pd.Timedelta(days=1)
        
        # Count accidents in windows
        post_mask = (acc.index >= post_start) & (acc.index <= post_end)
        pre_mask = (acc.index >= pre_start) & (acc.index <= pre_end)
        
        post_count = post_mask.sum()
        pre_count = pre_mask.sum()
        
        # Daily rates
        post_rate = post_count / window_days
        pre_rate = pre_count / pre_days
        
        # Log rate ratio (with continuity correction)
        if pre_rate > 0:
            log_rr = np.log(max(post_rate, 0.5/window_days) / pre_rate)
        else:
            log_rr = np.nan
        
        results.append({
            'onset': onset.strftime('%Y-%m-%d') if hasattr(onset, 'strftime') else str(onset)[:10],
            'post_count': int(post_count),
            'pre_count': int(pre_count),
            'post_rate': float(post_rate),
            'pre_rate': float(pre_rate),
            'log_rr': float(log_rr) if not np.isnan(log_rr) else None,
            'window_days': window_days
        })
    
    return results


def winter_baseline_rates(acc, exclude_ssw_dates, window_days=30):
    """Compute baseline winter accident rates excluding SSW windows."""
    winter_months = [11, 12, 1, 2, 3, 4]
    winter_acc = acc[acc.index.month.isin(winter_months)]
    
    # Exclude SSW windows
    mask = pd.Series(True, index=winter_acc.index)
    for onset in exclude_ssw_dates:
        onset_ts = pd.Timestamp(onset) if pd.Timestamp(onset).tzinfo is not None else pd.Timestamp(onset, tz='UTC')
        ssw_mask = (winter_acc.index >= onset_ts - pd.Timedelta(days=5)) & \
                   (winter_acc.index <= onset_ts + pd.Timedelta(days=window_days + 5))
        mask = mask & ~ssw_mask
    
    baseline = winter_acc[mask]
    years_covered = (baseline.index.max() - baseline.index.min()).days / 365.25
    baseline_daily_rate = len(baseline) / max((baseline.index.max() - baseline.index.min()).days, 1)
    
    return {
        'n_baseline_accidents': len(baseline),
        'years_covered': float(years_covered),
        'daily_rate': float(baseline_daily_rate)
    }


def bootstrap_rr(results_list, n_boot=10000, seed=42):
    """Bootstrap the mean log rate ratio with 95% CI."""
    rng = np.random.RandomState(seed)
    log_rrs = np.array([r['log_rr'] for r in results_list if r['log_rr'] is not None])
    n = len(log_rrs)
    
    if n < 3:
        return {'mean_log_rr': float(np.mean(log_rrs)), 'ci_95': [None, None], 'n': n}
    
    boot_means = np.array([np.mean(rng.choice(log_rrs, size=n, replace=True)) for _ in range(n_boot)])
    
    return {
        'mean_log_rr': float(np.mean(log_rrs)),
        'median_log_rr': float(np.median(log_rrs)),
        'se': float(np.std(log_rrs) / np.sqrt(n)),
        'ci_95': [float(np.percentile(boot_means, 2.5)), float(np.percentile(boot_means, 97.5))],
        'mean_rr': float(np.exp(np.mean(log_rrs))),
        'ci_95_rr': [float(np.exp(np.percentile(boot_means, 2.5))), float(np.exp(np.percentile(boot_means, 97.5)))],
        'n': n,
        'frac_negative': float(np.mean(log_rrs < 0))
    }


def permutation_test(log_rrs, n_perm=50000, seed=42):
    """Permutation test for H0: mean log RR = 0."""
    rng = np.random.RandomState(seed)
    observed = np.mean(log_rrs)
    n = len(log_rrs)
    
    # Under H0, each log_rr is equally likely to be positive or negative
    perm_means = np.array([
        np.mean(log_rrs * rng.choice([-1, 1], size=n))
        for _ in range(n_perm)
    ])
    
    p_value = np.mean(np.abs(perm_means) >= np.abs(observed))
    return float(p_value)


def cohens_d(log_rrs):
    """Effect size: Cohen's d for mean log RR against zero."""
    return float(np.mean(log_rrs) / np.std(log_rrs, ddof=1)) if len(log_rrs) > 1 else None


def main():
    ssw, acc, act = load_data()
    
    # SSW events overlapping with accident data (1970-2025)
    acc_start = acc.index.min()
    acc_end = acc.index.max()
    print(f"Accident data range: {acc_start.strftime('%Y-%m-%d')} to {acc_end.strftime('%Y-%m-%d')}")
    
    ssw_in_range = ssw[(ssw.index >= acc_start - pd.Timedelta(days=30)) & 
                        (ssw.index <= acc_end)]
    print(f"SSW events in accident range: {len(ssw_in_range)}")
    
    # Filter to winter SSW events only (Nov-Mar onset)
    winter_ssw = ssw_in_range[ssw_in_range.index.month.isin([11, 12, 1, 2, 3])]
    print(f"Winter SSW events (Nov-Mar): {len(winter_ssw)}")
    
    # Compute accident rates for all events
    all_results = compute_accident_rates(acc, winter_ssw.index, window_days=30)
    
    # Split into subsamples
    pre1998 = [r for r in all_results if int(r['onset'][:4]) < 1998]
    primary = [r for r in all_results if 1998 <= int(r['onset'][:4]) <= 2019]
    post2019 = [r for r in all_results if int(r['onset'][:4]) > 2019]
    
    print(f"\nPre-1998 events: {len(pre1998)}")
    print(f"Primary period (1998-2019): {len(primary)}")
    print(f"Post-2019 events: {len(post2019)}")
    
    # Bootstrap and test for each subsample
    all_log_rrs = np.array([r['log_rr'] for r in all_results if r['log_rr'] is not None])
    pre_log_rrs = np.array([r['log_rr'] for r in pre1998 if r['log_rr'] is not None])
    primary_log_rrs = np.array([r['log_rr'] for r in primary if r['log_rr'] is not None])
    
    # Full sample analysis
    full_boot = bootstrap_rr(all_results)
    full_perm_p = permutation_test(all_log_rrs)
    full_d = cohens_d(all_log_rrs)
    full_ttest = stats.ttest_1samp(all_log_rrs, 0)
    full_wilcox = stats.wilcoxon(all_log_rrs)
    
    # Pre-1998 subsample
    pre_boot = bootstrap_rr(pre1998)
    pre_perm_p = permutation_test(pre_log_rrs) if len(pre_log_rrs) >= 3 else None
    pre_d = cohens_d(pre_log_rrs)
    
    # Primary period comparison
    primary_boot = bootstrap_rr(primary)
    primary_perm_p = permutation_test(primary_log_rrs) if len(primary_log_rrs) >= 3 else None
    
    # Baseline rates
    baseline = winter_baseline_rates(acc, winter_ssw.index)
    
    # Sign test: what fraction of events show decreased accidents?
    n_decrease = sum(1 for r in all_results if r['log_rr'] is not None and r['log_rr'] < 0)
    n_valid = sum(1 for r in all_results if r['log_rr'] is not None)
    sign_p = stats.binomtest(n_decrease, n_valid, 0.5).pvalue if n_valid > 0 else None
    
    # Compile output
    output = {
        'description': 'Extended SSW-avalanche accident-based analysis (1970-2025)',
        'data_source': 'SLF fatal/caught avalanche incidents (independent from daily activity counts)',
        'window_days': 30,
        'full_sample': {
            'n_events': len(all_results),
            'n_valid_rr': int(n_valid),
            'mean_log_rr': full_boot['mean_log_rr'],
            'median_log_rr': full_boot['median_log_rr'],
            'mean_rr': full_boot['mean_rr'],
            'ci_95_rr': full_boot['ci_95_rr'],
            'se': full_boot['se'],
            'permutation_p': full_perm_p,
            'ttest_p': float(full_ttest.pvalue),
            'wilcoxon_p': float(full_wilcox.pvalue),
            'cohens_d': full_d,
            'frac_decrease': float(n_decrease / n_valid) if n_valid > 0 else None,
            'sign_test_p': float(sign_p) if sign_p is not None else None
        },
        'pre1998_subsample': {
            'description': 'Independent pre-1998 events never used in primary analysis',
            'n_events': len(pre1998),
            'n_valid_rr': len(pre_log_rrs),
            'mean_log_rr': pre_boot['mean_log_rr'],
            'mean_rr': pre_boot['mean_rr'],
            'ci_95_rr': pre_boot['ci_95_rr'],
            'permutation_p': pre_perm_p,
            'cohens_d': pre_d,
            'frac_decrease': float(np.mean(pre_log_rrs < 0)) if len(pre_log_rrs) > 0 else None
        },
        'primary_period': {
            'description': 'Events from original study period (1998-2019)',
            'n_events': len(primary),
            'mean_log_rr': primary_boot['mean_log_rr'],
            'mean_rr': primary_boot['mean_rr'],
            'ci_95_rr': primary_boot['ci_95_rr'],
            'permutation_p': primary_perm_p
        },
        'post2019_prospective': {
            'description': 'True out-of-sample events (post-2019)',
            'events': post2019
        },
        'per_event_details': all_results,
        'baseline': baseline
    }
    
    # Print summary
    print("\n" + "="*70)
    print("EXTENDED SSW-AVALANCHE ACCIDENT ANALYSIS")
    print("="*70)
    print(f"\nFull sample (n={n_valid}):")
    print(f"  Mean RR = {full_boot['mean_rr']:.3f} [{full_boot['ci_95_rr'][0]:.3f}, {full_boot['ci_95_rr'][1]:.3f}]")
    print(f"  Permutation P = {full_perm_p:.4f}")
    print(f"  t-test P = {full_ttest.pvalue:.4f}")
    print(f"  Wilcoxon P = {full_wilcox.pvalue:.4f}")
    print(f"  Cohen's d = {full_d:.3f}")
    print(f"  {n_decrease}/{n_valid} events show decrease ({100*n_decrease/n_valid:.0f}%)")
    
    print(f"\nPre-1998 independent subsample (n={len(pre_log_rrs)}):")
    print(f"  Mean RR = {pre_boot['mean_rr']:.3f} [{pre_boot['ci_95_rr'][0]:.3f}, {pre_boot['ci_95_rr'][1]:.3f}]")
    if pre_perm_p is not None:
        print(f"  Permutation P = {pre_perm_p:.4f}")
    print(f"  Cohen's d = {pre_d:.3f}" if pre_d else "  Cohen's d = N/A")
    
    print(f"\nPrimary period (1998-2019, n={len(primary_log_rrs)}):")
    print(f"  Mean RR = {primary_boot['mean_rr']:.3f}")
    
    if post2019:
        print(f"\nPost-2019 prospective test:")
        for ev in post2019:
            direction = "DECREASE" if ev['log_rr'] and ev['log_rr'] < 0 else "INCREASE"
            print(f"  {ev['onset']}: post={ev['post_count']}, pre={ev['pre_count']}, RR={np.exp(ev['log_rr']):.3f} ({direction})")
    
    # Save
    outpath = OUT / "32_extended_accident_ssw.json"
    with open(outpath, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nSaved to {outpath}")
    
    return output


if __name__ == '__main__':
    main()
