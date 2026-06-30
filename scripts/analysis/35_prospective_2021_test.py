#!/usr/bin/env python3
"""
35_prospective_2021_test.py

True out-of-sample prospective test using the January 2021 SSW event,
which occurred AFTER the primary study period (1998-2019).

Tests whether the SSW-avalanche association predicted by the primary analysis
holds for:
1. Swiss avalanche accidents (SLF, available through 2025)
2. NCEP stratospheric diagnostics (confirms SSW characteristics)
3. NCEP tropospheric response (Z500, SLP, u850)

This is the single most important analysis for addressing the "no prospective
validation" critique, as it uses data that was physically unavailable during
the primary analysis period.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "results"


def load_data():
    ssw = pd.read_parquet(ROOT / "data/processed/atmospheric/ssw_catalog.parquet")
    acc = pd.read_parquet(ROOT / "data/processed/cryosphere/slf_accidents.parquet")
    ncep_strat = pd.read_parquet(ROOT / "data/processed/atmospheric/ncep_stratosphere.parquet")
    ncep_trop = pd.read_parquet(ROOT / "data/processed/atmospheric/ncep_troposphere.parquet")
    
    dp_path = ROOT / "data/results/downward_propagation.json"
    with open(dp_path) as f:
        dp = json.load(f)
    
    return ssw, acc, ncep_strat, ncep_trop, dp


def analyze_2021_ssw(ncep_strat, ncep_trop, onset):
    """Characterize the Jan 2021 SSW event from NCEP data."""
    onset_ts = pd.Timestamp(onset) if pd.Timestamp(onset).tzinfo is not None else pd.Timestamp(onset, tz='UTC')
    
    pre = ncep_strat[(ncep_strat.index >= onset_ts - pd.Timedelta(days=30)) & 
                     (ncep_strat.index < onset_ts)]
    post = ncep_strat[(ncep_strat.index >= onset_ts) & 
                      (ncep_strat.index <= onset_ts + pd.Timedelta(days=30))]
    
    trop_post = ncep_trop[(ncep_trop.index >= onset_ts) & 
                          (ncep_trop.index <= onset_ts + pd.Timedelta(days=30))]
    trop_pre = ncep_trop[(ncep_trop.index >= onset_ts - pd.Timedelta(days=30)) & 
                         (ncep_trop.index < onset_ts)]
    
    # Stratospheric characterization
    u10_pre = pre['uwnd_ms_10hPa'].mean()
    u10_post = post['uwnd_ms_10hPa'].mean()
    u10_min = post['uwnd_ms_10hPa'].min()
    t10_pre = pre['air_K_10hPa'].mean()
    t10_post = post['air_K_10hPa'].mean()
    reversal_days = int((post['uwnd_ms_10hPa'] < 0).sum())
    
    # Climatological anomalies
    t10_clim = []
    for d in post.index:
        same_doy = ncep_strat[ncep_strat.index.dayofyear == d.dayofyear]
        t10_clim.append(same_doy['air_K_10hPa'].mean())
    t10_anom = post['air_K_10hPa'].values - np.array(t10_clim)
    
    # Tropospheric response
    z500_anom_vals = []
    for d in trop_post.index:
        same_doy = ncep_trop[ncep_trop.index.dayofyear == d.dayofyear]
        z500_anom_vals.append(trop_post.loc[d, 'hgt_500hPa_m'] - same_doy['hgt_500hPa_m'].mean())
    z500_anom = np.mean(z500_anom_vals)
    
    slp_anom_vals = []
    for d in trop_post.index:
        same_doy = ncep_trop[ncep_trop.index.dayofyear == d.dayofyear]
        slp_anom_vals.append(trop_post.loc[d, 'slp_Pa'] - same_doy['slp_Pa'].mean())
    slp_anom = np.mean(slp_anom_vals)
    
    return {
        'u10_pre_ms': float(u10_pre),
        'u10_post_ms': float(u10_post),
        'u10_decel_ms': float(u10_pre - u10_post),
        'u10_min_ms': float(u10_min),
        'reversal_days': reversal_days,
        't10_warming_K': float(t10_post - t10_pre),
        'cumulative_t10_anom': float(np.sum(t10_anom)),
        'peak_t10_anom': float(np.max(t10_anom)),
        'z500_anomaly_m': float(z500_anom),
        'slp_anomaly_Pa': float(slp_anom)
    }


def analyze_2021_accidents(acc, onset, window_days=30):
    """Analyze accident rates around the Jan 2021 SSW."""
    onset_ts = pd.Timestamp(onset) if pd.Timestamp(onset).tzinfo is not None else pd.Timestamp(onset, tz='UTC')
    
    post_start = onset_ts
    post_end = onset_ts + pd.Timedelta(days=window_days)
    pre_start = onset_ts - pd.Timedelta(days=window_days)
    pre_end = onset_ts - pd.Timedelta(days=1)
    
    post_acc = acc[(acc.index >= post_start) & (acc.index <= post_end)]
    pre_acc = acc[(acc.index >= pre_start) & (acc.index <= pre_end)]
    
    post_count = len(post_acc)
    pre_count = len(pre_acc)
    
    post_rate = post_count / window_days
    pre_rate = pre_count / window_days
    
    if pre_rate > 0:
        rr = post_rate / pre_rate
        log_rr = np.log(rr)
    else:
        rr = np.nan
        log_rr = np.nan
    
    # Activity breakdown
    post_activities = post_acc['activity'].value_counts().to_dict() if len(post_acc) > 0 else {}
    pre_activities = pre_acc['activity'].value_counts().to_dict() if len(pre_acc) > 0 else {}
    
    # Canton distribution  
    post_cantons = post_acc['canton'].value_counts().to_dict() if len(post_acc) > 0 else {}
    pre_cantons = pre_acc['canton'].value_counts().to_dict() if len(pre_acc) > 0 else {}
    
    # Fatalities
    post_dead = int(post_acc['number_dead'].sum()) if len(post_acc) > 0 else 0
    pre_dead = int(pre_acc['number_dead'].sum()) if len(pre_acc) > 0 else 0
    
    return {
        'window_days': window_days,
        'post_ssw_count': post_count,
        'pre_ssw_count': pre_count,
        'post_rate_per_day': float(post_rate),
        'pre_rate_per_day': float(pre_rate),
        'rate_ratio': float(rr) if not np.isnan(rr) else None,
        'log_rate_ratio': float(log_rr) if not np.isnan(log_rr) else None,
        'post_fatalities': post_dead,
        'pre_fatalities': pre_dead,
        'post_activities': {str(k): int(v) for k, v in post_activities.items()},
        'pre_activities': {str(k): int(v) for k, v in pre_activities.items()},
        'post_cantons': {str(k): int(v) for k, v in post_cantons.items()},
        'pre_cantons': {str(k): int(v) for k, v in pre_cantons.items()},
        'direction': 'DECREASE' if log_rr < 0 else 'INCREASE' if log_rr > 0 else 'NEUTRAL'
    }


def compare_to_primary_distribution(log_rr_2021, dp_events):
    """Compare 2021 result to the primary study distribution."""
    primary_log_rrs = np.array([e['log_rr'] for e in dp_events])
    
    # Where does 2021 fall in the primary distribution?
    percentile = float(stats.percentileofscore(primary_log_rrs, log_rr_2021))
    mean_primary = float(np.mean(primary_log_rrs))
    std_primary = float(np.std(primary_log_rrs))
    z_score = (log_rr_2021 - mean_primary) / std_primary if std_primary > 0 else np.nan
    
    # Is it within the bootstrap CI?
    rng = np.random.RandomState(42)
    boot_means = [np.mean(rng.choice(primary_log_rrs, size=len(primary_log_rrs), replace=True)) 
                  for _ in range(10000)]
    ci_lower = np.percentile(boot_means, 2.5)
    ci_upper = np.percentile(boot_means, 97.5)
    within_ci = ci_lower <= log_rr_2021 <= ci_upper
    
    return {
        'primary_mean_log_rr': mean_primary,
        'primary_std': std_primary,
        '2021_log_rr': float(log_rr_2021),
        '2021_percentile_in_primary': percentile,
        '2021_z_score': float(z_score) if not np.isnan(z_score) else None,
        'within_primary_95ci': bool(within_ci),
        'primary_boot_ci': [float(ci_lower), float(ci_upper)],
        'consistent_with_primary': percentile >= 5 and percentile <= 95
    }


def main():
    ssw, acc, ncep_strat, ncep_trop, dp = load_data()
    
    # The Jan 2021 SSW event
    onset_2021 = pd.Timestamp('2021-01-05', tz='UTC')
    print(f"Analyzing Jan 2021 SSW event (onset: {onset_2021.strftime('%Y-%m-%d')})")
    print(f"This is a TRUE OUT-OF-SAMPLE test — primary study ended May 2019")
    
    # 1. Stratospheric characterization
    strat_2021 = analyze_2021_ssw(ncep_strat, ncep_trop, onset_2021)
    print(f"\nStratospheric characterization:")
    print(f"  u10 deceleration: {strat_2021['u10_decel_ms']:.1f} m/s")
    print(f"  u10 minimum: {strat_2021['u10_min_ms']:.1f} m/s")
    print(f"  Reversal days: {strat_2021['reversal_days']}")
    print(f"  T10 warming: {strat_2021['t10_warming_K']:.1f} K")
    print(f"  Z500 anomaly: {strat_2021['z500_anomaly_m']:.1f} m")
    
    # 2. Accident analysis at multiple windows
    windows = [15, 30, 45]
    accident_results = {}
    for w in windows:
        result = analyze_2021_accidents(acc, onset_2021, window_days=w)
        accident_results[f'{w}day'] = result
        direction = result['direction']
        rr = result['rate_ratio']
        print(f"\n  {w}-day window: pre={result['pre_ssw_count']}, post={result['post_ssw_count']}, "
              f"RR={'N/A' if rr is None else f'{rr:.3f}'} ({direction})")
    
    # 3. Compare to primary distribution
    primary_result = accident_results['30day']
    comparison = None
    if primary_result['log_rate_ratio'] is not None:
        comparison = compare_to_primary_distribution(
            primary_result['log_rate_ratio'], dp['per_event_metrics'])
        print(f"\n  Percentile in primary distribution: {comparison['2021_percentile_in_primary']:.1f}%")
        print(f"  Consistent with primary findings: {comparison['consistent_with_primary']}")
    
    # 4. Also check any other post-2019 winter periods for broader context
    post2019_winters = [
        ('2019-20', '2019-12-01', '2020-03-31'),
        ('2020-21', '2020-12-01', '2021-03-31'),
        ('2021-22', '2021-12-01', '2022-03-31'),
        ('2022-23', '2022-12-01', '2023-03-31'),
        ('2023-24', '2023-12-01', '2024-03-31'),
        ('2024-25', '2024-12-01', '2025-03-31'),
    ]
    
    winter_context = []
    for name, start, end in post2019_winters:
        start_ts = pd.Timestamp(start, tz='UTC')
        end_ts = pd.Timestamp(end, tz='UTC')
        winter_acc = acc[(acc.index >= start_ts) & (acc.index <= end_ts)]
        n_days = (min(end_ts, acc.index.max()) - start_ts).days
        if n_days > 0:
            winter_context.append({
                'winter': name,
                'n_accidents': len(winter_acc),
                'n_days': n_days,
                'daily_rate': len(winter_acc) / n_days,
                'n_fatalities': int(winter_acc['number_dead'].sum()),
                'has_ssw': name == '2020-21'
            })
    
    print(f"\nPost-2019 winter context:")
    for wc in winter_context:
        ssw_flag = " *** SSW WINTER ***" if wc['has_ssw'] else ""
        print(f"  {wc['winter']}: {wc['n_accidents']} accidents in {wc['n_days']}d "
              f"({wc['daily_rate']:.2f}/day, {wc['n_fatalities']} deaths){ssw_flag}")
    
    output = {
        'description': 'Prospective out-of-sample test: Jan 2021 SSW event',
        'key_finding': 'True prospective validation using data unavailable during primary analysis',
        'ssw_onset': '2021-01-05',
        'primary_study_end': '2019-05-21',
        'stratospheric_characterization': strat_2021,
        'accident_analysis': accident_results,
        'comparison_to_primary': comparison,
        'post2019_winter_context': winter_context
    }
    
    print("\n" + "="*70)
    print("PROSPECTIVE 2021 TEST SUMMARY")
    print("="*70)
    rr30 = accident_results['30day']['rate_ratio']
    print(f"30-day accident RR: {'N/A' if rr30 is None else f'{rr30:.3f}'}")
    print(f"Direction: {accident_results['30day']['direction']}")
    if comparison:
        print(f"Consistent with primary: {comparison['consistent_with_primary']}")
    
    outpath = OUT / "35_prospective_2021_test.json"
    with open(outpath, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nSaved to {outpath}")


if __name__ == '__main__':
    main()
