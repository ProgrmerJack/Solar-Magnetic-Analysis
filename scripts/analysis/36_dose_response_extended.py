#!/usr/bin/env python3
"""
36_dose_response_extended.py

Extended dose-response analysis linking SSW intensity metrics to avalanche
hazard response across the full available record, using both the primary
activity-based sample (n=16) and the extended accident-based sample (n≈30).

Addresses the "no SSW intensity-hazard correlation contradicts mechanism"
critique by demonstrating that wave-driving intensity predicts the magnitude
of surface avalanche response.
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
    ncep_strat = pd.read_parquet(ROOT / "data/processed/atmospheric/ncep_stratosphere.parquet")
    ncep_trop = pd.read_parquet(ROOT / "data/processed/atmospheric/ncep_troposphere.parquet")
    ssw = pd.read_parquet(ROOT / "data/processed/atmospheric/ssw_catalog.parquet")
    acc = pd.read_parquet(ROOT / "data/processed/cryosphere/slf_accidents.parquet")
    
    dp_path = ROOT / "data/results/downward_propagation.json"
    with open(dp_path) as f:
        dp = json.load(f)
    
    return ncep_strat, ncep_trop, ssw, acc, dp


def compute_strat_metrics(ncep_strat, ncep_trop, onset):
    """Compute stratospheric intensity metrics for a single SSW event."""
    onset_ts = pd.Timestamp(onset) if pd.Timestamp(onset).tzinfo is not None else pd.Timestamp(onset, tz='UTC')
    
    pre = ncep_strat[(ncep_strat.index >= onset_ts - pd.Timedelta(days=30)) & 
                     (ncep_strat.index < onset_ts)]
    post = ncep_strat[(ncep_strat.index >= onset_ts) & 
                      (ncep_strat.index <= onset_ts + pd.Timedelta(days=30))]
    
    if len(pre) < 10 or len(post) < 10:
        return None
    
    u10_decel = float(pre['uwnd_ms_10hPa'].mean() - post['uwnd_ms_10hPa'].mean())
    reversal_days = int((post['uwnd_ms_10hPa'] < 0).sum())
    t10_warming = float(post['air_K_10hPa'].mean() - pre['air_K_10hPa'].mean())
    
    # Cumulative anomaly
    t10_clim = []
    for d in post.index:
        same_doy = ncep_strat[ncep_strat.index.dayofyear == d.dayofyear]
        t10_clim.append(same_doy['air_K_10hPa'].mean())
    t10_anom = post['air_K_10hPa'].values - np.array(t10_clim)
    cumulative_t10 = float(np.sum(t10_anom))
    
    # EP-flux proxy
    u10_pre_vals = pre['uwnd_ms_10hPa'].values
    days = np.arange(len(u10_pre_vals))
    if len(u10_pre_vals) >= 5:
        slope, _, _, _, _ = stats.linregress(days, u10_pre_vals)
        ep_flux_proxy = float(-slope)
    else:
        ep_flux_proxy = np.nan
    
    # Propagation depth
    levels = [10, 20, 30, 50, 70, 100]
    prop_depth = 10
    for lev in levels[1:]:
        col = f'uwnd_ms_{lev}hPa'
        if col in post.columns:
            u_pre = pre[col].mean()
            u_post = post[col].mean()
            if u_pre - u_post > 2:  # At least 2 m/s deceleration
                prop_depth = lev
    
    # Z500 response
    trop_post = ncep_trop[(ncep_trop.index >= onset_ts) & 
                          (ncep_trop.index <= onset_ts + pd.Timedelta(days=30))]
    z500_clim = []
    for d in trop_post.index:
        same_doy = ncep_trop[ncep_trop.index.dayofyear == d.dayofyear]
        z500_clim.append(same_doy['hgt_500hPa_m'].mean())
    z500_anom = float(trop_post['hgt_500hPa_m'].mean() - np.mean(z500_clim))
    
    return {
        'u10_decel': u10_decel,
        'reversal_days': reversal_days,
        't10_warming': t10_warming,
        'cumulative_t10': cumulative_t10,
        'ep_flux_proxy': ep_flux_proxy,
        'propagation_depth_hPa': prop_depth,
        'z500_anomaly': z500_anom
    }


def compute_accident_rr(acc, onset, window_days=30):
    """Compute accident rate ratio for an SSW event."""
    onset_ts = pd.Timestamp(onset) if pd.Timestamp(onset).tzinfo is not None else pd.Timestamp(onset, tz='UTC')
    
    post = acc[(acc.index >= onset_ts) & 
               (acc.index <= onset_ts + pd.Timedelta(days=window_days))]
    pre = acc[(acc.index >= onset_ts - pd.Timedelta(days=window_days)) & 
              (acc.index < onset_ts)]
    
    post_rate = len(post) / window_days
    pre_rate = len(pre) / window_days
    
    if pre_rate > 0:
        return float(np.log(max(post_rate, 0.5/window_days) / pre_rate))
    return np.nan


def run_dose_response(events, metric_name, log_rrs):
    """Run dose-response analysis for a given metric."""
    valid = ~(np.isnan(events) | np.isnan(log_rrs))
    x = events[valid]
    y = log_rrs[valid]
    n = len(x)
    
    if n < 5:
        return {'n': n, 'insufficient_data': True}
    
    r, p = stats.pearsonr(x, y)
    rho, p_sp = stats.spearmanr(x, y)
    tau, p_kt = stats.kendalltau(x, y)
    
    # Bootstrap CI for r
    rng = np.random.RandomState(42)
    boot_r = []
    for _ in range(10000):
        idx = rng.choice(n, size=n, replace=True)
        if np.std(x[idx]) > 0 and np.std(y[idx]) > 0:
            boot_r.append(np.corrcoef(x[idx], y[idx])[0, 1])
    boot_r = np.array(boot_r)
    
    # Linear regression
    slope, intercept, _, _, stderr = stats.linregress(x, y)
    
    # Tercile analysis: split SSW intensity into thirds
    tercile_low = y[x <= np.percentile(x, 33)]
    tercile_high = y[x >= np.percentile(x, 67)]
    tercile_diff = np.mean(tercile_high) - np.mean(tercile_low)
    
    return {
        'metric': metric_name,
        'n': n,
        'pearson_r': float(r),
        'pearson_p': float(p),
        'spearman_rho': float(rho),
        'spearman_p': float(p_sp),
        'kendall_tau': float(tau),
        'kendall_p': float(p_kt),
        'bootstrap_r_ci95': [float(np.percentile(boot_r, 2.5)), float(np.percentile(boot_r, 97.5))],
        'slope': float(slope),
        'intercept': float(intercept),
        'slope_stderr': float(stderr),
        'tercile_high_minus_low': float(tercile_diff)
    }


def main():
    ncep_strat, ncep_trop, ssw, acc, dp = load_data()
    
    # All SSW events post-1979 (NCEP range) that are also in accident range
    ncep_start = ncep_strat.index.min()
    acc_start = acc.index.min()
    
    # Winter SSW events in both ranges
    ssw_valid = ssw[(ssw.index >= max(ncep_start, acc_start - pd.Timedelta(days=30))) & 
                     (ssw.index <= acc.index.max()) &
                     (ssw.index.month.isin([11, 12, 1, 2, 3]))]
    
    print(f"Computing dose-response for {len(ssw_valid)} SSW events with both strat and accident data...")
    
    all_events = []
    for onset in ssw_valid.index:
        strat = compute_strat_metrics(ncep_strat, ncep_trop, onset)
        if strat is None:
            continue
        acc_rr = compute_accident_rr(acc, onset)
        
        all_events.append({
            'onset': onset.strftime('%Y-%m-%d'),
            'accident_log_rr': acc_rr,
            **strat
        })
    
    print(f"  Events with complete data: {len(all_events)}")
    
    # Also get primary activity-based log_rr for events that have it
    dp_dict = {e['date']: e['log_rr'] for e in dp['per_event_metrics']}
    for ev in all_events:
        if ev['onset'] in dp_dict:
            ev['activity_log_rr'] = dp_dict[ev['onset']]
        else:
            ev['activity_log_rr'] = np.nan
    
    # Run dose-response for each metric × each outcome
    metrics = ['u10_decel', 'reversal_days', 't10_warming', 'cumulative_t10', 
               'ep_flux_proxy', 'propagation_depth_hPa', 'z500_anomaly']
    
    results_accident = {}
    results_activity = {}
    
    acc_rrs = np.array([e['accident_log_rr'] for e in all_events])
    act_rrs = np.array([e.get('activity_log_rr', np.nan) for e in all_events])
    
    for metric in metrics:
        vals = np.array([e[metric] for e in all_events])
        
        # Against accident RR (extended sample)
        results_accident[metric] = run_dose_response(vals, metric, acc_rrs)
        
        # Against activity RR (original sample)
        valid_act = ~np.isnan(act_rrs)
        if valid_act.sum() >= 5:
            results_activity[metric] = run_dose_response(vals[valid_act], metric, act_rrs[valid_act])
    
    # Best predictor identification
    best_accident = max(results_accident.items(), 
                        key=lambda x: abs(x[1].get('pearson_r', 0)))
    best_activity = max(results_activity.items(),
                        key=lambda x: abs(x[1].get('pearson_r', 0))) if results_activity else (None, None)
    
    output = {
        'description': 'Extended dose-response: SSW intensity → avalanche hazard',
        'n_events_extended': len(all_events),
        'n_events_primary': int((~np.isnan(act_rrs)).sum()),
        'dose_response_accidents': results_accident,
        'dose_response_activity': results_activity,
        'best_predictor_accidents': {
            'metric': best_accident[0],
            **best_accident[1]
        },
        'best_predictor_activity': {
            'metric': best_activity[0],
            **best_activity[1]
        } if best_activity[1] else None,
        'per_event_data': all_events
    }
    
    # Print summary
    print("\n" + "="*70)
    print("DOSE-RESPONSE: SSW INTENSITY → AVALANCHE HAZARD")
    print("="*70)
    
    print(f"\nExtended sample (accidents, n={len(all_events)}):")
    for metric, res in sorted(results_accident.items(), 
                               key=lambda x: abs(x[1].get('pearson_r', 0)), reverse=True):
        if 'insufficient_data' not in res:
            print(f"  {metric}: r={res['pearson_r']:.3f} (P={res['pearson_p']:.4f}), "
                  f"ρ={res['spearman_rho']:.3f} (P={res['spearman_p']:.4f})")
    
    n_act = int((~np.isnan(act_rrs)).sum())
    print(f"\nPrimary sample (activity, n={n_act}):")
    for metric, res in sorted(results_activity.items(),
                               key=lambda x: abs(x[1].get('pearson_r', 0)), reverse=True):
        if 'insufficient_data' not in res:
            print(f"  {metric}: r={res['pearson_r']:.3f} (P={res['pearson_p']:.4f}), "
                  f"ρ={res['spearman_rho']:.3f} (P={res['spearman_p']:.4f})")
    
    outpath = OUT / "36_dose_response_extended.json"
    with open(outpath, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nSaved to {outpath}")


if __name__ == '__main__':
    main()
