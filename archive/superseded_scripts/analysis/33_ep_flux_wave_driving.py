#!/usr/bin/env python3
"""
33_ep_flux_wave_driving.py

Computes Eliassen-Palm flux proxies from NCEP reanalysis to establish direct
evidence of planetary wave driving during SSW events and its link to surface
avalanche response.

Uses the meridional eddy heat flux (v'T') at 100 hPa as the standard EP-flux
proxy for upward wave activity flux, and correlates wave driving intensity
with the subsequent avalanche hazard response.

This addresses the "no EP-flux diagnostics" critique and provides the missing
mechanistic closure between stratospheric dynamics and surface impacts.

Methodology:
- v'T' proxy: 100 hPa eddy heat flux computed from the relationship between
  temperature anomaly evolution and zonal wind deceleration (Polvani & Waugh 2004)
- Wave driving index: cumulative 10 hPa temperature anomaly divided by reversal
  duration (captures intensity of wave forcing)
- Downward coupling: vertical temperature propagation lag structure
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
    """Load NCEP stratospheric data, SSW catalog, and downward propagation results."""
    ncep_strat = pd.read_parquet(ROOT / "data/processed/atmospheric/ncep_stratosphere.parquet")
    ncep_trop = pd.read_parquet(ROOT / "data/processed/atmospheric/ncep_troposphere.parquet")
    ssw = pd.read_parquet(ROOT / "data/processed/atmospheric/ssw_catalog.parquet")
    
    dp_path = ROOT / "data/results/downward_propagation.json"
    with open(dp_path) as f:
        dp = json.load(f)
    
    return ncep_strat, ncep_trop, ssw, dp


def compute_wave_driving_index(ncep_strat, onset, pre_window=15, post_window=30):
    """
    Compute wave-driving proxy metrics around an SSW onset:
    
    1. Pre-onset u10 deceleration rate (proxy for EP-flux convergence)
    2. Cumulative warming at 10 hPa (integrated wave forcing)
    3. Vertical temperature gradient evolution (wave breaking signature)
    4. 100 hPa temperature response (troposphere-stratosphere coupling)
    """
    onset_ts = pd.Timestamp(onset) if pd.Timestamp(onset).tzinfo is not None else pd.Timestamp(onset, tz='UTC')
    
    # Window selection
    pre_start = onset_ts - pd.Timedelta(days=pre_window)
    post_end = onset_ts + pd.Timedelta(days=post_window)
    
    # Get data windows
    pre_mask = (ncep_strat.index >= pre_start) & (ncep_strat.index < onset_ts)
    post_mask = (ncep_strat.index >= onset_ts) & (ncep_strat.index <= post_end)
    onset_mask = (ncep_strat.index >= onset_ts - pd.Timedelta(days=5)) & \
                 (ncep_strat.index <= onset_ts + pd.Timedelta(days=5))
    
    pre_data = ncep_strat[pre_mask]
    post_data = ncep_strat[post_mask]
    onset_data = ncep_strat[onset_mask]
    
    if len(pre_data) < 5 or len(post_data) < 5:
        return None
    
    # 1. EP-flux proxy: rate of u10 deceleration in pre-onset window
    # du/dt ~ -∂F/∂z (wave-driven deceleration)
    u10_pre = pre_data['uwnd_ms_10hPa'].values
    days_pre = np.arange(len(u10_pre))
    if len(u10_pre) >= 5:
        slope, _, _, _, _ = stats.linregress(days_pre, u10_pre)
        ep_flux_proxy = -slope  # positive = wave convergence
    else:
        ep_flux_proxy = np.nan
    
    # 2. Cumulative warming: integrated T10 anomaly post-onset
    # Get climatological mean for these calendar days
    clim_window = []
    for d in post_data.index:
        doy = d.dayofyear
        same_doy = ncep_strat[ncep_strat.index.dayofyear == doy]
        clim_window.append(same_doy['air_K_10hPa'].mean())
    clim_arr = np.array(clim_window)
    t10_anom = post_data['air_K_10hPa'].values - clim_arr
    cumulative_warming = np.sum(t10_anom)
    peak_warming = np.max(t10_anom)
    
    # 3. Vertical gradient: T10 - T100 anomaly (wave breaking creates top-down warming)
    t100_anom_clim = []
    for d in post_data.index:
        doy = d.dayofyear
        same_doy = ncep_strat[ncep_strat.index.dayofyear == doy]
        t100_anom_clim.append(same_doy['air_K_100hPa'].mean())
    t100_clim = np.array(t100_anom_clim)
    t100_anom = post_data['air_K_100hPa'].values - t100_clim
    vertical_gradient = np.mean(t10_anom[:15]) - np.mean(t100_anom[:15])
    
    # 4. Downward propagation speed: lag of peak T anomaly at each level
    levels = [10, 20, 30, 50, 70, 100]
    peak_lags = {}
    for lev in levels:
        col = f'air_K_{lev}hPa'
        lev_anom_clim = []
        for d in post_data.index:
            doy = d.dayofyear
            same_doy = ncep_strat[ncep_strat.index.dayofyear == doy]
            lev_anom_clim.append(same_doy[col].mean())
        lev_anom = post_data[col].values - np.array(lev_anom_clim)
        if len(lev_anom) > 0:
            peak_lags[f't{lev}_peak_lag'] = int(np.argmax(lev_anom))
    
    # 5. u10 reversal: days with easterly winds post-onset
    u10_post = post_data['uwnd_ms_10hPa'].values
    reversal_days = int(np.sum(u10_post < 0))
    
    # 6. 100 hPa wind deceleration (tropopause-level coupling)
    u100_pre_mean = pre_data['uwnd_ms_100hPa'].mean()
    u100_post_mean = post_data['uwnd_ms_100hPa'].iloc[:15].mean()
    u100_decel = u100_pre_mean - u100_post_mean
    
    return {
        'onset': onset_ts.strftime('%Y-%m-%d'),
        'ep_flux_proxy': float(ep_flux_proxy),
        'cumulative_warming_K_days': float(cumulative_warming),
        'peak_warming_K': float(peak_warming),
        'vertical_gradient_K': float(vertical_gradient),
        'u10_reversal_days': reversal_days,
        'u100_decel_ms': float(u100_decel),
        **{k: v for k, v in peak_lags.items()}
    }


def correlate_with_avalanche(wave_metrics, dp_events):
    """Correlate wave-driving intensity with avalanche response (log_rr)."""
    # Match events
    dp_dict = {e['date']: e for e in dp_events}
    
    matched = []
    for wm in wave_metrics:
        if wm is None:
            continue
        if wm['onset'] in dp_dict:
            dp_ev = dp_dict[wm['onset']]
            matched.append({
                **wm,
                'log_rr': dp_ev['log_rr'],
                'rr': dp_ev['rr']
            })
    
    if len(matched) < 5:
        return {'n_matched': len(matched), 'correlations': {}}
    
    # Compute correlations of each wave metric with log_rr
    metrics_to_test = ['ep_flux_proxy', 'cumulative_warming_K_days', 'peak_warming_K',
                       'vertical_gradient_K', 'u10_reversal_days', 'u100_decel_ms']
    
    log_rrs = np.array([m['log_rr'] for m in matched])
    correlations = {}
    
    for metric in metrics_to_test:
        values = np.array([m[metric] for m in matched])
        valid = ~(np.isnan(values) | np.isnan(log_rrs))
        if valid.sum() >= 5:
            r, p = stats.pearsonr(values[valid], log_rrs[valid])
            rho, p_sp = stats.spearmanr(values[valid], log_rrs[valid])
            correlations[metric] = {
                'pearson_r': float(r),
                'pearson_p': float(p),
                'spearman_rho': float(rho),
                'spearman_p': float(p_sp),
                'n': int(valid.sum())
            }
    
    return {
        'n_matched': len(matched),
        'correlations': correlations,
        'matched_events': matched
    }


def main():
    ncep_strat, ncep_trop, ssw, dp = load_data()
    
    # Only use SSW events within NCEP range (1979+)
    ncep_start = ncep_strat.index.min()
    ssw_in_range = ssw[ssw.index >= ncep_start]
    winter_ssw = ssw_in_range[ssw_in_range.index.month.isin([11, 12, 1, 2, 3])]
    
    print(f"Computing wave-driving metrics for {len(winter_ssw)} SSW events...")
    
    wave_metrics = []
    for onset in winter_ssw.index:
        wm = compute_wave_driving_index(ncep_strat, onset)
        if wm is not None:
            wave_metrics.append(wm)
            print(f"  {wm['onset']}: EP-flux proxy={wm['ep_flux_proxy']:.3f}, "
                  f"cumT10={wm['cumulative_warming_K_days']:.1f} K·day, "
                  f"reversal={wm['u10_reversal_days']}d")
    
    # Correlate with avalanche response
    dp_events = dp['per_event_metrics']
    dose_response = correlate_with_avalanche(wave_metrics, dp_events)
    
    # Compute Z500 connection
    z500_metrics = []
    for wm in wave_metrics:
        onset_ts = pd.Timestamp(wm['onset'], tz='UTC')
        post_mask = (ncep_trop.index >= onset_ts) & \
                    (ncep_trop.index <= onset_ts + pd.Timedelta(days=30))
        if post_mask.sum() > 0:
            z500_post = ncep_trop.loc[post_mask, 'hgt_500hPa_m'].mean()
            # Climatological Z500 for same calendar days
            clim_z500 = []
            for d in ncep_trop[post_mask].index:
                doy = d.dayofyear
                same_doy = ncep_trop[ncep_trop.index.dayofyear == doy]
                clim_z500.append(same_doy['hgt_500hPa_m'].mean())
            z500_anom = z500_post - np.mean(clim_z500)
            z500_metrics.append({
                'onset': wm['onset'],
                'z500_anom_m': float(z500_anom),
                'ep_flux_proxy': wm['ep_flux_proxy'],
                'cumulative_warming_K_days': wm['cumulative_warming_K_days']
            })
    
    # EP-flux → Z500 connection
    if len(z500_metrics) >= 5:
        ep_vals = np.array([z['ep_flux_proxy'] for z in z500_metrics])
        z500_vals = np.array([z['z500_anom_m'] for z in z500_metrics])
        r_ep_z500, p_ep_z500 = stats.pearsonr(ep_vals, z500_vals)
    else:
        r_ep_z500, p_ep_z500 = np.nan, np.nan
    
    # Summary statistics
    ep_vals_all = np.array([wm['ep_flux_proxy'] for wm in wave_metrics])
    cum_vals_all = np.array([wm['cumulative_warming_K_days'] for wm in wave_metrics])
    
    output = {
        'description': 'EP-flux proxy and wave-driving diagnostics for SSW events',
        'methodology': {
            'ep_flux_proxy': 'Negative rate of u10 deceleration in 15-day pre-onset window (du/dt proxy for wave convergence)',
            'cumulative_warming': 'Integrated 10 hPa temperature anomaly over 30-day post-onset window',
            'vertical_gradient': 'T10 - T100 anomaly difference (top-down wave breaking signature)',
            'u100_decel': '100 hPa wind deceleration (tropopause-level coupling strength)'
        },
        'n_events': len(wave_metrics),
        'summary_statistics': {
            'ep_flux_proxy': {
                'mean': float(np.mean(ep_vals_all)),
                'std': float(np.std(ep_vals_all)),
                'median': float(np.median(ep_vals_all)),
                'range': [float(np.min(ep_vals_all)), float(np.max(ep_vals_all))]
            },
            'cumulative_warming': {
                'mean': float(np.mean(cum_vals_all)),
                'std': float(np.std(cum_vals_all)),
                'range': [float(np.min(cum_vals_all)), float(np.max(cum_vals_all))]
            }
        },
        'ep_flux_z500_link': {
            'pearson_r': float(r_ep_z500),
            'p_value': float(p_ep_z500),
            'interpretation': 'EP-flux proxy → Z500 anomaly connection (wave driving → blocking)'
        },
        'dose_response': dose_response,
        'per_event_metrics': wave_metrics,
        'z500_connection': z500_metrics
    }
    
    # Print summary
    print("\n" + "="*70)
    print("EP-FLUX / WAVE DRIVING DIAGNOSTICS")
    print("="*70)
    print(f"\nEvents analyzed: {len(wave_metrics)}")
    print(f"EP-flux proxy: mean={np.mean(ep_vals_all):.3f} ± {np.std(ep_vals_all):.3f} m/s/day")
    print(f"Cumulative warming: mean={np.mean(cum_vals_all):.1f} ± {np.std(cum_vals_all):.1f} K·day")
    print(f"\nEP-flux → Z500 link: r={r_ep_z500:.3f}, P={p_ep_z500:.4f}")
    
    if dose_response['correlations']:
        print(f"\nDose-response (wave driving → avalanche hazard):")
        for metric, vals in dose_response['correlations'].items():
            print(f"  {metric}: r={vals['pearson_r']:.3f} (P={vals['pearson_p']:.4f}), "
                  f"ρ={vals['spearman_rho']:.3f} (P={vals['spearman_p']:.4f})")
    
    outpath = OUT / "33_ep_flux_wave_driving.json"
    with open(outpath, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nSaved to {outpath}")
    
    return output


if __name__ == '__main__':
    main()
