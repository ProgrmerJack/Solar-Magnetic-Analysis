#!/usr/bin/env python3
"""Script 47: Download NCEP 100hPa eddy heat flux (v'T') for direct EP-flux computation.

Downloads temperature and meridional wind at 100hPa from NCEP-NCAR reanalysis
via OPeNDAP, computes zonal-mean poleward eddy heat flux [v'T'] at 45-75N,
and correlates with event-level avalanche response.

This addresses 4/5 reviewers' concern that the EP-flux proxy uses vortex 
deceleration rather than direct v'T' computation.
"""

import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
from scipy import stats
import json
import warnings
warnings.filterwarnings('ignore')

ROOT = Path("C:/Users/Jack0/Solar-Magnetic-Analysis")

# NCEP-NCAR Reanalysis OPeNDAP URLs
# Daily pressure-level data
NCEP_BASE = "https://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/pressure"

# Load SSW catalog
ssw_cat = pd.read_parquet(ROOT / "data/processed/atmospheric/ssw_catalog.parquet")

# Study period events (1999-2019)
activity = pd.read_parquet(ROOT / "data/processed/cryosphere/slf_activity.parquet")
activity_start = activity.index.min().tz_localize(None)
activity_end = activity.index.max().tz_localize(None)

ssw_events = []
for onset_date in ssw_cat.index:
    od = onset_date.tz_localize(None) if onset_date.tzinfo else onset_date
    if activity_start <= od <= activity_end:
        ssw_events.append(od)

print(f"SSW events in study period: {len(ssw_events)}")
print("Attempting to access NCEP 100hPa data via OPeNDAP...")

# Try to access NCEP data for each event
# We need temperature and v-wind at 100hPa, 45-75N latitude band
# NCEP data is organized by year

def compute_vt_flux_for_event(onset, pre_window=10):
    """Compute v'T' heat flux at 100hPa for an SSW event.
    
    Downloads a small subset via OPeNDAP for the pre-onset window.
    v'T' = [v*T] - [v]*[T] where [] is zonal mean and * is deviation.
    """
    year = onset.year
    
    # Time window: pre_window days before onset
    t_start = onset - pd.Timedelta(days=pre_window)
    t_end = onset
    
    try:
        # Open temperature dataset
        temp_url = f"{NCEP_BASE}/air.{year}.nc"
        vwnd_url = f"{NCEP_BASE}/vwnd.{year}.nc"
        
        # Open with OPeNDAP - select only 100hPa, 45-75N
        ds_t = xr.open_dataset(temp_url)
        ds_v = xr.open_dataset(vwnd_url)
        
        # Select 100hPa level, 45-75N, time window
        # NCEP levels are in hPa (or mbar)
        t_sel = ds_t['air'].sel(
            level=100, 
            lat=slice(75, 45),  # NCEP lat goes from 90 to -90
            time=slice(str(t_start.date()), str(t_end.date()))
        )
        v_sel = ds_v['vwnd'].sel(
            level=100,
            lat=slice(75, 45),
            time=slice(str(t_start.date()), str(t_end.date()))
        )
        
        # Compute eddy heat flux: v'T' = mean(v*T) - mean(v)*mean(T)
        # Zonal mean (average over longitude)
        vT = (v_sel * t_sel).mean(dim='lon')  # [v*T] zonal mean
        v_bar = v_sel.mean(dim='lon')  # [v]
        T_bar = t_sel.mean(dim='lon')  # [T]
        
        eddy_vT = vT - v_bar * T_bar  # v'T' = [vT] - [v][T]
        
        # Average over latitude band (45-75N) with cos(lat) weighting
        lats = eddy_vT.lat.values
        weights = np.cos(np.deg2rad(lats))
        weights = weights / weights.sum()
        
        # Time-mean eddy heat flux
        vT_mean = float((eddy_vT.mean(dim='time') * weights).sum())
        
        # Also compute daily time series for the window
        daily_vT = []
        for t in range(len(eddy_vT.time)):
            day_vT = float((eddy_vT.isel(time=t) * weights).sum())
            daily_vT.append(day_vT)
        
        ds_t.close()
        ds_v.close()
        
        return {
            'vT_mean': vT_mean,
            'vT_daily': daily_vT,
            'n_days': len(daily_vT),
            'success': True
        }
    except Exception as e:
        return {'success': False, 'error': str(e)}

# Process events - try one first to test connectivity
print(f"\nTesting OPeNDAP access with first event ({ssw_events[0].date()})...")
test = compute_vt_flux_for_event(ssw_events[0])

if test['success']:
    print(f"  SUCCESS! v'T' = {test['vT_mean']:.2f} K·m/s")
    
    # Process all events
    ep_flux_results = []
    for i, onset in enumerate(ssw_events):
        print(f"  Processing event {i+1}/{len(ssw_events)}: {onset.date()}...", end="")
        result = compute_vt_flux_for_event(onset)
        if result['success']:
            print(f" v'T' = {result['vT_mean']:.2f} K·m/s")
            ep_flux_results.append({
                'onset': str(onset.date()),
                'vT_mean': result['vT_mean'],
                'n_days': result['n_days'],
            })
        else:
            print(f" FAILED: {result['error'][:50]}")
    
    print(f"\n{len(ep_flux_results)}/{len(ssw_events)} events processed successfully")
    
    # Now correlate with event-level avalanche RR
    # Load the event-level RR data from mediation results
    med_path = ROOT / "data/results/46_formal_mediation.json"
    if med_path.exists():
        with open(med_path) as f:
            med = json.load(f)
        event_rr = {e['onset']: e['logRR'] for e in med['event_data']}
    else:
        # Compute RR from the primary analysis
        print("Mediation results not yet available, using placeholder")
        event_rr = {}
    
    # Match EP-flux to RR
    matched = []
    for ep in ep_flux_results:
        if ep['onset'] in event_rr:
            matched.append({
                'onset': ep['onset'],
                'vT_mean': ep['vT_mean'],
                'logRR': event_rr[ep['onset']],
            })
    
    if len(matched) >= 5:
        vT_arr = np.array([m['vT_mean'] for m in matched])
        logRR_arr = np.array([m['logRR'] for m in matched])
        
        # Correlation: stronger v'T' (more positive = more upward wave flux)
        # should correlate with more negative logRR (more suppression)
        rho, p_rho = stats.spearmanr(vT_arr, logRR_arr)
        r_pear, p_pear = stats.pearsonr(vT_arr, logRR_arr)
        
        print(f"\n=== Direct EP-flux dose-response ===")
        print(f"  n = {len(matched)} events")
        print(f"  Spearman ρ = {rho:.3f}, P = {p_rho:.4f}")
        print(f"  Pearson r = {r_pear:.3f}, P = {p_pear:.4f}")
        print(f"  v'T' range: [{min(vT_arr):.2f}, {max(vT_arr):.2f}] K·m/s")
    else:
        rho, p_rho, r_pear, p_pear = None, None, None, None
        print(f"Only {len(matched)} matched events, insufficient for correlation")
    
    # Save results
    output = {
        'n_events_processed': len(ep_flux_results),
        'n_events_total': len(ssw_events),
        'ep_flux_data': ep_flux_results,
        'matched_events': matched,
        'correlation': {
            'spearman_rho': float(rho) if rho is not None else None,
            'spearman_p': float(p_rho) if p_rho is not None else None,
            'pearson_r': float(r_pear) if r_pear is not None else None,
            'pearson_p': float(p_pear) if p_pear is not None else None,
            'n': len(matched),
        },
        'method': 'Direct v\'T\' at 100hPa, 45-75N, cos(lat)-weighted, 10-day pre-onset window',
        'source': 'NCEP-NCAR Reanalysis via OPeNDAP',
    }
    
    out_path = ROOT / "data/results/47_direct_ep_flux.json"
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nResults saved to {out_path}")

else:
    print(f"OPeNDAP access failed: {test.get('error', 'unknown')}")
    print("Falling back to alternative approach...")
    
    # Alternative: Compute from ERA5 if available, or use existing proxy
    # with better framing
    output = {
        'success': False,
        'error': test.get('error', 'OPeNDAP access failed'),
        'fallback': 'Use existing vortex-deceleration proxy with stronger framing as exploratory',
    }
    
    out_path = ROOT / "data/results/47_direct_ep_flux.json"
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nFallback results saved to {out_path}")
