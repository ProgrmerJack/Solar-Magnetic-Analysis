#!/usr/bin/env python3
"""Script 51: Cross-Hazard Regime-Redistribution Analysis

Shows that SSW-driven regime redistribution simultaneously shifts multiple
hazard-relevant weather variables, demonstrating the loaded-gun mechanism
is a GENERAL PRINCIPLE beyond avalanches.

Hazard domains tested:
  1. Cold-wave mortality risk (temperature deficit)
  2. Infrastructure/transport disruption (wind + snowfall)
  3. Rain-on-snow flood/landslide risk (rainfall fraction)
  4. Drought/fire risk (precipitation deficit in warm sector)
  5. Air quality (temperature inversion proxy)
"""
import json, warnings
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

warnings.filterwarnings('ignore')

ROOT = Path(__file__).resolve().parents[2]

# Load data
era5 = pd.read_parquet(ROOT / 'data/processed/era5_swiss_alps_extended.parquet')
ssw_cat = pd.read_parquet(ROOT / 'data/processed/atmospheric/ssw_catalog.parquet')
ncep = pd.read_parquet(ROOT / 'data/processed/atmospheric/ncep_troposphere.parquet')

# Align timezones
era5.index = pd.to_datetime(era5.index)
if era5.index.tz is None:
    era5.index = era5.index.tz_localize('UTC')
else:
    era5.index = era5.index.tz_convert('UTC')

ssw_dates = pd.to_datetime(ssw_cat.index).tz_convert('UTC') if ssw_cat.index.tz else pd.to_datetime(ssw_cat.index).tz_localize('UTC')

# Focus on DJF SSW events within ERA5 extended range (1998-2019)
mask = (ssw_dates >= era5.index.min()) & (ssw_dates <= era5.index.max())
ssw_in_range = ssw_dates[mask]
print(f"SSW events in ERA5 extended range: {len(ssw_in_range)}")

# Window: days -5 to +30 post-onset (core response window)
WINDOW_PRE = -5
WINDOW_POST = 30

def get_window_data(dates, df, pre=WINDOW_PRE, post=WINDOW_POST):
    """Get all data rows within window around event dates."""
    rows = []
    for d in dates:
        start = d + pd.Timedelta(days=pre)
        end = d + pd.Timedelta(days=post)
        chunk = df.loc[start:end].copy()
        chunk['event_date'] = d
        chunk['lag'] = (chunk.index - d).days
        rows.append(chunk)
    if rows:
        return pd.concat(rows)
    return pd.DataFrame()

def get_control_data(dates, df, pre=WINDOW_PRE, post=WINDOW_POST):
    """DOY-matched control from non-SSW winters."""
    rows = []
    for d in dates:
        doy_center = d.dayofyear
        for yr in df.index.year.unique():
            ctrl_date = pd.Timestamp(year=yr, month=d.month, day=d.day, tz='UTC')
            if ctrl_date == d:
                continue
            # Skip if this year had an SSW within ±60 days
            skip = False
            for sd in dates:
                if abs((ctrl_date - sd).days) < 60:
                    skip = True
                    break
            if skip:
                continue
            start = ctrl_date + pd.Timedelta(days=pre)
            end = ctrl_date + pd.Timedelta(days=post)
            chunk = df.loc[start:end].copy()
            if len(chunk) > 20:
                chunk['event_date'] = ctrl_date
                chunk['lag'] = (chunk.index - ctrl_date).days
                rows.append(chunk)
    if rows:
        return pd.concat(rows)
    return pd.DataFrame()

ssw_data = get_window_data(ssw_in_range, era5)
ctrl_data = get_control_data(ssw_in_range, era5)

print(f"SSW window days: {len(ssw_data)}")
print(f"Control window days: {len(ctrl_data)}")

# Define hazard-relevant metrics
results = {}

# 1. COLD-WAVE SEVERITY: Temperature deficit
ssw_temp = ssw_data['t2m_K'].values
ctrl_temp = ctrl_data['t2m_K'].values
temp_diff = ssw_temp.mean() - ctrl_temp.mean()
t_stat, p_val = stats.ttest_ind(ssw_temp, ctrl_temp, equal_var=False)
# Cold-wave threshold: days below 5th percentile of climatology
p5_temp = era5.loc[era5.index.month.isin([12,1,2,3]), 't2m_K'].quantile(0.05)
ssw_cold_frac = (ssw_temp < p5_temp).mean()
ctrl_cold_frac = (ctrl_temp < p5_temp).mean()
cold_rr = ssw_cold_frac / max(ctrl_cold_frac, 0.001)
results['cold_wave'] = {
    'hazard': 'Cold-wave mortality',
    'metric': 'Mean T2m anomaly (K)',
    'ssw_mean': round(float(ssw_temp.mean()), 2),
    'ctrl_mean': round(float(ctrl_temp.mean()), 2),
    'difference_K': round(float(temp_diff), 2),
    'p_value': round(float(p_val), 6),
    'extreme_cold_fraction_ssw': round(float(ssw_cold_frac), 4),
    'extreme_cold_fraction_ctrl': round(float(ctrl_cold_frac), 4),
    'cold_extreme_RR': round(float(cold_rr), 2),
    'direction': 'INCREASED' if temp_diff < 0 else 'decreased'
}

# 2. TOTAL PRECIPITATION: Drought/water resource risk
ssw_precip = ssw_data['tp_mm'].values
ctrl_precip = ctrl_data['tp_mm'].values
precip_ratio = ssw_precip.mean() / max(ctrl_precip.mean(), 0.001)
t_p, p_p = stats.ttest_ind(ssw_precip, ctrl_precip, equal_var=False)
results['precipitation'] = {
    'hazard': 'Drought / water resources',
    'metric': 'Daily precipitation (mm)',
    'ssw_mean': round(float(ssw_precip.mean()), 3),
    'ctrl_mean': round(float(ctrl_precip.mean()), 3),
    'ratio': round(float(precip_ratio), 3),
    'p_value': round(float(p_p), 6),
    'direction': 'DECREASED' if precip_ratio < 1 else 'increased'
}

# 3. SNOWFALL vs RAINFALL PARTITIONING: Rain-on-snow flood/landslide risk
ssw_sf = ssw_data['sf_mm'].values
ssw_rain = ssw_data['tp_mm'].values - ssw_data['sf_mm'].values
ctrl_sf = ctrl_data['sf_mm'].values
ctrl_rain = ctrl_data['tp_mm'].values - ctrl_data['sf_mm'].values
ssw_rain_frac = np.clip(ssw_rain, 0, None).sum() / max(ssw_data['tp_mm'].sum(), 0.001)
ctrl_rain_frac = np.clip(ctrl_rain, 0, None).sum() / max(ctrl_data['tp_mm'].sum(), 0.001)
results['rain_on_snow'] = {
    'hazard': 'Rain-on-snow flooding / landslides',
    'metric': 'Rainfall fraction of total precipitation',
    'ssw_rain_fraction': round(float(ssw_rain_frac), 4),
    'ctrl_rain_fraction': round(float(ctrl_rain_frac), 4),
    'ratio': round(float(ssw_rain_frac / max(ctrl_rain_frac, 0.001)), 3),
    'direction': 'DECREASED' if ssw_rain_frac < ctrl_rain_frac else 'increased'
}

# 4. WIND SPEED: Infrastructure/transport disruption
ssw_wind = ssw_data['wind_speed'].values
ctrl_wind = ctrl_data['wind_speed'].values
t_w, p_w = stats.ttest_ind(ssw_wind, ctrl_wind, equal_var=False)
# Extreme wind: days above 95th percentile
p95_wind = era5.loc[era5.index.month.isin([12,1,2,3]), 'wind_speed'].quantile(0.95)
ssw_extreme_wind = (ssw_wind > p95_wind).mean()
ctrl_extreme_wind = (ctrl_wind > p95_wind).mean()
results['wind_hazard'] = {
    'hazard': 'Wind-driven infrastructure damage',
    'metric': 'Mean wind speed (m/s)',
    'ssw_mean': round(float(ssw_wind.mean()), 3),
    'ctrl_mean': round(float(ctrl_wind.mean()), 3),
    'change_pct': round(float((ssw_wind.mean()/ctrl_wind.mean() - 1) * 100), 1),
    'p_value': round(float(p_w), 6),
    'extreme_wind_frac_ssw': round(float(ssw_extreme_wind), 4),
    'extreme_wind_frac_ctrl': round(float(ctrl_extreme_wind), 4),
    'direction': 'Variable (loading vs calm)'
}

# 5. SNOW DEPTH ANOMALY: Water resource / spring flood risk
ssw_sd = ssw_data['sd_m'].values
ctrl_sd = ctrl_data['sd_m'].values
t_sd, p_sd = stats.ttest_ind(ssw_sd, ctrl_sd, equal_var=False)
results['snow_loading'] = {
    'hazard': 'Spring snowmelt flooding',
    'metric': 'Snow depth (m)',
    'ssw_mean': round(float(ssw_sd.mean()), 4),
    'ctrl_mean': round(float(ctrl_sd.mean()), 4),
    'change_pct': round(float((ssw_sd.mean()/max(ctrl_sd.mean(), 0.001) - 1)*100), 1),
    'p_value': round(float(p_sd), 6),
    'direction': 'INCREASED' if ssw_sd.mean() > ctrl_sd.mean() else 'decreased'
}

# 6. TEMPERATURE INVERSION PROXY: Air quality risk
# Use: low wind + cold temperature = stagnation event
ssw_stag = ((ssw_data['wind_speed'] < era5['wind_speed'].quantile(0.25)) & 
            (ssw_data['t2m_K'] < era5['t2m_K'].quantile(0.25))).mean()
ctrl_stag = ((ctrl_data['wind_speed'] < era5['wind_speed'].quantile(0.25)) & 
             (ctrl_data['t2m_K'] < era5['t2m_K'].quantile(0.25))).mean()
results['air_quality'] = {
    'hazard': 'Air quality (cold stagnation events)',
    'metric': 'Fraction of cold+calm days',
    'ssw_fraction': round(float(ssw_stag), 4),
    'ctrl_fraction': round(float(ctrl_stag), 4),
    'ratio': round(float(ssw_stag / max(ctrl_stag, 0.001)), 2),
    'direction': 'INCREASED' if ssw_stag > ctrl_stag else 'decreased'
}

# Event-level analysis: per-SSW hazard vector
event_vectors = []
for d in ssw_in_range:
    start = d + pd.Timedelta(days=WINDOW_PRE)
    end = d + pd.Timedelta(days=WINDOW_POST)
    chunk = era5.loc[start:end]
    if len(chunk) < 20:
        continue
    
    # Get DOY-matched climatology
    doys = chunk['doy'].values if 'doy' in chunk.columns else chunk.index.dayofyear.values
    clim_mask = era5.index.month.isin([12,1,2,3])
    clim = era5[clim_mask]
    
    event_vectors.append({
        'event_date': str(d.date()),
        'temp_anom_K': round(float(chunk['t2m_K'].mean() - ctrl_temp.mean()), 2),
        'precip_ratio': round(float(chunk['tp_mm'].mean() / max(ctrl_precip.mean(), 0.001)), 3),
        'wind_change_pct': round(float((chunk['wind_speed'].mean()/max(ctrl_wind.mean(),0.001)-1)*100), 1),
        'snow_depth_change_pct': round(float((chunk['sd_m'].mean()/max(ctrl_sd.mean(),0.001)-1)*100), 1),
    })

# Consistency across events
n_cold = sum(1 for e in event_vectors if e['temp_anom_K'] < 0)
n_dry = sum(1 for e in event_vectors if e['precip_ratio'] < 1)
n_events = len(event_vectors)

# Multi-hazard simultaneity: how many hazard domains shift per event?
simultaneous_shifts = []
for e in event_vectors:
    n_hazards = 0
    if e['temp_anom_K'] < -0.5:
        n_hazards += 1  # cold wave
    if e['precip_ratio'] < 0.85:
        n_hazards += 1  # drought
    if e['wind_change_pct'] > 5:
        n_hazards += 1  # wind
    if e['snow_depth_change_pct'] > 2:
        n_hazards += 1  # snow loading
    simultaneous_shifts.append(n_hazards)

results['multi_hazard_summary'] = {
    'n_events': n_events,
    'n_colder_than_control': n_cold,
    'n_drier_than_control': n_dry,
    'cold_consistency': f"{n_cold}/{n_events}",
    'dry_consistency': f"{n_dry}/{n_events}",
    'mean_simultaneous_hazard_shifts': round(float(np.mean(simultaneous_shifts)), 2),
    'fraction_multi_hazard': round(float(np.mean([s >= 2 for s in simultaneous_shifts])), 3),
    'event_details': event_vectors
}

# Z500 analysis from NCEP (hemispheric blocking indicator)
ncep.index = pd.to_datetime(ncep.index)
if ncep.index.tz is None:
    ncep.index = ncep.index.tz_localize('UTC')

ssw_z500 = get_window_data(ssw_in_range, ncep)
ctrl_z500 = get_control_data(ssw_in_range, ncep)

if len(ssw_z500) > 0 and len(ctrl_z500) > 0:
    z500_ssw = ssw_z500['hgt_500hPa_m'].values
    z500_ctrl = ctrl_z500['hgt_500hPa_m'].values
    t_z, p_z = stats.ttest_ind(z500_ssw, z500_ctrl, equal_var=False)
    results['z500_blocking'] = {
        'hazard': 'Hemispheric blocking (multi-hazard driver)',
        'ssw_mean_gpm': round(float(z500_ssw.mean()), 1),
        'ctrl_mean_gpm': round(float(z500_ctrl.mean()), 1),
        'difference_gpm': round(float(z500_ssw.mean() - z500_ctrl.mean()), 1),
        'p_value': round(float(p_z), 6),
    }

# Print summary
print("\n" + "="*70)
print("CROSS-HAZARD REGIME-REDISTRIBUTION ANALYSIS")
print("="*70)
for key, val in results.items():
    if key in ('multi_hazard_summary', 'event_details'):
        continue
    print(f"\n{val.get('hazard', key)}:")
    for k, v in val.items():
        if k != 'hazard':
            print(f"  {k}: {v}")

mh = results['multi_hazard_summary']
print(f"\n{'='*70}")
print(f"MULTI-HAZARD SIMULTANEITY:")
print(f"  Events analysed: {mh['n_events']}")
print(f"  Colder than control: {mh['cold_consistency']}")
print(f"  Drier than control: {mh['dry_consistency']}")
print(f"  Mean simultaneous hazard shifts per event: {mh['mean_simultaneous_hazard_shifts']}")
print(f"  Fraction with ≥2 simultaneous shifts: {mh['fraction_multi_hazard']}")

# Save
out = ROOT / 'data/results/51_cross_hazard_regime.json'
with open(out, 'w') as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nSaved to {out}")
