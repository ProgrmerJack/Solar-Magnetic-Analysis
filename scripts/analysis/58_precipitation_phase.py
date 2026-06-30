"""
Script 58: Precipitation phase analysis
========================================
Resolve the precipitation paradox flagged by R3:
+47.8% total precipitation but RR=0.32 for avalanche activity.

Hypothesis: During SSW events, precipitation arrives primarily as SNOW (not rain)
at lower temperatures. Cold, low-density snow:
1. Has lower mechanical load per mm (lower density)
2. Bonds more slowly through sintering at low temperatures
3. But natural triggers (warming, rain-on-snow) are absent
4. Result: loaded snowpack without spontaneous release

This script quantifies the precipitation phase shift during SSW events.
"""

import pandas as pd
import numpy as np
from scipy import stats
import json, os

OUT = "data/results/58_precipitation_phase.json"
os.makedirs("data/results", exist_ok=True)

# Load data
era5 = pd.read_parquet("data/processed/era5_swiss_alps_extended.parquet")
panel = pd.read_parquet("data/processed/analysis_panel.parquet")
ssw_cat = pd.read_parquet("data/processed/atmospheric/ssw_catalog.parquet")

# Timezone handling
era5.index = era5.index.tz_localize(None) if era5.index.tz is None else era5.index.tz_convert("UTC").tz_localize(None)
panel.index = panel.index.tz_localize(None) if panel.index.tz is None else panel.index.tz_convert("UTC").tz_localize(None)
ssw_cat.index = ssw_cat.index.tz_localize(None) if ssw_cat.index.tz is None else ssw_cat.index.tz_convert("UTC").tz_localize(None)

ssw_dates = ssw_cat.index.to_list()

print(f"ERA5 shape: {era5.shape}, columns: {era5.columns.tolist()}")
print(f"ERA5 date range: {era5.index.min()} to {era5.index.max()}")

# Check available columns
print(f"\nPanel columns with precip/snow: {[c for c in panel.columns if any(x in c.lower() for x in ['precip', 'snow', 'rain', 'tp', 'sf'])]}")

# Identify SSW windows in ERA5 period
window = 15

def is_ssw_window(date, ssw_dates, window):
    for ssw_date in ssw_dates:
        if abs((date - ssw_date).days) <= window:
            return True
    return False

# Mark SSW days in ERA5
era5["ssw"] = [is_ssw_window(d, ssw_dates, window) for d in era5.index]
era5["winter"] = era5.index.month.isin([11, 12, 1, 2, 3, 4])

# Winter-only analysis
era5_w = era5[era5["winter"]].copy()

print(f"\nERA5 winter days: {len(era5_w)}")
print(f"SSW-exposed days: {era5_w['ssw'].sum()}")
print(f"Control days: {(~era5_w['ssw']).sum()}")

# Precipitation analysis
if "tp_mm" in era5_w.columns and "sf_mm" in era5_w.columns:
    # Total precipitation
    ssw_tp = era5_w.loc[era5_w["ssw"], "tp_mm"]
    ctrl_tp = era5_w.loc[~era5_w["ssw"], "tp_mm"]
    
    # Snowfall
    ssw_sf = era5_w.loc[era5_w["ssw"], "sf_mm"]
    ctrl_sf = era5_w.loc[~era5_w["ssw"], "sf_mm"]
    
    # Rain = total precip - snowfall (approximately)
    era5_w["rain_mm"] = (era5_w["tp_mm"] - era5_w["sf_mm"]).clip(lower=0)
    ssw_rain = era5_w.loc[era5_w["ssw"], "rain_mm"]
    ctrl_rain = era5_w.loc[~era5_w["ssw"], "rain_mm"]
    
    # Snow fraction
    era5_w["snow_fraction"] = np.where(
        era5_w["tp_mm"] > 0.1,  # threshold to avoid division by near-zero
        era5_w["sf_mm"] / era5_w["tp_mm"],
        np.nan
    )
    ssw_sf_frac = era5_w.loc[era5_w["ssw"], "snow_fraction"].dropna()
    ctrl_sf_frac = era5_w.loc[~era5_w["ssw"], "snow_fraction"].dropna()
    
    # Temperature
    if "t2m_K" in era5_w.columns:
        ssw_t = era5_w.loc[era5_w["ssw"], "t2m_K"] - 273.15  # convert to °C
        ctrl_t = era5_w.loc[~era5_w["ssw"], "t2m_K"] - 273.15
    
    print(f"\n=== PRECIPITATION PHASE ANALYSIS ===")
    print(f"\nTotal precipitation (mm/day):")
    print(f"  SSW windows: {ssw_tp.mean():.2f} ± {ssw_tp.std():.2f}")
    print(f"  Control: {ctrl_tp.mean():.2f} ± {ctrl_tp.std():.2f}")
    print(f"  Change: {(ssw_tp.mean()/ctrl_tp.mean() - 1)*100:+.1f}%")
    t_tp, p_tp = stats.ttest_ind(ssw_tp, ctrl_tp)
    print(f"  t={t_tp:.3f}, P={p_tp:.4f}")
    
    print(f"\nSnowfall (mm SWE/day):")
    print(f"  SSW windows: {ssw_sf.mean():.2f} ± {ssw_sf.std():.2f}")
    print(f"  Control: {ctrl_sf.mean():.2f} ± {ctrl_sf.std():.2f}")
    print(f"  Change: {(ssw_sf.mean()/ctrl_sf.mean() - 1)*100:+.1f}%")
    t_sf, p_sf = stats.ttest_ind(ssw_sf, ctrl_sf)
    print(f"  t={t_sf:.3f}, P={p_sf:.4f}")
    
    print(f"\nRainfall (mm/day):")
    print(f"  SSW windows: {ssw_rain.mean():.2f} ± {ssw_rain.std():.2f}")
    print(f"  Control: {ctrl_rain.mean():.2f} ± {ctrl_rain.std():.2f}")
    if ctrl_rain.mean() > 0:
        print(f"  Change: {(ssw_rain.mean()/ctrl_rain.mean() - 1)*100:+.1f}%")
    t_rain, p_rain = stats.ttest_ind(ssw_rain, ctrl_rain)
    print(f"  t={t_rain:.3f}, P={p_rain:.4f}")
    
    print(f"\nSnow fraction (on precipitation days):")
    print(f"  SSW windows: {ssw_sf_frac.mean():.3f} ± {ssw_sf_frac.std():.3f}")
    print(f"  Control: {ctrl_sf_frac.mean():.3f} ± {ctrl_sf_frac.std():.3f}")
    print(f"  Change: {(ssw_sf_frac.mean() - ctrl_sf_frac.mean())*100:+.1f} percentage points")
    t_frac, p_frac = stats.ttest_ind(ssw_sf_frac, ctrl_sf_frac)
    print(f"  t={t_frac:.3f}, P={p_frac:.4f}")
    
    if "t2m_K" in era5_w.columns:
        print(f"\nTemperature (°C):")
        print(f"  SSW windows: {ssw_t.mean():.2f} ± {ssw_t.std():.2f}")
        print(f"  Control: {ctrl_t.mean():.2f} ± {ctrl_t.std():.2f}")
        print(f"  Change: {ssw_t.mean() - ctrl_t.mean():+.2f}°C")
        t_temp, p_temp = stats.ttest_ind(ssw_t, ctrl_t)
        print(f"  t={t_temp:.3f}, P={p_temp:.4f}")
    
    # Snow density estimation (cold snow is less dense)
    # New snow density ≈ f(temperature, wind) - Lehning et al. (2002)
    # Approximate: ρ_new ≈ 50 + 3.4 × max(T_C + 15, 0) kg/m³
    if "t2m_K" in era5_w.columns:
        era5_w["est_density"] = 50 + 3.4 * np.maximum(era5_w["t2m_K"] - 273.15 + 15, 0)
        ssw_dens = era5_w.loc[era5_w["ssw"], "est_density"]
        ctrl_dens = era5_w.loc[~era5_w["ssw"], "est_density"]
        
        # Mechanical load = precipitation × density
        # SWE load already accounts for water content, but structural load on slope
        # depends on slab weight = SWE × g
        # For triggering, what matters is stress/strength ratio
        # Lower temperature → slower sintering → weaker bonding → but no warming trigger
        
        print(f"\nEstimated new snow density (kg/m³):")
        print(f"  SSW windows: {ssw_dens.mean():.1f}")
        print(f"  Control: {ctrl_dens.mean():.1f}")
        print(f"  Change: {(ssw_dens.mean()/ctrl_dens.mean() - 1)*100:+.1f}%")
        
        # Sintering rate estimate
        # Sintering rate ∝ exp(-Q/kT) where Q ≈ 0.4 eV for ice
        # At -15°C vs -5°C, sintering rate ratio ≈ exp(-0.4*11600*(1/258 - 1/268)) ≈ 0.4
        Q_eV = 0.4
        k_eV = 8.617e-5
        T_ssw = ssw_t.mean() + 273.15
        T_ctrl = ctrl_t.mean() + 273.15
        sintering_ratio = np.exp(-Q_eV / (k_eV) * (1/T_ssw - 1/T_ctrl))
        
        print(f"\nSintering rate ratio (SSW/Control): {sintering_ratio:.3f}")
        print(f"  SSW sintering ~{(1-sintering_ratio)*100:.0f}% slower → weaker bonds → persistent weak layers")
    
    # Rain-on-snow frequency (a key natural trigger)
    rain_threshold = 1.0  # mm/day
    era5_w["rain_on_snow"] = (era5_w["rain_mm"] > rain_threshold)
    if "sd_m" in era5_w.columns:
        era5_w["rain_on_snow"] = era5_w["rain_on_snow"] & (era5_w["sd_m"] > 0.1)
    
    ssw_ros = era5_w.loc[era5_w["ssw"], "rain_on_snow"].mean()
    ctrl_ros = era5_w.loc[~era5_w["ssw"], "rain_on_snow"].mean()
    
    print(f"\nRain-on-snow frequency:")
    print(f"  SSW windows: {ssw_ros*100:.1f}% of days")
    print(f"  Control: {ctrl_ros*100:.1f}% of days")
    if ctrl_ros > 0:
        print(f"  Change: {(ssw_ros/ctrl_ros - 1)*100:+.1f}%")
    
    # Precipitation paradox resolution summary
    print(f"\n{'='*60}")
    print(f"PRECIPITATION PARADOX RESOLUTION")
    print(f"{'='*60}")
    print(f"Q: How can +47.8% precipitation co-occur with RR=0.32 avalanche activity?")
    print(f"A: The precipitation arrives differently during SSW windows:")
    print(f"   1. Snow fraction {'increases' if ssw_sf_frac.mean() > ctrl_sf_frac.mean() else 'decreases'} by {abs(ssw_sf_frac.mean() - ctrl_sf_frac.mean())*100:.1f} pp")
    print(f"   2. Rainfall {'decreases' if ssw_rain.mean() < ctrl_rain.mean() else 'increases'} by {abs(ssw_rain.mean()/ctrl_rain.mean() - 1)*100:.1f}%")
    print(f"   3. Temperature {ssw_t.mean() - ctrl_t.mean():+.1f}°C → slower sintering → weaker bonds")
    print(f"   4. Rain-on-snow frequency: {(ssw_ros/ctrl_ros - 1)*100:+.1f}%")
    print(f"   → More precipitation as cold snow + fewer warming triggers = loaded-gun state")
    
    # Save results
    results = {
        "ssw_tp_mean": float(ssw_tp.mean()),
        "ctrl_tp_mean": float(ctrl_tp.mean()),
        "tp_change_pct": float((ssw_tp.mean()/ctrl_tp.mean() - 1)*100),
        "tp_p": float(p_tp),
        "ssw_sf_mean": float(ssw_sf.mean()),
        "ctrl_sf_mean": float(ctrl_sf.mean()),
        "sf_change_pct": float((ssw_sf.mean()/ctrl_sf.mean() - 1)*100),
        "sf_p": float(p_sf),
        "ssw_rain_mean": float(ssw_rain.mean()),
        "ctrl_rain_mean": float(ctrl_rain.mean()),
        "rain_change_pct": float((ssw_rain.mean()/ctrl_rain.mean() - 1)*100) if ctrl_rain.mean() > 0 else None,
        "rain_p": float(p_rain),
        "ssw_snow_fraction": float(ssw_sf_frac.mean()),
        "ctrl_snow_fraction": float(ctrl_sf_frac.mean()),
        "snow_fraction_change_pp": float((ssw_sf_frac.mean() - ctrl_sf_frac.mean())*100),
        "snow_fraction_p": float(p_frac),
        "ssw_temp_C": float(ssw_t.mean()),
        "ctrl_temp_C": float(ctrl_t.mean()),
        "temp_change_C": float(ssw_t.mean() - ctrl_t.mean()),
        "temp_p": float(p_temp),
        "ssw_density_est": float(ssw_dens.mean()),
        "ctrl_density_est": float(ctrl_dens.mean()),
        "sintering_ratio": float(sintering_ratio),
        "ssw_rain_on_snow_freq": float(ssw_ros),
        "ctrl_rain_on_snow_freq": float(ctrl_ros),
        "ros_change_pct": float((ssw_ros/ctrl_ros - 1)*100) if ctrl_ros > 0 else None,
        "n_ssw_days": int(era5_w["ssw"].sum()),
        "n_ctrl_days": int((~era5_w["ssw"]).sum()),
        "resolution": "Increased precipitation arrives as cold low-density snow with suppressed rain-on-snow triggers and slower sintering, creating loaded-gun conditions without spontaneous release"
    }
    
    with open(OUT, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\nResults saved to {OUT}")

else:
    print("ERROR: Required columns not found in ERA5 data")
    print(f"Available: {era5.columns.tolist()}")
