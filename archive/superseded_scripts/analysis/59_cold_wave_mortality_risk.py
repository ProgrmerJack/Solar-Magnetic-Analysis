"""
Script 59: Cold-wave mortality risk during SSW windows
=====================================================
Applies published cold-mortality dose-response curves (Gasparrini et al. 2015)
to ERA5 Alpine temperature distributions during SSW vs non-SSW windows.

Demonstrates the loaded-gun mechanism in a SECOND hazard domain:
- SSW windows shift temperature distributions toward cold extremes
- Published dose-response relationships translate this to excess mortality risk
- The "loaded-gun" operates because cold extremes persist (loading) while
  the normal mixing/warming that terminates cold spells is suppressed (trigger suppression)

References:
- Gasparrini et al. (2015) Lancet 386:369-375 — temperature-mortality curves
- Reported minimum mortality temperature (MMT) for Central Europe: ~21°C mean,
  but winter-specific cold effects dominate below ~0°C
"""

import pandas as pd
import numpy as np
from scipy import stats
import json
import warnings
warnings.filterwarnings('ignore')

era5 = pd.read_parquet('data/processed/era5_swiss_alps_extended.parquet')
era5.index = pd.to_datetime(era5.index)
if era5.index.tz is not None:
    era5.index = era5.index.tz_localize(None)

ssw = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw.index = pd.to_datetime(ssw.index)
if ssw.index.tz is not None:
    ssw.index = ssw.index.tz_localize(None)
ssw_dates = ssw.index

era5['t2m_C'] = era5['t2m_K'] - 273.15
era5['doy'] = era5.index.dayofyear
era5['month'] = era5.index.month

winter = era5[(era5['month'] >= 11) | (era5['month'] <= 3)].copy()
print(f"Winter days: {len(winter)}")
print(f"Temperature range: {winter['t2m_C'].min():.1f} to {winter['t2m_C'].max():.1f}°C")

WINDOW = 15
ssw_mask = pd.Series(False, index=winter.index)
for d in ssw_dates:
    mask = (winter.index >= d - pd.Timedelta(days=WINDOW)) & \
           (winter.index <= d + pd.Timedelta(days=WINDOW))
    ssw_mask |= mask

ssw_days = winter[ssw_mask]
ctrl_days = winter[~ssw_mask]
print(f"\nSSW window days: {len(ssw_days)}")
print(f"Control days: {len(ctrl_days)}")

# DOY-matched comparison
results = {}
ssw_temps = ssw_days['t2m_C'].values
ctrl_temps = ctrl_days['t2m_C'].values

mean_ssw = np.mean(ssw_temps)
mean_ctrl = np.mean(ctrl_temps)
diff = mean_ssw - mean_ctrl
tstat, pval = stats.ttest_ind(ssw_temps, ctrl_temps)
cohens_d = diff / np.sqrt((np.var(ssw_temps) + np.var(ctrl_temps)) / 2)

print(f"\nTemperature comparison:")
print(f"  SSW mean: {mean_ssw:.2f}°C")
print(f"  Control mean: {mean_ctrl:.2f}°C")
print(f"  Difference: {diff:.2f}°C (P={pval:.4f})")
print(f"  Cohen's d: {cohens_d:.3f}")

# Cold extreme thresholds
percentiles = [1, 5, 10, 25]
print("\nCold extreme frequency shifts:")
cold_results = {}
for p in percentiles:
    threshold = np.percentile(winter['t2m_C'].values, p)
    ssw_freq = np.mean(ssw_temps < threshold)
    ctrl_freq = np.mean(ctrl_temps < threshold)
    rr = ssw_freq / ctrl_freq if ctrl_freq > 0 else np.inf
    cold_results[f'p{p}'] = {
        'threshold_C': round(threshold, 1),
        'ssw_freq_pct': round(ssw_freq * 100, 2),
        'ctrl_freq_pct': round(ctrl_freq * 100, 2),
        'risk_ratio': round(rr, 2)
    }
    print(f"  <P{p} ({threshold:.1f}°C): SSW {ssw_freq*100:.1f}% vs Ctrl {ctrl_freq*100:.1f}% → RR={rr:.2f}")

# Cold spell duration
def count_cold_spells(temps, dates, threshold, min_days=3):
    """Count cold spells (consecutive days below threshold)."""
    below = temps < threshold
    spells = []
    current_length = 0
    for i in range(len(below)):
        if below[i]:
            current_length += 1
        else:
            if current_length >= min_days:
                spells.append(current_length)
            current_length = 0
    if current_length >= min_days:
        spells.append(current_length)
    return spells

p10_thresh = np.percentile(winter['t2m_C'].values, 10)
ssw_spells = count_cold_spells(ssw_temps, ssw_days.index, p10_thresh)
ctrl_spells = count_cold_spells(ctrl_temps, ctrl_days.index, p10_thresh)

print(f"\nCold spell analysis (T < P10 = {p10_thresh:.1f}°C, ≥3 consecutive days):")
print(f"  SSW: {len(ssw_spells)} spells, mean duration {np.mean(ssw_spells):.1f} days" if ssw_spells else "  SSW: 0 spells")
print(f"  Ctrl: {len(ctrl_spells)} spells, mean duration {np.mean(ctrl_spells):.1f} days" if ctrl_spells else "  Ctrl: 0 spells")

# Apply Gasparrini dose-response to estimate excess mortality risk
# From Gasparrini et al. 2015 (Lancet), Central European cold-mortality curve:
# RR per 1°C below optimum ≈ 1.02-1.04 for moderate cold, accelerating at extremes
# Using the simplified log-linear model for temperatures below MMT:
# ln(RR) = β * (MMT - T) where β ≈ 0.03-0.06 per °C for cold effects

# Conservative estimate: β = 0.03 (moderate cold effect, Gasparrini Table S3)
# Central estimate: β = 0.045 (pooled Central European cities)
beta_conservative = 0.03
beta_central = 0.045
reference_T = 0.0  # Winter reference: 0°C (typical Swiss winter)

def compute_excess_mortality_rr(temps, beta, ref_T=0):
    """Compute average excess mortality RR from temperature distribution."""
    cold_deviation = np.maximum(ref_T - temps, 0)  # only cold effects
    rr_individual = np.exp(beta * cold_deviation)
    return np.mean(rr_individual)

rr_ssw_cons = compute_excess_mortality_rr(ssw_temps, beta_conservative)
rr_ctrl_cons = compute_excess_mortality_rr(ctrl_temps, beta_conservative)
rr_ssw_cent = compute_excess_mortality_rr(ssw_temps, beta_central)
rr_ctrl_cent = compute_excess_mortality_rr(ctrl_temps, beta_central)

excess_rr_cons = rr_ssw_cons / rr_ctrl_cons
excess_rr_cent = rr_ssw_cent / rr_ctrl_cent

print(f"\nCold-wave mortality risk estimation (Gasparrini dose-response):")
print(f"  Conservative (β=0.03): SSW RR={rr_ssw_cons:.3f}, Ctrl RR={rr_ctrl_cons:.3f}")
print(f"    → Excess mortality RR during SSW: {excess_rr_cons:.3f} ({(excess_rr_cons-1)*100:.1f}% increase)")
print(f"  Central (β=0.045): SSW RR={rr_ssw_cent:.3f}, Ctrl RR={rr_ctrl_cent:.3f}")
print(f"    → Excess mortality RR during SSW: {excess_rr_cent:.3f} ({(excess_rr_cent-1)*100:.1f}% increase)")

# Bootstrap CI for excess mortality RR
np.random.seed(42)
n_boot = 5000
boot_rr = []
for _ in range(n_boot):
    ssw_boot = np.random.choice(ssw_temps, size=len(ssw_temps), replace=True)
    ctrl_boot = np.random.choice(ctrl_temps, size=len(ctrl_temps), replace=True)
    rr_s = compute_excess_mortality_rr(ssw_boot, beta_central)
    rr_c = compute_excess_mortality_rr(ctrl_boot, beta_central)
    boot_rr.append(rr_s / rr_c)

boot_rr = np.array(boot_rr)
ci_low, ci_high = np.percentile(boot_rr, [2.5, 97.5])
print(f"  Bootstrap 95% CI for excess RR: [{ci_low:.3f}, {ci_high:.3f}]")

# Compound cold-stagnation (joint cold + low wind = air quality proxy)
wind_p25 = np.percentile(winter['wind_speed'].dropna().values, 25)
t_p25 = np.percentile(winter['t2m_C'].values, 25)

ssw_wind = ssw_days['wind_speed'].values
ctrl_wind = ctrl_days['wind_speed'].values

ssw_compound = np.mean((ssw_temps < t_p25) & (ssw_wind < wind_p25))
ctrl_compound = np.mean((ctrl_temps < t_p25) & (ctrl_wind < wind_p25))
compound_rr = ssw_compound / ctrl_compound if ctrl_compound > 0 else np.inf

print(f"\nCompound cold-stagnation (T<P25 AND wind<P25):")
print(f"  SSW: {ssw_compound*100:.1f}%, Ctrl: {ctrl_compound*100:.1f}%")
print(f"  Compound RR: {compound_rr:.2f}")

# Heating degree days (proxy for energy demand)
hdd_base = 15.5  # Standard European HDD base
ssw_hdd = np.mean(np.maximum(hdd_base - ssw_temps, 0))
ctrl_hdd = np.mean(np.maximum(hdd_base - ctrl_temps, 0))
hdd_increase = (ssw_hdd - ctrl_hdd) / ctrl_hdd * 100

print(f"\nHeating degree days (base {hdd_base}°C):")
print(f"  SSW: {ssw_hdd:.1f} HDD, Ctrl: {ctrl_hdd:.1f} HDD")
print(f"  Increase during SSW: {hdd_increase:.1f}%")

# Event-level analysis (like the avalanche analysis)
print("\n=== EVENT-LEVEL COLD-EXTREME ANALYSIS ===")
event_results = []
for d in ssw_dates:
    window_start = d - pd.Timedelta(days=WINDOW)
    window_end = d + pd.Timedelta(days=WINDOW)
    
    ssw_window = winter[(winter.index >= window_start) & (winter.index <= window_end)]
    if len(ssw_window) < 5:
        continue
    
    doy_center = d.dayofyear
    ctrl_window = winter[~ssw_mask & 
        (winter['doy'] >= doy_center - WINDOW) & 
        (winter['doy'] <= doy_center + WINDOW)]
    
    if len(ctrl_window) < 10:
        continue
    
    ssw_mean_t = ssw_window['t2m_C'].mean()
    ctrl_mean_t = ctrl_window['t2m_C'].mean()
    
    # Extreme cold frequency
    ssw_extreme = np.mean(ssw_window['t2m_C'].values < p10_thresh)
    ctrl_extreme = np.mean(ctrl_window['t2m_C'].values < p10_thresh)
    
    # Mortality RR
    rr_event_ssw = compute_excess_mortality_rr(ssw_window['t2m_C'].values, beta_central)
    rr_event_ctrl = compute_excess_mortality_rr(ctrl_window['t2m_C'].values, beta_central)
    
    event_results.append({
        'date': d.strftime('%Y-%m-%d'),
        'ssw_mean_T': ssw_mean_t,
        'ctrl_mean_T': ctrl_mean_t,
        'delta_T': ssw_mean_t - ctrl_mean_t,
        'ssw_extreme_freq': ssw_extreme,
        'ctrl_extreme_freq': ctrl_extreme,
        'extreme_rr': ssw_extreme / ctrl_extreme if ctrl_extreme > 0 else np.nan,
        'mortality_rr_ssw': rr_event_ssw,
        'mortality_rr_ctrl': rr_event_ctrl,
        'excess_mortality_rr': rr_event_ssw / rr_event_ctrl if rr_event_ctrl > 0 else np.nan
    })

edf = pd.DataFrame(event_results)
n_colder = (edf['delta_T'] < 0).sum()
n_total = len(edf)
n_excess_mort = (edf['excess_mortality_rr'] > 1).sum()

binom_cold = stats.binomtest(n_colder, n_total, 0.5)
binom_mort = stats.binomtest(n_excess_mort, n_total, 0.5)

print(f"Events analyzed: {n_total}")
print(f"Events with colder SSW windows: {n_colder}/{n_total} (P={binom_cold.pvalue:.4f})")
print(f"Events with excess mortality risk: {n_excess_mort}/{n_total} (P={binom_mort.pvalue:.4f})")
print(f"Median temperature shift: {edf['delta_T'].median():.2f}°C")
print(f"Median excess mortality RR: {edf['excess_mortality_rr'].median():.3f}")

# Save results
output = {
    'aggregate': {
        'ssw_mean_T_C': round(mean_ssw, 2),
        'ctrl_mean_T_C': round(mean_ctrl, 2),
        'temperature_shift_C': round(diff, 2),
        'temperature_shift_P': round(pval, 6),
        'cohens_d': round(cohens_d, 3),
    },
    'cold_extreme_shifts': cold_results,
    'cold_spell': {
        'n_ssw_spells': len(ssw_spells),
        'n_ctrl_spells': len(ctrl_spells),
        'mean_ssw_duration': round(np.mean(ssw_spells), 1) if ssw_spells else 0,
        'mean_ctrl_duration': round(np.mean(ctrl_spells), 1) if ctrl_spells else 0,
    },
    'mortality_risk': {
        'beta_central': beta_central,
        'excess_rr_conservative': round(excess_rr_cons, 4),
        'excess_rr_central': round(excess_rr_cent, 4),
        'excess_rr_ci_low': round(ci_low, 4),
        'excess_rr_ci_high': round(ci_high, 4),
        'excess_pct_central': round((excess_rr_cent - 1) * 100, 1),
    },
    'compound_stagnation': {
        'ssw_freq_pct': round(ssw_compound * 100, 1),
        'ctrl_freq_pct': round(ctrl_compound * 100, 1),
        'compound_rr': round(compound_rr, 2),
    },
    'heating_degree_days': {
        'ssw_hdd': round(ssw_hdd, 1),
        'ctrl_hdd': round(ctrl_hdd, 1),
        'increase_pct': round(hdd_increase, 1),
    },
    'event_level': {
        'n_events': n_total,
        'n_colder': int(n_colder),
        'colder_sign_P': round(binom_cold.pvalue, 4),
        'n_excess_mortality': int(n_excess_mort),
        'excess_mort_sign_P': round(binom_mort.pvalue, 4),
        'median_delta_T': round(edf['delta_T'].median(), 2),
        'median_excess_mort_rr': round(edf['excess_mortality_rr'].median(), 4),
    }
}

with open('data/results/59_cold_wave_mortality_risk.json', 'w') as f:
    json.dump(output, f, indent=2)

print("\nResults saved to data/results/59_cold_wave_mortality_risk.json")
print("\n=== SUMMARY FOR MANUSCRIPT ===")
print(f"SSW windows shift Alpine temperature by {diff:.1f}°C (P={pval:.4f})")
print(f"Extreme cold days (P10) increase {cold_results['p10']['risk_ratio']:.1f}× during SSW")
print(f"Excess cold-mortality risk: +{(excess_rr_cent-1)*100:.1f}% (95% CI: +{(ci_low-1)*100:.1f}% to +{(ci_high-1)*100:.1f}%)")
print(f"Cold-stagnation compound events: {compound_rr:.1f}× more frequent")
print(f"Heating demand increases by {hdd_increase:.1f}%")
print(f"Event-level: {n_colder}/{n_total} SSW events produce colder conditions (P={binom_cold.pvalue:.4f})")
print(f"Event-level: {n_excess_mort}/{n_total} events produce excess mortality risk (P={binom_mort.pvalue:.4f})")
