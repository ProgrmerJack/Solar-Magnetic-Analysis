"""
Sintering timescale calculation for the temperature anomaly.
Estimates how +0.81 K warming affects sintering rate using
the McClung & Schaerer (2006) temperature-dependent grain bonding model.

Key physics: Sintering rate depends on ice surface diffusion,
which follows Arrhenius kinetics: rate ∝ exp(-Ea / (k*T))
For ice, Ea ≈ 0.6 eV (activation energy for surface diffusion)
"""
import numpy as np

# Physical constants
kB = 8.617e-5  # eV/K (Boltzmann constant)
Ea = 0.6  # eV, activation energy for ice surface diffusion (Hobbs, 1974)

# Temperature range (typical deep-winter Swiss Alps)
T_base_C = -5.0  # baseline temperature (°C)
T_base_K = T_base_C + 273.15  # convert to Kelvin
delta_T = 0.81  # our observed SSW anomaly (K)
T_ssw_K = T_base_K + delta_T

# Arrhenius ratio: how much faster is sintering at T_ssw vs T_base?
rate_ratio = np.exp(-Ea/kB * (1/T_ssw_K - 1/T_base_K))

print("=== SINTERING TIMESCALE CALCULATION ===")
print(f"Baseline temperature: {T_base_C}°C ({T_base_K:.2f} K)")
print(f"SSW anomaly: +{delta_T} K")
print(f"SSW temperature: {T_base_C + delta_T}°C ({T_ssw_K:.2f} K)")
print(f"Activation energy (ice surface diffusion): {Ea} eV")
print(f"\nArrhenius rate enhancement: {rate_ratio:.3f}x")
print(f"  → Sintering proceeds {(rate_ratio-1)*100:.1f}% faster during SSW")
print(f"  → 15-day cumulative effect: {15*(rate_ratio-1)*100:.0f}% more sintering")

# Sensitivity to baseline temperature
print("\n=== SENSITIVITY TO BASELINE TEMPERATURE ===")
for T_base in [-10, -7, -5, -3, -1, 0]:
    T_K = T_base + 273.15
    T_ssw = T_K + delta_T
    ratio = np.exp(-Ea/kB * (1/T_ssw - 1/T_K))
    print(f"  T_base={T_base:+3d}°C: rate enhancement = {ratio:.3f}x ({(ratio-1)*100:.1f}%)")

# Sensitivity to anomaly magnitude
print("\n=== SENSITIVITY TO ANOMALY MAGNITUDE ===")
for dT in [0.45, 0.81, 0.92]:
    T_ssw = T_base_K + dT
    ratio = np.exp(-Ea/kB * (1/T_ssw - 1/T_base_K))
    label = ""
    if dT == 0.45: label = " (pre-SSW)"
    elif dT == 0.81: label = " (full window)"
    elif dT == 0.92: label = " (post-SSW peak)"
    print(f"  ΔT={dT:+.2f}K{label}: rate = {ratio:.3f}x ({(ratio-1)*100:.1f}%)")

# McClung bond-growth timescale
print("\n=== BOND GROWTH TIMESCALE (McClung & Schaerer 2006) ===")
# Bond radius grows as r ∝ t^(1/5) for viscous sintering (Hobbs & Mason 1964)
# Time to double bond radius: t2 = t1 * 2^5 = 32*t1
# At -5°C, characteristic sintering time for new snow ≈ 12-48 hours
# The rate enhancement means the same bond growth is reached faster
tau_base_hr = 24  # characteristic sintering time at baseline (hours, mid-range)
tau_ssw_hr = tau_base_hr / rate_ratio
print(f"Characteristic sintering time at {T_base_C}°C: ~{tau_base_hr}h")
print(f"Characteristic sintering time at {T_base_C+delta_T}°C: ~{tau_ssw_hr:.1f}h")
print(f"Time savings: {tau_base_hr - tau_ssw_hr:.1f}h per cycle")
print(f"\nOver 15 days with continuous warming:")
print(f"  Baseline: equivalent to {15} sintering cycles")
print(f"  SSW: equivalent to {15*rate_ratio:.1f} effective sintering cycles")
print(f"  Net advantage: {15*(rate_ratio-1):.1f} additional effective days of sintering")

# Clausius-Clapeyron: vapor pressure enhancement
print("\n=== VAPOR PRESSURE ENHANCEMENT (Clausius-Clapeyron) ===")
L_sub = 2.83e6  # J/kg, latent heat of sublimation
Rv = 461.5  # J/(kg*K), water vapor gas constant
# Saturation vapor pressure over ice: es(T) = es0 * exp(L/Rv * (1/T0 - 1/T))
es_ratio = np.exp(L_sub/Rv * (1/T_base_K - 1/T_ssw_K))
print(f"Vapor pressure enhancement: {es_ratio:.4f}x ({(es_ratio-1)*100:.2f}%)")
print(f"This drives faster vapor transport through snowpack pore space")
print(f"  → Enhanced hand-to-hand vapor transfer (Yosida, 1955)")
print(f"  → Accelerated grain rounding and equilibrium metamorphism")
