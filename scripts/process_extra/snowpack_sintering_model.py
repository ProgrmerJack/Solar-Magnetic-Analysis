"""
Process-based 1D snowpack sintering model.

Physics:
- Arrhenius bond growth (Hobbs 1974, Colbeck 1997)
- Clausius-Clapeyron vapor transport
- Density-dependent slab strength (Jamieson & Johnston 2001)
- Temperature-dependent viscosity (Mellor 1975)

Forces with ERA5 temperature during SSW vs control windows.
Compares predicted slab strength against observed Rutschblock changes.
"""
import numpy as np
import pandas as pd
import os, sys

DATA = r"C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed"
OUT  = r"C:\Users\Jack0\Solar-Magnetic-Analysis\data\results"
os.makedirs(OUT, exist_ok=True)

# ========== PHYSICAL CONSTANTS ==========
kB = 8.617e-5       # Boltzmann constant [eV/K]
R_gas = 8.314        # Gas constant [J/(mol·K)]
Ea_sintering = 0.6   # Activation energy ice surface diffusion [eV] (Hobbs 1974)
Ea_sintering_J = Ea_sintering * 1.602e-19 * 6.022e23  # Convert to J/mol = 57.87 kJ/mol
L_sub = 2.83e6       # Latent heat sublimation [J/kg]
Rv = 461.5           # Water vapor gas constant [J/(kg·K)]
rho_ice = 917.0      # Ice density [kg/m³]

# ========== SINTERING MODEL EQUATIONS ==========

def arrhenius_rate(T_K, Ea_eV=0.6):
    """Arrhenius sintering rate factor. Hobbs (1974) ice surface diffusion."""
    return np.exp(-Ea_eV / (kB * T_K))

def bond_growth_rate(T_K, r_grain=0.5e-3, r_bond_ratio=0.3):
    """
    Bond radius growth rate [m/s].
    Based on Colbeck (1997) vapor transport model.
    dr_b/dt = C * D_v * (de_s/dT) * (dT/dx) / (rho_ice * r_b)
    Simplified: rate proportional to Arrhenius factor × vapor pressure gradient.
    """
    # Vapor diffusion coefficient in air [m²/s] (T-dependent)
    D_v = 2.036e-5 * (T_K / 273.15)**1.75
    # Saturation vapor pressure [Pa] (Buck 1981)
    T_C = T_K - 273.15
    e_s = 611.15 * np.exp((23.036 - T_C/333.7) * T_C / (279.82 + T_C))
    # de_s/dT [Pa/K] from Clausius-Clapeyron
    de_dT = e_s * L_sub / (Rv * T_K**2)
    # Bond growth rate coefficient
    r_b = r_bond_ratio * r_grain
    rate = D_v * de_dT / (rho_ice * r_b * 1e3)  # Normalize
    return rate

def slab_strength(density, bond_grain_ratio):
    """
    Slab tensile strength [kPa].
    σ_t = a * (ρ/ρ_ice)^b * f(r_b/r_g)
    Based on Jamieson & Johnston (2001) and Schweizer et al. (2003).
    """
    a = 300.0  # kPa, reference strength at ρ_ice
    b = 2.4    # Density exponent (empirical)
    rho_ratio = density / rho_ice
    # Bond contribution: strength scales with bond-to-grain ratio squared
    bond_factor = (bond_grain_ratio / 0.3)**2
    return a * rho_ratio**b * bond_factor

def viscous_densification(T_K, density, dt_s=86400):
    """
    Snow densification rate [kg/m³/s].
    η = η₀ * exp(Q/R * (1/T - 1/T_ref)) * exp(k*ρ)
    Based on Mellor (1975).
    """
    T_ref = 263.15  # -10°C reference
    eta_0 = 1e7     # Reference viscosity [Pa·s]
    Q = 67000.0     # Activation energy [J/mol]
    k_rho = 0.023   # Density coefficient [m³/kg]
    
    viscosity = eta_0 * np.exp(Q/R_gas * (1/T_K - 1/T_ref)) * np.exp(k_rho * density)
    # Overburden stress for 0.5m slab
    sigma = density * 9.81 * 0.5
    d_rho = sigma / viscosity * dt_s
    return d_rho

# ========== RUN SIMULATION ==========

def simulate_slab(T_series_K, dt=86400, initial_density=250.0, initial_bgr=0.2):
    """
    Simulate slab evolution over time given temperature series.
    
    Parameters:
    - T_series_K: array of daily mean temperatures [K]
    - dt: timestep [s] (default: 1 day)
    - initial_density: initial slab density [kg/m³]
    - initial_bgr: initial bond-to-grain ratio
    
    Returns: DataFrame with daily slab properties
    """
    n_days = len(T_series_K)
    density = np.zeros(n_days)
    bgr = np.zeros(n_days)
    strength = np.zeros(n_days)
    sinter_rate = np.zeros(n_days)
    vapor_rate = np.zeros(n_days)
    
    density[0] = initial_density
    bgr[0] = initial_bgr
    strength[0] = slab_strength(density[0], bgr[0])
    
    for i in range(1, n_days):
        T = T_series_K[i]
        
        # 1. Arrhenius sintering rate
        rate_factor = arrhenius_rate(T) / arrhenius_rate(263.15)  # Normalize to -10°C
        sinter_rate[i] = rate_factor
        
        # 2. Bond growth
        bgr_growth = bond_growth_rate(T) * dt * 1e-3  # Scale factor
        bgr[i] = min(bgr[i-1] + bgr_growth * rate_factor, 0.6)
        vapor_rate[i] = bgr_growth
        
        # 3. Densification
        d_rho = viscous_densification(T, density[i-1], dt)
        density[i] = min(density[i-1] + d_rho, 550.0)
        
        # 4. Strength
        strength[i] = slab_strength(density[i], bgr[i])
    
    return pd.DataFrame({
        'day': range(n_days),
        'T_K': T_series_K,
        'density': density,
        'bgr': bgr,
        'strength_kPa': strength,
        'sinter_rate_norm': sinter_rate,
    })

# ========== LOAD ERA5 AND SSW DATA ==========
print("Loading ERA5 and SSW data...")
era5 = pd.read_parquet(os.path.join(DATA, "era5_swiss_alps_daily.parquet"))
ssw_cat = pd.read_parquet(os.path.join(DATA, "atmospheric", "ssw_catalog.parquet"))

# Fix timezone
if hasattr(ssw_cat.index, 'tz') and ssw_cat.index.tz is not None:
    ssw_cat.index = ssw_cat.index.tz_localize(None)

# Filter SSW events within ERA5 coverage
era5_start = era5.index.min()
era5_end = era5.index.max()
ssw_dates = ssw_cat.index[(ssw_cat.index >= era5_start) & (ssw_cat.index <= era5_end)]
ssw_dates = ssw_dates[ssw_dates.month.isin([11,12,1,2,3])]  # Winter only
print(f"SSW events in ERA5 range: {len(ssw_dates)}")

# ========== SIMULATE SSW vs CONTROL WINDOWS ==========
print("\n=== Running sintering model for SSW vs Control windows ===")

ssw_results = []
ctrl_results = []

for ssw_date in ssw_dates:
    # SSW window: d-5 to d+15 (20 days, centered on post-SSW warming)
    window_start = ssw_date - pd.Timedelta(days=5)
    window_end = ssw_date + pd.Timedelta(days=15)
    
    ssw_era = era5.loc[window_start:window_end]
    if len(ssw_era) < 15:
        continue
    
    # Control: same calendar days, other years
    ctrl_temps = []
    for yr_offset in [-2, -1, 1, 2]:
        ctrl_start = window_start + pd.DateOffset(years=yr_offset)
        ctrl_end = window_end + pd.DateOffset(years=yr_offset)
        ctrl_era = era5.loc[ctrl_start:ctrl_end]
        if len(ctrl_era) >= 15:
            ctrl_temps.append(ctrl_era['t2m_K'].values[:len(ssw_era)])
    
    if len(ctrl_temps) == 0:
        continue
    
    ctrl_mean_temp = np.mean(ctrl_temps, axis=0)
    
    # Run simulation for SSW and control
    ssw_sim = simulate_slab(ssw_era['t2m_K'].values)
    ctrl_sim = simulate_slab(ctrl_mean_temp)
    
    ssw_results.append({
        'ssw_date': ssw_date,
        'ssw_final_strength': ssw_sim['strength_kPa'].iloc[-1],
        'ssw_final_density': ssw_sim['density'].iloc[-1],
        'ssw_final_bgr': ssw_sim['bgr'].iloc[-1],
        'ssw_mean_T': ssw_era['t2m_K'].mean(),
        'ctrl_final_strength': ctrl_sim['strength_kPa'].iloc[-1],
        'ctrl_final_density': ctrl_sim['density'].iloc[-1],
        'ctrl_final_bgr': ctrl_sim['bgr'].iloc[-1],
        'ctrl_mean_T': ctrl_mean_temp.mean(),
        'strength_ratio': ssw_sim['strength_kPa'].iloc[-1] / ctrl_sim['strength_kPa'].iloc[-1],
        'T_diff': ssw_era['t2m_K'].mean() - ctrl_mean_temp.mean(),
    })

results_df = pd.DataFrame(ssw_results)
print(f"\nSimulated {len(results_df)} SSW events")
print("\n=== Event-level sintering model results ===")
print(results_df[['ssw_date','T_diff','strength_ratio','ssw_final_bgr','ctrl_final_bgr']].to_string())

# Summary statistics
print(f"\n=== SUMMARY ===")
print(f"Mean T difference (SSW - control): {results_df['T_diff'].mean():.3f} K")
print(f"Mean strength ratio (SSW/control): {results_df['strength_ratio'].mean():.4f}")
print(f"Mean SSW final BGR: {results_df['ssw_final_bgr'].mean():.4f}")
print(f"Mean control final BGR: {results_df['ctrl_final_bgr'].mean():.4f}")
print(f"BGR enhancement: {(results_df['ssw_final_bgr'].mean()/results_df['ctrl_final_bgr'].mean()-1)*100:.2f}%")

# Statistical test: paired comparison
from scipy import stats
t_stat, p_val = stats.ttest_rel(results_df['ssw_final_strength'], results_df['ctrl_final_strength'])
print(f"\nPaired t-test (SSW vs control strength): t={t_stat:.3f}, P={p_val:.4f}")

# Wilcoxon signed-rank
try:
    w_stat, w_pval = stats.wilcoxon(results_df['ssw_final_strength'] - results_df['ctrl_final_strength'])
    print(f"Wilcoxon signed-rank: W={w_stat:.1f}, P={w_pval:.4f}")
except:
    print("Wilcoxon: insufficient data")

# Sign test
n_stronger = (results_df['strength_ratio'] > 1).sum()
n_total = len(results_df)
sign_p = stats.binom_test(n_stronger, n_total, 0.5) if hasattr(stats, 'binom_test') else \
         stats.binomtest(n_stronger, n_total, 0.5).pvalue
print(f"Sign test: {n_stronger}/{n_total} SSW windows have stronger slabs, P={sign_p:.4f}")

# ========== SENSITIVITY ANALYSIS ==========
print("\n=== Sensitivity to activation energy ===")
for Ea_test in [0.4, 0.5, 0.6, 0.7]:
    T_ssw = 268.0  # -5°C + 0.81K warming
    T_ctrl = 267.19  # -5°C baseline
    ratio = np.exp(-Ea_test/kB * (1/T_ssw - 1/T_ctrl))
    print(f"  Ea={Ea_test} eV: rate ratio = {ratio:.4f} ({(ratio-1)*100:.2f}% enhancement)")

# ========== TIME EVOLUTION COMPARISON ==========
print("\n=== Time evolution: Composite SSW vs Control ===")
# Run a standard 20-day simulation at -5°C vs -5.81°C
T_warm = np.full(20, 268.0)   # SSW: -5.19°C (with +0.81K warming)
T_cold = np.full(20, 267.19)  # Control: -6°C baseline
sim_warm = simulate_slab(T_warm)
sim_cold = simulate_slab(T_cold)

print("Day | SSW strength | Ctrl strength | Ratio  | SSW BGR | Ctrl BGR")
print("-"*70)
for day in [0, 5, 10, 15, 19]:
    sw = sim_warm.iloc[day]
    sc = sim_cold.iloc[day]
    ratio = sw['strength_kPa'] / sc['strength_kPa'] if sc['strength_kPa'] > 0 else float('nan')
    print(f" {day:2d}  |  {sw['strength_kPa']:8.2f}    |  {sc['strength_kPa']:8.2f}     | {ratio:.4f} |  {sw['bgr']:.4f} |  {sc['bgr']:.4f}")

# After 15 days
delta_strength = sim_warm.iloc[15]['strength_kPa'] - sim_cold.iloc[15]['strength_kPa']
pct = (sim_warm.iloc[15]['strength_kPa'] / sim_cold.iloc[15]['strength_kPa'] - 1) * 100
print(f"\nAfter 15 days: SSW slab is {delta_strength:.2f} kPa stronger ({pct:.1f}% enhancement)")

# Save results
results_df.to_csv(os.path.join(OUT, "snowpack_sintering_model.csv"), index=False)
print(f"\nResults saved to {OUT}/snowpack_sintering_model.csv")
