"""
R77: Height-time composite of stratospheric deceleration propagation.

Creates a Hovmöller-style composite showing zonal wind tendency (du/dt)
at 4 pressure levels (10, 30, 50, 100 hPa) as a function of lag from
SSW onset, demonstrating the downward propagation signature.

Input:  data/results/r57_epflux_diagnostics.json
Output: data/results/r77_height_time_composite.json
"""
import json
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "results" / "r57_epflux_diagnostics.json"
OUTPUT = ROOT / "data" / "results" / "r77_height_time_composite.json"

with open(INPUT) as f:
    data = json.load(f)

ct = data["composite_timeseries"]
dp = data["downward_propagation_speed"]

levels = [10, 30, 50, 100]
level_keys = {10: "dUdt_10", 30: "dUdt_30", 50: "dUdt_50", 100: "dUdt_100"}
u_keys = {10: "U10", 30: "U30", 50: "U50", 100: "U100"}

days = ct["dUdt_10"]["days"]

composite = {}
for lev in levels:
    key = level_keys[lev]
    series = ct[key]
    mean_vals = series["mean"]
    std_vals = series["std"]
    
    # Find peak deceleration day
    peak_idx = int(np.argmin(mean_vals))
    peak_day = days[peak_idx]
    peak_val = mean_vals[peak_idx]
    
    # t-statistic at peak (n=16 events)
    n = 16
    t_stat = abs(peak_val) / (std_vals[peak_idx] / np.sqrt(n))
    
    composite[f"{lev}hPa"] = {
        "days": days,
        "mean_dUdt": mean_vals,
        "std_dUdt": std_vals,
        "peak_deceleration_day": peak_day,
        "peak_dUdt_ms_per_day": round(peak_val, 3),
        "t_statistic_at_peak": round(t_stat, 2),
        "n_events": n,
    }

# Add zonal wind composites
for lev in levels:
    key = u_keys[lev]
    if key in ct:
        series = ct[key]
        composite[f"{lev}hPa"]["mean_U"] = series["mean"]
        composite[f"{lev}hPa"]["std_U"] = series["std"]

# Downward propagation summary
prop_summary = {}
for lev_str, vals in dp["composite"].items():
    prop_summary[lev_str] = {
        "mean_lag_days": vals["mean_lag_days"],
        "median_lag_days": vals["median_lag_days"],
        "std_lag_days": vals["std_lag_days"],
    }

# Compute height-time matrix for the Hovmöller description
hovmoller = {
    "pressure_levels_hPa": levels,
    "lag_days": days,
    "matrix_dUdt_mean": [],
    "peak_deceleration_by_level": {},
}

for lev in levels:
    key = level_keys[lev]
    hovmoller["matrix_dUdt_mean"].append(ct[key]["mean"])
    lev_data = composite[f"{lev}hPa"]
    hovmoller["peak_deceleration_by_level"][f"{lev}hPa"] = {
        "peak_day": lev_data["peak_deceleration_day"],
        "peak_magnitude": lev_data["peak_dUdt_ms_per_day"],
        "t_stat": lev_data["t_statistic_at_peak"],
    }

# Propagation speed calculation
prop_10_to_100 = dp["composite"]["100hPa"]["mean_lag_days"]
height_diff_km = 16 - 31  # approximate heights: 10hPa~31km, 100hPa~16km
speed_km_per_day = abs(height_diff_km) / max(prop_10_to_100, 0.5)

result = {
    "metadata": {
        "description": "Height-time composite of SSW deceleration propagation",
        "n_events": 16,
        "pressure_levels": levels,
        "lag_range_days": [days[0], days[-1]],
    },
    "composite_by_level": composite,
    "propagation_summary": prop_summary,
    "hovmoller": hovmoller,
    "propagation_speed": {
        "10_to_100hPa_mean_lag_days": prop_10_to_100,
        "approximate_descent_speed_km_per_day": round(speed_km_per_day, 1),
    },
}

with open(OUTPUT, "w") as f:
    json.dump(result, f, indent=2)

print(f"Height-time composite saved to {OUTPUT}")
print(f"\nPeak deceleration by level:")
for lev in levels:
    d = composite[f"{lev}hPa"]
    print(f"  {lev:>4} hPa: day {d['peak_deceleration_day']:+3d}, "
          f"du/dt = {d['peak_dUdt_ms_per_day']:+.3f} m/s/d, "
          f"t = {d['t_statistic_at_peak']:.1f}")
print(f"\nPropagation 10→100 hPa: {prop_10_to_100:.1f} days")
