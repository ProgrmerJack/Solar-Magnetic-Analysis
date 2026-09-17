"""
67_regime_occupancy_decomposition.py
Quantifies per-event regime-occupancy shifts and decomposes the SSW-avalanche
association into between-regime (frequency redistribution) and within-regime
(rate change) components for each of the 16 SSW events.

Output: data/results/r67_regime_occupancy_decomposition.json
"""
import pandas as pd
import numpy as np
import json
from pathlib import Path
from scipy import stats

ROOT = Path(r"C:\Users\Jack0\Solar-Magnetic-Analysis")
DATA = ROOT / "data"
OUT = DATA / "results" / "r67_regime_occupancy_decomposition.json"

# Load SNOWPACK weather data
weather = pd.read_csv(DATA / "cryosphere" / "envidat" / "weather_snowpack_danger.csv")
weather["date"] = pd.to_datetime(weather["datum"])

# Load daily avalanche counts (semicolon-delimited)
daily = pd.read_csv(DATA / "cryosphere" / "davos_avalanches" / "daily_activity.csv", sep=";")
daily["date"] = pd.to_datetime(daily.iloc[:, 0])

# Load SSW events
events = pd.read_csv(DATA / "results" / "ssw_event_catalog.csv")
events["onset"] = pd.to_datetime(events["date"])

# Define weather regimes using station-averaged daily TA and precipitation
station_daily = weather.groupby("date").agg(
    TA_mean=("TA", "mean"),
    precip=("MS_Rain", "sum"),
    HN24=("HN24", "mean"),
    VW=("VW", "mean"),
).reset_index()

# Merge with daily counts
if "dry.NATURAL.SIZE_1234" in daily.columns:
    count_col = "dry.NATURAL.SIZE_1234"
elif "AAI.dry.NATURAL" in daily.columns:
    count_col = "AAI.dry.NATURAL"
elif "dry_natural" in daily.columns:
    count_col = "dry_natural"
elif "count" in daily.columns:
    count_col = "count"
else:
    num_cols = daily.select_dtypes(include=[np.number]).columns.tolist()
    print(f"Available numeric columns: {num_cols}")
    count_col = num_cols[0] if num_cols else None

print(f"Using count column: {count_col}")
print(f"Daily data columns: {daily.columns.tolist()}")

merged = station_daily.merge(daily[["date", count_col]], on="date", how="inner") if count_col else None
if merged is None:
    raise ValueError("No count column found")

# Define four regimes based on temperature and precipitation
# Cold-dry, Cold-wet, Warm-dry, Warm-wet using median splits
ta_median = merged["TA_mean"].median()
hn_median = merged["HN24"].median()

def classify_regime(row):
    cold = row["TA_mean"] < ta_median
    snowy = row["HN24"] > hn_median
    if cold and not snowy:
        return "cold_dry"
    elif cold and snowy:
        return "cold_wet"
    elif not cold and not snowy:
        return "warm_dry"
    else:
        return "warm_wet"

merged["regime"] = merged.apply(classify_regime, axis=1)

# Define SSW windows and control windows
HALF_WIN = 15  # ±15 days

ssw_winters = set()
for _, ev in events.iterrows():
    winter_year = ev["onset"].year if ev["onset"].month >= 10 else ev["onset"].year - 1
    ssw_winters.add(winter_year)

merged["winter_year"] = merged["date"].apply(
    lambda d: d.year if d.month >= 10 else d.year - 1
)
merged["doy"] = merged["date"].dt.dayofyear

# Tag SSW window days
merged["in_ssw_window"] = False
merged["event_id"] = None
for i, ev in events.iterrows():
    onset = ev["onset"]
    start = onset - pd.Timedelta(days=HALF_WIN)
    end = onset + pd.Timedelta(days=HALF_WIN)
    mask = (merged["date"] >= start) & (merged["date"] <= end)
    merged.loc[mask, "in_ssw_window"] = True
    merged.loc[mask, "event_id"] = i

# Control: same DOY range from non-SSW winters
merged["is_ssw_winter"] = merged["winter_year"].isin(ssw_winters)

results = {"per_event": [], "aggregate": {}}
regimes = ["cold_dry", "cold_wet", "warm_dry", "warm_wet"]

# Per-event regime occupancy
for i, ev in events.iterrows():
    onset = ev["onset"]
    start_doy = (onset - pd.Timedelta(days=HALF_WIN)).dayofyear
    end_doy = (onset + pd.Timedelta(days=HALF_WIN)).dayofyear
    
    # SSW window days
    ssw_mask = merged["event_id"] == i
    ssw_days = merged[ssw_mask]
    
    # Control: same DOY from non-SSW winters
    if start_doy <= end_doy:
        ctrl_mask = (~merged["is_ssw_winter"]) & (merged["doy"] >= start_doy) & (merged["doy"] <= end_doy)
    else:
        ctrl_mask = (~merged["is_ssw_winter"]) & ((merged["doy"] >= start_doy) | (merged["doy"] <= end_doy))
    ctrl_days = merged[ctrl_mask]
    
    if len(ssw_days) == 0 or len(ctrl_days) == 0:
        continue
    
    ev_result = {
        "event_date": str(ev["date"]),
        "n_ssw_days": int(len(ssw_days)),
        "n_ctrl_days": int(len(ctrl_days)),
        "ssw_regime_fractions": {},
        "ctrl_regime_fractions": {},
        "regime_shift": {},
        "ssw_count_mean": float(ssw_days[count_col].mean()),
        "ctrl_count_mean": float(ctrl_days[count_col].mean()),
        "event_rr": float(ev["rr"]),
    }
    
    for r in regimes:
        ssw_frac = (ssw_days["regime"] == r).mean()
        ctrl_frac = (ctrl_days["regime"] == r).mean()
        ev_result["ssw_regime_fractions"][r] = round(float(ssw_frac), 4)
        ev_result["ctrl_regime_fractions"][r] = round(float(ctrl_frac), 4)
        ev_result["regime_shift"][r] = round(float(ssw_frac - ctrl_frac), 4)
    
    results["per_event"].append(ev_result)

# Aggregate regime occupancy
all_ssw = merged[merged["in_ssw_window"]]
all_ctrl = merged[~merged["is_ssw_winter"]]

agg = {}
for r in regimes:
    ssw_frac = (all_ssw["regime"] == r).mean()
    ctrl_frac = (all_ctrl["regime"] == r).mean()
    
    # Within-regime avalanche rates
    ssw_rate = all_ssw[all_ssw["regime"] == r][count_col].mean() if (all_ssw["regime"] == r).sum() > 0 else 0
    ctrl_rate = all_ctrl[all_ctrl["regime"] == r][count_col].mean() if (all_ctrl["regime"] == r).sum() > 0 else 0
    within_rr = ssw_rate / ctrl_rate if ctrl_rate > 0 else float('nan')
    
    agg[r] = {
        "ssw_fraction": round(float(ssw_frac), 4),
        "ctrl_fraction": round(float(ctrl_frac), 4),
        "shift": round(float(ssw_frac - ctrl_frac), 4),
        "ratio": round(float(ssw_frac / ctrl_frac), 3) if ctrl_frac > 0 else None,
        "ssw_daily_rate": round(float(ssw_rate), 4),
        "ctrl_daily_rate": round(float(ctrl_rate), 4),
        "within_regime_rr": round(float(within_rr), 4) if not np.isnan(within_rr) else None,
    }

results["aggregate"]["regime_occupancy"] = agg

# Decomposition: between-regime vs within-regime
# Expected RR from occupancy shift alone (holding rates constant at control levels)
ctrl_overall_rate = all_ctrl[count_col].mean()
counterfactual_rate = sum(
    agg[r]["ssw_fraction"] * agg[r]["ctrl_daily_rate"] for r in regimes
)
occupancy_only_rr = counterfactual_rate / ctrl_overall_rate if ctrl_overall_rate > 0 else float('nan')

# Expected RR from rate change alone (holding occupancy constant at control)
rate_change_rate = sum(
    agg[r]["ctrl_fraction"] * agg[r]["ssw_daily_rate"] for r in regimes
)
rate_only_rr = rate_change_rate / ctrl_overall_rate if ctrl_overall_rate > 0 else float('nan')

# Actual aggregate
actual_ssw_rate = all_ssw[count_col].mean()
actual_rr = actual_ssw_rate / ctrl_overall_rate if ctrl_overall_rate > 0 else float('nan')

results["aggregate"]["decomposition"] = {
    "actual_rr": round(float(actual_rr), 4),
    "occupancy_shift_only_rr": round(float(occupancy_only_rr), 4),
    "rate_change_only_rr": round(float(rate_only_rr), 4),
    "pct_explained_by_occupancy": round(
        float((1 - occupancy_only_rr) / (1 - actual_rr) * 100), 1
    ) if actual_rr != 1 else None,
    "pct_explained_by_rate_change": round(
        float((1 - rate_only_rr) / (1 - actual_rr) * 100), 1
    ) if actual_rr != 1 else None,
}

# Per-event: count how many events show cold-dry increase
n_cold_dry_increase = sum(
    1 for ev in results["per_event"] if ev["regime_shift"].get("cold_dry", 0) > 0
)
n_warm_wet_decrease = sum(
    1 for ev in results["per_event"] if ev["regime_shift"].get("warm_wet", 0) < 0
)

results["aggregate"]["event_consistency"] = {
    "n_events": len(results["per_event"]),
    "n_cold_dry_increase": n_cold_dry_increase,
    "n_warm_wet_decrease": n_warm_wet_decrease,
    "cold_dry_sign_test_p": round(float(
        1 - stats.binom.cdf(n_cold_dry_increase - 1, len(results["per_event"]), 0.5)
    ), 6),
    "warm_wet_sign_test_p": round(float(
        1 - stats.binom.cdf(n_warm_wet_decrease - 1, len(results["per_event"]), 0.5)
    ), 6),
}

# Correlation: event-level cold-dry shift vs avalanche RR
shifts = [ev["regime_shift"]["cold_dry"] for ev in results["per_event"]]
rrs = [ev["event_rr"] for ev in results["per_event"]]
if len(shifts) >= 5:
    r_val, p_val = stats.spearmanr(shifts, rrs)
    results["aggregate"]["cold_dry_shift_vs_rr"] = {
        "spearman_r": round(float(r_val), 4),
        "p_value": round(float(p_val), 6),
    }

OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w") as f:
    json.dump(results, f, indent=2, default=str)

print(f"\nSaved to {OUT}")
print(f"\n=== AGGREGATE REGIME OCCUPANCY ===")
for r in regimes:
    a = agg[r]
    print(f"  {r:12s}: SSW={a['ssw_fraction']:.3f}  Ctrl={a['ctrl_fraction']:.3f}  "
          f"Shift={a['shift']:+.3f}  Ratio={a['ratio']}  "
          f"Within-RR={a['within_regime_rr']}")

print(f"\n=== DECOMPOSITION ===")
d = results["aggregate"]["decomposition"]
print(f"  Actual RR: {d['actual_rr']:.4f}")
print(f"  Occupancy-shift-only RR: {d['occupancy_shift_only_rr']:.4f}")
print(f"  Rate-change-only RR: {d['rate_change_only_rr']:.4f}")
print(f"  % explained by occupancy: {d['pct_explained_by_occupancy']}%")
print(f"  % explained by rate change: {d['pct_explained_by_rate_change']}%")

print(f"\n=== EVENT CONSISTENCY ===")
ec = results["aggregate"]["event_consistency"]
print(f"  {ec['n_cold_dry_increase']}/{ec['n_events']} events show cold-dry increase (P={ec['cold_dry_sign_test_p']})")
print(f"  {ec['n_warm_wet_decrease']}/{ec['n_events']} events show warm-wet decrease (P={ec['warm_wet_sign_test_p']})")

if "cold_dry_shift_vs_rr" in results["aggregate"]:
    c = results["aggregate"]["cold_dry_shift_vs_rr"]
    print(f"\n=== COLD-DRY SHIFT vs AVALANCHE RR ===")
    print(f"  Spearman r = {c['spearman_r']}, P = {c['p_value']}")
