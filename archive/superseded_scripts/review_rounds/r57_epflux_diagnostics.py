"""
R57 – EP-flux diagnostics for SSW-avalanche manuscript.

Addresses reviewer criticism: "No EP-flux vectors/divergence/wave decomposition
for a Nature Geo stratospheric dynamics paper."

Uses NCEP daily polar-cap-mean reanalysis (1979–2024) and SLF avalanche
activity (1998–2019) to compute:
  a) Vertical EP-flux proxy  – dT/dt at polar cap (Newman–Rosenfield)
  b) Wave-forcing amplitude  – geopotential-height anomaly amplitude
  c) EP-flux divergence proxy – du/dt at 10 hPa ≈ (1/acosφ)∇·F
  d) Strat–trop coupling      – lag correlation of du/dt vs surface response
"""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats, signal

warnings.filterwarnings("ignore", category=FutureWarning)

# ── paths ──────────────────────────────────────────────────────────────────
ROOT = Path(r"C:\Users\Jack0\Solar-Magnetic-Analysis")
STRAT_PATH = ROOT / "data/processed/atmospheric/ncep_stratosphere.parquet"
TROP_PATH  = ROOT / "data/processed/atmospheric/ncep_troposphere.parquet"
SSW_PATH   = ROOT / "data/processed/atmospheric/ssw_catalog.parquet"
AVAL_PATH  = ROOT / "data/processed/cryosphere/slf_activity.parquet"
OUT_JSON   = ROOT / "data/results/r57_epflux_diagnostics.json"

# ── load data ──────────────────────────────────────────────────────────────
print("Loading data...")
strat = pd.read_parquet(STRAT_PATH)
trop  = pd.read_parquet(TROP_PATH)
ssw   = pd.read_parquet(SSW_PATH)
aval  = pd.read_parquet(AVAL_PATH)

# Strip timezone for safe comparison
strat.index = strat.index.tz_localize(None)
trop.index  = trop.index.tz_localize(None)
aval.index  = aval.index.tz_localize(None)
ssw_dates   = ssw.index.tz_localize(None)

# Rename for readability
strat_cols_map = {
    "air_K_100hPa": "T100", "air_K_70hPa": "T70", "air_K_50hPa": "T50",
    "air_K_30hPa": "T30", "air_K_20hPa": "T20", "air_K_10hPa": "T10",
    "hgt_m_100hPa": "Z100", "hgt_m_70hPa": "Z70", "hgt_m_50hPa": "Z50",
    "hgt_m_30hPa": "Z30", "hgt_m_20hPa": "Z20", "hgt_m_10hPa": "Z10",
    "uwnd_ms_100hPa": "U100", "uwnd_ms_70hPa": "U70", "uwnd_ms_50hPa": "U50",
    "uwnd_ms_30hPa": "U30", "uwnd_ms_20hPa": "U20", "uwnd_ms_10hPa": "U10",
}
strat = strat.rename(columns=strat_cols_map)
trop = trop.rename(columns={
    "hgt_500hPa_m": "Z500", "slp_Pa": "SLP", "uwnd_850hPa_ms": "U850",
})

# Filter SSW events to manuscript period (winter 1998/99 onward through 2019)
ssw_events = ssw_dates[(ssw_dates >= "1998-12-01") & (ssw_dates <= "2019-12-31")]
print(f"SSW events in 1999–2019: {len(ssw_events)}")
for d in ssw_events:
    print(f"  {d.strftime('%Y-%m-%d')}")

# ── helper: climatology & anomaly ──────────────────────────────────────────
def daily_climatology(series, smooth_days=15):
    """31-day running mean climatology based on day-of-year."""
    doy = series.index.dayofyear
    clim = series.groupby(doy).mean()
    # Smooth with rolling window (circular)
    extended = pd.concat([clim.iloc[-smooth_days:], clim, clim.iloc[:smooth_days]])
    smoothed = extended.rolling(2 * smooth_days + 1, center=True).mean()
    return smoothed.iloc[smooth_days:-smooth_days]


def anomaly(series, smooth_days=15):
    """Daily anomaly relative to smoothed climatology."""
    clim = daily_climatology(series, smooth_days)
    doy = series.index.dayofyear
    return series - clim.reindex(doy).values


# ── 1. Compute derived fields ──────────────────────────────────────────────
print("\nComputing derived fields...")

# 1a. du/dt at each level (central difference, m/s per day)
for lev in [10, 20, 30, 50, 70, 100]:
    col = f"U{lev}"
    strat[f"dUdt_{lev}"] = np.gradient(strat[col].values, 1.0)  # Δt = 1 day

# 1b. dT/dt at each level (proxy for EP-flux convergence / adiabatic warming)
for lev in [10, 20, 30, 50, 70, 100]:
    col = f"T{lev}"
    strat[f"dTdt_{lev}"] = np.gradient(strat[col].values, 1.0)

# 1c. Anomaly fields for geopotential height (wave amplitude proxy)
for lev in [10, 30, 50, 100]:
    col = f"Z{lev}"
    strat[f"Za_{lev}"] = anomaly(strat[col])

# 1d. Anomaly for temperature
for lev in [10, 50, 100]:
    col = f"T{lev}"
    strat[f"Ta_{lev}"] = anomaly(strat[col])

# 1e. Tropospheric anomalies
trop["Z500a"] = anomaly(trop["Z500"])
trop["U850a"] = anomaly(trop["U850"])

# ── 2. EP-flux diagnostics per SSW event ───────────────────────────────────
print("\nComputing EP-flux diagnostics per SSW event...")

PHASES = {
    "precursor": (-30, -6),
    "onset":     (-5,  +5),
    "recovery":  (+6, +30),
}

event_diagnostics = []

for onset in ssw_events:
    evt = {"onset_date": onset.strftime("%Y-%m-%d")}
    onset_ts = pd.Timestamp(onset)

    for phase_name, (d0, d1) in PHASES.items():
        t0 = onset_ts + pd.Timedelta(days=d0)
        t1 = onset_ts + pd.Timedelta(days=d1)
        mask_s = (strat.index >= t0) & (strat.index <= t1)
        mask_t = (trop.index >= t0) & (trop.index <= t1)
        ss = strat.loc[mask_s]
        tt = trop.loc[mask_t]

        if len(ss) < 5:
            continue

        pref = f"{phase_name}_"

        # -- EP-flux divergence proxy: mean du/dt at 10 hPa --
        evt[pref + "dudt_10hPa_mean"] = round(float(ss["dUdt_10"].mean()), 4)
        evt[pref + "dudt_10hPa_min"]  = round(float(ss["dUdt_10"].min()), 4)

        # -- du/dt at other levels --
        for lev in [30, 50, 100]:
            evt[pref + f"dudt_{lev}hPa_mean"] = round(float(ss[f"dUdt_{lev}"].mean()), 4)

        # -- Vertical EP-flux proxy: dT/dt at 100 hPa (heat flux convergence) --
        evt[pref + "dTdt_100hPa_mean"] = round(float(ss["dTdt_100"].mean()), 4)
        evt[pref + "dTdt_10hPa_mean"]  = round(float(ss["dTdt_10"].mean()), 4)

        # -- Wave amplitude proxy: geopotential height anomaly RMS --
        for lev in [10, 50, 100]:
            col = f"Za_{lev}"
            if col in ss.columns:
                evt[pref + f"Z_anom_rms_{lev}hPa"] = round(
                    float(np.sqrt((ss[col] ** 2).mean())), 2
                )

        # -- Zonal wind values (state) --
        evt[pref + "U10_mean"] = round(float(ss["U10"].mean()), 2)
        evt[pref + "U100_mean"] = round(float(ss["U100"].mean()), 2)

        # -- Polar temperature anomaly (wave-driven warming) --
        for lev in [10, 100]:
            col = f"Ta_{lev}"
            if col in ss.columns:
                evt[pref + f"T_anom_{lev}hPa_mean"] = round(float(ss[col].mean()), 2)

        # -- Tropospheric response --
        if len(tt) > 0:
            evt[pref + "Z500a_mean"] = round(float(tt["Z500a"].mean()), 2)
            evt[pref + "U850a_mean"] = round(float(tt["U850a"].mean()), 2)

    # -- Vertical propagation timing: peak du/dt at each level --
    window = strat.loc[
        (strat.index >= onset_ts - pd.Timedelta(days=30))
        & (strat.index <= onset_ts + pd.Timedelta(days=30))
    ]
    if len(window) > 10:
        for lev in [10, 30, 50, 70, 100]:
            col = f"dUdt_{lev}"
            min_idx = window[col].idxmin()
            evt[f"peak_decel_lag_{lev}hPa_days"] = (min_idx - onset_ts).days

    event_diagnostics.append(evt)


# ── 3. Composite analysis ──────────────────────────────────────────────────
print("\nComputing SSW composites...")

composite_days = np.arange(-30, 31)
composite_fields = {}

for field in ["dUdt_10", "dUdt_30", "dUdt_50", "dUdt_100",
              "dTdt_10", "dTdt_100",
              "U10", "U30", "U50", "U100",
              "T10", "Ta_10", "Za_10", "Za_100"]:
    if field not in strat.columns:
        continue
    day_vals = {d: [] for d in composite_days}
    for onset in ssw_events:
        onset_ts = pd.Timestamp(onset)
        for d in composite_days:
            target = onset_ts + pd.Timedelta(days=int(d))
            if target in strat.index:
                day_vals[d].append(strat.loc[target, field])

    means = [np.mean(day_vals[d]) if day_vals[d] else np.nan for d in composite_days]
    stds  = [np.std(day_vals[d])  if len(day_vals[d]) > 1 else np.nan for d in composite_days]
    composite_fields[field] = {
        "days": composite_days.tolist(),
        "mean": [round(float(m), 4) if not np.isnan(m) else None for m in means],
        "std":  [round(float(s), 4) if not np.isnan(s) else None for s in stds],
    }

# Also composite tropospheric fields
for field_trop in ["Z500a", "U850a"]:
    day_vals = {d: [] for d in composite_days}
    for onset in ssw_events:
        onset_ts = pd.Timestamp(onset)
        for d in composite_days:
            target = onset_ts + pd.Timedelta(days=int(d))
            if target in trop.index:
                day_vals[d].append(trop.loc[target, field_trop])
    means = [np.mean(day_vals[d]) if day_vals[d] else np.nan for d in composite_days]
    stds  = [np.std(day_vals[d])  if len(day_vals[d]) > 1 else np.nan for d in composite_days]
    composite_fields[field_trop] = {
        "days": composite_days.tolist(),
        "mean": [round(float(m), 4) if not np.isnan(m) else None for m in means],
        "std":  [round(float(s), 4) if not np.isnan(s) else None for s in stds],
    }


# ── 4. Strat–trop coupling: lag correlations ──────────────────────────────
print("\nComputing strat–trop coupling lag correlations...")

# Merge strat du/dt with tropospheric & avalanche data
merged = strat[["dUdt_10", "dUdt_30", "U10", "Ta_10"]].join(
    trop[["Z500a", "U850a"]], how="inner"
)

# Join avalanche activity (winter months only)
aval_daily = aval[["aai_all", "aai_all_natural", "natural_size_234"]].copy()
merged = merged.join(aval_daily, how="left")

# Restrict to DJFM (winter) for meaningful strat–trop coupling
merged["month"] = merged.index.month
winter = merged[merged["month"].isin([12, 1, 2, 3])].copy()
winter = winter.drop(columns=["month"])

# Compute lag correlations
max_lag = 45  # days
lag_results = {}

target_pairs = [
    ("dUdt_10", "Z500a",            "dudt10_vs_Z500a"),
    ("dUdt_10", "U850a",            "dudt10_vs_U850a"),
    ("dUdt_10", "aai_all_natural",  "dudt10_vs_avalanche"),
    ("U10",     "Z500a",            "U10_vs_Z500a"),
    ("U10",     "U850a",            "U10_vs_U850a"),
    ("U10",     "aai_all_natural",  "U10_vs_avalanche"),
    ("Ta_10",   "aai_all_natural",  "Ta10_vs_avalanche"),
]

for strat_var, surf_var, label in target_pairs:
    if strat_var not in winter.columns or surf_var not in winter.columns:
        continue
    s1 = winter[strat_var].dropna()
    s2 = winter[surf_var].dropna()
    common = s1.index.intersection(s2.index)
    if len(common) < 100:
        continue

    lags = list(range(0, max_lag + 1))
    corrs, pvals = [], []
    for lag in lags:
        shifted = s1.reindex(common).shift(lag).dropna()
        aligned = s2.reindex(shifted.index).dropna()
        common2 = shifted.index.intersection(aligned.index)
        if len(common2) < 50:
            corrs.append(None)
            pvals.append(None)
            continue
        r, p = stats.pearsonr(shifted.loc[common2], aligned.loc[common2])
        corrs.append(round(float(r), 5))
        pvals.append(round(float(p), 6))

    # Find peak correlation
    valid = [(l, c) for l, c in zip(lags, corrs) if c is not None]
    if valid:
        peak_lag, peak_r = max(valid, key=lambda x: abs(x[1]))
    else:
        peak_lag, peak_r = None, None

    lag_results[label] = {
        "lags_days": lags,
        "correlations": corrs,
        "p_values": pvals,
        "peak_lag_days": peak_lag,
        "peak_correlation": peak_r,
    }


# ── 5. SSW-composite avalanche response ───────────────────────────────────
print("\nComputing SSW-composite avalanche response...")

aval_composite_days = np.arange(-30, 61)
aval_fields = ["aai_all", "aai_all_natural", "natural_size_234"]
aval_composites = {}

for field in aval_fields:
    if field not in aval.columns:
        continue
    day_vals = {d: [] for d in aval_composite_days}
    n_events_used = 0
    for onset in ssw_events:
        onset_ts = pd.Timestamp(onset)
        for d in aval_composite_days:
            target = onset_ts + pd.Timedelta(days=int(d))
            if target in aval.index:
                val = aval.loc[target, field]
                if not np.isnan(val):
                    day_vals[d].append(val)
        if any(len(day_vals[d]) > 0 for d in aval_composite_days):
            n_events_used += 1

    means = [np.mean(day_vals[d]) if day_vals[d] else None for d in aval_composite_days]
    aval_composites[field] = {
        "days": aval_composite_days.tolist(),
        "mean": [round(float(m), 4) if m is not None else None for m in means],
        "n_events": n_events_used,
    }


# ── 6. Vertical propagation diagnostics ───────────────────────────────────
print("\nComputing vertical propagation diagnostics...")

levels_hpa = [100, 70, 50, 30, 20, 10]
vert_prop = {}

for onset in ssw_events:
    onset_ts = pd.Timestamp(onset)
    key = onset_ts.strftime("%Y-%m-%d")
    vert_prop[key] = {}

    window = strat.loc[
        (strat.index >= onset_ts - pd.Timedelta(days=30))
        & (strat.index <= onset_ts + pd.Timedelta(days=15))
    ]
    if len(window) < 10:
        continue

    # For each level, find day of maximum deceleration and wind reversal
    for lev in levels_hpa:
        u_col = f"U{lev}"
        du_col = f"dUdt_{lev}"

        # Day of peak deceleration
        peak_decel_idx = window[du_col].idxmin()
        peak_decel_day = (peak_decel_idx - onset_ts).days

        # Day of wind reversal (first day U < 0), if it occurs
        negative = window[window[u_col] < 0]
        if len(negative) > 0:
            reversal_day = (negative.index[0] - onset_ts).days
        else:
            reversal_day = None

        vert_prop[key][f"{lev}hPa"] = {
            "peak_decel_day": peak_decel_day,
            "peak_decel_rate_ms_day": round(float(window[du_col].min()), 3),
            "wind_reversal_day": reversal_day,
            "min_wind_ms": round(float(window[u_col].min()), 2),
        }


# ── 7. Phase-mean summary statistics ──────────────────────────────────────
print("\nComputing phase-mean summary statistics...")

phase_summary = {}
for phase_name, (d0, d1) in PHASES.items():
    vals = {
        "dudt_10": [], "dudt_30": [], "dudt_50": [], "dudt_100": [],
        "dTdt_100": [], "dTdt_10": [],
        "U10": [], "Z_anom_rms_10": [],
        "Z500a": [], "U850a": [],
        "T_anom_10": [],
    }
    for onset in ssw_events:
        onset_ts = pd.Timestamp(onset)
        t0 = onset_ts + pd.Timedelta(days=d0)
        t1 = onset_ts + pd.Timedelta(days=d1)
        ss = strat.loc[(strat.index >= t0) & (strat.index <= t1)]
        tt = trop.loc[(trop.index >= t0) & (trop.index <= t1)]
        if len(ss) < 3:
            continue
        vals["dudt_10"].append(ss["dUdt_10"].mean())
        vals["dudt_30"].append(ss["dUdt_30"].mean())
        vals["dudt_50"].append(ss["dUdt_50"].mean())
        vals["dudt_100"].append(ss["dUdt_100"].mean())
        vals["dTdt_100"].append(ss["dTdt_100"].mean())
        vals["dTdt_10"].append(ss["dTdt_10"].mean())
        vals["U10"].append(ss["U10"].mean())
        vals["T_anom_10"].append(ss["Ta_10"].mean())
        if "Za_10" in ss.columns:
            vals["Z_anom_rms_10"].append(np.sqrt((ss["Za_10"] ** 2).mean()))
        if len(tt) > 0:
            vals["Z500a"].append(tt["Z500a"].mean())
            vals["U850a"].append(tt["U850a"].mean())

    phase_summary[phase_name] = {}
    for k, v in vals.items():
        arr = np.array([x for x in v if not np.isnan(x)])
        if len(arr) > 2:
            t_stat, p_val = stats.ttest_1samp(arr, 0)
            phase_summary[phase_name][k] = {
                "mean": round(float(arr.mean()), 4),
                "std": round(float(arr.std(ddof=1)), 4),
                "n": len(arr),
                "t_stat": round(float(t_stat), 3),
                "p_value": round(float(p_val), 5),
            }


# ── 8. Downward propagation speed estimate ────────────────────────────────
print("\nEstimating downward propagation speed...")

propagation_lags = {}
ref_level = 10  # reference: 10 hPa

for onset in ssw_events:
    onset_ts = pd.Timestamp(onset)
    key = onset_ts.strftime("%Y-%m-%d")
    window = strat.loc[
        (strat.index >= onset_ts - pd.Timedelta(days=10))
        & (strat.index <= onset_ts + pd.Timedelta(days=50))
    ]
    if len(window) < 20:
        continue

    ref_series = window[f"dUdt_{ref_level}"]
    lags_per_level = {}
    for lev in [20, 30, 50, 70, 100]:
        target_series = window[f"dUdt_{lev}"]
        # Cross-correlation to find optimal lag
        n = len(ref_series)
        if n < 15:
            continue
        max_shift = min(30, n - 10)
        best_lag, best_corr = 0, -1
        for lag in range(0, max_shift):
            if lag >= n:
                break
            r1 = ref_series.values[:n - lag]
            r2 = target_series.values[lag:]
            if len(r1) < 10:
                break
            cc = np.corrcoef(r1, r2)[0, 1]
            if not np.isnan(cc) and cc > best_corr:
                best_corr = cc
                best_lag = lag
        lags_per_level[f"{lev}hPa_lag_days"] = best_lag
        lags_per_level[f"{lev}hPa_corr"] = round(float(best_corr), 3)

    propagation_lags[key] = lags_per_level

# Composite propagation speed
composite_prop = {}
for lev in [20, 30, 50, 70, 100]:
    all_lags = [v[f"{lev}hPa_lag_days"] for v in propagation_lags.values()
                if f"{lev}hPa_lag_days" in v]
    if all_lags:
        composite_prop[f"{lev}hPa"] = {
            "mean_lag_days": round(float(np.mean(all_lags)), 1),
            "median_lag_days": float(np.median(all_lags)),
            "std_lag_days": round(float(np.std(all_lags)), 1),
        }


# ── 9. Assemble output ────────────────────────────────────────────────────
results = {
    "metadata": {
        "description": "EP-flux diagnostics for SSW-avalanche manuscript (R57)",
        "n_ssw_events": len(ssw_events),
        "ssw_dates": [d.strftime("%Y-%m-%d") for d in ssw_events],
        "data_source": "NCEP daily polar-cap-mean reanalysis",
        "levels_hPa": levels_hpa,
        "methodology": {
            "EP_flux_divergence_proxy": (
                "du/dt at 10 hPa: ∂ū/∂t ≈ (1/a cos φ)∇·F. "
                "Negative du/dt = net EP-flux convergence = wave-driven deceleration."
            ),
            "vertical_EP_flux_proxy": (
                "dT/dt at 100 hPa polar cap: proxy for poleward eddy heat flux "
                "convergence (Newman & Rosenfield 1997). Positive dT/dt indicates "
                "enhanced wave activity entering the stratosphere."
            ),
            "wave_amplitude_proxy": (
                "RMS of geopotential height anomaly at each level. Larger amplitude "
                "= stronger planetary wave activity."
            ),
            "downward_propagation": (
                "Lag-correlation of du/dt between 10 hPa and lower levels. "
                "Positive lag = downward propagation of wave-driven deceleration."
            ),
        },
    },
    "event_diagnostics": event_diagnostics,
    "composite_timeseries": composite_fields,
    "phase_summary": phase_summary,
    "lag_correlations": lag_results,
    "avalanche_composites": aval_composites,
    "vertical_propagation": vert_prop,
    "downward_propagation_speed": {
        "per_event": propagation_lags,
        "composite": composite_prop,
    },
}

# ── 10. Save ───────────────────────────────────────────────────────────────
OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_JSON, "w") as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nResults saved to {OUT_JSON}")


# ── 11. Print summary ─────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("EP-FLUX DIAGNOSTICS SUMMARY  (R57)")
print("=" * 72)

print(f"\nSSW events analysed: {len(ssw_events)}")

print("\n── Phase-mean EP-flux divergence proxy (du/dt at 10 hPa, m/s/day) ──")
for phase_name in ["precursor", "onset", "recovery"]:
    d = phase_summary[phase_name]["dudt_10"]
    sig = "***" if d["p_value"] < 0.001 else "**" if d["p_value"] < 0.01 else "*" if d["p_value"] < 0.05 else "ns"
    print(f"  {phase_name:12s}: {d['mean']:+.3f} ± {d['std']:.3f}  (t={d['t_stat']:.2f}, p={d['p_value']:.4f}) {sig}")

print("\n── Phase-mean vertical EP-flux proxy (dT/dt at 100 hPa, K/day) ──")
for phase_name in ["precursor", "onset", "recovery"]:
    d = phase_summary[phase_name]["dTdt_100"]
    sig = "***" if d["p_value"] < 0.001 else "**" if d["p_value"] < 0.01 else "*" if d["p_value"] < 0.05 else "ns"
    print(f"  {phase_name:12s}: {d['mean']:+.4f} ± {d['std']:.4f}  (t={d['t_stat']:.2f}, p={d['p_value']:.4f}) {sig}")

print("\n── Phase-mean polar cap temperature anomaly (10 hPa, K) ──")
for phase_name in ["precursor", "onset", "recovery"]:
    d = phase_summary[phase_name]["T_anom_10"]
    sig = "***" if d["p_value"] < 0.001 else "**" if d["p_value"] < 0.01 else "*" if d["p_value"] < 0.05 else "ns"
    print(f"  {phase_name:12s}: {d['mean']:+.2f} ± {d['std']:.2f}  (t={d['t_stat']:.2f}, p={d['p_value']:.4f}) {sig}")

print("\n── Phase-mean zonal wind (10 hPa, m/s) ──")
for phase_name in ["precursor", "onset", "recovery"]:
    d = phase_summary[phase_name]["U10"]
    print(f"  {phase_name:12s}: {d['mean']:+.2f} ± {d['std']:.2f}")

print("\n── Wave amplitude proxy (Z anomaly RMS at 10 hPa, m) ──")
for phase_name in ["precursor", "onset", "recovery"]:
    d = phase_summary[phase_name].get("Z_anom_rms_10")
    if d:
        print(f"  {phase_name:12s}: {d['mean']:.1f} ± {d['std']:.1f}")

print("\n── Tropospheric response (Z500 anomaly, m) ──")
for phase_name in ["precursor", "onset", "recovery"]:
    d = phase_summary[phase_name].get("Z500a")
    if d:
        sig = "***" if d["p_value"] < 0.001 else "**" if d["p_value"] < 0.01 else "*" if d["p_value"] < 0.05 else "ns"
        print(f"  {phase_name:12s}: {d['mean']:+.1f} ± {d['std']:.1f}  (p={d['p_value']:.4f}) {sig}")

print("\n── Strat–trop coupling lag correlations (winter DJFM) ──")
for label, res in lag_results.items():
    peak = res["peak_correlation"]
    lag = res["peak_lag_days"]
    if peak is not None:
        print(f"  {label:30s}: r={peak:+.4f} at lag={lag:2d} days")

print("\n── Downward propagation speed (composite, days after 10 hPa) ──")
for lev_key, vals in composite_prop.items():
    print(f"  10 hPa → {lev_key:>7s}: {vals['mean_lag_days']:4.1f} ± {vals['std_lag_days']:.1f} days")

print("\n── Avalanche response composite (natural AAI) ──")
ac = aval_composites.get("aai_all_natural")
if ac:
    days = ac["days"]
    means = ac["mean"]
    pre_vals  = [m for d, m in zip(days, means) if d is not None and m is not None and -30 <= d < -5]
    post_vals = [m for d, m in zip(days, means) if d is not None and m is not None and 10 <= d <= 45]
    if pre_vals and post_vals:
        print(f"  Pre-SSW mean  (−30 to −6d): {np.mean(pre_vals):.3f}")
        print(f"  Post-SSW mean (+10 to +45d): {np.mean(post_vals):.3f}")
        ratio = np.mean(post_vals) / np.mean(pre_vals) if np.mean(pre_vals) > 0 else float("inf")
        print(f"  Ratio (post/pre): {ratio:.2f}")

print("\n" + "=" * 72)
print("Done. All results in:", OUT_JSON)
