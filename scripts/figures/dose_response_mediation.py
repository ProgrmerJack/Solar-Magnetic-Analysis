#!/usr/bin/env python3
"""
Dose-response mediation figure: stratospheric vs tropospheric predictors
of event-level avalanche response.

Generates Extended Data Fig 3: two-panel scatter showing that event-level
avalanche sensitivity concentrates at the tropospheric (Z500) level
rather than the stratospheric level.

Input:  data/processed/analysis_panel_v2.parquet (master daily panel)
        data/processed/atmospheric/ncep_stratosphere.parquet (50 hPa T)
Output: data/figures/extended_dose_response.pdf
        data/results/event_level_dose_response.csv
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "data" / "results"
FIGURES = ROOT / "data" / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

# SSW central dates (Butler et al. 2017 catalog)
SSW_DATES = pd.to_datetime([
    "1998-12-15", "1999-02-26", "2001-02-11", "2001-12-30",
    "2002-02-17", "2003-01-18", "2004-01-05", "2006-01-21",
    "2007-02-24", "2008-02-22", "2009-01-24", "2010-02-09",
    "2012-01-11", "2013-01-07", "2018-02-12", "2019-01-01",
])

# Load master panel
df = pd.read_parquet(ROOT / "data" / "processed" / "analysis_panel_v2.parquet")
df.index = pd.to_datetime(df.index.tz_localize(None) if df.index.tz else df.index)

# Non-SSW winter reference days
non_ssw = df[(df["is_winter"] == 1) & (df["ssw_within_15d"] == 0)]

# Compute event-level statistics from raw daily data
events = []
for sd in SSW_DATES:
    window = df[
        (df.index >= sd - pd.Timedelta(days=15))
        & (df.index <= sd + pd.Timedelta(days=15))
    ]
    obs = window["dry_natural_size_1234"].sum()

    # DOY-matched expected count (±3 day bandwidth)
    exp = 0
    for offset in range(31):
        d = (sd - pd.Timedelta(days=15) + pd.Timedelta(days=offset)).timetuple().tm_yday
        ref = non_ssw[
            (non_ssw.index.dayofyear >= d - 3) & (non_ssw.index.dayofyear <= d + 3)
        ]
        if len(ref) > 0:
            exp += ref["dry_natural_size_1234"].mean()

    rr = obs / exp if exp > 0 else np.nan
    log_rr = np.log(rr) if rr and rr > 0 else np.nan

    z500 = window["ncep_z500_nh"].mean()
    t50 = window["ncep_t_50hpa"].mean()
    slp = window["ncep_slp_nh"].mean()

    events.append({
        "date": sd.strftime("%Y-%m-%d"),
        "rr": rr, "log_rr": log_rr,
        "z500_nh_m": z500, "t50_K": t50, "slp_Pa": slp,
    })

edf = pd.DataFrame(events)
valid = edf.dropna(subset=["log_rr", "z500_nh_m", "t50_K"])

# Save event-level data for reproducibility
valid.to_csv(RESULTS / "event_level_dose_response.csv", index=False)

z500 = valid["z500_nh_m"].values
t50 = valid["t50_K"].values
log_rr = valid["log_rr"].values
n = len(valid)

# Regressions
slope_z, intercept_z, r_z, p_z, se_z = stats.linregress(z500, log_rr)
slope_t, intercept_t, r_t, p_t, se_t = stats.linregress(t50, log_rr)

# --- Figure ---
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)

# Panel A: Stratospheric T50 vs log(RR) — null relationship
ax = axes[0]
ax.scatter(t50, log_rr, c="#d62728", s=60, alpha=0.7,
           edgecolors="k", linewidths=0.5, zorder=3)
x_fit = np.linspace(t50.min() - 2, t50.max() + 2, 100)
ax.plot(x_fit, slope_t * x_fit + intercept_t, "--",
        color="#d62728", alpha=0.5, linewidth=1.5)
# Proper confidence band for Panel A too
resid_t = log_rr - (slope_t * t50 + intercept_t)
s_resid_t = np.sqrt(np.sum(resid_t**2) / (n - 2))
t50_bar = t50.mean()
Sxx_t = np.sum((t50 - t50_bar)**2)
se_mean_t = s_resid_t * np.sqrt(1.0 / n + (x_fit - t50_bar)**2 / Sxx_t)
t_crit_a = stats.t.ppf(0.975, df=n - 2)
ax.fill_between(x_fit,
                slope_t * x_fit + intercept_t - t_crit_a * se_mean_t,
                slope_t * x_fit + intercept_t + t_crit_a * se_mean_t,
                alpha=0.10, color="#d62728")
ax.axhline(0, color="grey", linestyle=":", linewidth=0.8)
ax.set_xlabel("Mean 50 hPa temperature during SSW window (K)", fontsize=10)
ax.set_ylabel("log(Rate Ratio)", fontsize=10)
ax.set_title("a  Stratospheric intensity", fontsize=11,
             fontweight="bold", loc="left")
ax.text(0.05, 0.95,
        f"$r$ = {r_t:.2f}, $R^2$ = {r_t**2:.3f}\n$P$ = {p_t:.2f}, $n$ = {n}",
        transform=ax.transAxes, fontsize=9, verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

# Panel B: Z500 vs log(RR) — significant relationship
ax = axes[1]
ax.scatter(z500, log_rr, c="#1f77b4", s=60, alpha=0.7,
           edgecolors="k", linewidths=0.5, zorder=3)
x_fit = np.linspace(z500.min() - 5, z500.max() + 5, 100)
y_fit = slope_z * x_fit + intercept_z
ax.plot(x_fit, y_fit, "-", color="#1f77b4", alpha=0.7, linewidth=2)
# Proper confidence band for mean response (narrowest at x̄, widens toward extremes)
resid = log_rr - (slope_z * z500 + intercept_z)
s_resid = np.sqrt(np.sum(resid**2) / (n - 2))
x_bar = z500.mean()
Sxx = np.sum((z500 - x_bar)**2)
se_mean = s_resid * np.sqrt(1.0 / n + (x_fit - x_bar)**2 / Sxx)
t_crit = stats.t.ppf(0.975, df=n - 2)
ax.fill_between(x_fit, y_fit - t_crit * se_mean, y_fit + t_crit * se_mean,
                alpha=0.15, color="#1f77b4")
ax.axhline(0, color="grey", linestyle=":", linewidth=0.8)
ax.set_xlabel("NH-mean Z500 during SSW window (m)", fontsize=10)
ax.set_title("b  Tropospheric circulation (Z500)", fontsize=11,
             fontweight="bold", loc="left")
ax.text(0.05, 0.95,
        f"$r$ = {r_z:.2f}, $R^2$ = {r_z**2:.3f}\n$P$ = {p_z:.3f}, $n$ = {n}",
        transform=ax.transAxes, fontsize=9, verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

fig.suptitle(
    "Event-level avalanche response tracks tropospheric, not stratospheric, intensity",
    fontsize=11, fontweight="bold", y=1.02)
plt.tight_layout()
fig.savefig(FIGURES / "extended_dose_response.pdf", bbox_inches="tight", dpi=300)
fig.savefig(FIGURES / "extended_dose_response.png", bbox_inches="tight", dpi=150)
plt.close()

print("Figure saved to data/figures/extended_dose_response.pdf")
print(f"\nPanel A (Stratospheric T50): r={r_t:.3f}, R²={r_t**2:.4f}, P={p_t:.3f}")
print(f"Panel B (Tropospheric Z500): r={r_z:.3f}, R²={r_z**2:.4f}, P={p_z:.3f}")
print(f"\nConclusion: Event-level sensitivity concentrated at tropospheric level")
print(f"  Z500 explains {r_z**2*100:.1f}% of variance vs {r_t**2*100:.1f}% for strat T")
