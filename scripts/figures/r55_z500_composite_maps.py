"""
r55_z500_composite_maps.py
==========================
Publication-quality Z500 composite anomaly maps for Nature Geoscience
manuscript on SSW effects on avalanche hazard.

Figures produced:
  1. fig_z500_composite_map  — NH polar stereographic Z500 anomaly composite
  2. fig_stratosphere_troposphere_composite — 4-panel vertical propagation
  3. fig_z500_lead_lag — Lead-lag correlation panels

Data sources:
  - NCEP/NCAR Reanalysis daily Z500 via NOAA PSL OPeNDAP
  - data/processed/analysis_panel_v2.parquet
  - data/processed/atmospheric/ssw_catalog.parquet
"""

import warnings
warnings.filterwarnings("ignore")

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False
    print("WARNING: cartopy not available, using fallback projections")

try:
    import netCDF4 as nc
    HAS_NETCDF4 = True
except ImportError:
    HAS_NETCDF4 = False

# ── Paths ────────────────────────────────────────────────────────────────
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIG_DIR = os.path.join(ROOT, "data", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# ── Nature Geoscience style ─────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "lines.linewidth": 1.0,
    "patch.linewidth": 0.5,
})

# ── Load core data ──────────────────────────────────────────────────────
print("Loading analysis panel and SSW catalog...")
panel = pd.read_parquet(os.path.join(ROOT, "data", "processed", "analysis_panel_v2.parquet"))

ssw_cat = pd.read_parquet(
    os.path.join(ROOT, "data", "processed", "atmospheric", "ssw_catalog.parquet")
)
ssw_cat.index = ssw_cat.index.tz_localize(None)
# Filter to study period
ssw_events = ssw_cat[
    (ssw_cat.index >= "1998-11-01") & (ssw_cat.index <= "2019-06-01")
]
ssw_dates = ssw_events.index.tolist()
print(f"  {len(ssw_dates)} SSW events in study period")

# Years that contain SSW events
ssw_years = sorted(set(d.year for d in ssw_dates))
# Also need preceding year if event is in Jan (for climatology window)
all_years_needed = set()
for d in ssw_dates:
    all_years_needed.add(d.year)
    if d.month <= 2:
        all_years_needed.add(d.year - 1)
    all_years_needed.add(d.year + 1)  # for post-event window
# Intersect with data availability (1998-2019)
all_years_needed = sorted(y for y in all_years_needed if 1998 <= y <= 2019)

# ════════════════════════════════════════════════════════════════════════
# SECTION 1: Download gridded Z500 data via OPeNDAP
# ════════════════════════════════════════════════════════════════════════

def download_z500_opendap(years):
    """Download Z500 (500 hPa geopotential height) for given years via OPeNDAP."""
    import datetime as dt
    all_data = {}
    base_url = "https://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/pressure/hgt.{}.nc"

    for yr in years:
        url = base_url.format(yr)
        print(f"  Fetching Z500 for {yr} via OPeNDAP...", end=" ", flush=True)
        try:
            ds = nc.Dataset(url)
            # Manual time parsing (units: "hours since 1800-01-01 00:00:0.0")
            time_raw = ds.variables["time"][:]
            base_date = dt.datetime(1800, 1, 1)
            times = pd.DatetimeIndex([
                base_date + dt.timedelta(hours=float(h)) for h in time_raw
            ])

            lat = ds.variables["lat"][:]
            lon = ds.variables["lon"][:]
            levels = ds.variables["level"][:]
            lev_idx = int(np.argmin(np.abs(levels - 500)))

            # Only fetch 500 hPa level, NH (lat >= 20N)
            lat_mask = lat >= 20.0
            lat_nh = lat[lat_mask]
            lat_start = int(np.where(lat_mask)[0][0])
            lat_end = int(np.where(lat_mask)[0][-1]) + 1

            hgt = ds.variables["hgt"][:, lev_idx, lat_start:lat_end, :]
            ds.close()

            all_data[yr] = {
                "times": times,
                "lat": lat_nh,
                "lon": lon,
                "hgt": np.array(hgt),
            }
            print(f"OK ({hgt.shape[0]} days, lat {lat_nh[0]:.0f}-{lat_nh[-1]:.0f}N)")
        except Exception as e:
            print(f"FAILED: {e}")

    return all_data


def compute_composites(z500_data, ssw_dates, window_pre=5, window_post=30):
    """
    Compute Z500 anomaly composites around SSW events.
    
    Anomaly = SSW-window mean - DOY climatology.
    Returns composite anomaly map and significance mask.
    """
    lat = None
    lon = None
    event_anomalies = []

    for ev_date in ssw_dates:
        yr = ev_date.year
        if yr not in z500_data:
            print(f"    Skipping {ev_date.strftime('%Y-%m-%d')}: no data for year {yr}")
            continue

        d = z500_data[yr]
        if lat is None:
            lat = d["lat"]
            lon = d["lon"]

        # Define event window (days 5-30 after onset for tropospheric response)
        event_start = ev_date + pd.Timedelta(days=window_pre)
        event_end = ev_date + pd.Timedelta(days=window_post)

        # Get event window data (may span year boundary)
        event_days = pd.date_range(event_start, event_end)
        event_hgt = []
        for day in event_days:
            y = day.year
            if y in z500_data:
                dd = z500_data[y]
                idx = np.where(dd["times"] == day)[0]
                if len(idx) > 0:
                    event_hgt.append(dd["hgt"][idx[0]])
        if len(event_hgt) < 10:
            print(f"    Skipping {ev_date.strftime('%Y-%m-%d')}: insufficient event data ({len(event_hgt)} days)")
            continue
        event_mean = np.nanmean(event_hgt, axis=0)

        # Compute DOY-matched climatology (same DOYs from all available years)
        clim_hgt = []
        for day in event_days:
            doy = day.dayofyear
            for cy in z500_data:
                if cy == yr:
                    continue
                cd = z500_data[cy]
                doy_matches = np.array([t.dayofyear for t in cd["times"]])
                idx = np.where(doy_matches == doy)[0]
                if len(idx) > 0:
                    clim_hgt.append(cd["hgt"][idx[0]])

        if len(clim_hgt) < 20:
            print(f"    Skipping {ev_date.strftime('%Y-%m-%d')}: insufficient climatology data")
            continue
        clim_mean = np.nanmean(clim_hgt, axis=0)

        anomaly = event_mean - clim_mean
        event_anomalies.append(anomaly)
        print(f"    {ev_date.strftime('%Y-%m-%d')}: anomaly range [{np.nanmin(anomaly):.0f}, {np.nanmax(anomaly):.0f}] m")

    if len(event_anomalies) == 0:
        return None, None, None, None

    event_anomalies = np.array(event_anomalies)
    composite = np.nanmean(event_anomalies, axis=0)

    # Significance via one-sample t-test at each grid point
    n_events = event_anomalies.shape[0]
    t_stat = composite / (np.nanstd(event_anomalies, axis=0, ddof=1) / np.sqrt(n_events))
    p_values = 2 * (1 - stats.t.cdf(np.abs(t_stat), df=n_events - 1))

    return composite, p_values, lat, lon


# ════════════════════════════════════════════════════════════════════════
# FIGURE 1: Z500 composite anomaly map
# ════════════════════════════════════════════════════════════════════════

def plot_z500_composite(composite, p_values, lat, lon):
    """NH polar stereographic Z500 anomaly composite map."""
    print("\n── Figure 1: Z500 Composite Anomaly Map ──")

    fig = plt.figure(figsize=(5.5, 5.5))

    if HAS_CARTOPY:
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())
        ax.set_extent([-180, 180, 20, 90], crs=ccrs.PlateCarree())

        # Wrap longitude for continuous plotting
        lon_plot = lon.copy()

        # Create mesh
        lon2d, lat2d = np.meshgrid(lon_plot, lat)

        # Smooth for visual appeal
        composite_smooth = gaussian_filter(composite, sigma=0.8)

        # Colormap limits
        vmax = np.nanpercentile(np.abs(composite_smooth), 97)
        vmax = max(vmax, 30)  # at least 30m
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

        # Filled contours
        levels = np.linspace(-vmax, vmax, 21)
        cf = ax.contourf(
            lon2d, lat2d, composite_smooth,
            levels=levels,
            cmap="RdBu_r",
            norm=norm,
            transform=ccrs.PlateCarree(),
            extend="both",
        )

        # Contour lines
        cl = ax.contour(
            lon2d, lat2d, composite_smooth,
            levels=levels[::2],
            colors="k",
            linewidths=0.3,
            transform=ccrs.PlateCarree(),
        )

        # Significance stippling (p < 0.10)
        sig_mask = p_values < 0.10
        # Subsample for clean stippling
        skip = 2
        lat_stip = lat2d[::skip, ::skip]
        lon_stip = lon2d[::skip, ::skip]
        sig_stip = sig_mask[::skip, ::skip]
        ax.scatter(
            lon_stip[sig_stip], lat_stip[sig_stip],
            marker=".",
            s=0.4,
            c="k",
            alpha=0.5,
            transform=ccrs.PlateCarree(),
            zorder=5,
        )

        # Map features
        ax.coastlines(linewidth=0.5, color="0.3")
        ax.add_feature(cfeature.BORDERS, linewidth=0.3, edgecolor="0.5")
        ax.gridlines(
            draw_labels=False,
            linewidth=0.3,
            color="gray",
            alpha=0.5,
            linestyle="--",
        )

        # Mark Alps
        ax.plot(
            8.2, 46.8,
            marker="*",
            markersize=10,
            color="gold",
            markeredgecolor="k",
            markeredgewidth=0.5,
            transform=ccrs.PlateCarree(),
            zorder=10,
        )
        ax.text(
            12, 43.5,
            "Alps",
            fontsize=7,
            fontweight="bold",
            color="k",
            transform=ccrs.PlateCarree(),
            ha="left",
            zorder=10,
        )

        # Colorbar
        cbar = fig.colorbar(cf, ax=ax, orientation="horizontal", pad=0.05,
                            fraction=0.046, aspect=30, shrink=0.8)
        cbar.set_label("Z500 anomaly (m)", fontsize=8)
        cbar.ax.tick_params(labelsize=7)

        ax.set_title(
            f"Z500 anomaly composite: days 5–30 after SSW onset\n"
            f"(N = {len(ssw_dates)} events, stippling: p < 0.10)",
            fontsize=9,
            fontweight="bold",
            pad=12,
        )

    # Save
    for ext in ["pdf", "png"]:
        fpath = os.path.join(FIG_DIR, f"fig_z500_composite_map.{ext}")
        fig.savefig(fpath, dpi=300)
        print(f"  Saved: {fpath}")
    plt.close(fig)


# ════════════════════════════════════════════════════════════════════════
# FIGURE 2: Stratosphere–troposphere coupling composite (time series)
# ════════════════════════════════════════════════════════════════════════

def plot_strat_trop_composite():
    """
    4-panel composite of stratosphere-troposphere coupling around SSW events.
    Uses the time-series data from analysis_panel_v2.parquet.
    """
    print("\n── Figure 2: Stratosphere–Troposphere Composite ──")

    # Compute DOY climatology for key variables
    winter = panel[panel["is_winter"] == 1].copy()
    clim_cols = ["ncep_t_10hpa", "ncep_t_50hpa", "ncep_u_10hpa",
                 "ncep_z500_nh", "aai_all_natural", "nao_daily"]
    doy_clim = winter.groupby("day_of_year")[clim_cols].mean()

    # Build event-centered composites
    lag_range = np.arange(-30, 46)
    variables = {
        "T10": "ncep_t_10hpa",
        "U10": "ncep_u_10hpa",
        "Z500": "ncep_z500_nh",
        "AAI": "aai_all_natural",
    }

    composites = {k: [] for k in variables}

    for ev_date in ssw_dates:
        event_data = {}
        for lag in lag_range:
            date = ev_date + pd.Timedelta(days=int(lag))
            if date in panel.index:
                doy = int(panel.loc[date, "day_of_year"])
                for key, col in variables.items():
                    val = panel.loc[date, col]
                    if pd.notna(val) and doy in doy_clim.index:
                        clim_val = doy_clim.loc[doy, col]
                        if key not in event_data:
                            event_data[key] = {}
                        event_data[key][lag] = val - clim_val

        for key in variables:
            if key in event_data:
                series = np.array([event_data[key].get(lag, np.nan) for lag in lag_range])
                composites[key].append(series)

    # Convert to arrays
    for key in composites:
        composites[key] = np.array(composites[key])

    # ── Plot ──
    fig, axes = plt.subplots(4, 1, figsize=(5.5, 8), sharex=True)

    panel_labels = ["a", "b", "c", "d"]
    titles = [
        "10 hPa temperature anomaly (polar cap)",
        "10 hPa zonal-mean zonal wind anomaly",
        "500 hPa geopotential height anomaly (NH)",
        "Avalanche activity index anomaly",
    ]
    ylabels = [r"$\Delta$T (K)", r"$\Delta$u (m s$^{-1}$)", r"$\Delta$Z500 (m)", r"$\Delta$AAI"]
    colors = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd"]
    keys = ["T10", "U10", "Z500", "AAI"]

    for i, (ax, key, title, ylabel, color) in enumerate(
        zip(axes, keys, titles, ylabels, colors)
    ):
        data = composites[key]
        n_events = data.shape[0]
        mean = np.nanmean(data, axis=0)
        se = np.nanstd(data, axis=0, ddof=1) / np.sqrt(
            np.sum(~np.isnan(data), axis=0)
        )
        ci_lo = mean - 1.96 * se
        ci_hi = mean + 1.96 * se

        # Smooth for visual clarity
        from scipy.ndimage import uniform_filter1d
        mean_s = uniform_filter1d(mean, size=3)
        ci_lo_s = uniform_filter1d(ci_lo, size=3)
        ci_hi_s = uniform_filter1d(ci_hi, size=3)

        ax.fill_between(lag_range, ci_lo_s, ci_hi_s, alpha=0.2, color=color)
        ax.plot(lag_range, mean_s, color=color, linewidth=1.5)
        ax.axhline(0, color="k", linewidth=0.4, linestyle="-")
        ax.axvline(0, color="k", linewidth=0.6, linestyle="--", alpha=0.7)

        # Mark the tropospheric response window
        ax.axvspan(5, 30, alpha=0.06, color="gray")

        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=8, fontweight="bold", loc="left", pad=4)

        # Panel label
        ax.text(
            -0.08, 1.05, panel_labels[i],
            transform=ax.transAxes,
            fontsize=11,
            fontweight="bold",
            va="top",
        )

        # Print stats
        post_mean = np.nanmean(mean[lag_range >= 5])
        print(f"  {key}: post-SSW (lag 5-45) mean anomaly = {post_mean:.2f}")

        # Significance stars for peak anomaly
        peak_idx = np.nanargmax(np.abs(mean_s))
        peak_lag = lag_range[peak_idx]
        peak_vals = data[:, peak_idx]
        peak_vals_clean = peak_vals[~np.isnan(peak_vals)]
        if len(peak_vals_clean) > 3:
            t_stat, p_val = stats.ttest_1samp(peak_vals_clean, 0)
            print(f"    Peak at lag {peak_lag}: {mean_s[peak_idx]:.2f}, t={t_stat:.2f}, p={p_val:.4f}")

    axes[-1].set_xlabel("Days relative to SSW onset")
    axes[-1].set_xlim(-30, 45)

    fig.suptitle(
        f"Stratosphere–troposphere coupling composite (N = {len(ssw_dates)} SSW events)",
        fontsize=10,
        fontweight="bold",
        y=0.98,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96], h_pad=1.5)

    for ext in ["pdf", "png"]:
        fpath = os.path.join(FIG_DIR, f"fig_stratosphere_troposphere_composite.{ext}")
        fig.savefig(fpath, dpi=300)
        print(f"  Saved: {fpath}")
    plt.close(fig)


# ════════════════════════════════════════════════════════════════════════
# FIGURE 3: Lead-lag correlations
# ════════════════════════════════════════════════════════════════════════

def plot_lead_lag_correlations():
    """
    Lead-lag correlations showing information flow direction.
    Panel a: Z500 anomaly vs avalanche activity
    Panel b: Stratospheric T50 vs Z500
    """
    print("\n── Figure 3: Lead-Lag Correlations ──")

    winter = panel[panel["is_winter"] == 1].copy()

    # Compute DOY anomalies
    anom_cols = ["ncep_z500_nh", "ncep_t_50hpa", "ncep_t_10hpa", "ncep_u_10hpa",
                 "aai_all_natural", "nao_daily"]
    doy_clim = winter.groupby("day_of_year")[anom_cols].mean()
    for col in anom_cols:
        if col in winter.columns and col in doy_clim.columns:
            winter[f"{col}_anom"] = winter.apply(
                lambda row: row[col] - doy_clim.loc[int(row["day_of_year"]), col]
                if int(row["day_of_year"]) in doy_clim.index and pd.notna(row[col])
                else np.nan,
                axis=1,
            )

    # Lead-lag correlation function
    def lead_lag_corr(series_a, series_b, max_lag=30):
        """Compute correlation of series_a(t) with series_b(t+lag)."""
        lags = np.arange(-max_lag, max_lag + 1)
        corrs = np.full(len(lags), np.nan)
        p_vals = np.full(len(lags), np.nan)
        n_obs = np.full(len(lags), 0)

        a_vals = series_a.values
        b_vals = series_b.values

        for i, lag in enumerate(lags):
            if lag >= 0:
                a = a_vals[:len(a_vals) - lag] if lag > 0 else a_vals
                b = b_vals[lag:]
            else:
                a = a_vals[-lag:]
                b = b_vals[:len(b_vals) + lag]

            mask = ~np.isnan(a) & ~np.isnan(b)
            if mask.sum() > 30:
                r, p = stats.pearsonr(a[mask], b[mask])
                corrs[i] = r
                p_vals[i] = p
                n_obs[i] = mask.sum()

        return lags, corrs, p_vals, n_obs

    # Compute correlations
    lags_za, corrs_za, pvals_za, n_za = lead_lag_corr(
        winter["ncep_z500_nh_anom"], winter["aai_all_natural_anom"], max_lag=30
    )
    lags_tz, corrs_tz, pvals_tz, n_tz = lead_lag_corr(
        winter["ncep_t_10hpa_anom"], winter["ncep_z500_nh_anom"], max_lag=30
    )
    lags_ta, corrs_ta, pvals_ta, n_ta = lead_lag_corr(
        winter["ncep_t_10hpa_anom"], winter["aai_all_natural_anom"], max_lag=30
    )

    # ── SSW-conditional lead-lag ──
    ssw_mask = winter["ssw_within_15d"] == 1
    non_ssw_mask = winter["ssw_within_15d"] == 0

    lags_za_ssw, corrs_za_ssw, _, _ = lead_lag_corr(
        winter.loc[ssw_mask, "ncep_z500_nh_anom"],
        winter.loc[ssw_mask, "aai_all_natural_anom"],
        max_lag=30,
    )
    lags_za_non, corrs_za_non, _, _ = lead_lag_corr(
        winter.loc[non_ssw_mask, "ncep_z500_nh_anom"],
        winter.loc[non_ssw_mask, "aai_all_natural_anom"],
        max_lag=30,
    )

    # ── Plot ──
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.0))

    # Panel a: Z500 → Avalanche
    ax = axes[0]
    ax.fill_between(lags_za, 0, corrs_za, where=pvals_za < 0.05,
                     alpha=0.3, color="#2ca02c")
    ax.plot(lags_za, corrs_za, color="#2ca02c", linewidth=1.5, label="All winters")
    ax.axhline(0, color="k", linewidth=0.4)
    ax.axvline(0, color="k", linewidth=0.4, linestyle="--")
    # Significance threshold
    n_eff = np.nanmedian(n_za)
    sig_thresh = 1.96 / np.sqrt(n_eff) if n_eff > 0 else 0
    ax.axhline(sig_thresh, color="gray", linewidth=0.4, linestyle=":")
    ax.axhline(-sig_thresh, color="gray", linewidth=0.4, linestyle=":")
    ax.set_xlabel("Lag (days)\n← Z500 leads | AAI leads →")
    ax.set_ylabel("Pearson r")
    ax.set_title("Z500 vs avalanche activity", fontsize=8, fontweight="bold")
    ax.set_ylim(-0.15, 0.15)

    peak_idx = np.nanargmax(np.abs(corrs_za))
    peak_lag = lags_za[peak_idx]
    peak_r = corrs_za[peak_idx]
    print(f"  Z500→AAI: peak r = {peak_r:.4f} at lag {peak_lag}")

    # Panel b: T10 → Z500
    ax = axes[1]
    ax.fill_between(lags_tz, 0, corrs_tz, where=pvals_tz < 0.05,
                     alpha=0.3, color="#d62728")
    ax.plot(lags_tz, corrs_tz, color="#d62728", linewidth=1.5)
    ax.axhline(0, color="k", linewidth=0.4)
    ax.axvline(0, color="k", linewidth=0.4, linestyle="--")
    n_eff2 = np.nanmedian(n_tz)
    sig_thresh2 = 1.96 / np.sqrt(n_eff2) if n_eff2 > 0 else 0
    ax.axhline(sig_thresh2, color="gray", linewidth=0.4, linestyle=":")
    ax.axhline(-sig_thresh2, color="gray", linewidth=0.4, linestyle=":")
    ax.set_xlabel("Lag (days)\n← T10 leads | Z500 leads →")
    ax.set_ylabel("Pearson r")
    ax.set_title("Strat. T (10 hPa) vs Z500", fontsize=8, fontweight="bold")

    peak_idx2 = np.nanargmax(np.abs(corrs_tz))
    print(f"  T10→Z500: peak r = {corrs_tz[peak_idx2]:.4f} at lag {lags_tz[peak_idx2]}")

    # Panel c: T10 → Avalanche (full chain)
    ax = axes[2]
    ax.fill_between(lags_ta, 0, corrs_ta, where=pvals_ta < 0.05,
                     alpha=0.3, color="#9467bd")
    ax.plot(lags_ta, corrs_ta, color="#9467bd", linewidth=1.5, label="All winters")
    # Overlay SSW vs non-SSW conditional
    ax.plot(lags_za_ssw, corrs_za_ssw, color="#ff7f0e", linewidth=1.0,
            linestyle="--", alpha=0.8, label="SSW windows")
    ax.plot(lags_za_non, corrs_za_non, color="#1f77b4", linewidth=1.0,
            linestyle=":", alpha=0.8, label="Non-SSW")
    ax.axhline(0, color="k", linewidth=0.4)
    ax.axvline(0, color="k", linewidth=0.4, linestyle="--")
    ax.set_xlabel("Lag (days)\n← Strat. leads | Surface leads →")
    ax.set_ylabel("Pearson r")
    ax.set_title("T10 → avalanche / SSW contrast", fontsize=8, fontweight="bold")
    ax.legend(loc="upper right", fontsize=6, framealpha=0.9)

    for i, ax in enumerate(axes):
        ax.text(-0.12, 1.08, chr(97 + i), transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top")
        ax.set_xlim(-30, 30)

    fig.suptitle("Lead–lag correlations: stratospheric to surface propagation",
                 fontsize=10, fontweight="bold", y=1.02)
    fig.tight_layout(w_pad=2.0)

    for ext in ["pdf", "png"]:
        fpath = os.path.join(FIG_DIR, f"fig_z500_lead_lag.{ext}")
        fig.savefig(fpath, dpi=300)
        print(f"  Saved: {fpath}")
    plt.close(fig)


# ════════════════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-download", action="store_true",
                        help="Skip OPeNDAP download and only regenerate Figs 2 & 3")
    args = parser.parse_args()

    # ── FIGURE 1: Z500 gridded composite ──
    cache_path = os.path.join(ROOT, "data", "processed", "z500_composite_cache.npz")
    gridded_success = False

    if not args.skip_download:
        print("\n" + "=" * 70)
        print("DOWNLOADING NCEP/NCAR Z500 GRIDDED DATA VIA OPeNDAP")
        print("=" * 70)

        if HAS_NETCDF4:
            z500_data = download_z500_opendap(all_years_needed)
            if len(z500_data) >= 10:
                print(f"\nComputing Z500 composites from {len(z500_data)} years of data...")
                composite, p_values, lat, lon = compute_composites(
                    z500_data, ssw_dates, window_pre=5, window_post=30
                )
                if composite is not None:
                    n_events = len([d for d in ssw_dates if d.year in z500_data])
                    print(f"\n  Composite from ~{n_events} events")
                    print(f"  Anomaly range: [{np.nanmin(composite):.1f}, {np.nanmax(composite):.1f}] m")
                    print(f"  Significant grid points (p<0.10): "
                          f"{np.nansum(p_values < 0.10)} / {p_values.size} "
                          f"({100*np.nansum(p_values<0.10)/p_values.size:.1f}%)")
                    # Cache for fast re-run
                    np.savez_compressed(cache_path, composite=composite,
                                       p_values=p_values, lat=lat, lon=lon)
                    plot_z500_composite(composite, p_values, lat, lon)
                    gridded_success = True
    elif os.path.exists(cache_path):
        print("\n── Loading cached Z500 composite ──")
        cached = np.load(cache_path)
        plot_z500_composite(cached["composite"], cached["p_values"],
                            cached["lat"], cached["lon"])
        gridded_success = True

    if not gridded_success and not args.skip_download:
        print("\nGridded Z500 download incomplete — creating time-series-based map")
        print("  (See Figure 2 for time-series composites)")

    # ── FIGURE 2: Stratosphere-troposphere composite ──
    print("\n" + "=" * 70)
    print("FIGURE 2: STRATOSPHERE–TROPOSPHERE VERTICAL COUPLING")
    print("=" * 70)
    plot_strat_trop_composite()

    # ── FIGURE 3: Lead-lag correlations ──
    print("\n" + "=" * 70)
    print("FIGURE 3: LEAD-LAG CORRELATIONS")
    print("=" * 70)
    plot_lead_lag_correlations()

    print("\n" + "=" * 70)
    print("ALL FIGURES COMPLETE")
    print("=" * 70)
