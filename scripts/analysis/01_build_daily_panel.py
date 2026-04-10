"""
01_build_daily_panel.py — Construct the unified daily analysis panel
=====================================================================
Merges all data sources into a single daily panel for NDJFM winters.
Output: data/processed/analysis_panel.parquet
"""
import sys
from pathlib import Path
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _analysis_utils import PROCESSED, LOG, save_analysis

# ═══════════════════════════════════════════════════════════════════════════════
# LOAD RAW DATASETS
# ═══════════════════════════════════════════════════════════════════════════════

def load_kp() -> pd.DataFrame:
    """Load 3-hourly Kp and compute daily max Kp / max Ap."""
    kp = pd.read_parquet(PROCESSED / "geomagnetic" / "kp_index.parquet")
    kp.index = kp.index.tz_localize(None) if kp.index.tz else kp.index
    daily = kp.resample("D").agg({"kp": "max", "ap": "max", "ap_daily": "first"})
    daily.columns = ["kp_max", "ap_max", "ap_daily"]
    return daily


def load_dst() -> pd.DataFrame:
    """Load hourly Dst and compute daily min Dst (most disturbed)."""
    dst = pd.read_parquet(PROCESSED / "geomagnetic" / "dst_index.parquet")
    dst.index = dst.index.tz_localize(None) if dst.index.tz else dst.index
    daily = dst[["dst"]].resample("D").agg({"dst": "min"})
    daily.columns = ["dst_min"]
    return daily


def load_mls_ozone() -> pd.DataFrame:
    """Load MLS polar ozone at multiple pressure levels."""
    o3 = pd.read_parquet(PROCESSED / "atmospheric" / "mls_ozone_polar.parquet")
    o3.index = o3.index.tz_localize(None) if o3.index.tz else o3.index
    o3.columns = [f"mls_o3_{c}" for c in o3.columns]
    return o3


def load_mls_temp() -> pd.DataFrame:
    """Load MLS polar temperature at multiple pressure levels."""
    t = pd.read_parquet(PROCESSED / "atmospheric" / "mls_temperature_polar.parquet")
    t.index = t.index.tz_localize(None) if t.index.tz else t.index
    t.columns = [f"mls_t_{c}" for c in t.columns]
    return t


def load_ncep_strat() -> pd.DataFrame:
    """Load NCEP stratospheric data and return wide format with ncep_t/u/z_XXhpa columns.

    Handles both the legacy long format (source_file, level_hpa, polar_cap_mean)
    and the rebuilt wide format (air_K_10hPa, hgt_m_10hPa, uwnd_ms_10hPa, …).
    """
    ns = pd.read_parquet(PROCESSED / "atmospheric" / "ncep_stratosphere.parquet")
    ns.index = ns.index.tz_localize(None) if ns.index.tz else ns.index

    # ── New wide format produced by fix_ncep_redownload.py ───────────────────
    # Columns like: air_K_10hPa, air_K_20hPa, hgt_m_10hPa, uwnd_ms_10hPa …
    if "air_K_10hPa" in ns.columns or any(c.startswith("air_K_") for c in ns.columns):
        rename = {}
        for col in ns.columns:
            # air_K_10hPa → ncep_t_10hpa
            # hgt_m_10hPa → ncep_z_10hpa
            # uwnd_ms_10hPa → ncep_u_10hpa
            if col.startswith("air_K_"):
                lev = col.split("_")[-1].lower()          # "10hpa"
                rename[col] = f"ncep_t_{lev}"
            elif col.startswith("hgt_m_"):
                lev = col.split("_")[-1].lower()
                rename[col] = f"ncep_z_{lev}"
            elif col.startswith("uwnd_ms_"):
                lev = col.split("_")[-1].lower()
                rename[col] = f"ncep_u_{lev}"
        wide = ns.rename(columns=rename)
        wide.index.name = "date"
        return wide

    # ── Legacy long format (source_file, level_hpa, polar_cap_mean) ──────────
    def var_from_file(f):
        f = str(f).lower()
        if "air" in f:
            return "t"
        elif "uwnd" in f:
            return "u"
        elif "hgt" in f:
            return "z"
        return "unknown"

    ns["var"] = ns["source_file"].apply(var_from_file)
    ns["col"] = "ncep_" + ns["var"] + "_" + ns["level_hpa"].astype(int).astype(str) + "hpa"
    wide = ns.pivot_table(index=ns.index, columns="col", values="polar_cap_mean", aggfunc="first")
    wide.index.name = "date"
    return wide


def load_ncep_trop() -> pd.DataFrame:
    """Load NCEP tropospheric NH means.

    Handles both the legacy long format (source_file, nh_mean) and the
    current wide format (hgt_500hPa_m, slp_Pa, uwnd_850hPa_ms).
    """
    nt = pd.read_parquet(PROCESSED / "atmospheric" / "ncep_troposphere.parquet")
    nt.index = nt.index.tz_localize(None) if nt.index.tz else nt.index

    # ── New wide format ───────────────────────────────────────────────────────
    if "hgt_500hPa_m" in nt.columns:
        rename = {
            "hgt_500hPa_m":   "ncep_z500_nh",
            "slp_Pa":         "ncep_slp_nh",
            "uwnd_850hPa_ms": "ncep_u850_nh",
        }
        wide = nt.rename(columns={k: v for k, v in rename.items() if k in nt.columns})
        wide.index.name = "date"
        return wide

    # ── Legacy long format (source_file, nh_mean) ────────────────────────────
    def var_from_file(f):
        f = str(f).lower()
        if "hgt_500" in f:
            return "ncep_z500_nh"
        elif "slp" in f:
            return "ncep_slp_nh"
        elif "uwnd_850" in f:
            return "ncep_u850_nh"
        return "ncep_trop_other"

    nt["col"] = nt["source_file"].apply(var_from_file)
    wide = nt.pivot_table(index=nt.index, columns="col", values="nh_mean", aggfunc="first")
    wide.index.name = "date"
    return wide


def load_climate_indices() -> pd.DataFrame:
    """Load climate indices (monthly) and forward-fill to daily."""
    ci = pd.read_parquet(PROCESSED / "atmospheric" / "climate_indices.parquet")
    ci.index = ci.index.tz_localize(None) if ci.index.tz else ci.index

    # Select key indices
    cols = ["nao_daily", "nao_monthly", "qbo_u50", "qbo_u30_cpc",
            "mei_v2_bimonthly", "pdo_monthly", "pna_monthly", "amo_monthly"]
    cols = [c for c in cols if c in ci.columns]
    ci = ci[cols]

    # Resample to daily and forward-fill monthly values
    daily = ci.resample("D").first().ffill(limit=35)
    return daily


def load_slf_activity() -> pd.DataFrame:
    """Load SLF daily avalanche activity — the PRIMARY endpoint."""
    act = pd.read_parquet(PROCESSED / "cryosphere" / "slf_activity.parquet")
    act.index = act.index.tz_localize(None) if act.index.tz else act.index
    act.index.name = "date"

    # Key columns for analysis
    key_cols = [
        "aai_all", "aai_all_natural", "aai_all_human",
        "aai_all_dry", "aai_all_wet",
        "dry_natural_size_1234", "wet_natural_size_1234",
        "natural_size_234", "natural_size_1234",
        "size_1234", "max_size", "max_danger_corr",
    ]
    key_cols = [c for c in key_cols if c in act.columns]
    return act[key_cols]


def load_slf_accidents_daily() -> pd.DataFrame:
    """Load SLF accidents and aggregate to daily counts (negative control)."""
    acc = pd.read_parquet(PROCESSED / "cryosphere" / "slf_accidents.parquet")
    acc.index = acc.index.tz_localize(None) if acc.index.tz else acc.index
    acc.index.name = "date"
    daily = acc.groupby(acc.index.date).size().rename("accident_count")
    daily.index = pd.DatetimeIndex(daily.index, name="date")
    return daily.to_frame()


def load_flares_daily() -> pd.DataFrame:
    """Load flare catalog and aggregate to daily max class / count."""
    fl = pd.read_parquet(PROCESSED / "solar" / "flares.parquet")
    fl.index = fl.index.tz_localize(None) if fl.index.tz else fl.index

    # Map flare class to numeric scale
    def class_to_numeric(ct):
        if pd.isna(ct) or not isinstance(ct, str):
            return 0
        ct = ct.strip()
        cls = ct[0].upper()
        try:
            mag = float(ct[1:])
        except (ValueError, IndexError):
            mag = 1.0
        scale = {"A": 1, "B": 2, "C": 3, "M": 4, "X": 5}.get(cls, 0)
        return scale + mag / 10.0

    fl["class_numeric"] = fl["classType"].apply(class_to_numeric)
    daily = fl.resample("D").agg(
        flare_count=("classType", "count"),
        flare_max_class=("class_numeric", "max"),
    )
    return daily


def load_omni_daily() -> pd.DataFrame:
    """Load OMNI hourly solar wind and compute daily summaries."""
    om = pd.read_parquet(PROCESSED / "solar" / "omni_hourly.parquet")
    om.index = om.index.tz_localize(None) if om.index.tz else om.index

    # Replace fill values with NaN
    for col in om.columns:
        if om[col].dtype in [np.float64, np.float32]:
            om.loc[om[col] > 9999, col] = np.nan

    # Select key solar wind parameters
    cols_agg = {}
    if "Bz_GSM" in om.columns:
        cols_agg["sw_bz_min"] = ("Bz_GSM", "min")
    if "flow_speed" in om.columns:
        cols_agg["sw_speed_max"] = ("flow_speed", "max")
    elif "V" in om.columns:
        cols_agg["sw_speed_max"] = ("V", "max")
    if "proton_density" in om.columns:
        cols_agg["sw_density_max"] = ("proton_density", "max")

    if cols_agg:
        daily = om.resample("D").agg(**cols_agg)
    else:
        daily = pd.DataFrame(index=om.resample("D").first().index)
    return daily


def load_ssw_events() -> pd.DataFrame:
    """Load SSW catalog as binary daily indicator."""
    ssw = pd.read_parquet(PROCESSED / "atmospheric" / "ssw_catalog.parquet")
    ssw.index = ssw.index.tz_localize(None) if ssw.index.tz else ssw.index
    ssw["ssw_event"] = 1
    return ssw[["ssw_event"]]


# ═══════════════════════════════════════════════════════════════════════════════
# BUILD PANEL
# ═══════════════════════════════════════════════════════════════════════════════

def build_panel() -> pd.DataFrame:
    """Merge all data into a single daily panel."""
    LOG.info("Loading datasets...")

    kp = load_kp()
    dst = load_dst()
    mls_o3 = load_mls_ozone()
    mls_t = load_mls_temp()
    ncep_s = load_ncep_strat()
    ncep_t = load_ncep_trop()
    climate = load_climate_indices()
    slf = load_slf_activity()
    accidents = load_slf_accidents_daily()
    flares = load_flares_daily()
    omni = load_omni_daily()
    ssw = load_ssw_events()

    LOG.info("Data loaded. Building panel...")

    # Create daily date index spanning the full SLF activity period
    dates = pd.date_range("1998-11-01", "2019-06-01", freq="D", name="date")
    panel = pd.DataFrame(index=dates)

    # Merge all datasets
    for df, name in [
        (kp, "kp"), (dst, "dst"), (mls_o3, "mls_o3"),
        (mls_t, "mls_t"), (ncep_s, "ncep_strat"), (ncep_t, "ncep_trop"),
        (climate, "climate"), (slf, "slf"), (accidents, "accidents"),
        (flares, "flares"), (omni, "omni"),
    ]:
        before = len(panel.columns)
        panel = panel.join(df, how="left")
        LOG.info("  Merged %s: +%d cols", name, len(panel.columns) - before)

    # SSW: create proximity indicator (within N days of SSW onset)
    ssw_dates = ssw.index.values
    panel["ssw_within_15d"] = 0
    for sd in ssw_dates:
        mask = (panel.index >= sd - pd.Timedelta(days=15)) & \
               (panel.index <= sd + pd.Timedelta(days=15))
        panel.loc[mask, "ssw_within_15d"] = 1

    # ─── Derived columns ──────────────────────────────────────────────
    panel["month"] = panel.index.month
    panel["day_of_year"] = panel.index.dayofyear

    # Winter assignment: Nov-Mar → winter labeled by the January year
    def assign_winter(dt):
        if dt.month >= 11:
            return f"{dt.year}/{dt.year+1}"
        elif dt.month <= 3:
            return f"{dt.year-1}/{dt.year}"
        return None

    panel["winter_id"] = [assign_winter(d) for d in panel.index]

    # Day of season: days since Nov 1 of that winter
    def day_of_season(dt):
        if dt.month >= 11:
            return (dt - pd.Timestamp(dt.year, 11, 1)).days
        elif dt.month <= 3:
            return (dt - pd.Timestamp(dt.year - 1, 11, 1)).days
        return np.nan

    panel["day_of_season"] = [day_of_season(d) for d in panel.index]
    panel["day_of_season_sq"] = panel["day_of_season"] ** 2

    # Season flag
    panel["is_winter"] = panel["month"].isin([11, 12, 1, 2, 3]).astype(int)
    panel["is_summer"] = panel["month"].isin([5, 6, 7, 8, 9, 10]).astype(int)

    # ─── Geomagnetic event definition ─────────────────────────────────
    panel["geo_event_raw"] = (
        (panel["kp_max"] >= 5.0) | (panel["dst_min"] <= -50.0)
    ).astype(int)

    # Decluster: 10-day washout
    panel["geo_event"] = 0
    last_event = pd.Timestamp("1900-01-01")
    for i, (dt, row) in enumerate(panel.iterrows()):
        if row["geo_event_raw"] == 1 and (dt - last_event).days >= 10:
            panel.iloc[i, panel.columns.get_loc("geo_event")] = 1
            last_event = dt

    # Post-event windows
    event_dates = panel.index[panel["geo_event"] == 1]
    panel["post_event_1_3d"] = 0    # Fast pathway
    panel["post_event_5_21d"] = 0   # Stratospheric pathway
    panel["post_event_0_30d"] = 0   # Full window

    for ed in event_dates:
        m1 = (panel.index > ed + pd.Timedelta(days=0)) & \
             (panel.index <= ed + pd.Timedelta(days=3))
        m2 = (panel.index >= ed + pd.Timedelta(days=5)) & \
             (panel.index <= ed + pd.Timedelta(days=21))
        m3 = (panel.index >= ed) & \
             (panel.index <= ed + pd.Timedelta(days=30))
        panel.loc[m1, "post_event_1_3d"] = 1
        panel.loc[m2, "post_event_5_21d"] = 1
        panel.loc[m3, "post_event_0_30d"] = 1

    # ─── Summary ──────────────────────────────────────────────────────
    n_events = panel["geo_event"].sum()
    n_winter_days = panel["is_winter"].sum()
    n_with_slf = panel.loc[panel["is_winter"] == 1, "aai_all_natural"].notna().sum()

    LOG.info("Panel built: %d rows, %d cols", len(panel), len(panel.columns))
    LOG.info("Geomagnetic events (declustered): %d", n_events)
    LOG.info("Winter days: %d, with SLF data: %d", n_winter_days, n_with_slf)

    return panel


def main():
    import logging
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")

    panel = build_panel()

    # Save
    out = PROCESSED / "analysis_panel.parquet"
    save_analysis(panel, out, {
        "title": "Daily analysis panel for solar-avalanche study",
        "created_by": "01_build_daily_panel.py",
    })

    # Print summary statistics
    winter = panel[panel["is_winter"] == 1]
    print(f"\n{'='*60}")
    print("PANEL SUMMARY (Winter NDJFM only)")
    print(f"{'='*60}")
    print(f"Total winter days: {len(winter)}")
    print(f"Winters: {winter['winter_id'].nunique()}")
    print(f"Date range: {winter.index.min().date()} to {winter.index.max().date()}")
    print(f"\nGeomagnetic events in winter: {winter['geo_event'].sum()}")
    print(f"Days in post-event 5-21d window: {winter['post_event_5_21d'].sum()}")
    print(f"Days in fast 1-3d window: {winter['post_event_1_3d'].sum()}")
    print(f"\nSLF natural activity coverage: "
          f"{winter['aai_all_natural'].notna().sum()} / {len(winter)} "
          f"({100*winter['aai_all_natural'].notna().mean():.1f}%)")
    print(f"Mean daily natural AAI: {winter['aai_all_natural'].mean():.2f}")
    print(f"Max daily natural AAI: {winter['aai_all_natural'].max():.1f}")
    print(f"\nColumns: {len(panel.columns)}")


if __name__ == "__main__":
    main()
