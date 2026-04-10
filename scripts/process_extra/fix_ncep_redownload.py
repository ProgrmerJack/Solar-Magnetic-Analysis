"""
fix_ncep_redownload.py
Re-downloads NCEP/NCAR Reanalysis I stratospheric polar-cap data from PSL THREDDS
one year at a time and rebuilds ncep_stratosphere.parquet with correct values.

Also reshapes ncep_troposphere.parquet from long (mixed nh_mean) to proper wide format.

Run from repo root:  python scripts/fix_ncep_redownload.py

ROOT CAUSE of all-zeros bug:
  Requesting all N time steps at once via OPeNDAP returns zeros silently.
  FIX: fetch in 30-day batches; each batch returns correct data in ~3 seconds.
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).parent / "process"))
from _utils import PROCESSED_ROOT, save_parquet, LOG, setup_logging

THREDDS = "https://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/pressure"

ATM_OUT = PROCESSED_ROOT / "atmospheric"
POLAR_MIN_LAT = 60.0
YEARS = range(1979, 2025)
BATCH_DAYS = 30   # max safe OPeNDAP time-slice; >~365 returns all zeros


# ---------------------------------------------------------------------------
# Stratosphere: air temperature, geopotential height, zonal wind
# ---------------------------------------------------------------------------
STRAT_VARS = {
    "air":  "air_K",
    "hgt":  "hgt_m",
    "uwnd": "uwnd_ms",
}


def download_strat_year(var: str, year: int) -> pd.DataFrame | None:
    """
    Download one year via OPeNDAP in BATCH_DAYS-day chunks.
    Requesting all time steps at once silently returns zeros (OPeNDAP server limit).
    """
    url = f"{THREDDS}/{var}.{year}.nc"
    try:
        ds = xr.open_dataset(url, engine="netcdf4")
        da = ds[var]

        lats = da["lat"].values
        levs = da["level"].values
        n_times = da.sizes["time"]
        times = pd.DatetimeIndex(pd.to_datetime(da["time"].values), tz="UTC")

        strat_idx = np.where(levs <= 100.0)[0]
        polar_idx = np.where(lats >= POLAR_MIN_LAT)[0]
        if len(strat_idx) == 0 or len(polar_idx) == 0:
            LOG.warning("  %s %d: no strat/polar grid points", var, year)
            ds.close()
            return None

        strat_levs = levs[strat_idx]
        polar_lats = lats[polar_idx]
        weights = np.cos(np.deg2rad(polar_lats))

        # Fetch in BATCH_DAYS chunks to avoid the all-zeros OPeNDAP bug
        all_polar_means = []
        for t_start in range(0, n_times, BATCH_DAYS):
            t_end = min(t_start + BATCH_DAYS, n_times)
            chunk = da.isel(
                time=slice(t_start, t_end),
                level=strat_idx.tolist(),
                lat=polar_idx.tolist(),
            ).values   # shape: (chunk_days, n_strat_levels, n_polar_lats, n_lon)

            lon_mean = chunk.mean(axis=3)                              # (days, lev, lat)
            polar_mean = np.average(lon_mean, axis=2, weights=weights)  # (days, lev)
            all_polar_means.append(polar_mean)

        ds.close()

        polar_mean_year = np.concatenate(all_polar_means, axis=0)    # (n_times, n_lev)

        col_prefix = STRAT_VARS[var]
        rows = {f"{col_prefix}_{int(lv)}hPa": polar_mean_year[:, li].tolist()
                for li, lv in enumerate(strat_levs)}
        rows["time"] = times
        return pd.DataFrame(rows).set_index("time")

    except Exception as exc:
        LOG.warning("  %s %d: %s", var, year, exc)
        return None


def rebuild_ncep_stratosphere() -> None:
    out = ATM_OUT / "ncep_stratosphere.parquet"

    # Collect one DataFrame per year across all 3 variables
    year_dfs: dict[int, list[pd.DataFrame]] = {}

    for var in STRAT_VARS:
        LOG.info("Downloading NCEP strat var=%s (1979-2024)...", var)
        for year in YEARS:
            df = download_strat_year(var, year)
            if df is not None:
                year_dfs.setdefault(year, []).append(df)
                LOG.info("  %s %d: %d rows, %d cols", var, year, len(df), len(df.columns))
            else:
                LOG.warning("  %s %d: no data", var, year)
            time.sleep(0.3)   # be polite to PSL server

    if not year_dfs:
        LOG.error("No stratosphere data downloaded — aborting")
        return

    # Merge variables for each year by joining on time index
    all_years: list[pd.DataFrame] = []
    for year in sorted(year_dfs):
        frames = year_dfs[year]
        merged = frames[0]
        for frame in frames[1:]:
            merged = merged.join(frame, how="outer")
        all_years.append(merged)

    full = pd.concat(all_years).sort_index()

    meta = {
        "title": "NCEP/NCAR Reanalysis I — Stratospheric Polar-Cap Daily Means (≥60°N)",
        "source": "NOAA/NCEP-NCAR Reanalysis I via PSL HTTP fileServer (direct NC download)",
        "references": "https://doi.org/10.1175/1520-0477(1996)077<0437:TNYRP>2.0.CO;2",
        "units": "air_K=Kelvin; hgt_m=geopotential meters; uwnd_ms=m/s; polar cap >=60N cos-lat weighted mean",
        "levels": "10,20,30,50,70,100 hPa",
        "time_range": f"{full.index.min()} / {full.index.max()}",
    }
    # Delete the all-zero stub before writing
    if out.exists():
        out.unlink()
    save_parquet(full, out, meta)
    LOG.info("Rebuilt ncep_stratosphere.parquet: %d rows × %d cols", len(full), len(full.columns))
    LOG.info("Sample (first row):\n%s", full.iloc[0].to_dict())


# ---------------------------------------------------------------------------
# Troposphere: reshape from mixed long → proper wide format
# ---------------------------------------------------------------------------
def fix_ncep_troposphere() -> None:
    out = ATM_OUT / "ncep_troposphere.parquet"
    if not out.exists():
        LOG.warning("ncep_troposphere.parquet not found — skipping reshape")
        return

    df = pd.read_parquet(out)
    if "source_file" not in df.columns:
        LOG.info("ncep_troposphere.parquet already in expected format — skipping")
        return

    # Map source file names → clean column names and units
    src_map = {
        "hgt_500hPa_NH_1979_2024.csv":  ("hgt_500hPa_m",   "geopotential meters"),
        "slp_NH_1979_2024.csv":         ("slp_Pa",          "Pa"),
        "uwnd_850hPa_NH_1979_2024.csv": ("uwnd_850hPa_ms",  "m s-1"),
    }

    pivoted_parts = []
    for src_file, (col_name, _) in src_map.items():
        subset = df[df["source_file"] == src_file][["nh_mean"]].rename(columns={"nh_mean": col_name})
        pivoted_parts.append(subset)

    if not pivoted_parts:
        LOG.warning("No recognised source files in ncep_troposphere.parquet")
        return

    wide = pivoted_parts[0].join(pivoted_parts[1:], how="outer").sort_index()
    wide.index.name = "time"

    meta = {
        "title": "NCEP/NCAR Reanalysis I — Tropospheric NH Daily Means (30–90°N)",
        "source": "NOAA/NCEP-NCAR Reanalysis I via PSL THREDDS OPeNDAP",
        "references": "https://doi.org/10.1175/1520-0477(1996)077<0437:TNYRP>2.0.CO;2",
        "units": "hgt_500hPa_m=geopotential meters; slp_Pa=Pascals; uwnd_850hPa_ms=m/s",
        "region": "NH area-weighted mean 30–90°N",
    }
    if out.exists():
        out.unlink()
    save_parquet(wide, out, meta)
    LOG.info("Reshaped ncep_troposphere.parquet: %d rows × %d cols", len(wide), len(wide.columns))
    LOG.info("Columns: %s", list(wide.columns))
    LOG.info("Sample:\n%s", wide.head(2))


if __name__ == "__main__":
    setup_logging()
    LOG.info("=== NCEP data fix — stratosphere re-download + troposphere reshape ===")
    ATM_OUT.mkdir(parents=True, exist_ok=True)

    LOG.info("--- Step 1: Fix troposphere reshape (quick) ---")
    fix_ncep_troposphere()

    LOG.info("--- Step 2: Re-download stratosphere from THREDDS ---")
    rebuild_ncep_stratosphere()

    LOG.info("=== Done ===")
