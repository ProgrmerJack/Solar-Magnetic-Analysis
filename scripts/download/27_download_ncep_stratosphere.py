"""
NCEP/NCAR Reanalysis — Stratospheric Polar Cap Data
Downloads daily-mean temperature (air), zonal wind (uwnd), and geopotential
height (hgt) at stratospheric pressure levels (10–100 hPa) averaged over the
polar cap (60–90 °N).

Data are free via PSL THREDDS OPeNDAP — no registration required.
Uses xarray with OPeNDAP for lazy remote access + efficient subsetting.
Output: merged multi-year CSV per variable with columns: date, level_hPa, polar_cap_mean.

NCEP pressure levels: [1000,925,850,700,600,500,400,300,250,200,150,100,70,50,30,20,10] hPa
Stratospheric levels (10–100 hPa): [100, 70, 50, 30, 20, 10]
10 hPa is the canonical SSW detection level (zonal wind reversal at 60°N, 10 hPa).
"""
import sys
import time
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger

logger = get_logger("ncep_strat")
OUT = DATA_DIR / "atmospheric" / "ncep_stratosphere"
OUT.mkdir(parents=True, exist_ok=True)

THREDDS = "https://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/pressure"

VARIABLES = ["air", "uwnd", "hgt"]
YEARS = range(1979, 2025)
STRAT_MAX_HPA = 100.0   # include 100 hPa and above (lower pressure = higher altitude)
POLAR_MIN_LAT = 60.0    # polar cap: 60–90 °N


def process_year(var: str, year: int) -> pd.DataFrame | None:
    """Open one year via OPeNDAP, subset stratosphere + polar cap, return daily means."""
    url = f"{THREDDS}/{var}.{year}.nc"
    try:
        ds = xr.open_dataset(url, engine="netcdf4")
        da = ds[var]
        # Subset stratospheric levels (pressure <= 100 hPa in NCEP convention)
        strat = da.sel(level=da.level[da.level <= STRAT_MAX_HPA])
        # Subset polar cap (lat 60–90 N; NCEP lat goes 90→-90)
        polar = strat.sel(lat=strat.lat[strat.lat >= POLAR_MIN_LAT])
        # Area-weighted mean over lat/lon (cos-lat weighting)
        weights = np.cos(np.deg2rad(polar.lat))
        polar_mean = polar.weighted(weights).mean(dim=["lat", "lon"]).compute()
        ds.close()
    except Exception as exc:
        logger.warning(f"  {var} {year}: {exc}")
        return None

    # Convert to long-format DataFrame
    times = pd.to_datetime(polar_mean.time.values)
    levels = polar_mean.level.values.astype(int)
    records = []
    for li, lev in enumerate(levels):
        for ti, t in enumerate(times):
            records.append({
                "date": t.date(),
                "level_hPa": int(lev),
                "polar_cap_mean": float(polar_mean.values[ti, li]),
            })
    return pd.DataFrame(records)


logger.info("=== NCEP/NCAR Stratospheric Polar-Cap Download ===")
logger.info(f"Variables: {VARIABLES}  |  Years: {YEARS[0]}–{YEARS[-1]}")
logger.info(f"Output: {OUT}")

for var in VARIABLES:
    merged_file = OUT / f"{var}_strat_polarcap_1979_2024.csv"
    if merged_file.exists():
        logger.info(f"  {var}: merged file already exists, skipping")
        continue

    var_out = OUT / var
    var_out.mkdir(exist_ok=True)
    all_dfs = []

    for year in YEARS:
        out_file = var_out / f"{var}_{year}.csv"
        if out_file.exists():
            all_dfs.append(pd.read_csv(out_file, parse_dates=["date"]))
            continue

        logger.info(f"  {var} {year} …")
        df = process_year(var, year)
        if df is not None and len(df) > 0:
            df.to_csv(out_file, index=False)
            all_dfs.append(df)
            logger.info(f"    {year}: {len(df):,} records ({df['level_hPa'].nunique()} levels)")
        else:
            logger.warning(f"    {year}: no data returned")
        time.sleep(0.3)

    if all_dfs:
        merged = pd.concat(all_dfs, ignore_index=True).sort_values(["date", "level_hPa"])
        merged["date"] = pd.to_datetime(merged["date"])
        merged.to_csv(merged_file, index=False)
        logger.info(f"  {var} merged → {merged_file.name} ({len(merged):,} records)")

logger.info("=== NCEP stratospheric download complete ===")
