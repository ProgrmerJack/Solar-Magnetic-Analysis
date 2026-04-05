"""
NCEP/NCAR Reanalysis — Tropospheric Variables for Atmospheric Bridge Analysis
Downloads daily-mean 500 hPa geopotential height and sea-level pressure (SLP)
for the Northern Hemisphere mid-latitudes and polar cap (30–90 °N).

These data characterize blocking patterns, planetary wave activity, and the
tropospheric response to stratospheric sudden warmings — the final link in the
solar → stratosphere → troposphere → snowpack mechanism chain.

Run AFTER 27_download_ncep_stratosphere.py completes.
"""
import sys
import time
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger

logger = get_logger("ncep_tropo")
OUT = DATA_DIR / "atmospheric" / "ncep_troposphere"
OUT.mkdir(parents=True, exist_ok=True)

THREDDS = "https://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/pressure"
SLP_THREDDS = "https://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/surface"

YEARS = range(1979, 2025)
NH_MIN_LAT = 30.0   # Northern Hemisphere region: 30–90 °N

CONFIG = {
    # Pressure-level variables at 500 hPa (blocking patterns / planetary waves)
    "hgt_500hPa": {
        "var": "hgt", "level": 500.0, "thredds": THREDDS,
        "desc": "Geopotential height at 500 hPa (NH 30-90N)"
    },
    # SLP is in the surface directory
    "slp": {
        "var": "slp", "level": None, "thredds": SLP_THREDDS,
        "desc": "Sea-level pressure NH daily mean"
    },
    # Zonal wind at 850 hPa (surface circulation, warm air advection)
    "uwnd_850hPa": {
        "var": "uwnd", "level": 850.0, "thredds": THREDDS,
        "desc": "Zonal wind at 850 hPa (NH 30-90N)"
    },
}


def process_year(cfg: dict, year: int) -> pd.DataFrame | None:
    """Download one year, subset NH region, return daily spatial mean."""
    url = f"{cfg['thredds']}/{cfg['var']}.{year}.nc"
    try:
        ds = xr.open_dataset(url, engine="netcdf4")
        da = ds[cfg["var"]]

        if cfg["level"] is not None:
            da = da.sel(level=cfg["level"], method="nearest")

        # Subset NH (30–90 N)
        da_nh = da.sel(lat=da.lat[da.lat >= NH_MIN_LAT])
        weights = np.cos(np.deg2rad(da_nh.lat))
        nh_mean = da_nh.weighted(weights).mean(dim=["lat", "lon"]).compute()
        ds.close()
    except Exception as exc:
        logger.warning(f"  {cfg['var']} {year}: {exc}")
        return None

    times = pd.to_datetime(nh_mean.time.values)
    vals = nh_mean.values
    df = pd.DataFrame({"date": [t.date() for t in times], "nh_mean": vals.tolist()})
    return df


logger.info("=== NCEP/NCAR Tropospheric NH Download ===")
logger.info(f"Output: {OUT}")

for name, cfg in CONFIG.items():
    merged_file = OUT / f"{name}_NH_1979_2024.csv"
    if merged_file.exists():
        logger.info(f"  {name}: merged file already exists, skipping")
        continue

    var_out = OUT / name
    var_out.mkdir(exist_ok=True)
    all_dfs = []

    for year in YEARS:
        out_file = var_out / f"{name}_{year}.csv"
        if out_file.exists():
            all_dfs.append(pd.read_csv(out_file, parse_dates=["date"]))
            continue

        logger.info(f"  {name} {year} …")
        df = process_year(cfg, year)
        if df is not None and len(df) > 0:
            df.to_csv(out_file, index=False)
            all_dfs.append(df)
            logger.info(f"    {year}: {len(df)} records")
        else:
            logger.warning(f"    {year}: no data")
        time.sleep(0.3)

    if all_dfs:
        merged = pd.concat(all_dfs, ignore_index=True)
        merged["date"] = pd.to_datetime(merged["date"])
        merged = merged.sort_values("date")
        merged.to_csv(merged_file, index=False)
        logger.info(f"  {name} merged → {merged_file.name} ({len(merged):,} records)")

logger.info("=== NCEP tropospheric download complete ===")
