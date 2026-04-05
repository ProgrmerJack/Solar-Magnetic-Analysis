"""
seNorge Norwegian Snow Water Equivalent (SWE) — Area-Averaged Time Series
Downloads annual SWE from MET Norway's seNorge gridded product (1957–2025)
via OPeNDAP. Extracts national and latitudinal zone means rather than full
grids (200 MB/year), giving a compact daily time series for analysis.

Data: MET Norway seNorge snow model, 1km grid, daily, 1957–2025
License: Norwegian Licence for Open Government Data (NLOD) / CC-BY 4.0
Attribution: Norwegian Meteorological Institute (MET Norway)
Thredds: https://thredds.met.no/thredds/catalog/senorge/seNorge_snow/swe/
"""
import sys
import time
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger

logger = get_logger("senorge_swe")
OUT = DATA_DIR / "cryosphere" / "ngi_norway" / "senorge_swe"
OUT.mkdir(parents=True, exist_ok=True)

DODSBASE = "https://thredds.met.no/thredds/dodsC/senorge/seNorge_snow/swe"
YEARS = range(1957, 2026)

# Norway latitudinal zones (for avalanche region analysis)
# seNorge UTM Zone 33 → lat/lon available as coords
LAT_ZONES = {
    "south": (57.0, 63.0),   # Sørlandet, Vestlandet, Telemark (main ski/avalanche terrain)
    "central": (63.0, 67.0), # Trøndelag, Møre og Romsdal
    "north": (67.0, 72.0),   # Nordland, Troms, Finnmark
}


def process_year(year: int) -> pd.DataFrame | None:
    """Extract daily national mean SWE via OPeNDAP with stride sampling.
    
    Samples every 10th pixel (stride=10) to reduce download ~100x while
    still accurately representing the national/zonal mean snowpack state.
    Grid: 1550 y × 1195 x → 155 y × 120 x sampled pixels per timestep.
    """
    url = f"{DODSBASE}/swe_{year}.nc"
    try:
        # Open without loading data (lazy)
        ds = xr.open_dataset(url, engine="netcdf4")
        swe = ds["snow_water_equivalent"]
        lat = ds["lat"]

        # Stride-sample every 20th pixel (reduces download by ~400x; ~18s/year vs 300s full)
        stride = 20
        swe_sub = swe[:, ::stride, ::stride]   # (time, y_sub, x_sub)
        lat_sub = lat[::stride, ::stride]       # (y_sub, x_sub)

        # Fetch subsampled data in one go (~2 MB vs 200 MB per year)
        swe_vals = swe_sub.values               # triggers OPeNDAP fetch
        lat_vals = lat_sub.values
        times = pd.to_datetime(ds["time"].values)

        nat_mean = np.nanmean(swe_vals, axis=(1, 2))

        zone_vals = {}
        for zone, (lat_min, lat_max) in LAT_ZONES.items():
            mask = (lat_vals >= lat_min) & (lat_vals <= lat_max)
            if mask.any():
                zone_vals[zone] = np.nanmean(
                    np.where(mask[np.newaxis, :, :], swe_vals, np.nan),
                    axis=(1, 2)
                )
            else:
                zone_vals[zone] = np.full(len(times), np.nan)

        ds.close()

        records = []
        for i, t in enumerate(times):
            records.append({
                "date": t.date(),
                "swe_national_mean_mm": float(nat_mean[i]),
                "swe_south_mean_mm": float(zone_vals["south"][i]),
                "swe_central_mean_mm": float(zone_vals["central"][i]),
                "swe_north_mean_mm": float(zone_vals["north"][i]),
            })
        return pd.DataFrame(records)

    except Exception as exc:
        logger.warning(f"  {year}: {exc}")
        return None


logger.info("=== seNorge Norwegian SWE Daily Time Series ===")
logger.info(f"Years: {YEARS[0]}–{YEARS[-1]}  |  Output: {OUT}")

all_dfs = []
merged_file = OUT / "norway_swe_daily_1957_2025.csv"

if merged_file.exists():
    logger.info(f"Merged file exists: {merged_file.name}")
else:
    for year in YEARS:
        out_file = OUT / f"swe_norway_{year}.csv"
        if out_file.exists():
            all_dfs.append(pd.read_csv(out_file, parse_dates=["date"]))
            continue

        logger.info(f"  {year} …")
        df = process_year(year)
        if df is not None and len(df) > 0:
            df.to_csv(out_file, index=False)
            all_dfs.append(df)
            logger.info(f"    {year}: {len(df)} days, "
                       f"max SWE={df['swe_national_mean_mm'].max():.1f} mm")
        else:
            logger.warning(f"    {year}: no data")
        time.sleep(0.5)

    if all_dfs:
        merged = pd.concat(all_dfs, ignore_index=True)
        merged["date"] = pd.to_datetime(merged["date"])
        merged = merged.sort_values("date").reset_index(drop=True)
        merged.to_csv(merged_file, index=False)
        logger.info(f"\nMerged {len(merged):,} records → {merged_file.name}")

logger.info("=== seNorge SWE download complete ===")
