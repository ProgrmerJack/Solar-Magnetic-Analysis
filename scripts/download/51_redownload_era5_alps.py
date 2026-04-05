"""
Re-download ERA5 Swiss Alps surface data for the analysis period (2004-2025).
Uses CDS API. Downloads year-by-year to data/atmospheric/era5/swiss_alps/.
Only downloads files not already on disk.
"""
import cdsapi
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parents[1]))

OUT = Path(__file__).parents[2] / "data" / "atmospheric" / "era5"
ALPS_DIR = OUT / "swiss_alps"
ALPS_DIR.mkdir(parents=True, exist_ok=True)

c = cdsapi.Client(
    url="https://cds.climate.copernicus.eu/api",
    key="a9746760-8d00-43f6-95f2-6d1f9a0d8b47",
    quiet=False,
)

ALL_MONTHS = [f"{m:02d}" for m in range(1, 13)]
ALL_DAYS   = [f"{d:02d}" for d in range(1, 32)]

ALPS_PARAMS = {
    "product_type": "reanalysis",
    "variable": [
        "2m_temperature", "total_precipitation",
        "snowfall", "snow_depth",
        "10m_u_component_of_wind", "10m_v_component_of_wind",
        "surface_pressure", "mean_sea_level_pressure",
        "2m_dewpoint_temperature",
    ],
    "time": ["00:00", "06:00", "12:00", "18:00"],
    "area": [48, 5, 44, 11],   # N/W/S/E  (Swiss Alps: 44-48N, 5-11E)
    "data_format": "netcdf",
    "download_format": "unarchived",
}

# Download 2004-2025 only (matching Aura/MLS + avalanche analysis period)
YEARS = list(range(2004, 2026))

print(f"ERA5 Swiss Alps re-download: {len(YEARS)} years ({YEARS[0]}-{YEARS[-1]})")
print(f"Output: {ALPS_DIR}")
total_written = 0

for yr in YEARS:
    out_path = ALPS_DIR / f"era5_swiss_alps_{yr}.nc"
    if out_path.exists() and out_path.stat().st_size > 5000:
        print(f"  SKIP {yr}: already on disk ({out_path.stat().st_size/1e6:.1f} MB)")
        continue
    print(f"  Submitting {yr} ...")
    try:
        req = {**ALPS_PARAMS, "year": str(yr), "month": ALL_MONTHS, "day": ALL_DAYS}
        c.retrieve("reanalysis-era5-single-levels", req, str(out_path))
        sz = out_path.stat().st_size / 1e6
        print(f"  ✓ {yr}: {sz:.1f} MB")
        total_written += 1
    except Exception as e:
        print(f"  FAIL {yr}: {e}")

print(f"\nDone. Downloaded {total_written} new year files.")
