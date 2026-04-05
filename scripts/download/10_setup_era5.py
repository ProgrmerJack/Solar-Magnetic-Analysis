"""
Script 10 — ERA5 Reanalysis Setup
Configures Copernicus CDS API access for ERA5 downloads.

ERA5 (ECMWF Reanalysis 5th generation):
  - Coverage: 1940–present (global, hourly)
  - Resolution: 0.25° (~31 km), 37 pressure levels + surface
  - Key variables: Temperature, wind (u/v), geopotential, vorticity,
    specific humidity, total precipitation, snow depth, sea-level pressure
  - Access: https://cds.climate.copernicus.eu/

REGISTRATION REQUIRED (free):
  1. Create account at https://cds.climate.copernicus.eu/user/register
  2. Accept terms for ERA5 datasets
  3. Get your API key from https://cds.climate.copernicus.eu/profile
"""
import sys
import os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger, write_instructions

logger = get_logger("10_era5_setup")
OUT = DATA_DIR / "atmospheric" / "era5"
OUT.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------- #
# CDS API configuration template                                               #
# --------------------------------------------------------------------------- #
CDS_RC_TEMPLATE = """url: https://cds.climate.copernicus.eu/api/v2
key: {UID}:{API_KEY}
verify: 1
"""

cds_rc_path = Path.home() / ".cdsapirc"
if not cds_rc_path.exists():
    logger.info("  Creating template ~/.cdsapirc — FILL IN YOUR CREDENTIALS")
    cds_rc_path.write_text(
        "url: https://cds.climate.copernicus.eu/api/v2\n"
        "key: YOUR_UID:YOUR_API_KEY\n"
        "verify: 1\n",
        encoding="utf-8"
    )
else:
    logger.info(f"  ~/.cdsapirc already exists")

# --------------------------------------------------------------------------- #
# ERA5 download script (run AFTER setting up .cdsapirc)                       #
# --------------------------------------------------------------------------- #
ERA5_DOWNLOAD_SCRIPT = '''\
"""ERA5 Download Script — run after configuring ~/.cdsapirc"""
import cdsapi
from pathlib import Path

OUT = Path(__file__).parents[2] / "data" / "atmospheric" / "era5"
OUT.mkdir(parents=True, exist_ok=True)

c = cdsapi.Client()

# ── 1. Polar stratospheric temperature (60-90N, 1-100 hPa, daily) ──────────
# Key for SSW detection and EPP-NOx analysis
c.retrieve(
    "reanalysis-era5-pressure-levels",
    {
        "product_type": "reanalysis",
        "variable": ["temperature", "u_component_of_wind",
                     "v_component_of_wind", "geopotential"],
        "pressure_level": ["1", "2", "3", "5", "7", "10", "20", "30",
                           "50", "70", "100", "150", "200", "250", "300"],
        "year":  [str(y) for y in range(1979, 2026)],
        "month": [f"{m:02d}" for m in range(1, 13)],
        "day":   [f"{d:02d}" for d in range(1, 32)],
        "time":  ["00:00", "12:00"],
        "area":  [90, -180, 60, 180],  # N/W/S/E — Arctic polar cap
        "format": "netcdf",
    },
    str(OUT / "era5_polar_strat_60N90N_1979_2025.nc")
)

# ── 2. North Atlantic / European troposphere (for teleconnection analysis) ──
c.retrieve(
    "reanalysis-era5-pressure-levels",
    {
        "product_type": "reanalysis",
        "variable": ["geopotential", "temperature",
                     "u_component_of_wind", "v_component_of_wind"],
        "pressure_level": ["200", "300", "500", "700", "850", "925", "1000"],
        "year":  [str(y) for y in range(1979, 2026)],
        "month": [f"{m:02d}" for m in range(1, 13)],
        "day":   [f"{d:02d}" for d in range(1, 32)],
        "time":  ["00:00", "12:00"],
        "area":  [75, -30, 35, 50],   # Europe + North Atlantic
        "format": "netcdf",
    },
    str(OUT / "era5_europe_trop_1979_2025.nc")
)

# ── 3. Alpine surface (Switzerland + Alps, for avalanche-weather link) ──────
c.retrieve(
    "reanalysis-era5-single-levels",
    {
        "product_type": "reanalysis",
        "variable": ["2m_temperature", "total_precipitation",
                     "snowfall", "snow_depth",
                     "10m_u_component_of_wind", "10m_v_component_of_wind",
                     "surface_pressure", "mean_sea_level_pressure"],
        "year":  [str(y) for y in range(1959, 2026)],
        "month": [f"{m:02d}" for m in range(1, 13)],
        "day":   [f"{d:02d}" for d in range(1, 32)],
        "time":  ["00:00", "06:00", "12:00", "18:00"],
        "area":  [48, 5, 44, 11],    # Swiss Alps bounding box
        "format": "netcdf",
    },
    str(OUT / "era5_swiss_alps_surface_1959_2025.nc")
)

# ── 4. Colorado Rockies surface (for CAIC cross-validation) ─────────────────
c.retrieve(
    "reanalysis-era5-single-levels",
    {
        "product_type": "reanalysis",
        "variable": ["2m_temperature", "total_precipitation",
                     "snowfall", "snow_depth"],
        "year":  [str(y) for y in range(1959, 2026)],
        "month": [f"{m:02d}" for m in range(1, 13)],
        "day":   [f"{d:02d}" for d in range(1, 32)],
        "time":  ["00:00", "06:00", "12:00", "18:00"],
        "area":  [42, -109, 36, -104],  # Colorado Rockies
        "format": "netcdf",
    },
    str(OUT / "era5_colorado_rockies_surface_1959_2025.nc")
)

print("ERA5 downloads complete.")
'''

era5_script_path = OUT / "download_era5.py"
if not era5_script_path.exists():
    era5_script_path.write_text(ERA5_DOWNLOAD_SCRIPT, encoding="utf-8")
    logger.info(f"  ✓  ERA5 download script written → {era5_script_path}")

write_instructions(
    OUT,
    "ERA5 Reanalysis — Setup & Download Instructions",
    f"""
STEP 1: Register at https://cds.climate.copernicus.eu/user/register (free)
STEP 2: Accept terms for the datasets at:
          https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels
          https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels
STEP 3: Get your API key at https://cds.climate.copernicus.eu/profile
STEP 4: Edit ~/.cdsapirc and replace YOUR_UID:YOUR_API_KEY with your credentials
STEP 5: pip install cdsapi
STEP 6: python {era5_script_path}

WHAT WILL BE DOWNLOADED (~100-500 GB total, downloads run on CDS servers):
  • Polar stratospheric temperature/wind/GPH (60-90N, 1-300 hPa, 1979-2025)
    → For SSW detection and polar vortex analysis
  • European troposphere (1979-2025)
    → For blocking pattern and teleconnection analysis  
  • Swiss Alps surface (snowfall, temp, precip, 1959-2025)
    → For avalanche meteorological baseline
  • Colorado Rockies surface (1959-2025)
    → For CAIC cross-validation

NOTE: ERA5 is the gold-standard dataset. MERRA-2 is a good alternative for
faster access (NASA Earthdata, see 11_setup_merra2.py).
"""
)

logger.info("  Run 'pip install cdsapi' then configure ~/.cdsapirc to enable ERA5 download.")
logger.info("=== Script 10 complete (setup only — credentials required) ===")
