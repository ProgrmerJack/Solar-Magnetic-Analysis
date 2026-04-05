"""
ERA5 Download — year-by-year, small domains only to stay within CDS cost limits.

Scope (reduced to pass CDS free-tier quotas):
  1. Swiss Alps surface   (44-48N, 5-11E)  1959-2025 — avalanche met baseline
  2. Norwegian mountains  (58-72N, 5-30E)  1959-2025 — Varsom baseline
  3. Colorado Rockies     (36-42N,-109--104W) 1959-2025 — CAIC baseline
  4. European trop 500hPa (35-75N,-30-50E)  1979-2025 — blocking (MONTHLY MEANS)

Polar stratosphere is NOT requested here — NCEP/NCAR polar-cap daily data
(already downloaded) covers SSW detection needs at higher temporal resolution.

CDS limit: ~120,000 fields/request AND a per-day compute quota.
Each surface year = 9 vars × 4 times × 365 days = ~13,140 fields × small domain → OK.
Monthly means = 4 vars × 5 levels × 12 months × small-ish domain → OK.
"""
import cdsapi
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parents[1]))
from utils import get_logger

LOG = get_logger("38_era5")
OUT = Path(__file__).parents[2] / "data" / "atmospheric" / "era5"
OUT.mkdir(parents=True, exist_ok=True)

c = cdsapi.Client()

ALL_MONTHS = [f"{m:02d}" for m in range(1, 13)]
ALL_DAYS   = [f"{d:02d}" for d in range(1, 32)]

def retrieve_year(dataset, params_base, year, out_path):
    if out_path.exists() and out_path.stat().st_size > 5000:
        LOG.info("    skip: %s", out_path.name)
        return
    LOG.info("    submitting %d …", year)
    req = {**params_base, "year": str(year), "month": ALL_MONTHS, "day": ALL_DAYS}
    c.retrieve(dataset, req, str(out_path))
    LOG.info("    ✓ %s", out_path.name)

# ─── 1. Swiss Alps surface (1959-2025, 44-48N 5-11E = ~16×24 grid) ────────────
LOG.info("=== ERA5: Swiss Alps surface (1959-2025) ===")
alps_dir = OUT / "swiss_alps"
alps_dir.mkdir(exist_ok=True)
alps_params = {
    "product_type": "reanalysis",
    "variable": [
        "2m_temperature", "total_precipitation",
        "snowfall", "snow_depth",
        "10m_u_component_of_wind", "10m_v_component_of_wind",
        "surface_pressure", "mean_sea_level_pressure",
        "2m_dewpoint_temperature",
    ],
    "time": ["00:00", "06:00", "12:00", "18:00"],
    "area": [48, 5, 44, 11],        # N/W/S/E
    "data_format": "netcdf",
    "download_format": "unarchived",
}
for yr in range(1959, 2026):
    retrieve_year("reanalysis-era5-single-levels", alps_params, yr,
                  alps_dir / f"era5_swiss_alps_{yr}.nc")

# ─── 2. Norwegian mountains surface (1959-2025, 58-72N 5-30E = ~56×100 grid) ──
LOG.info("=== ERA5: Norwegian mountains surface (1959-2025) ===")
norway_dir = OUT / "norway"
norway_dir.mkdir(exist_ok=True)
norway_params = {
    "product_type": "reanalysis",
    "variable": [
        "2m_temperature", "total_precipitation",
        "snowfall", "snow_depth",
        "10m_u_component_of_wind", "10m_v_component_of_wind",
        "surface_pressure", "2m_dewpoint_temperature",
    ],
    "time": ["00:00", "06:00", "12:00", "18:00"],
    "area": [72, 5, 58, 30],
    "data_format": "netcdf",
    "download_format": "unarchived",
}
for yr in range(1959, 2026):
    retrieve_year("reanalysis-era5-single-levels", norway_params, yr,
                  norway_dir / f"era5_norway_{yr}.nc")

# ─── 3. Colorado Rockies surface (1959-2025, 36-42N -109--104W = ~24×20 grid) ─
LOG.info("=== ERA5: Colorado Rockies surface (1959-2025) ===")
colorado_dir = OUT / "colorado"
colorado_dir.mkdir(exist_ok=True)
colorado_params = {
    "product_type": "reanalysis",
    "variable": [
        "2m_temperature", "total_precipitation",
        "snowfall", "snow_depth",
        "10m_u_component_of_wind", "10m_v_component_of_wind",
        "surface_pressure",
    ],
    "time": ["00:00", "06:00", "12:00", "18:00"],
    "area": [42, -109, 36, -104],
    "data_format": "netcdf",
    "download_format": "unarchived",
}
for yr in range(1959, 2026):
    retrieve_year("reanalysis-era5-single-levels", colorado_params, yr,
                  colorado_dir / f"era5_colorado_{yr}.nc")

# ─── 4. European troposphere MONTHLY MEANS (1979-2025) — for blocking index ──
# Monthly means are tiny and stay well within quota.
LOG.info("=== ERA5: European trop monthly means (1979-2025) ===")
europe_dir = OUT / "europe_trop_monthly"
europe_dir.mkdir(exist_ok=True)
for yr in range(1979, 2026):
    out_path = europe_dir / f"era5_europe_trop_monthly_{yr}.nc"
    if out_path.exists() and out_path.stat().st_size > 5000:
        LOG.info("    skip: %s", out_path.name)
        continue
    LOG.info("    submitting Europe monthly %d …", yr)
    c.retrieve(
        "reanalysis-era5-pressure-levels-monthly-means",
        {
            "product_type": "monthly_averaged_reanalysis",
            "variable": ["geopotential", "temperature",
                         "u_component_of_wind", "v_component_of_wind"],
            "pressure_level": ["200", "300", "500", "700", "850"],
            "year": str(yr),
            "month": ALL_MONTHS,
            "time": "00:00",
            "area": [75, -30, 35, 50],
            "data_format": "netcdf",
            "download_format": "unarchived",
        },
        str(out_path)
    )
    LOG.info("    ✓ %s", out_path.name)

LOG.info("=== All ERA5 requests complete ===")
