"""
ERA5 Polar Stratosphere Monthly Means (1979-2025)

Downloads polar cap (60-90N) monthly mean pressure-level data for the
stratospheric mediation test (EPP -> NOx -> polar vortex -> surface).

Variables:
  T          = temperature (K)             -- SSW warm stratosphere signal
  U          = zonal wind (m/s)            -- polar vortex strength
  V          = meridional wind (m/s)       -- wave-driven circulation
  geopotential = height (m^2/s^2)          -- 10hPa Z for SSW detection
  ozone_mass_mixing_ratio = O3             -- ozone depletion diagnostic

Pressure levels: 1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100 hPa (stratosphere)
Area: 60-90N, 0-360E (full polar cap longitude)
Time: monthly means, 1979-2025 (covers all 38 Butler et al. SSW events)

Output: data/atmospheric/era5/polar_strat_monthly/era5_polar_strat_YYYY.nc
Estimated size: ~2-5 MB/year × 47 years = ~100-250 MB total
"""
import cdsapi
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parents[1]))
from utils import get_logger

LOG = get_logger("47_era5_polar_strat")
OUT = Path(__file__).parents[2] / "data" / "atmospheric" / "era5" / "polar_strat_monthly"
OUT.mkdir(parents=True, exist_ok=True)

c = cdsapi.Client()

ALL_MONTHS = [f"{m:02d}" for m in range(1, 13)]

LOG.info("=== ERA5 Polar Stratosphere Monthly Means (1979-2025) ===")
LOG.info("Area: 60-90N full longitude, levels: 1-100 hPa")

for yr in range(1979, 2026):
    out_path = OUT / f"era5_polar_strat_{yr}.nc"
    if out_path.exists() and out_path.stat().st_size > 5000:
        LOG.info("  skip: %s", out_path.name)
        continue
    LOG.info("  submitting %d ...", yr)
    c.retrieve(
        "reanalysis-era5-pressure-levels-monthly-means",
        {
            "product_type": "monthly_averaged_reanalysis",
            "variable": [
                "temperature",
                "u_component_of_wind",
                "v_component_of_wind",
                "geopotential",
                "ozone_mass_mixing_ratio",
            ],
            "pressure_level": ["1", "2", "3", "5", "7", "10",
                                "20", "30", "50", "70", "100"],
            "year": str(yr),
            "month": ALL_MONTHS,
            "time": "00:00",
            "area": [90, 0, 60, 360],   # N/W/S/E
            "data_format": "netcdf",
            "download_format": "unarchived",
        },
        str(out_path)
    )
    size_kb = out_path.stat().st_size / 1024 if out_path.exists() else 0
    LOG.info("  ✓ %s (%.0f KB)", out_path.name, size_kb)

LOG.info("=== ERA5 polar stratosphere complete ===")
