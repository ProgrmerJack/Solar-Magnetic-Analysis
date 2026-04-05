"""
Script 02 — OMNIWeb Combined Solar Wind & Geomagnetic Indices
Downloads the NASA GSFC OMNIWeb datasets from SPDF:
  • OMNI2 low-resolution (hourly, 1963 – present)  ~55 MB
  • OMNI high-resolution 5-min (1995 – present)     yearly ASCII files
  • OMNI high-resolution 1-min (1995 – present)     yearly ASCII files

OMNI contains: solar wind (speed, density, Bz, proton temp), Kp, Dst, AE/AL/AU,
plasma beta, Alfvénic Mach, sunspot number, F10.7 and more — all in one place.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger, download_file

logger = get_logger("02_omni")
OUT = DATA_DIR / "atmospheric" / "omni_solar_wind"
OUT.mkdir(parents=True, exist_ok=True)

SPDF = "https://spdf.gsfc.nasa.gov/pub/data/omni"

# --------------------------------------------------------------------------- #
# 1. OMNI2 low-res yearly files (hourly, 1963-present) — correct URL scheme  #
# --------------------------------------------------------------------------- #
logger.info("=== OMNI2 Low-resolution (hourly, 1963–present) ===")
LR_OUT = OUT / "low_res"
LR_OUT.mkdir(exist_ok=True)

# Format description first
download_file(
    f"{SPDF}/low_res_omni/omni2.text",
    LR_OUT / "omni2_format.txt",
    desc="OMNI2 format description"
)

# Yearly files (1963–present) — confirmed path: /pub/data/omni/low_res_omni/omni2_YYYY.dat
for year in range(1963, 2026):
    download_file(
        f"{SPDF}/low_res_omni/omni2_{year}.dat",
        LR_OUT / f"omni2_{year}.dat",
        desc=f"OMNI2 hourly {year}"
    )

# --------------------------------------------------------------------------- #
# 2. OMNI 5-min high-resolution (1995–present), yearly ASCII                  #
# --------------------------------------------------------------------------- #
logger.info("=== OMNI High-resolution 5-min (1995–2025) ===")
HR5_OUT = OUT / "high_res_5min"
HR5_OUT.mkdir(exist_ok=True)

for year in range(1995, 2026):
    url = f"{SPDF}/high_res_omni/omni_5min{year}.asc"
    download_file(url, HR5_OUT / f"omni_5min_{year}.asc", desc=f"OMNI 5min {year}")

# --------------------------------------------------------------------------- #
# 3. OMNI 1-min high-resolution (1995–present), yearly ASCII                  #
# --------------------------------------------------------------------------- #
logger.info("=== OMNI High-resolution 1-min (1995–2025) ===")
HR1_OUT = OUT / "high_res_1min"
HR1_OUT.mkdir(exist_ok=True)

for year in range(1995, 2026):
    url = f"{SPDF}/high_res_omni/omni_min{year}.asc"
    download_file(url, HR1_OUT / f"omni_1min_{year}.asc", desc=f"OMNI 1min {year}")

# --------------------------------------------------------------------------- #
# 4. README / format documentation                                             #
# --------------------------------------------------------------------------- #
logger.info("=== OMNI documentation ===")
docs = {
    "high_res_omni_readme.txt": f"{SPDF}/high_res_omni/hro_format.txt",
    "omni2_readme.txt":         f"{SPDF}/low_res_omni/omni2.text",
}
for fname, url in docs.items():
    download_file(url, OUT / fname, desc=f"OMNI doc {fname}")

logger.info("=== Script 02 complete ===")
