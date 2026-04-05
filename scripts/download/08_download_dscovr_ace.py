"""
Script 08 — NOAA DSCOVR / ACE Solar Wind Data
Downloads real-time and archived solar wind data from DSCOVR (2016-present)
and ACE archive via NOAA NCEI and NASA CDAWeb.

DSCOVR (Deep Space Climate Observatory, L1 Lagrange point):
  - Operational: 2015-present
  - Instruments: WIND, MAGNETOMETER (FC, NISTAR)
  - Key data: solar wind speed, proton density, temperature, Bz (GSM)

ACE (Advanced Composition Explorer, L1):
  - Operational: 1997-present
  - Key data: solar wind, energetic particles, magnetic field

WIND spacecraft: 1994-present (NASA)
"""
import sys
import re
from pathlib import Path
from datetime import date, timedelta
import requests
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger, download_file, write_instructions

logger = get_logger("08_dscovr_ace")
OUT = DATA_DIR / "solar" / "ace_wind_dscovr"
OUT.mkdir(parents=True, exist_ok=True)

session = requests.Session()
session.headers.update({"User-Agent": "Solar-Magnetic-Analysis/1.0"})


# --------------------------------------------------------------------------- #
# 1. NOAA SWPC Real-time and recent DSCOVR data                               #
# --------------------------------------------------------------------------- #
logger.info("=== NOAA SWPC DSCOVR Real-time Data ===")
swpc_rt = {
    "dscovr_plasma_7d.json":     "https://services.swpc.noaa.gov/products/solar-wind/plasma-7-day.json",
    "dscovr_mag_7d.json":        "https://services.swpc.noaa.gov/products/solar-wind/mag-7-day.json",
    "dscovr_plasma_1d.json":     "https://services.swpc.noaa.gov/products/solar-wind/plasma-1-day.json",
    "dscovr_mag_1d.json":        "https://services.swpc.noaa.gov/products/solar-wind/mag-1-day.json",
    "ace_swepam_1d.json":        "https://services.swpc.noaa.gov/json/ace/swepam/ace-swepam-1-day.json",
    "ace_swepam_6d.json":        "https://services.swpc.noaa.gov/json/ace/swepam/ace-swepam-6-hour.json",
    "ace_mag_1d.json":           "https://services.swpc.noaa.gov/json/ace/mag/ace-magnetometer-1-day.json",
    "ace_epam_1d.json":          "https://services.swpc.noaa.gov/json/ace/epam/ace-epam-1-day.json",
    "ace_sis_1d.json":           "https://services.swpc.noaa.gov/json/ace/sis/ace-sis-1-day.json",
}
for fname, url in swpc_rt.items():
    download_file(url, OUT / "realtime" / fname, desc=f"DSCOVR/ACE {fname}", session=session)


# --------------------------------------------------------------------------- #
# 2. NOAA NCEI DSCOVR archive (netCDF files, 2016–present)                    #
# --------------------------------------------------------------------------- #
logger.info("=== NOAA NCEI DSCOVR Archive ===")
NCEI_DSCOVR = "https://www.ngdc.noaa.gov/dscovr/data"
dscovr_archive = OUT / "ncei_archive"
dscovr_archive.mkdir(exist_ok=True)

# DSCOVR archive by year/month
for year in range(2016, date.today().year + 1):
    for month in range(1, 13):
        if year == date.today().year and month > date.today().month:
            break
        yyyy, mm = str(year), f"{month:02d}"
        dir_url = f"{NCEI_DSCOVR}/{yyyy}/{mm}/"
        try:
            r = session.get(dir_url, timeout=30)
            r.raise_for_status()
            files = re.findall(r'href="(oe_[^"]+\.nc\.gz)"', r.text)
            for fname in files:
                dest = dscovr_archive / yyyy / mm / fname
                if not dest.exists():
                    download_file(f"{dir_url}{fname}", dest,
                                  desc=f"DSCOVR {yyyy}-{mm} {fname}",
                                  session=session)
        except Exception as exc:
            logger.debug(f"  DSCOVR {yyyy}/{mm}: {exc}")


# --------------------------------------------------------------------------- #
# 3. ACE Real Time Solar Wind (RTSW) archive via NOAA SWPC                    #
# --------------------------------------------------------------------------- #
logger.info("=== ACE RTSW Archive (NOAA SWPC) ===")
ACE_ARCHIVE = "https://sohoftp.nascom.nasa.gov/sdb/ace/daily"
ace_archive = OUT / "ace_archive"
ace_archive.mkdir(exist_ok=True)

for year in range(1998, date.today().year + 1):
    for month in range(1, 13):
        if year == date.today().year and month > date.today().month:
            break
        yyyy, mm = str(year), f"{month:02d}"
        # Try multiple naming conventions for ACE archive
        for fname_pattern in [
            f"{yyyy}{mm}01_ace_swepam_1m.txt",
            f"{yyyy}{mm}_ace_swepam_1m.txt",
        ]:
            url = f"{ACE_ARCHIVE}/{year}/{fname_pattern}"
            dest = ace_archive / str(year) / fname_pattern
            if download_file(url, dest, desc=f"ACE RTSW {yyyy}-{mm}", session=session):
                break


# --------------------------------------------------------------------------- #
# 4. CDAWeb instructions for full ACE/WIND/DSCOVR CDF files                   #
# --------------------------------------------------------------------------- #
write_instructions(
    OUT / "cdaweb_setup",
    "CDAWeb NASA — ACE / WIND / DSCOVR Full Archive",
    """
NASA CDAWeb (https://cdaweb.gsfc.nasa.gov/) hosts the full scientific archive
for ACE, WIND, and DSCOVR in CDF format.

AUTOMATED DOWNLOAD via CDAWeb REST API (no registration required):

  # List available datasets
  curl "https://cdaweb.gsfc.nasa.gov/WS/cdasws/1/dataviews/sp_phys/datasets"

  # Download ACE SWEPAM solar wind (1 hour averages):
  Base: https://cdaweb.gsfc.nasa.gov/pub/data/ace/swepam/level_2_cdaweb/swe_h2/
  Files: ac_h2_swe_YYYYMMDD_v0N.cdf

  # Download DSCOVR FC solar wind:
  Base: https://cdaweb.gsfc.nasa.gov/pub/data/dscovr/h1/
  Files: dscovr_h1_fc_YYYYMMDD_v01.cdf

  # Download WIND SWE solar wind:
  Base: https://cdaweb.gsfc.nasa.gov/pub/data/wind/swe/swe_h1/
  Files: wi_h1_swe_YYYYMMDD_vNN.cdf

Python access via hapi-client:
  pip install hapiclient
  from hapiclient import hapi
  data, meta = hapi('https://cdaweb.gsfc.nasa.gov/hapi',
                    'AC_H2_SWE',
                    'Np,Vp,Tpr',
                    '1998-01-01', '2026-01-01')
"""
)

logger.info("=== Script 08 complete ===")
