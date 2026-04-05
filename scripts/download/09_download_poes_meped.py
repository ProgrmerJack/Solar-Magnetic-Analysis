"""
Script 09 — NOAA POES / MEPED Energetic Particle Precipitation
Downloads NOAA POES (Polar Orbiting Environmental Satellites) particle data.

POES instruments: MEPED (Medium Energy Proton and Electron Detector)
                  TED (Total Energy Detector)
Provides: Hemispheric Power Index (HPI), energetic particle fluxes
Coverage: 1978–present (TIROS-N through NOAA-19, MetOp-A/B/C)

NOAA NGDC archive: https://www.ngdc.noaa.gov/stp/satellite/poes/
Also available via: https://www.ncei.noaa.gov/products/space-weather/
"""
import sys
import re
from pathlib import Path
from datetime import date
import requests
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger, download_file, write_instructions

logger = get_logger("09_poes_meped")
OUT = DATA_DIR / "atmospheric" / "poes_meped"
OUT.mkdir(parents=True, exist_ok=True)

session = requests.Session()
session.headers.update({"User-Agent": "Solar-Magnetic-Analysis/1.0"})


# --------------------------------------------------------------------------- #
# 1. Hemispheric Power Index (HPI) — key product for EPP quantification       #
# --------------------------------------------------------------------------- #
logger.info("=== NOAA SWPC Hemispheric Power Index ===")
SWPC = "https://services.swpc.noaa.gov"
hpi_files = {
    "hemispheric_power_7d.json":      f"{SWPC}/json/hemispheric_power/hemispheric-power-7-day.json",
    "hemispheric_power_nowcast.json": f"{SWPC}/products/nowcast-hemispherical-power.json",
}
for fname, url in hpi_files.items():
    download_file(url, OUT / fname, desc=f"HPI {fname}", session=session)


# --------------------------------------------------------------------------- #
# 2. NOAA NCEI POES archive (SEM/MEPED files, ASCII)                          #
# --------------------------------------------------------------------------- #
logger.info("=== NOAA NCEI POES SEM Archive ===")
POES_BASE = "https://www.ngdc.noaa.gov/stp/satellite/poes/data/accessory"

# Satellite list (active / historical POES)
POES_SATS = {
    "noaa15": (1998, 2024),
    "noaa16": (2000, 2014),
    "noaa17": (2002, 2013),
    "noaa18": (2005, date.today().year),
    "noaa19": (2009, date.today().year),
    "metop01": (2006, date.today().year),  # MetOp-B
    "metop02": (2012, date.today().year),  # MetOp-A
    "metop03": (2018, date.today().year),  # MetOp-C
}

# Try to list POES archive structure
poes_ncei = OUT / "ncei"
poes_ncei.mkdir(exist_ok=True)

for sat, (start_yr, end_yr) in POES_SATS.items():
    sat_url = f"https://www.ngdc.noaa.gov/stp/satellite/poes/data/{sat}/"
    sat_out = poes_ncei / sat
    sat_out.mkdir(exist_ok=True)
    try:
        r = session.get(sat_url, timeout=30)
        r.raise_for_status()
        # Find yearly subdirectories
        years = re.findall(r'href="(\d{4})/"', r.text)
        for yr in years:
            if not (start_yr <= int(yr) <= end_yr):
                continue
            yr_url = f"{sat_url}{yr}/"
            yr_out = sat_out / yr
            yr_out.mkdir(exist_ok=True)
            try:
                r2 = session.get(yr_url, timeout=30)
                r2.raise_for_status()
                files = re.findall(r'href="([^"]+\.(txt|asc|gz))"', r2.text)
                for fmatch in files[:12]:  # limit to avoid huge downloads
                    fname = fmatch[0]
                    download_file(f"{yr_url}{fname}", yr_out / fname,
                                  desc=f"POES {sat} {yr} {fname}", session=session)
            except Exception as exc:
                logger.debug(f"  POES {sat}/{yr}: {exc}")
    except Exception as exc:
        logger.debug(f"  POES satellite {sat}: {exc}")


# --------------------------------------------------------------------------- #
# 3. GFZ Potsdam HPI archive (cross-calibrated hemispheric power)             #
# --------------------------------------------------------------------------- #
logger.info("=== GFZ / WDC Hemispheric Power ===")
gfz_hpi_urls = {
    "noaa_hpi_ovation.txt": "https://www.ngdc.noaa.gov/stp/ovation_prime/ovation_prime_realtime.txt",
}
for fname, url in gfz_hpi_urls.items():
    download_file(url, OUT / fname, desc=fname, session=session)


# --------------------------------------------------------------------------- #
# 4. CDAWeb POES data setup instructions                                       #
# --------------------------------------------------------------------------- #
write_instructions(
    OUT / "cdaweb_meped",
    "CDAWeb — NOAA POES MEPED Full Science Archive",
    """
Full MEPED science-grade CDF files are available from NASA CDAWeb:
  https://cdaweb.gsfc.nasa.gov/

DATASET IDs for POES MEPED:
  NOAA15_H1_SEM  — NOAA-15 SEM particle fluxes (hourly)
  NOAA16_H1_SEM  — NOAA-16 SEM particle fluxes
  NOAA17_H1_SEM  — NOAA-17
  NOAA18_H1_SEM  — NOAA-18
  NOAA19_H1_SEM  — NOAA-19

Also available via the SPASE Virtual Observatory:
  https://hapi-server.org/servers/

Python download:
  pip install hapiclient
  from hapiclient import hapi
  data, meta = hapi('https://cdaweb.gsfc.nasa.gov/hapi',
                    'NOAA15_H1_SEM', 'E1_90', '1998-06-01', '2020-01-01')

For the Hemispheric Power Index specifically:
  Dataset: POES_HPI
  Or use: https://www.ngdc.noaa.gov/stp/satellite/poes/dataaccess.html
"""
)

logger.info("=== Script 09 complete ===")
