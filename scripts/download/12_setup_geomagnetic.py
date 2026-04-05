"""
Script 12 — SuperMAG & INTERMAGNET Geomagnetic Network Data
Sets up access for 500+ ground magnetometer station networks.

SuperMAG (supermag.jhuapl.edu):
  - 500+ stations worldwide, 1970–present
  - Key indices: SML, SMU, SME (auroral electrojet)
  - REGISTRATION REQUIRED: https://supermag.jhuapl.edu/info/?page=services

INTERMAGNET (intermagnet.org):
  - ~150 geomagnetic observatories, 1991–present
  - 1-minute and 1-second definitive data
  - REGISTRATION: free, contact individual IMOs or use portal

Also includes WDC Kyoto Dst / AE index downloads.
"""
import sys
import json
from pathlib import Path
import requests
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger, download_file, write_instructions

logger = get_logger("12_supermag_intermagnet")

session = requests.Session()
session.headers.update({"User-Agent": "Solar-Magnetic-Analysis/1.0"})

# --------------------------------------------------------------------------- #
# 1. WDC Kyoto Dst Index (publicly accessible web service)                    #
# --------------------------------------------------------------------------- #
logger.info("=== WDC Kyoto Dst Index ===")
DST_OUT = DATA_DIR / "geomagnetic" / "dst_ae_index"
DST_OUT.mkdir(parents=True, exist_ok=True)

# WDC Kyoto provides Dst via a service URL
# Annual Dst files: quick-look (provisional + definitive)
kyoto_dst_base = "https://wdc.kugi.kyoto-u.ac.jp"

# Try direct file access for yearly Dst
for year in range(1957, 2026):
    fname = f"dst_{year}.dat"
    # Try multiple URL patterns WDC Kyoto uses
    urls_to_try = [
        f"{kyoto_dst_base}/dst_final/{year}/dst{year:02d}.for.request",
        f"https://wdc.kugi.kyoto-u.ac.jp/dstdir/dst_YYYY.for.request".replace("YYYY", str(year)),
    ]
    for url in urls_to_try:
        if download_file(url, DST_OUT / fname, desc=f"Dst {year}", session=session):
            break


# --------------------------------------------------------------------------- #
# 2. Kp / Dst / AE from NOAA NCEI (alternative source)                       #
# --------------------------------------------------------------------------- #
logger.info("=== NOAA NCEI Geomagnetic Indices ===")
NCEI_GEO = "https://www.ngdc.noaa.gov/geomag-web/calculators/calculateIGRF"
noaa_geo = {
    "noaa_kp_ap_monthly.json": "https://services.swpc.noaa.gov/json/planetary_k_index_1m.json",
    "noaa_geomag_storm_cat.json": "https://services.swpc.noaa.gov/products/noaa-geomagnetic-storm-forecast-27-day.json",
    "noaa_geomag_alerts.json": "https://services.swpc.noaa.gov/products/alerts.json",
    "noaa_geomag_3day.json": "https://services.swpc.noaa.gov/products/geomagnetic-forecast.json",
}
geo_out = DATA_DIR / "geomagnetic" / "dst_ae_index"
for fname, url in noaa_geo.items():
    download_file(url, geo_out / fname, desc=fname, session=session)


# --------------------------------------------------------------------------- #
# 3. SuperMAG instructions                                                     #
# --------------------------------------------------------------------------- #
write_instructions(
    DATA_DIR / "geomagnetic" / "supermag",
    "SuperMAG 500+ Station Geomagnetic Network",
    """
SuperMAG provides 1-minute resolution data from 500+ ground magnetometer stations,
with cross-calibrated, rotated, and baseline-subtracted perturbation vectors.
Includes: SML, SMU, SME indices (global auroral electrojet estimates).

REGISTRATION: https://supermag.jhuapl.edu/info/?page=services (free, academic)

PYTHON ACCESS (after registration):
  import pyspedas  # pip install pyspedas
  # OR use SuperMAG web service directly:
  
  import requests
  r = requests.get("https://supermag.jhuapl.edu/mag/lib/php/", params={
      "user":    "YOUR_USERNAME",
      "start":   "2003-11-20T00:00:00",
      "extent":  "86400",          # seconds (1 day)
      "stations": "OTT,BOU,ABK",   # station codes
      "fmt":     "json",
      "elements": "NEDS",
  })

KEY PRODUCTS TO REQUEST:
  1. SME/SML/SMU indices (1970–present):
     https://supermag.jhuapl.edu/indices/?layers=SME
  2. Geomagnetic storm list:
     https://supermag.jhuapl.edu/storm/?storm=LIST
  3. Station data for polar cap stations (ABK, TRO, HOR, etc.)
     for particle precipitation proxy

ALTERNATIVE (free, no registration): 
  WDC for Geomagnetism Edinburgh (BGS):
  https://wdc.bgs.ac.uk/catalog/master.html
"""
)

# --------------------------------------------------------------------------- #
# 4. INTERMAGNET instructions                                                  #
# --------------------------------------------------------------------------- #
write_instructions(
    DATA_DIR / "geomagnetic" / "intermagnet",
    "INTERMAGNET Global Geomagnetic Observatory Network",
    """
INTERMAGNET provides definitive 1-minute geomagnetic data from ~150 observatories.

ACCESS: https://intermagnet.org/data-donnee/download-eng.php

KEY OBSERVATORIES for this research:
  ABK  — Abisko, Sweden (68°N, near Aurora zone)
  TRO  — Tromsø, Norway (69°N)  
  HRN  — Hornsund, Svalbard (77°N)
  SOD  — Sodankylä, Finland (67°N)
  BOU  — Boulder, Colorado (40°N — near CAIC region)
  FRN  — Fresno, California (37°N)
  NEW  — Newport, Washington (48°N)

PYTHON ACCESS:
  pip install intermagnet  # if available
  # OR use the INTERMAGNET web service:
  import requests
  url = "https://imag-data.bgs.ac.uk/GIN_V1/GINServices"
  params = {
      "Request":         "GetData",
      "observatoryIaga": "ABK",
      "startdate":       "2003-10-29",
      "enddate":         "2003-11-01",
      "samplesPerDay":   "1440",   # 1-minute data
      "dataStandard":    "IAGA2002",
      "publicationState": "definitive",
  }
  r = requests.get(url, params=params)

BGS WDC Edinburgh also provides INTERMAGNET data (no registration):
  https://wdc.bgs.ac.uk/
  
For bulk access, contact: info@intermagnet.org
"""
)

# --------------------------------------------------------------------------- #
# 5. BGS WDC Edinburgh — public geomagnetic data                              #
# --------------------------------------------------------------------------- #
logger.info("=== BGS WDC Edinburgh — Public Geomagnetic Data ===")
bgs_out = DATA_DIR / "geomagnetic" / "dst_ae_index" / "bgs_wdc"
bgs_out.mkdir(exist_ok=True)

# BGS provides some indices directly
bgs_urls = {
    "bgs_geomag_indices_info.json": "https://wdc.bgs.ac.uk/catalog/master.html",
}
for fname, url in bgs_urls.items():
    download_file(url, bgs_out / fname, desc=fname, session=session)


# --------------------------------------------------------------------------- #
# 6. Kyoto Geomagnetic indices web service (AE, AU, AL)                       #
# --------------------------------------------------------------------------- #
logger.info("=== WDC Kyoto AE Index ===")
ae_out = DST_OUT / "ae_index"
ae_out.mkdir(exist_ok=True)

# AE index from OMNI is the most accessible path (already in Script 02)
# Direct Kyoto access requires form submission, but we can try API
kyoto_ae_urls = {
    "ae_index_info.txt": "https://wdc.kugi.kyoto-u.ac.jp/ae_realtime/",
}
for fname, url in kyoto_ae_urls.items():
    download_file(url, ae_out / fname, desc=fname, session=session)

write_instructions(
    ae_out,
    "WDC Kyoto AE/AL/AU Auroral Electrojet Indices",
    """
AE (Auroral Electrojet) index measures geomagnetic activity at auroral latitudes.
Best proxy for energetic particle precipitation intensity.

PRIMARY SOURCE (no registration): OMNI dataset
  → Already downloaded via Script 02 (data/atmospheric/omni_solar_wind/)
  → AE is column in OMNI hourly files (1975–present)

DIRECT WDC KYOTO ACCESS:
  Provisional: https://wdc.kugi.kyoto-u.ac.jp/ae_provisional/
  Final:       https://wdc.kugi.kyoto-u.ac.jp/ae_final/
  
  Form-based download for monthly files.
  OR contact wdc-geo@kugi.kyoto-u.ac.jp for bulk access.

ALTERNATIVE: SuperMAG SME index (better spatial coverage, modern reprocessing)
"""
)

logger.info("=== Script 12 complete ===")
