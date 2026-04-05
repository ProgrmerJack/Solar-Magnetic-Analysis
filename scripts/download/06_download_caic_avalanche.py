"""
Script 06 — CAIC Colorado Avalanche Information Center
Downloads avalanche accident statistics from the Colorado Avalanche Information Center.

Primary data:
  • CAIC accident CSV (all accidents since 1950)
  • CAIC annual reports (PDF)
  • Utah Avalanche Center (UAC) API — machine-readable accident data
  • NWAC (Northwest Avalanche Center) data

URL: https://avalanche.state.co.us/accidents/statistics-and-reporting/
"""
import sys
import json
from pathlib import Path
from datetime import date
import requests
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger, download_file, write_instructions

logger = get_logger("06_caic")
OUT = DATA_DIR / "cryosphere" / "caic"
OUT.mkdir(parents=True, exist_ok=True)

session = requests.Session()
session.headers.update({
    "User-Agent": "Mozilla/5.0 (compatible; Solar-Magnetic-Analysis/1.0)"
})


# --------------------------------------------------------------------------- #
# 1. CAIC direct accident data endpoints                                       #
# --------------------------------------------------------------------------- #
logger.info("=== CAIC Avalanche Accident Data ===")

caic_endpoints = {
    # CAIC accident data - try multiple known endpoints
    "caic_accidents.csv": "https://avalanche.state.co.us/caic/media/doc/accidents.csv",
    "caic_accidents_v2.csv": "https://avalanche.state.co.us/accidents/accidents.csv",
    "caic_statistics.json": "https://avalanche.state.co.us/api/accidents",
}

for fname, url in caic_endpoints.items():
    download_file(url, OUT / fname, desc=f"CAIC {fname}", session=session)


# --------------------------------------------------------------------------- #
# 2. Utah Avalanche Center (UAC) — public API, machine-readable                #
# --------------------------------------------------------------------------- #
logger.info("=== Utah Avalanche Center (UAC) Accident Data ===")
uac_out = OUT / "uac"
uac_out.mkdir(exist_ok=True)

uac_endpoints = {
    "uac_accidents_all.json": "https://utahavalanchecenter.org/avalanche_accidents.json",
    "uac_avalanches.json":    "https://utahavalanchecenter.org/avalanches.json",
}
for fname, url in uac_endpoints.items():
    download_file(url, OUT / fname, desc=f"UAC {fname}", session=session)


# --------------------------------------------------------------------------- #
# 3. American Avalanche Association — public accident data                     #
# --------------------------------------------------------------------------- #
logger.info("=== American Avalanche Association Data ===")
aaa_endpoints = {
    "aaa_accidents.json": "https://www.americanavalancheassociation.org/wp-json/wp/v2/avalanche-accidents?per_page=100",
    "aaa_stats.json":     "https://avalanche.org/avalanche-accidents/",
}
for fname, url in aaa_endpoints.items():
    download_file(url, OUT / fname, desc=f"AAA {fname}", session=session)


# --------------------------------------------------------------------------- #
# 4. NWAC (Northwest Avalanche Center) — Pacific Northwest data               #
# --------------------------------------------------------------------------- #
logger.info("=== Northwest Avalanche Center (NWAC) ===")
nwac_out = OUT / "nwac"
nwac_out.mkdir(exist_ok=True)

download_file(
    "https://nwac.us/avalanche-data/accident-observations.json",
    nwac_out / "nwac_accidents.json",
    desc="NWAC accidents",
    session=session
)


# --------------------------------------------------------------------------- #
# 5. European Avalanche Warning Services (EAWS) data                          #
# --------------------------------------------------------------------------- #
logger.info("=== EAWS Avalanche Bulletin Archive ===")
eaws_out = DATA_DIR / "cryosphere" / "eaws"
eaws_out.mkdir(parents=True, exist_ok=True)

eaws_endpoints = {
    "eaws_regions.geojson": "https://regions.avalanches.org/micro-regions.geojson",
    "eaws_latest_bulletins.json": "https://api.avalanches.org/v1/bulletins?lang=en",
}
for fname, url in eaws_endpoints.items():
    download_file(url, eaws_out / fname, desc=f"EAWS {fname}", session=session)


# --------------------------------------------------------------------------- #
# 6. Write detailed instructions for SLF Davos (requires direct contact)      #
# --------------------------------------------------------------------------- #
slf_instructions = """
The WSL Institute for Snow and Avalanche Research SLF (Davos, Switzerland) 
maintains the world's most comprehensive avalanche database (1970–present).

TO REQUEST DATA ACCESS:
  1. Visit: https://www.slf.ch/en/avalanche-bulletin-and-snow-situation/archive.html
  2. Contact: info@slf.ch
  3. Mention you are conducting academic research on avalanche-solar forcing coupling

KEY CONTACTS (check latest publications from SLF for current personnel):
  - Dr. Jürg Schweizer (Head of Snow Avalanches):  juerg.schweizer@slf.ch
  - Dr. Christoph Marty (Climate-snow links):       christoph.marty@slf.ch
  - Data Access:                                    data@slf.ch

WHAT TO REQUEST:
  1. Daily avalanche activity data (1970–present) for the Swiss Alps
     - Minimum: activity level (number of avalanches per day/region)
     - Ideal: size class, aspect, elevation, triggering type
  2. Snowpack stability data (ECT / RB test results)
  3. If possible: runout distance or release area estimates for power-law analysis

TEMPLATE EMAIL:
-----------------------------------------------------------------
Subject: Academic data request — Swiss avalanche database for solar forcing research

Dear Dr. Schweizer / Data Team,

I am conducting academic research investigating whether self-organized criticality
statistics in snow avalanche activity are modulated by solar magnetic activity and
stratospheric sudden warming events. This project aims to test the hypothesis that
solar energetic particle precipitation modulates polar vortex stability and, through
atmospheric teleconnections, periodically drives alpine snowpacks toward criticality.

To perform the power-law (SOC) analysis described in Clauset et al. (2009), I would
require the Swiss avalanche catalog (1970–present) with:
  - Event date, region/aspect, size class (SLF 1–5 scale)
  - If available: release area or runout estimate

I am happy to sign a data sharing agreement, provide institutional affiliation, and
include SLF co-authorship for any publication using this data.

Kind regards,
[Your name and institutional affiliation]
-----------------------------------------------------------------

ALTERNATIVE PUBLIC DATA SOURCES:
  • SLF Avalanche Bulletin Archive (bulletins only):
    https://www.slf.ch/en/avalanche-bulletin-and-snow-situation/archive.html
  • WSL Open Research Data:
    https://www.envidat.ch/#/metadata/avalanche-accidents-switzerland
  • EnviDat (Swiss data repository, may have SLF datasets):
    https://www.envidat.ch/
"""
write_instructions(
    DATA_DIR / "cryosphere" / "slf_avalanche",
    "WSL/SLF Davos Avalanche Database — Data Access Instructions",
    slf_instructions
)

# --------------------------------------------------------------------------- #
# 7. Try EnviDat (Swiss environmental data portal) — may have SLF data        #
# --------------------------------------------------------------------------- #
logger.info("=== EnviDat Swiss Environmental Data Portal ===")
slf_open = DATA_DIR / "cryosphere" / "slf_avalanche" / "envidat"
slf_open.mkdir(exist_ok=True)

envidat_urls = {
    "envidat_search_avalanche.json":
        "https://www.envidat.ch/api/3/action/package_search?q=avalanche&rows=50",
    "envidat_slf_packages.json":
        "https://www.envidat.ch/api/3/action/package_search?q=SLF+avalanche&rows=20",
}
for fname, url in envidat_urls.items():
    download_file(url, slf_open / fname, desc=f"EnviDat {fname}", session=session)


logger.info("=== Script 06 complete ===")
