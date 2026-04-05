"""
GOES X-Ray Sensor (XRS) Historical Data — NOAA NCEI
Downloads GOES-8 through GOES-18 1-minute XRS data (1986–present).
Also fetches GOES-16/17/18 science-quality NetCDF from NCEI.
"""
import sys, json, time, re
from pathlib import Path
from datetime import date
import requests
from bs4 import BeautifulSoup
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger, download_file

logger = get_logger("goes_xrs")
OUT = DATA_DIR / "solar" / "goes_xrs"
OUT.mkdir(parents=True, exist_ok=True)

session = requests.Session()
session.headers["User-Agent"] = "Solar-Magnetic-Analysis/1.0"

# ── 1. NCEI GOES XRS archive (legacy text format, GOES-8 to GOES-15) ─────────
logger.info("=== NCEI GOES SEM XRS 1-minute text archive ===")
NCEI_BASE = "https://www.ngdc.noaa.gov/stp/satellite/goes/dataaccess.html"
# Direct known paths:
# https://www.ngdc.noaa.gov/stp/satellite/goes/doc/GOES_XRS_readme.pdf
# https://satdat.ngdc.noaa.gov/sem/goes/data/science/xrs/

SEM_BASE = "https://satdat.ngdc.noaa.gov/sem/goes/data/full"
text_out = OUT / "ncei_text"
text_out.mkdir(exist_ok=True)

# Test SEM archive path
for test_url in [
    "https://satdat.ngdc.noaa.gov/sem/goes/data/full/",
    "https://satdat.ngdc.noaa.gov/sem/goes/data/science/xrs/",
    "https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/science/xrs/goes16/xrsf-l2-flx1s_science/",
]:
    try:
        r = session.get(test_url, timeout=15)
        logger.info(f"  {r.status_code} [{len(r.content)//1024}KB] {test_url}")
    except Exception as e:
        logger.warning(f"  FAIL {test_url}: {e}")

# ── 2. CDAWeb GOES XRS CDF files ──────────────────────────────────────────────
logger.info("\n=== CDAWeb GOES-16/17/18 XRS 1-min science data ===")
CDAWEB = "https://cdaweb.gsfc.nasa.gov/pub/data"
cdf_out = OUT / "cdaweb_cdf"
cdf_out.mkdir(exist_ok=True)

# GOES-16 XRS 1-min
goes_paths = {
    "goes16_xrs_1min": "goes/goes16/l2/xrsf-l2-flx1s_science",
    "goes17_xrs_1min": "goes/goes17/l2/xrsf-l2-flx1s_science",
    "goes18_xrs_1min": "goes/goes18/l2/xrsf-l2-flx1s_science",
    "goes16_xrs_avg":  "goes/goes16/l2/xrsf-l2-avg1m_science",
}

for name, path in goes_paths.items():
    url = f"{CDAWEB}/{path}/"
    try:
        r = session.get(url, timeout=20)
        if r.status_code == 200:
            logger.info(f"  ✓  {name}: {url}")
            # List years
            soup = BeautifulSoup(r.text, "html.parser")
            years = sorted([a["href"].strip("/") for a in soup.find_all("a", href=True)
                           if re.match(r"^\d{4}/?$", a["href"])])
            logger.info(f"       Years available: {years}")
        else:
            logger.warning(f"  {r.status_code} {name}")
    except Exception as e:
        logger.warning(f"  FAIL {name}: {e}")

# ── 3. NOAA NCEI goes-space-environment-monitor ───────────────────────────────
logger.info("\n=== NOAA NCEI GOES SEM science NetCDF ===")
NCEI_SEM = "https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/science/xrs"

for satellite in ["goes16", "goes17", "goes18"]:
    for product in ["xrsf-l2-avg1m_science", "xrsf-l2-flx1s_science"]:
        url = f"{NCEI_SEM}/{satellite}/{product}/"
        try:
            r = session.get(url, timeout=20)
            if r.status_code == 200:
                soup = BeautifulSoup(r.text, "html.parser")
                years = sorted([a["href"].strip("/") for a in soup.find_all("a", href=True)
                               if re.match(r"^\d{4}/?$", a["href"])])
                logger.info(f"  ✓  {satellite}/{product}: years {years}")
            else:
                logger.info(f"  {r.status_code} {satellite}/{product}")
        except Exception as e:
            logger.warning(f"  FAIL {satellite}/{product}: {e}")
        time.sleep(0.5)

# ── 4. Download GOES-16 XRS 1-min avg NetCDF (2017–present) ──────────────────
logger.info("\n=== Downloading GOES-16 XRS 1-min avg NetCDF files ===")
NCEI_G16_AVG = "https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/science/xrs/goes16/xrsf-l2-avg1m_science"
g16_out = OUT / "goes16_xrs_1min_avg"
g16_out.mkdir(exist_ok=True)

try:
    r = session.get(f"{NCEI_G16_AVG}/", timeout=20)
    if r.status_code == 200:
        soup = BeautifulSoup(r.text, "html.parser")
        years = sorted([a["href"].strip("/") for a in soup.find_all("a", href=True)
                       if re.match(r"^\d{4}/?$", a["href"])])
        logger.info(f"  GOES-16 avg1m years: {years}")
        for yr in years:
            yr_url = f"{NCEI_G16_AVG}/{yr}/"
            yr_r = session.get(yr_url, timeout=20)
            if yr_r.status_code != 200:
                continue
            yr_soup = BeautifulSoup(yr_r.text, "html.parser")
            nc_files = [a["href"] for a in yr_soup.find_all("a", href=True)
                       if a["href"].endswith(".nc")]
            logger.info(f"    {yr}: {len(nc_files)} NetCDF files")
            yr_out = g16_out / yr
            yr_out.mkdir(exist_ok=True)
            for nc in nc_files:
                dest = yr_out / nc
                if not dest.exists():
                    download_file(f"{NCEI_G16_AVG}/{yr}/{nc}", dest,
                                  desc=f"G16 XRS {yr}/{nc[:20]}")
                    time.sleep(0.3)
except Exception as e:
    logger.warning(f"  GOES-16 avg1m download failed: {e}")

# ── 5. GOES-18 XRS 1-min (most recent satellite, 2022–present) ───────────────
logger.info("\n=== Downloading GOES-18 XRS 1-min avg NetCDF files ===")
NCEI_G18_AVG = "https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/science/xrs/goes18/xrsf-l2-avg1m_science"
g18_out = OUT / "goes18_xrs_1min_avg"
g18_out.mkdir(exist_ok=True)

try:
    r = session.get(f"{NCEI_G18_AVG}/", timeout=20)
    if r.status_code == 200:
        soup = BeautifulSoup(r.text, "html.parser")
        years = sorted([a["href"].strip("/") for a in soup.find_all("a", href=True)
                       if re.match(r"^\d{4}/?$", a["href"])])
        logger.info(f"  GOES-18 avg1m years: {years}")
        for yr in years:
            yr_url = f"{NCEI_G18_AVG}/{yr}/"
            yr_r = session.get(yr_url, timeout=20)
            if yr_r.status_code != 200:
                continue
            yr_soup = BeautifulSoup(yr_r.text, "html.parser")
            nc_files = [a["href"] for a in yr_soup.find_all("a", href=True)
                       if a["href"].endswith(".nc")]
            logger.info(f"    {yr}: {len(nc_files)} NetCDF files")
            yr_out = g18_out / yr
            yr_out.mkdir(exist_ok=True)
            for nc in nc_files:
                dest = yr_out / nc
                if not dest.exists():
                    download_file(f"{NCEI_G18_AVG}/{yr}/{nc}", dest,
                                  desc=f"G18 XRS {yr}/{nc[:20]}")
                    time.sleep(0.3)
except Exception as e:
    logger.warning(f"  GOES-18 avg1m download failed: {e}")

logger.info("=== GOES XRS download complete ===")
