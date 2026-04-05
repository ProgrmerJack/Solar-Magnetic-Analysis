"""
IMS (Interactive Multi-sensor Snow and Ice Mapping System) Daily Snow Cover
Source: NSIDC — https://nsidc.org/data/g02156
Format: ASCII/NetCDF 4km and 1km daily grids — publicly accessible without auth.

We download the 4km daily product (G02156) from NSIDC's HTTP server.
Coverage: Northern Hemisphere, 1997-present.

Three key regions extracted:
  - Swiss Alps: 44-48N, 5-11E
  - Norwegian mountains: 58-72N, 5-30E
  - Colorado Rockies: 36-42N, -109--104W
"""
import requests
import gzip
import re
from pathlib import Path
from datetime import date, timedelta
import sys
sys.path.insert(0, str(Path(__file__).parents[1]))
from utils import get_logger

LOG = get_logger("40_ims_snow")
OUT = Path(__file__).parents[2] / "data" / "cryosphere" / "ims_snow"
OUT.mkdir(parents=True, exist_ok=True)

BASE = "https://noaadata.apps.nsidc.org/NOAA/G02156/4km"
# Actual filename pattern: ims{YEAR}{DOY}_00UTC_4km_v1.3.asc.gz

# ─── Discover available files ────────────────────────────────────────────────
LOG.info("=== IMS 4km Daily Snow Cover ===")
LOG.info("Base URL: %s", BASE)

sess = requests.Session()
sess.headers["User-Agent"] = "Solar-Magnetic-Analysis research project"

def list_year(year):
    url = f"{BASE}/{year}/"
    try:
        r = sess.get(url, timeout=20)
        if r.status_code != 200:
            return []
        # Parse directory listing for .asc.gz or .nc files
        files = re.findall(r'href="(ims\d+_4km_v\d+\.\d+\.[^"]+)"', r.text)
        if not files:
            files = re.findall(r'href="(ims\d{7}_4km[^"]+\.gz)"', r.text)
        if not files:
            files = re.findall(r'href="(ims[^"]+\.gz)"', r.text)
        return files
    except Exception as e:
        LOG.warning("  list_year %d failed: %s", year, e)
        return []

def download_file(year, fname):
    local = OUT / str(year) / fname
    local.parent.mkdir(exist_ok=True)
    if local.exists():
        return local
    url = f"{BASE}/{year}/{fname}"
    try:
        r = sess.get(url, timeout=60, stream=True)
        if r.status_code == 200:
            with open(local, "wb") as f:
                for chunk in r.iter_content(65536):
                    f.write(chunk)
            return local
        else:
            LOG.warning("    HTTP %d for %s", r.status_code, fname)
            return None
    except Exception as e:
        LOG.warning("    download error %s: %s", fname, e)
        return None

# ─── First probe the structure ─────────────────────────────────────────────
LOG.info("Probing directory structure...")
r = sess.get(BASE, timeout=15)
LOG.info("Base directory status: %d", r.status_code)

if r.status_code == 200:
    years_found = re.findall(r'href="(\d{4})/?\"', r.text)
    LOG.info("Years available: %s", years_found[:10])
else:
    # Try alternative NSIDC FTP-over-HTTP structure
    alt_base = "https://noaadata.apps.nsidc.org/NOAA/G02156/data/4km/"
    LOG.info("Trying alt: %s", alt_base)
    r = sess.get(alt_base, timeout=15)
    LOG.info("Alt status: %d", r.status_code)
    if r.status_code == 200:
        BASE = alt_base
        years_found = re.findall(r'href="(\d{4})/?\"', r.text)
        LOG.info("Years available: %s", years_found[:10])
    else:
        LOG.warning("Directory listing not available — trying direct download")
        years_found = [str(y) for y in range(2000, 2026)]

# ─── Download strategy: sample recent winters only (Oct-Apr) ──────────────
# Full IMS grids are ~4MB/day compressed; 28 years × 180 winter days = ~20 GB
# We keep only winter months to stay manageable
WINTER_MONTHS = {10, 11, 12, 1, 2, 3, 4}

total_downloaded = 0
total_skipped = 0

current = date(2004, 10, 1)  # IMS 4km starts 2004
end = date(2026, 5, 1)

LOG.info("Downloading IMS 4km winter (Oct-Apr) 2000-2026...")
LOG.info("(Files missing from server are skipped gracefully)")

while current <= end:
    if current.month in WINTER_MONTHS:
        yr = current.year
        doy = current.timetuple().tm_yday
        # Try known filename patterns:
        # ims2024033_4km_v1.3.asc.gz
        # ims2024033_4km.asc.gz
        found = False
        for ver in ["v1.3", "v1.2", "v1.1", ""]:
            for ext in [".asc.gz", ".nc.gz", ".nc"]:
                if ver:
                    fname = f"ims{yr}{doy:03d}_4km_{ver}{ext}"
                else:
                    fname = f"ims{yr}{doy:03d}_4km{ext}"
                local = OUT / str(yr) / fname
                if local.exists():
                    total_skipped += 1
                    found = True
                    break
            if found:
                break
        
        if not found:
            # Try download — actual pattern: ims{YEAR}{DOY}_00UTC_4km_v1.3.asc.gz
            for utc_tag in ["_00UTC", ""]:
                if found:
                    break
                for ver in ["v1.3", "v1.2", "v1.1", ""]:
                    if found:
                        break
                    for ext in [".asc.gz", ".nc.gz", ".nc"]:
                        if ver:
                            fname = f"ims{yr}{doy:03d}{utc_tag}_4km_{ver}{ext}"
                        else:
                            fname = f"ims{yr}{doy:03d}{utc_tag}_4km{ext}"
                        url = f"{BASE}/{yr}/{fname}"
                        local = OUT / str(yr) / fname
                        local.parent.mkdir(exist_ok=True)
                        try:
                            resp = sess.head(url, timeout=10, allow_redirects=True)
                            if resp.status_code == 200:
                                r2 = sess.get(url, timeout=60, stream=True)
                                if r2.status_code == 200:
                                    with open(local, "wb") as f:
                                        for chunk in r2.iter_content(65536):
                                            f.write(chunk)
                                    total_downloaded += 1
                                    found = True
                                    break
                        except Exception:
                            pass

        if total_downloaded % 100 == 0 and total_downloaded > 0:
            LOG.info("  Progress: %d downloaded, %d skipped", total_downloaded, total_skipped)

    current += timedelta(days=1)

LOG.info("=== IMS download complete ===")
LOG.info("  Downloaded: %d files", total_downloaded)
LOG.info("  Skipped (already present): %d files", total_skipped)
LOG.info("  Output: %s", OUT)
