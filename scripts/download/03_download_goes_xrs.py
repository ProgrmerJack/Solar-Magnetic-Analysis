"""
Script 03 — NOAA GOES X-ray Sensor (XRS) Data
Downloads GOES-R series (GOES-16 + GOES-18) XRS 1-minute averages
and 1-second science-grade data from NOAA NCEI.

GOES-16:  operational 2017-present  (GOES East until 2022, then overlap)
GOES-18:  operational 2022-present  (current GOES East)

Files are daily NetCDF-4 at:
  https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/
  goes/{satellite}/l2/data/{product}/YYYY/MM/
"""
import sys
import re
from pathlib import Path
from datetime import date, timedelta
import requests
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger, download_file

logger = get_logger("03_goes_xrs")

NCEI = "https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes"

# Products we want
PRODUCTS = {
    "xrsf-l2-flx1s_science": "1-second XRS science grade",
    "xrsf-l2-avg1m_science":  "1-minute XRS science grade",
}

# Satellites and their operational windows
SATELLITES = {
    "goes16": (date(2017, 2, 1),  date.today()),
    "goes18": (date(2022, 6, 1),  date.today()),
}

def list_ncei_dir(url: str, session: requests.Session) -> list[str]:
    """Return filenames listed in an NCEI THREDDS/Apache directory."""
    try:
        r = session.get(url, timeout=30, headers={"User-Agent": "Solar-Magnetic-Analysis/1.0"})
        r.raise_for_status()
        # find .nc links
        return re.findall(r'href="([^"]+\.nc)"', r.text)
    except Exception as exc:
        logger.warning(f"  Could not list {url}: {exc}")
        return []

def download_goes_product(sat: str, product: str, start: date, end: date):
    out_root = DATA_DIR / "solar" / "goes_xrs" / sat / product
    out_root.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    current = date(start.year, start.month, 1)
    end_month = date(end.year, end.month, 1)

    downloaded, skipped, failed = 0, 0, 0

    while current <= end_month:
        yyyy = current.strftime("%Y")
        mm   = current.strftime("%m")
        dir_url = f"{NCEI}/{sat}/l2/data/{product}/{yyyy}/{mm}/"
        files   = list_ncei_dir(dir_url, session)

        if not files:
            logger.debug(f"  No files found: {dir_url}")
        else:
            for fname in files:
                dest = out_root / yyyy / mm / fname
                if dest.exists():
                    skipped += 1
                    continue
                dest.parent.mkdir(parents=True, exist_ok=True)
                ok = download_file(f"{dir_url}{fname}", dest,
                                   desc=f"{sat}/{product} {yyyy}-{mm} {fname}",
                                   session=session)
                if ok:
                    downloaded += 1
                else:
                    failed += 1

        # advance one month
        if current.month == 12:
            current = date(current.year + 1, 1, 1)
        else:
            current = date(current.year, current.month + 1, 1)

    logger.info(f"  {sat}/{product}: {downloaded} new, {skipped} skipped, {failed} failed")


if __name__ == "__main__":
    for sat, (start, end) in SATELLITES.items():
        for product, desc in PRODUCTS.items():
            logger.info(f"=== {sat} / {product} ({desc}) ===")
            download_goes_product(sat, product, start, end)

    logger.info("=== Script 03 complete ===")
