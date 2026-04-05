"""
Download GOES-R era (16/17/18) XRS 1-min averaged data from NCEI data.ngdc.noaa.gov.
This extends the GOES 8-15 archive (script 22) to cover Solar Cycle 25 (2017–present).
Source: https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/

GOES-16: primary 2017–2024 (East)
GOES-17: backup 2018–2023 (West, decommissioned)
GOES-18: primary 2022–present (East, replaced GOES-16 as primary)

Product: xrsf-l2-avg1m_science (1-minute averaged X-ray flux)
Output: data/solar/goes_xrs/goes{16,17,18}/
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from download.utils import get_logger, download_file

import re
import time
import requests

logger = get_logger("goes_r_xrs")
BASE_DATA = Path(__file__).resolve().parents[2] / "data" / "solar" / "goes_xrs"
BASE_DATA.mkdir(parents=True, exist_ok=True)

NGDC = "https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes"

SATELLITES = {
    "goes16": range(2017, 2025),
    "goes18": range(2022, 2027),
}
PRODUCT = "xrsf-l2-avg1m_science"

def list_months(sat: str, year: int, session: requests.Session) -> list[str]:
    url = f"{NGDC}/{sat}/l2/data/{PRODUCT}/{year}/"
    try:
        r = session.get(url, timeout=20)
        if r.status_code == 404:
            return []
        r.raise_for_status()
        return re.findall(r'href="(\d{2}/)"', r.text)
    except Exception as e:
        logger.warning(f"  {sat}/{year}: {e}")
        return []

def list_files(sat: str, year: int, month: str, session: requests.Session) -> list[str]:
    url = f"{NGDC}/{sat}/l2/data/{PRODUCT}/{year}/{month}"
    try:
        r = session.get(url, timeout=20)
        if r.status_code == 404:
            return []
        r.raise_for_status()
        return re.findall(r'href="(sci_xrsf[^"]+\.nc)"', r.text)
    except Exception as e:
        logger.warning(f"  {sat}/{year}/{month}: {e}")
        return []

def main():
    logger.info("=== GOES-R XRS 1-min Averaged Data Download ===")

    session = requests.Session()
    session.headers["User-Agent"] = "solar-research/1.0 (academic)"

    total_new = 0
    total_skip = 0

    for sat, years in SATELLITES.items():
        sat_dir = BASE_DATA / sat
        sat_dir.mkdir(exist_ok=True)
        logger.info(f"\n  Satellite: {sat.upper()}")

        for year in years:
            months = list_months(sat, year, session)
            if not months:
                logger.info(f"    {year}: no data")
                continue

            year_new = 0
            for month in months:
                files = list_files(sat, year, month, session)
                for fname in files:
                    dest = sat_dir / fname
                    if dest.exists() and dest.stat().st_size > 10_000:
                        total_skip += 1
                        continue
                    url = f"{NGDC}/{sat}/l2/data/{PRODUCT}/{year}/{month}{fname}"
                    ok = download_file(url, dest, desc=f"{sat}/{year}/{month}{fname[:30]}", retries=3)
                    if ok:
                        total_new += 1
                        year_new += 1
                    time.sleep(0.1)

            logger.info(f"    {year}: {year_new} new files downloaded")

    logger.info(f"\n✓ GOES-R XRS download complete: {total_new} new, {total_skip} already present")

if __name__ == "__main__":
    main()
