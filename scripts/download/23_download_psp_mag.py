"""
Download Parker Solar Probe (PSP) FIELDS L2 magnetometer data (RTN, 4 sa/cycle)
for SOC analysis of solar coronal magnetic field fluctuations.
Source: SPDF https://spdf.gsfc.nasa.gov/pub/data/psp/fields/l2/mag_rtn_4_per_cycle/
Coverage: 2018–present (PSP launched Aug 2018)
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from download.utils import get_logger, download_file

import re
import time
import requests

logger = get_logger("psp_mag")
BASE = Path(__file__).resolve().parents[2] / "data" / "solar" / "psp_mag"
BASE.mkdir(parents=True, exist_ok=True)

SPDF_BASE = "https://spdf.gsfc.nasa.gov/pub/data/psp/fields/l2/mag_rtn_4_per_cycle"
YEARS = range(2018, 2027)

def list_files(year: int) -> list[str]:
    url = f"{SPDF_BASE}/{year}/"
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 404:
                return []
            r.raise_for_status()
            return re.findall(r'href="(psp_fld_l2_mag_rtn[^"]+\.cdf)"', r.text)
        except Exception as e:
            logger.warning(f"  Listing {year} attempt {attempt+1}/3: {e}")
            time.sleep(5)
    return []

def main():
    logger.info("=== Parker Solar Probe MAG RTN 4sa/cycle download ===")
    logger.info(f"Output: {BASE}")

    total_downloaded = 0
    total_skipped = 0

    for year in YEARS:
        year_dir = BASE / str(year)
        year_dir.mkdir(exist_ok=True)

        files = list_files(year)
        if not files:
            logger.info(f"  {year}: no files found (likely not yet available)")
            continue

        logger.info(f"  {year}: {len(files)} files")
        for fname in files:
            dest = year_dir / fname
            if dest.exists() and dest.stat().st_size > 1000:
                total_skipped += 1
                continue
            url = f"{SPDF_BASE}/{year}/{fname}"
            ok = download_file(url, dest, desc=f"PSP {year}/{fname[:35]}", retries=3)
            if ok:
                total_downloaded += 1
            time.sleep(0.3)  # be polite to SPDF

        logger.info(f"    {year} done: {total_downloaded} downloaded, {total_skipped} already present")

    logger.info(f"\n✓ PSP MAG download complete: {total_downloaded} new files, {total_skipped} pre-existing")

if __name__ == "__main__":
    main()
