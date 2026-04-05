"""
Download Kyoto AE/AL/AU index (1957–present) from WDC Kyoto.
Mirrors the working approach used for Dst (script 15).
Note: OMNI2 already contains AE from 1963; this gives standalone 1957+ and cross-validation.
Paths:
  ae_final/YYYYMM/     1957–~2020  (definitive)
  ae_provisional/YYYYMM/  recent provisional
  ae_realtime/YYYYMM/   near-real-time
Output: data/geomagnetic/dst_ae_index/kyoto_ae_1957_present.csv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from download.utils import get_logger

import requests
import re
import time
import pandas as pd
from datetime import datetime, date
from bs4 import BeautifulSoup

logger = get_logger("kyoto_ae")
OUT_DIR = Path(__file__).resolve().parents[2] / "data" / "geomagnetic" / "dst_ae_index"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_FILE = OUT_DIR / "kyoto_ae_1957_present.csv"

BASE_URL = "https://wdc.kugi.kyoto-u.ac.jp"

# Define which path tier to try per year-month
def get_url(year: int, month: int) -> list[str]:
    ym = f"{year:04d}{month:02d}"
    today = date.today()
    yr_mo = date(year, month, 1)
    # Rough tiers (same as Dst approach)
    # definitive lags ~2 years
    cutoff_final = date(today.year - 2, today.month, 1)
    cutoff_prov  = date(today.year - 1, today.month, 1)
    candidates = []
    if yr_mo < cutoff_final:
        candidates.append(f"{BASE_URL}/ae_final/{ym}/index.html")
    if yr_mo < cutoff_prov:
        candidates.append(f"{BASE_URL}/ae_provisional/{ym}/index.html")
    candidates.append(f"{BASE_URL}/ae_realtime/{ym}/index.html")
    return candidates

def fetch_ae_month(year: int, month: int, session: requests.Session) -> list[dict]:
    """Fetch AE hourly data for one month from WDC Kyoto HTML page."""
    urls = get_url(year, month)
    for url in urls:
        for attempt in range(3):
            try:
                r = session.get(url, timeout=20)
                if r.status_code == 404:
                    break  # try next tier
                r.raise_for_status()
                return parse_ae_html(r.text, year, month, url)
            except Exception as e:
                if attempt < 2:
                    time.sleep(2)
                else:
                    logger.debug(f"  {year}-{month:02d} [{url}]: {e}")
    return []

def parse_ae_html(html: str, year: int, month: int, source_url: str) -> list[dict]:
    """
    Parse Kyoto AE HTML page.  The data appears as a <pre> block with lines like:
      AE   2023  JAN    1   XXX  XXX  XXX  XXX ... (hourly values)
    or as formatted text with 24 hourly values per day.
    """
    soup = BeautifulSoup(html, "html.parser")
    pre = soup.find("pre")
    if not pre:
        return []

    text = pre.get_text()
    records = []

    # Kyoto AE pages typically format as:
    # Day-line with 24 hourly values after the date prefix
    # Pattern depends on sub-page type; try a few approaches:

    # Format 1: "YYYYMM  DD  H00  H01 ... H23"
    # Format 2: Lines like "AE 2020 JAN  1   23  18  15 ..."
    # Use a generic approach: find all integers after date tokens

    lines = text.split("\n")
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("Day"):
            continue
        # Match lines with year, possibly month abbreviation, then day, then numbers
        # Example: " 1   23  18  15  12  10   8  14  25  35  42  18  16 ..."
        # or: "AE 2020 JAN  1  23  18 ..."
        # Try to find day number followed by 24 hourly values
        nums = re.findall(r"-?\d+", line)
        if len(nums) < 6:
            continue

        # Try to identify day of month
        # If line contains a recognizable date structure
        # Approach: look for lines where first significant number is 1-31 (day)
        # preceded by month abbrev or year
        day = None
        start_idx = 0

        # Check if line starts with month abbrev + day
        m = re.match(r"[A-Za-z]+\s+(\d{1,2})\s+([\d\s]+)", line)
        if m:
            day = int(m.group(1))
            val_str = m.group(2)
            hourly = re.findall(r"-?\d+", val_str)
        else:
            # Check for lines like "  1   xxx xxx xxx ..."
            m2 = re.match(r"\s*(\d{1,2})\s+((?:\s*-?\d+){24,})", line)
            if m2:
                day = int(m2.group(1))
                hourly = re.findall(r"-?\d+", m2.group(2))
            else:
                continue

        if day is None or day < 1 or day > 31:
            continue
        if len(hourly) < 24:
            continue

        for hour in range(24):
            val = int(hourly[hour])
            if abs(val) > 9990:  # fill value
                val = None
            try:
                dt = datetime(year, month, day, hour)
                records.append({"datetime_utc": dt.isoformat() + "Z",
                                "ae_nT": val})
            except ValueError:
                pass

    return records

def main():
    logger.info("=== Kyoto AE Index Download ===")

    if OUT_FILE.exists():
        existing = pd.read_csv(OUT_FILE)
        logger.info(f"  Resuming from existing {len(existing):,} records")
    else:
        existing = pd.DataFrame()

    session = requests.Session()
    session.headers["User-Agent"] = "solar-research/1.0 (academic)"

    all_records = []
    total_months = 0
    empty_months = 0

    start_year = 1957
    end_year = datetime.now().year
    end_month = datetime.now().month

    for year in range(start_year, end_year + 1):
        m_end = 12 if year < end_year else end_month
        year_records = 0
        for month in range(1, m_end + 1):
            recs = fetch_ae_month(year, month, session)
            if recs:
                all_records.extend(recs)
                year_records += len(recs)
                total_months += 1
            else:
                empty_months += 1
            time.sleep(0.2)
        if year_records > 0:
            logger.info(f"  {year}: {year_records} hourly records")
        else:
            logger.debug(f"  {year}: no data")

    if all_records:
        df = pd.DataFrame(all_records)
        df["ae_nT"] = pd.to_numeric(df["ae_nT"], errors="coerce")
        df = df.dropna(subset=["ae_nT"])
        df.sort_values("datetime_utc", inplace=True)
        df.drop_duplicates(subset=["datetime_utc"], inplace=True)
        df.to_csv(OUT_FILE, index=False)
        logger.info(f"\n✓  Kyoto AE → {OUT_FILE.name}")
        logger.info(f"   {len(df):,} records | {total_months} months OK | {empty_months} empty")
    else:
        logger.warning("No AE records retrieved — check Kyoto URL structure")
        # Save instructions
        instr = OUT_DIR / "AE_INDEX_INSTRUCTIONS.txt"
        instr.write_text(
            "Kyoto AE index is available at:\n"
            "  https://wdc.kugi.kyoto-u.ac.jp/ae_realtime/YYYYMM/index.html\n"
            "  https://wdc.kugi.kyoto-u.ac.jp/ae_provisional/YYYYMM/index.html\n"
            "  https://wdc.kugi.kyoto-u.ac.jp/ae_final/YYYYMM/index.html\n\n"
            "Note: AE data 1963-2025 is also embedded in OMNI2 (column 42).\n"
            "Use data/atmospheric/omni_solar_wind/low_res/omni2_YYYY.dat column 42.\n"
        )
        logger.info(f"  Instructions written to {instr.name}")

if __name__ == "__main__":
    main()
