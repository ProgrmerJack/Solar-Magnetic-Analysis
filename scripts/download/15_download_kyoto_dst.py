"""
Kyoto Dst Historical Download — HTML scraper
WDC Kyoto provides Dst via monthly HTML pages; this script parses them.
Final (1957-2020), Provisional (2021-2024), Realtime (2025-present)
"""
import re
import csv
import sys
import time
from pathlib import Path
from datetime import date
import requests
from bs4 import BeautifulSoup
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger

logger = get_logger("dst_kyoto")
OUT = DATA_DIR / "geomagnetic" / "dst_ae_index" / "kyoto_dst_historical"
OUT.mkdir(parents=True, exist_ok=True)

session = requests.Session()
session.headers.update({"User-Agent": "Solar-Magnetic-Analysis/1.0"})

BASE = "https://wdc.kugi.kyoto-u.ac.jp"
TIERS = [
    ("dst_final",       1957, 2020),
    ("dst_provisional", 2021, 2024),
    ("dst_realtime",    2025, date.today().year),
]

all_rows = []  # (datetime_str, dst_nT)

def parse_dst_page(html: str, year: int, month: int) -> list[tuple[str, int]]:
    """Extract hourly Dst values from a WDC Kyoto monthly page."""
    soup = BeautifulSoup(html, "lxml")
    rows = []
    # Kyoto pages use a PRE block or TD cells with the hourly values
    # Format: 24 values per day, each value in nT
    pre = soup.find("pre")
    if pre:
        text = pre.get_text()
        # Each line: "DD HH:00  val val val..."
        for line in text.splitlines():
            m = re.match(r"\s*(\d{1,2})\s+(.*)", line)
            if m:
                day = int(m.group(1))
                vals = re.findall(r"-?\d+", m.group(2))
                for hour, val in enumerate(vals[:24]):
                    rows.append((f"{year:04d}-{month:02d}-{day:02d}T{hour:02d}:00:00Z",
                                 int(val)))
    # Fallback: look for table cells with numeric values
    if not rows:
        tds = [td.get_text(strip=True) for td in soup.find_all("td")
               if re.match(r"-?\d{1,4}$", td.get_text(strip=True))]
        # Arrange as 24-hour blocks per day
        day = 1
        for i in range(0, len(tds), 24):
            block = tds[i:i+24]
            for hour, val in enumerate(block):
                rows.append((f"{year:04d}-{month:02d}-{day:02d}T{hour:02d}:00:00Z",
                             int(val)))
            day += 1
    return rows


for tier, start_yr, end_yr in TIERS:
    for year in range(start_yr, end_yr + 1):
        for month in range(1, 13):
            if year == date.today().year and month > date.today().month:
                break
            ym = f"{year:04d}{month:02d}"
            out_file = OUT / f"dst_{tier}_{ym}.csv"
            if out_file.exists():
                continue

            url = f"{BASE}/{tier}/{ym}/index.html"
            try:
                r = session.get(url, timeout=30)
                if r.status_code == 404:
                    continue
                r.raise_for_status()
                rows = parse_dst_page(r.text, year, month)
                if rows:
                    with open(out_file, "w", newline="") as f:
                        w = csv.writer(f)
                        w.writerow(["datetime_utc", "dst_nT"])
                        w.writerows(rows)
                    logger.info(f"  ✓  {tier}/{ym}: {len(rows)} hourly values")
                    all_rows.extend(rows)
                time.sleep(0.3)
            except Exception as exc:
                logger.debug(f"  {tier}/{ym}: {exc}")

# Merge into single CSV
if all_rows:
    all_rows.sort(key=lambda x: x[0])
    merged = OUT.parent / "kyoto_dst_1957_present.csv"
    with open(merged, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["datetime_utc", "dst_nT"])
        w.writerows(all_rows)
    logger.info(f"  ✓  Merged Dst: {len(all_rows)} records → {merged}")

logger.info("Done.")
