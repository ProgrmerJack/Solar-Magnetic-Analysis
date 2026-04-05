"""
GSFC GOES X-ray Solar Flare Event Catalog 1975–2026
All flare classes A/B/C/M/X — from NASA Goddard Space Flight Center
Source: https://hesperia.gsfc.nasa.gov/goes/goes_event_listings/

This is the most comprehensive GOES-era flare catalog, spanning 5+ solar cycles.
Used in Step 1 (SOC analysis: power-law fitting on flare energies).
"""
import requests
import re
import pandas as pd
from pathlib import Path
from datetime import datetime
import sys

sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger

logger = get_logger("gsfc_goes_flare")
OUT = DATA_DIR / "solar" / "flare_catalog" / "gsfc_complete"
OUT.mkdir(parents=True, exist_ok=True)

BASE = "https://hesperia.gsfc.nasa.gov/goes/goes_event_listings"
YEARS = range(1975, 2027)

MONTH_MAP = {m: i+1 for i, m in enumerate(
    ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
)}


def parse_line(line: str) -> dict | None:
    """Parse one event line: ' 3-Jan-2020 21:02   21:03   21:04   A1.0   N22E30   12756'"""
    line = line.strip()
    if not line or not re.match(r'\d+', line):
        return None
    # Format: D-Mon-YYYY HH:MM  HH:MM  HH:MM  CLASS  [POS]  [AR]
    parts = line.split()
    if len(parts) < 5:
        return None
    try:
        date_str, start, peak, end, cls = parts[0], parts[1], parts[2], parts[3], parts[4]
        pos = parts[5] if len(parts) > 5 and re.match(r'[NS]\d', parts[5]) else None
        ar = parts[-1] if len(parts) > 5 and re.match(r'^\d{5}$', parts[-1]) else None

        d, mon, yr = date_str.split("-")
        month_num = MONTH_MAP[mon]
        date = datetime(int(yr), month_num, int(d))
        start_dt = f"{yr}-{month_num:02d}-{int(d):02d}T{start}:00"
        peak_dt = f"{yr}-{month_num:02d}-{int(d):02d}T{peak}:00"

        # Compute log10 peak flux proxy: class letter + number → W/m²
        # A=1e-8, B=1e-7, C=1e-6, M=1e-5, X=1e-4
        flux_base = {"A": 1e-8, "B": 1e-7, "C": 1e-6, "M": 1e-5, "X": 1e-4}
        cls_letter = cls[0].upper()
        cls_number = float(cls[1:]) if len(cls) > 1 else 1.0
        peak_flux = flux_base.get(cls_letter, 1e-8) * cls_number

        return {
            "date": date.date(),
            "start_time_utc": start_dt,
            "peak_time_utc": peak_dt,
            "class": cls,
            "class_letter": cls_letter,
            "position": pos,
            "active_region": ar,
            "peak_flux_Wm2": peak_flux,
        }
    except Exception:
        return None


logger.info("=== GSFC GOES Solar Flare Event Catalog 1975–2026 ===")
all_records = []

for year in YEARS:
    outfile = OUT / f"flares_{year}.txt"
    url = f"{BASE}/goes_xray_event_list_{year}.txt"
    try:
        if not outfile.exists():
            r = requests.get(url, timeout=30)
            r.raise_for_status()
            outfile.write_text(r.text, encoding="utf-8")
        else:
            text = outfile.read_text(encoding="utf-8")
            r_text = text
            r = type("R", (), {"text": r_text})()  # mock

        lines = outfile.read_text().splitlines()
        year_records = [parse_line(l) for l in lines]
        year_records = [r for r in year_records if r is not None]
        all_records.extend(year_records)
        logger.info(f"  {year}: {len(year_records):4d} flares")
    except Exception as exc:
        logger.warning(f"  {year}: {exc}")

df = pd.DataFrame(all_records)
df["date"] = pd.to_datetime(df["date"])
df = df.sort_values("date").reset_index(drop=True)

out_csv = DATA_DIR / "solar" / "flare_catalog" / "gsfc_goes_flares_1975_2026.csv"
df.to_csv(out_csv, index=False)
logger.info(f"\nTotal flares: {len(df):,}")
logger.info(f"  A: {(df.class_letter=='A').sum():,}  B: {(df.class_letter=='B').sum():,}  "
           f"C: {(df.class_letter=='C').sum():,}  M: {(df.class_letter=='M').sum():,}  "
           f"X: {(df.class_letter=='X').sum():,}")
logger.info(f"Saved → {out_csv.name}")
logger.info("=== Done ===")
