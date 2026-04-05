"""
Utah Avalanche Center (UAC) Avalanche Incident Database
All seasons from earliest available (~1980s) to present.
Source: https://utahavalanchecenter.org/api/v2/avalanches?season=YYYY-YYYY

Fields: date, aspect, elevation, avalancheType, trigger, depth, width,
        caught/carried/buried/killed/injured counts, slopeAngle,
        avalancheProblem, weakLayer, lat/lng (if available)
"""
import json
import time
import requests
import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger

logger = get_logger("uac_avalanche")
OUT = DATA_DIR / "cryosphere" / "uac_utah"
OUT.mkdir(parents=True, exist_ok=True)

BASE = "https://utahavalanchecenter.org/api/v2/avalanches"
HDR = {"User-Agent": "research-project-solar-avalanche@academic"}

# Determine seasons: run from 1980-1981 up to current
import datetime
current_year = datetime.date.today().year
seasons = [f"{y}-{y+1}" for y in range(1980, current_year + 1)]


def fetch_season(season: str) -> list | None:
    url = f"{BASE}?season={season}"
    try:
        r = requests.get(url, headers=HDR, timeout=30)
        r.raise_for_status()
        return r.json()
    except Exception as exc:
        logger.warning(f"  {season}: {exc}")
        return None


logger.info("=== Utah Avalanche Center Incident Database ===")
all_records = []

for season in seasons:
    out_file = OUT / f"uac_{season}.json"
    if out_file.exists():
        with open(out_file) as f:
            data = json.load(f)
    else:
        data = fetch_season(season)
        if data is None or len(data) == 0:
            logger.warning(f"  {season}: no data — skipping")
            continue
        with open(out_file, "w") as f:
            json.dump(data, f)

    all_records.extend(data)
    logger.info(f"  {season}: {len(data):4d} events")
    time.sleep(0.3)

if all_records:
    df = pd.DataFrame(all_records)
    # Parse date (format: YYYYMMDD)
    df["date"] = pd.to_datetime(df["date"], format="%Y%m%d", errors="coerce")
    df = df.sort_values("date").reset_index(drop=True)
    out_csv = OUT / "uac_avalanche_incidents_all.csv"
    df.to_csv(out_csv, index=False)
    logger.info(f"\nTotal events: {len(df):,}")
    logger.info(f"Date range: {df.date.min().date()} to {df.date.max().date()}")
    logger.info(f"Saved → {out_csv.name}")

logger.info("=== Done ===")
