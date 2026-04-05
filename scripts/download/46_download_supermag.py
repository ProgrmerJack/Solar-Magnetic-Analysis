"""
SuperMAG geomagnetic indices and station data download.

SuperMAG provides:
  1. SMR index (ring current, better Dst) 1970-present
  2. SML/SMU/SME indices (auroral, better AE) 1970-present
  3. Individual station data (500+ stations) for SEA analysis

Registration required at: https://supermag.jhuapl.edu/
  - Go to the site, click "Sign In / Register"
  - Enter your email as the logon
  - Check email for confirmation

Set SUPERMAG_LOGON environment variable or update LOGON below.
"""
import os
import time
import gzip
import json
import requests
import pandas as pd
from pathlib import Path
from datetime import datetime, timedelta
import sys
sys.path.insert(0, str(Path(__file__).parents[1]))
from utils import get_logger

LOG = get_logger("46_supermag")

# ── Configure logon ──────────────────────────────────────────────────────────
LOGON = os.environ.get("SUPERMAG_LOGON", "andersonmark578")
BASE_URL = "https://supermag.jhuapl.edu/services"

OUT = Path(__file__).parents[2] / "data" / "geomagnetic" / "supermag"
OUT.mkdir(parents=True, exist_ok=True)

# ── Verify logon ─────────────────────────────────────────────────────────────
LOG.info("=== SuperMAG download (logon: %s) ===", LOGON)

def test_logon():
    r = requests.get(
        f"{BASE_URL}/data-api.php",
        params={"fmt": "json", "logon": LOGON, "start": "2010-01-01T00:00:00",
                "extent": 60, "station": "ABK", "delta": "baseline"},
        timeout=20
    )
    return r.status_code == 200 and "Invalid username" not in r.text

if not test_logon():
    LOG.error("SuperMAG logon '%s' is invalid.", LOGON)
    LOG.error("Register at https://supermag.jhuapl.edu/ then set:")
    LOG.error("  $env:SUPERMAG_LOGON='your@email.com'")
    LOG.error("  python scripts/download/46_download_supermag.py")
    sys.exit(1)

LOG.info("✓ Logon validated")

# ── Helper ───────────────────────────────────────────────────────────────────
def api_get(endpoint, params, retries=3, delay=5):
    for attempt in range(retries):
        try:
            r = requests.get(f"{BASE_URL}/{endpoint}", params={**params, "logon": LOGON, "fmt": "json"},
                             timeout=60)
            if r.status_code == 200 and not r.text.startswith("ERROR"):
                return r
            LOG.warning("  attempt %d/%d: %s", attempt + 1, retries, r.text[:80])
        except Exception as exc:
            LOG.warning("  attempt %d/%d exc: %s", attempt + 1, retries, exc)
        time.sleep(delay)
    return None

# ── 1. Indices: SMR, SML, SMU, SME — daily means 1970–2025 ───────────────────
LOG.info("=== SMR / SML / SMU / SME indices (daily) ===")

indices_out = OUT / "indices"
indices_out.mkdir(exist_ok=True)

INDICES_START = datetime(1970, 1, 1)
INDICES_END = datetime(2025, 12, 31)

chunk_days = 30  # 30-day request chunks to stay within API limits

rows = []
current = INDICES_START
while current <= INDICES_END:
    end_chunk = min(current + timedelta(days=chunk_days), INDICES_END)
    extent_secs = int((end_chunk - current).total_seconds())

    r = api_get("indices.php", {
        "start": current.strftime("%Y-%m-%dT%H:%M:%S"),
        "extent": extent_secs,
        "indices": "SMR,SML,SMU,SME",
    })
    if r:
        try:
            data = r.json()
            if isinstance(data, list):
                rows.extend(data)
            elif isinstance(data, dict) and "data" in data:
                rows.extend(data["data"])
        except Exception:
            LOG.warning("  JSON parse error for %s", current.strftime("%Y-%m"))

    if current.month % 3 == 0 and current.day == 1:
        LOG.info("  progress: %s, rows so far: %d", current.strftime("%Y-%m"), len(rows))

    current = end_chunk + timedelta(seconds=1)
    time.sleep(0.5)

if rows:
    df = pd.DataFrame(rows)
    out_file = indices_out / "supermag_indices_1970_2025.csv"
    df.to_csv(out_file, index=False)
    LOG.info("✓ Indices saved: %d rows → %s", len(df), out_file)
else:
    LOG.warning("No index data retrieved — check logon registration status")

# ── 2. Station inventory ─────────────────────────────────────────────────────
LOG.info("=== Station inventory ===")
r_inv = api_get("inventory.php", {
    "start": "2010-01-01T00:00:00",
    "extent": 86400,
})
if r_inv:
    try:
        inv = r_inv.json()
        inv_file = OUT / "station_inventory.json"
        with open(inv_file, "w") as f:
            json.dump(inv, f, indent=2)
        n_stations = len(inv) if isinstance(inv, list) else len(inv.get("stations", []))
        LOG.info("✓ Inventory saved: %d stations → %s", n_stations, inv_file)
    except Exception as exc:
        LOG.warning("  Inventory parse error: %s", exc)

# ── 3. Polar station data for EPP-geomagnetic coupling analysis ───────────────
#    Download daily-resolution data for key high-latitude stations (60-90°N)
#    that are most sensitive to energetic particle precipitation effects.
LOG.info("=== Polar station daily data (geomagnetic coupling) ===")

POLAR_STATIONS = [
    "ABK",   # Abisko, Sweden          (68.4°N)
    "TRO",   # Tromsø, Norway           (69.7°N)
    "BJN",   # Bjørnøya, Norway         (74.5°N)
    "LYR",   # Longyearbyen, Svalbard   (78.2°N)
    "NAL",   # Ny-Ålesund, Svalbard     (78.9°N)
    "THL",   # Thule, Greenland         (77.5°N)
    "GDH",   # Qeqertarsuaq, Greenland  (69.2°N)
    "BRW",   # Barrow (Utqiaġvik), AK   (71.3°N)
    "CMO",   # College, Alaska          (64.9°N)
    "YKC",   # Yellowknife, Canada      (62.5°N)
    "MEA",   # Meanook, Canada          (54.6°N)
    "OTT",   # Ottawa, Canada           (45.4°N — mid-lat reference)
]

station_out = OUT / "polar_stations"
station_out.mkdir(exist_ok=True)

STATION_START = datetime(1970, 1, 1)
STATION_END = datetime(2025, 12, 31)

for station in POLAR_STATIONS:
    out_file = station_out / f"{station}_1970_2025.csv.gz"
    if out_file.exists():
        LOG.info("  %s: already exists, skipping", station)
        continue

    LOG.info("  Downloading station %s ...", station)
    station_rows = []
    current = STATION_START
    while current <= STATION_END:
        end_chunk = min(current + timedelta(days=30), STATION_END)
        extent_secs = int((end_chunk - current).total_seconds())

        r = api_get("data-api.php", {
            "start": current.strftime("%Y-%m-%dT%H:%M:%S"),
            "extent": extent_secs,
            "station": station,
            "delta": "baseline",
            "flagging": "1",
        })
        if r:
            try:
                data = r.json()
                if isinstance(data, list):
                    station_rows.extend(data)
                elif isinstance(data, dict) and "data" in data:
                    station_rows.extend(data["data"])
            except Exception:
                pass
        current = end_chunk + timedelta(seconds=1)
        time.sleep(0.3)

    if station_rows:
        df = pd.DataFrame(station_rows)
        with gzip.open(out_file, "wt") as f:
            df.to_csv(f, index=False)
        LOG.info("  ✓ %s: %d records → %s", station, len(df), out_file.name)
    else:
        LOG.warning("  %s: no data retrieved", station)

LOG.info("=== SuperMAG download complete ===")
LOG.info("Output: %s", OUT)
