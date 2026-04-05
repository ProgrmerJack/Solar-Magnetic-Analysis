"""
Script 07 — NRCS SNOTEL Snowpack Telemetry Network
Downloads snowpack data from the USDA NRCS SNOTEL network
(900+ stations across the western United States, ~1980–present).

API: https://wcc.sc.egov.usda.gov/awdbRestApi/swagger-ui/index.html

Downloads:
  • Full station metadata list
  • Snow Water Equivalent (SWE) daily data for all active stations
  • Snow depth daily data
  • Precipitation data
  • Air temperature
"""
import sys
import json
import time
from pathlib import Path
from datetime import date, timedelta
import requests
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger, download_file

logger = get_logger("07_snotel")
OUT = DATA_DIR / "cryosphere" / "snotel"
OUT.mkdir(parents=True, exist_ok=True)

BASE = "https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1"
session = requests.Session()
session.headers.update({"User-Agent": "Solar-Magnetic-Analysis/1.0"})


def api_get(endpoint: str, params: dict = None) -> dict | list | None:
    try:
        r = session.get(f"{BASE}/{endpoint}", params=params, timeout=60)
        r.raise_for_status()
        return r.json()
    except Exception as exc:
        logger.warning(f"  API error {endpoint}: {exc}")
        return None


# --------------------------------------------------------------------------- #
# 1. Station metadata (all SNOTEL stations)                                   #
# --------------------------------------------------------------------------- #
logger.info("=== SNOTEL Station Metadata ===")
meta_file = OUT / "snotel_stations.json"
if not meta_file.exists():
    stations = api_get("stations", params={
        "networkCds": "SNTL",
        "activeOnly": "false",
        "returnObjects": "true",
    })
    if stations:
        meta_file.write_text(json.dumps(stations, indent=2), encoding="utf-8")
        logger.info(f"  ✓  {len(stations)} stations metadata saved")
    else:
        logger.warning("  Could not retrieve station metadata")
        stations = []
else:
    stations = json.loads(meta_file.read_text(encoding="utf-8"))
    logger.info(f"  SKIP  snotel_stations.json  ({len(stations)} stations, already exists)")


# --------------------------------------------------------------------------- #
# 2. Bulk SWE / Snow Depth / Temperature data per state                       #
# --------------------------------------------------------------------------- #
logger.info("=== SNOTEL Bulk Data Download ===")

ELEMENTS = {
    "WTEQ": "Snow Water Equivalent (in)",
    "SNWD": "Snow Depth (in)",
    "PREC": "Accumulated Precipitation (in)",
    "TMAX": "Max Air Temperature (degF)",
    "TMIN": "Min Air Temperature (degF)",
    "TAVG": "Mean Air Temperature (degF)",
}

# Western US states with significant SNOTEL coverage
STATES = ["AK", "AZ", "CA", "CO", "ID", "MT", "NM", "NV", "OR", "UT", "WA", "WY"]

bulk_out = OUT / "bulk_data"
bulk_out.mkdir(exist_ok=True)

start_date = "1980-01-01"
end_date   = date.today().strftime("%Y-%m-%d")

for state in STATES:
    state_out = bulk_out / state
    state_out.mkdir(exist_ok=True)
    for elem_cd, elem_name in ELEMENTS.items():
        out_file = state_out / f"{state}_{elem_cd}_daily.json"
        if out_file.exists():
            logger.info(f"  SKIP  {state}/{elem_cd}")
            continue
        data = api_get("data", params={
            "stationTriplets": f"*:{state}:SNTL",
            "elementCd":       elem_cd,
            "beginDate":       start_date,
            "endDate":         end_date,
            "duration":        "DAILY",
            "getFlags":        "true",
        })
        if data:
            out_file.write_text(json.dumps(data, indent=2), encoding="utf-8")
            logger.info(f"  ✓  {state}/{elem_cd} ({elem_name})")
        else:
            logger.warning(f"  ✗  {state}/{elem_cd} — no data returned")
        time.sleep(0.5)  # polite rate limiting


# --------------------------------------------------------------------------- #
# 3. Station-level metadata with coordinates (for GIS mapping)                #
# --------------------------------------------------------------------------- #
logger.info("=== SNOTEL Station Details with Coordinates ===")
if stations:
    station_details = []
    detail_file = OUT / "snotel_station_details.json"
    if not detail_file.exists():
        for i, sta in enumerate(stations[:50]):  # sample first 50 for testing
            triplet = sta.get("stationTriplet", "")
            if not triplet:
                continue
            detail = api_get(f"stations/{triplet}", params={"includeFlags": "false"})
            if detail:
                station_details.append(detail)
            if i % 10 == 0:
                logger.info(f"  {i}/{len(stations)} stations detailed...")
            time.sleep(0.2)
        detail_file.write_text(json.dumps(station_details, indent=2), encoding="utf-8")
        logger.info(f"  ✓  Station details: {len(station_details)} records")


# --------------------------------------------------------------------------- #
# 4. Current season snapshot (all states, SWE)                                #
# --------------------------------------------------------------------------- #
logger.info("=== Current Season SWE Snapshot ===")
current_year = date.today().year
current_season_start = f"{current_year - 1}-10-01"

snapshot_data = api_get("data", params={
    "stationTriplets": "*:*:SNTL",
    "elementCd":       "WTEQ",
    "beginDate":       current_season_start,
    "endDate":         end_date,
    "duration":        "DAILY",
})
if snapshot_data:
    snap_file = OUT / f"snotel_swe_season_{current_year}.json"
    snap_file.write_text(json.dumps(snapshot_data, indent=2), encoding="utf-8")
    logger.info(f"  ✓  Current season SWE snapshot saved")


logger.info("=== Script 07 complete ===")
