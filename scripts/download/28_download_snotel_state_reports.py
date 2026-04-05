"""
SNOTEL State-Level SWE Reports — Corrected per-station approach
Downloads individual SNOTEL station daily data via NRCS Report Generator CSV format
(simpler flat CSV vs. AWDB JSON) for key high-elevation stations in avalanche terrain.

Focuses on stations above 9,000 ft in CO, WA, MT, ID, UT, WY — prime avalanche terrain.
Uses the customSingleStationReport endpoint which returns clean daily CSVs.
"""
import sys
import time
import json
import requests
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger

logger = get_logger("snotel_states")
OUT = DATA_DIR / "cryosphere" / "snotel" / "state_reports"
OUT.mkdir(parents=True, exist_ok=True)

session = requests.Session()
session.headers.update({"User-Agent": "Solar-Magnetic-Analysis/1.0 (research)"})

STATES = ["CO", "WA", "MT", "ID", "UT", "WY", "OR", "CA"]
MIN_ELEV_FT = {  # Minimum elevation filter per state (ft) — focus on high alpine terrain
    "CO": 9000, "WA": 4500, "MT": 5500, "ID": 5500,
    "UT": 8000, "WY": 8000, "OR": 4500, "CA": 7000,
}
ELEMENTS = "WTEQ::value,SNWD::value,PREC::value,TMAX::value,TMIN::value"
REPORT_BASE = "https://wcc.sc.egov.usda.gov/reportGenerator/view_csv/customSingleStationReport/daily/start_of_period"

logger.info("=== SNOTEL High-Elevation Station Reports ===")

for state in STATES:
    state_dir = OUT / state
    state_dir.mkdir(exist_ok=True)
    min_elev = MIN_ELEV_FT.get(state, 7000)

    # Get station list for this state
    r = session.get("https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1/stations",
                    params={"stateCode": state, "networkCd": "SNTL",
                            "returnForecastPointData": "false"}, timeout=30)
    if r.status_code != 200:
        logger.warning(f"  {state}: station list HTTP {r.status_code}")
        continue
    stations = [x for x in json.loads(r.text)
                if x.get("stateCode") == state and x.get("networkCode") == "SNTL"
                and (x.get("elevation") or 0) >= min_elev]
    logger.info(f"  {state}: {len(stations)} stations above {min_elev} ft")

    for sta in stations:
        triplet = sta["stationTriplet"]
        sta_id = triplet.split(":")[0]
        out_file = state_dir / f"{sta_id}_daily.csv"
        if out_file.exists():
            continue

        url = f"{REPORT_BASE}/{triplet}/POR_BEGIN,POR_END/{ELEMENTS}"
        try:
            r = session.get(url, timeout=90)
            if r.status_code == 200 and len(r.content) > 500:
                out_file.write_bytes(r.content)
        except Exception as exc:
            logger.debug(f"  {triplet}: {exc}")
        time.sleep(0.2)

    saved = len(list(state_dir.glob("*.csv")))
    logger.info(f"    {state}: {saved} CSV files saved")

logger.info("=== SNOTEL state reports complete ===")

