"""
SNOTEL Bulk Download — Fixed version
Uses station-by-station queries with correct NRCS AWDB REST API format.
Gets SWE + snow depth + precipitation for key high-elevation stations.
"""
import sys
import json
import time
from pathlib import Path
from datetime import date
import requests
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger

logger = get_logger("snotel_fixed")
OUT = DATA_DIR / "cryosphere" / "snotel"
OUT.mkdir(parents=True, exist_ok=True)

BASE  = "https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1"
session = requests.Session()
session.headers.update({"User-Agent": "Solar-Magnetic-Analysis/1.0"})

# Load station list saved by script 07
stations_file = OUT / "snotel_stations.json"
if not stations_file.exists():
    logger.error("Run 07_download_snotel.py first to get station list")
    raise SystemExit(1)

stations = json.loads(stations_file.read_text())
logger.info(f"Loaded {len(stations)} SNOTEL stations")

# ── Key elements ─────────────────────────────────────────────────────────────
ELEMENTS = ["WTEQ", "SNWD", "PREC", "TMAX", "TMIN", "TAVG"]
START = "1980-01-01"
END   = str(date.today())

# ── Filter to stations with long records in key mountain regions ──────────────
# Focus on Colorado, Utah, Idaho, Montana, Washington, Wyoming, Oregon (priority)
PRIORITY_STATES = {"CO", "UT", "ID", "MT", "WA", "WY", "OR", "CA", "NV", "AK"}

priority = [s for s in stations
            if s.get("stationTriplet", "").split(":")[1] in PRIORITY_STATES]
logger.info(f"Priority stations: {len(priority)}")

# ── Download per-station data (individual requests) ───────────────────────────
bulk_out = OUT / "station_data"
bulk_out.mkdir(exist_ok=True)

errors = 0
for i, sta in enumerate(priority):
    triplet = sta.get("stationTriplet", "")
    if not triplet:
        continue

    parts = triplet.split(":")
    state = parts[1] if len(parts) > 1 else "XX"
    sta_id = parts[0]

    sta_out = bulk_out / state / sta_id
    done_file = sta_out / "_complete.flag"
    if done_file.exists():
        continue

    sta_out.mkdir(parents=True, exist_ok=True)
    success = True

    for elem in ELEMENTS:
        out_file = sta_out / f"{elem}_daily.json"
        if out_file.exists():
            continue
        try:
            r = session.get(f"{BASE}/data", params={
                "stationTriplets": triplet,
                "elements":        elem,
                "beginDate":       START,
                "endDate":         END,
                "duration":        "DAILY",
                "getFlags":        "false",
            }, timeout=45)
            if r.status_code == 200:
                data = r.json()
                if data:
                    out_file.write_text(json.dumps(data, indent=2), encoding="utf-8")
            elif r.status_code != 404:
                success = False
            time.sleep(0.15)
        except Exception as exc:
            logger.debug(f"  {triplet}/{elem}: {exc}")
            success = False
            errors += 1

    if success:
        done_file.touch()

    if (i + 1) % 50 == 0:
        logger.info(f"  Progress: {i+1}/{len(priority)} stations ({errors} errors)")

logger.info(f"SNOTEL station download complete. Errors: {errors}")

# ── Also download the NRCS Report Generator bulk CSV for state-level SWE ─────
logger.info("=== NRCS Report Generator — State SWE CSVs ===")
REPORT_BASE = "https://wcc.sc.egov.usda.gov/reportGenerator/view_csv"
report_out = OUT / "state_reports"
report_out.mkdir(exist_ok=True)

for state in sorted(PRIORITY_STATES):
    out_file = report_out / f"{state}_SNTL_SWE_daily_all.csv"
    if out_file.exists():
        continue
    # Report Generator URL for all SNOTEL stations in a state, SWE, full period
    url = (f"{REPORT_BASE}/customMultiTimeSeriesGroupByStationReport/daily/"
           f"start_of_period/{state},SNTL/POR_BEGIN,POR_END/WTEQ::value,SNWD::value,"
           f"PREC::value,TMAX::value,TMIN::value")
    try:
        r = session.get(url, timeout=120)
        if r.status_code == 200 and len(r.content) > 200:
            out_file.write_bytes(r.content)
            logger.info(f"  ✓  {state} SWE report ({len(r.content)//1024} KB)")
        else:
            logger.warning(f"  {state}: status {r.status_code}, {len(r.content)} bytes")
    except Exception as exc:
        logger.warning(f"  {state}: {exc}")
    time.sleep(1)

logger.info("=== SNOTEL fixed script complete ===")
