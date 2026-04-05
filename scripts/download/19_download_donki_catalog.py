"""
DONKI Flare Catalog — NASA CCMC (complete, 1995–present)
Uses the DONKI REST API: https://kauai.ccmc.gsfc.nasa.gov/DONKI/WS/get/FLR
Returns JSON with classification, timing, location, active region number.
Also downloads SWPC edited_events.json (recent 60-day cleaned event list).
"""
import sys, json, time
from pathlib import Path
from datetime import date, timedelta
import requests
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger

logger = get_logger("donki_flares")
OUT = DATA_DIR / "solar" / "flare_catalog"
OUT.mkdir(parents=True, exist_ok=True)

DONKI = "https://kauai.ccmc.gsfc.nasa.gov/DONKI/WS/get"
session = requests.Session()
session.headers["User-Agent"] = "Solar-Magnetic-Analysis/1.0"
session.headers["Accept"]     = "application/json"

def get_json(url, desc=""):
    for attempt in range(4):
        try:
            r = session.get(url, timeout=60)
            r.raise_for_status()
            return r.json()
        except Exception as exc:
            logger.warning(f"  Attempt {attempt+1}/4 ({desc}): {exc}")
            time.sleep(8 * (attempt + 1))
    return None

# ── DONKI FLR (Solar Flares) — monthly chunks ─────────────────────────────────
logger.info("=== DONKI Solar Flare Catalog (1995–present) ===")
flares_dir = OUT / "donki_monthly"
flares_dir.mkdir(exist_ok=True)

all_flares = []
start = date(1995, 1, 1)
today = date.today()
current = start

while current <= today:
    end_m = (current.replace(day=28) + timedelta(days=4)).replace(day=1) - timedelta(days=1)
    end_m = min(end_m, today)
    label = current.strftime("%Y-%m")
    out_file = flares_dir / f"flares_{label}.json"

    if not out_file.exists():
        url = (f"{DONKI}/FLR?"
               f"startDate={current.strftime('%Y-%m-%d')}"
               f"&endDate={end_m.strftime('%Y-%m-%d')}")
        data = get_json(url, desc=label)
        if data is not None:
            out_file.write_text(json.dumps(data, indent=2), encoding="utf-8")
            if data:
                logger.info(f"  {label}: {len(data)} flares")
            all_flares.extend(data or [])
        time.sleep(0.8)
    else:
        try:
            month_data = json.loads(out_file.read_text())
            all_flares.extend(month_data or [])
        except Exception:
            pass

    # Advance to first day of next month
    current = (current.replace(day=28) + timedelta(days=4)).replace(day=1)

# Save merged catalog
if all_flares:
    merged = OUT / "donki_flares_1995_present.json"
    merged.write_text(json.dumps(all_flares, indent=2), encoding="utf-8")
    logger.info(f"  ✓  Merged: {len(all_flares)} flares → {merged.name} "
                f"({merged.stat().st_size//1024} KB)")

# ── DONKI CME catalog too ────────────────────────────────────────────────────
logger.info("=== DONKI CME Catalog (1995–present) ===")
cme_dir = OUT / "donki_cme_monthly"
cme_dir.mkdir(exist_ok=True)
all_cme = []
current = date(1995, 1, 1)

while current <= today:
    end_m = (current.replace(day=28) + timedelta(days=4)).replace(day=1) - timedelta(days=1)
    end_m = min(end_m, today)
    label = current.strftime("%Y-%m")
    out_file = cme_dir / f"cme_{label}.json"

    if not out_file.exists():
        url = (f"{DONKI}/CME?"
               f"startDate={current.strftime('%Y-%m-%d')}"
               f"&endDate={end_m.strftime('%Y-%m-%d')}")
        data = get_json(url, desc=f"CME {label}")
        if data is not None:
            out_file.write_text(json.dumps(data, indent=2), encoding="utf-8")
            if data:
                logger.info(f"  CME {label}: {len(data)} events")
            all_cme.extend(data or [])
        time.sleep(0.8)
    else:
        try:
            all_cme.extend(json.loads(out_file.read_text()) or [])
        except Exception:
            pass

    current = (current.replace(day=28) + timedelta(days=4)).replace(day=1)

if all_cme:
    merged_cme = OUT / "donki_cme_1995_present.json"
    merged_cme.write_text(json.dumps(all_cme, indent=2), encoding="utf-8")
    logger.info(f"  ✓  CME merged: {len(all_cme)} events → {merged_cme.name}")

# ── DONKI Geomagnetic Storm catalog ─────────────────────────────────────────
logger.info("=== DONKI GST (Geomagnetic Storm) Catalog ===")
gst_url = f"{DONKI}/GST?startDate=1995-01-01&endDate={today}"
gst = get_json(gst_url, desc="GST all")
if gst:
    (OUT / "donki_gst_1995_present.json").write_text(json.dumps(gst, indent=2))
    logger.info(f"  ✓  {len(gst)} geomagnetic storms")

# ── SWPC edited_events (recent 60-day cleaned catalog) ──────────────────────
logger.info("=== SWPC edited_events.json (recent 60-day) ===")
r = session.get("https://services.swpc.noaa.gov/json/edited_events.json", timeout=30)
if r.status_code == 200:
    (OUT / "swpc_edited_events.json").write_bytes(r.content)
    logger.info(f"  ✓  edited_events.json: {len(r.content)//1024} KB")

logger.info("=== DONKI flare/CME/GST catalog download complete ===")
