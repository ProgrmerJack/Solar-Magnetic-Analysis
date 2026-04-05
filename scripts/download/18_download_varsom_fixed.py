"""
Varsom (NVE) Avalanche Warnings — Correct API v6.3.0
Uses Archive/Warning/All endpoint: all 46 regions, full detail, 2012-present.
Downloads by avalanche season (Oct–May) to keep file sizes manageable.
"""
import sys, json, time
from pathlib import Path
from datetime import date
import requests
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger

logger = get_logger("varsom_fixed")
BASE  = "https://api01.nve.no/hydrology/forecast/avalanche/v6.3.0/api"
OUT   = DATA_DIR / "cryosphere" / "ngi_norway" / "varsom"
OUT.mkdir(parents=True, exist_ok=True)

session = requests.Session()
session.headers["User-Agent"] = "Solar-Magnetic-Analysis/1.0"
session.headers["Accept"]     = "application/json"

def get_json(url, desc=""):
    for attempt in range(3):
        try:
            r = session.get(url, timeout=120)
            r.raise_for_status()
            return r.json()
        except Exception as exc:
            logger.warning(f"  Attempt {attempt+1}/3 failed ({desc}): {exc}")
            time.sleep(5 * (attempt + 1))
    return None

# ── Season definitions (Oct–May) ──────────────────────────────────────────────
def season_windows(start_year=2012, end_year=None):
    """Yields (label, start_date, end_date) tuples per season."""
    if end_year is None:
        end_year = date.today().year
    for yr in range(start_year, end_year + 1):
        yield f"{yr}-{yr+1}", f"{yr}-10-01", f"{yr+1}-05-31"

logger.info("=== Varsom Archive/Warning/All (2012–present) ===")
seasons_dir = OUT / "seasons"
seasons_dir.mkdir(exist_ok=True)

for label, start, end in season_windows():
    out_file = seasons_dir / f"warnings_{label}.json"
    if out_file.exists() and out_file.stat().st_size > 1000:
        logger.info(f"  ✓  {label} already saved ({out_file.stat().st_size//1024} KB)")
        continue

    # Archive endpoint: all regions, English (langkey=2), JSON format
    url = f"{BASE}/Archive/Warning/All/2/{start}/{end}/json"
    logger.info(f"  Downloading {label} ({start} → {end}) …")
    data = get_json(url, desc=label)
    if data is None:
        logger.error(f"  ✗  {label} failed")
        continue

    out_file.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    kb = out_file.stat().st_size // 1024
    logger.info(f"  ✓  {label}: {len(data)} warnings, {kb} KB")
    time.sleep(2)  # polite delay between seasons

# ── Also grab per-region simple warnings for quick analysis ──────────────────
logger.info("=== Per-region simple warnings (all regions, 2012–present) ===")
simple_dir = OUT / "simple"
simple_dir.mkdir(exist_ok=True)

# Load region list
regions_file = OUT / "regions.json"
if not regions_file.exists():
    logger.info("  Fetching region list …")
    regions = get_json(f"{BASE}/Region/")
    if regions:
        regions_file.write_text(json.dumps(regions, ensure_ascii=False, indent=2), encoding="utf-8")
else:
    regions = json.loads(regions_file.read_text())

if regions:
    region_ids = [r["Id"] for r in regions]
    logger.info(f"  {len(region_ids)} regions")
    for rid in region_ids:
        out_file = simple_dir / f"region_{rid}_simple_2012_present.json"
        if out_file.exists() and out_file.stat().st_size > 500:
            continue
        url = f"{BASE}/AvalancheWarningByRegion/Simple/{rid}/2/2012-10-01/{date.today()}"
        data = get_json(url, desc=f"region {rid}")
        if data:
            out_file.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
            logger.info(f"  ✓  Region {rid}: {len(data)} warnings")
        time.sleep(0.5)

logger.info("=== Varsom download complete ===")
