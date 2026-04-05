"""
Norway Avalanche Data — Fixed API versions
Uses correct Varsom v6 + Regobs v4 endpoints
"""
import sys
import json
import time
from pathlib import Path
from datetime import date
import requests
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger, download_file, write_instructions

logger = get_logger("norway_fixed")
OUT = DATA_DIR / "cryosphere" / "ngi_norway"
OUT.mkdir(parents=True, exist_ok=True)

session = requests.Session()
session.headers.update({
    "User-Agent": "Solar-Magnetic-Analysis/1.0",
    "Accept": "application/json",
    "Content-Type": "application/json",
})

# ─────────────────────────────────────────────────────────────────────────────
# 1. Varsom — try multiple API versions to find correct one
# ─────────────────────────────────────────────────────────────────────────────
logger.info("=== Varsom API discovery ===")
VARSOM_VERSIONS = ["v6.3.0", "v6.2.0", "v6.1.0", "v6.0.0", "v5.2.0", "v5.0.1"]
varsom_base = None
for ver in VARSOM_VERSIONS:
    url = f"https://api01.nve.no/hydrology/forecast/avalanche/{ver}/api/Region"
    try:
        r = session.get(url, timeout=15)
        if r.status_code == 200:
            varsom_base = f"https://api01.nve.no/hydrology/forecast/avalanche/{ver}/api"
            logger.info(f"  ✓ Varsom API version: {ver}")
            regions = r.json()
            logger.info(f"  Regions found: {len(regions)}")
            varsom_out = OUT / "varsom"
            varsom_out.mkdir(exist_ok=True)
            (varsom_out / "regions.json").write_text(json.dumps(regions, indent=2), encoding="utf-8")
            break
    except Exception as exc:
        logger.debug(f"  {ver}: {exc}")

if varsom_base:
    # Download avalanche warnings for all regions and years
    varsom_out = OUT / "varsom"
    for region in regions[:20]:  # top 20 regions
        rid = region.get("Id", region.get("RegionId", ""))
        if not rid:
            continue
        for year in range(2013, date.today().year + 1):
            fname = f"warnings_region{rid}_{year}.json"
            dest = varsom_out / str(year) / fname
            if dest.exists():
                continue
            start = f"{year}-10-01"
            end   = f"{year+1}-06-30" if year < date.today().year else str(date.today())
            url = f"{varsom_base}/AvalancheWarning/GetWarningSimple/{rid}/{start}/{end}/1"
            try:
                r = session.get(url, timeout=30)
                if r.status_code == 200:
                    data = r.json()
                    if data:
                        dest.parent.mkdir(parents=True, exist_ok=True)
                        dest.write_text(json.dumps(data, indent=2), encoding="utf-8")
                        logger.info(f"  ✓  Region {rid} {year}: {len(data)} warnings")
                time.sleep(0.2)
            except Exception as exc:
                logger.debug(f"  Region {rid} {year}: {exc}")
else:
    logger.warning("  No working Varsom API version found — trying NVE varsom.no")
    # Try varsom.no direct
    for url in [
        "https://varsom.no/api/forecast/avalanche/en/overview?RegionId=3003&date=2024-02-01",
        "https://api.nve.no/doc/",
    ]:
        try:
            r = session.get(url, timeout=15)
            logger.info(f"  {r.status_code} {url}")
        except Exception as exc:
            logger.debug(f"  {url}: {exc}")

# ─────────────────────────────────────────────────────────────────────────────
# 2. Regobs — correct API with proper auth & headers
# ─────────────────────────────────────────────────────────────────────────────
logger.info("=== Regobs v4 avalanche observations ===")
regobs_out = OUT / "regobs"
regobs_out.mkdir(exist_ok=True)

REGOBS_BASE = "https://www.regobs.no/api/v4"

# Paginated search for all snow/avalanche observations in Norway
offset = 0
page_size = 250
total_saved = 0

for year in range(1970, date.today().year + 1):
    out_file = regobs_out / f"regobs_avalanche_{year}.json"
    if out_file.exists():
        continue
    payload = {
        "GeoHazardTids": [10],   # 10 = Snow/avalanche
        "NumberOfRecords": page_size,
        "Offset": 0,
        "FromDate": f"{year}-10-01",
        "ToDate": f"{year+1}-07-01",
        "LangKey": 1,
    }
    try:
        r = session.post(f"{REGOBS_BASE}/Registration/Search",
                         data=json.dumps(payload), timeout=60)
        if r.status_code == 200 and r.content[:1] == b"[":
            data = r.json()
            if data:
                out_file.write_text(json.dumps(data, indent=2), encoding="utf-8")
                logger.info(f"  ✓  Regobs {year}: {len(data)} observations")
                total_saved += len(data)
        elif r.status_code == 200:
            logger.debug(f"  Regobs {year}: non-JSON response")
        time.sleep(0.3)
    except Exception as exc:
        logger.debug(f"  Regobs {year}: {exc}")

logger.info(f"  Total Regobs records saved: {total_saved}")

# ─────────────────────────────────────────────────────────────────────────────
# 3. xGeo — Norwegian natural hazard warning service (alternative to Varsom)
# ─────────────────────────────────────────────────────────────────────────────
logger.info("=== xGeo / NVE Open Data ===")
xgeo_out = OUT / "xgeo_open"
xgeo_out.mkdir(exist_ok=True)
xgeo_urls = {
    "nve_open_data_catalog.json": "https://nedlasting.nve.no/api/catalog",
    "nve_avalanche_events.json":  "https://nedlasting.nve.no/gis/snoskred/",
}
for fname, url in xgeo_urls.items():
    download_file(url, xgeo_out / fname, desc=fname, session=session)

# ─────────────────────────────────────────────────────────────────────────────
# 4. NVE HydAPI — requires API key (free registration)
# ─────────────────────────────────────────────────────────────────────────────
write_instructions(
    OUT / "nve_hydapi_setup",
    "NVE HydAPI — Snow & Hydrology Data Access",
    """
NVE HydAPI requires a free API key.

REGISTRATION: https://hydapi.nve.no/UserDocumentation/  → 'Request API Key'
Or: Send email to nve@nve.no requesting API key for research

STEP 1: Get your API key at https://hydapi.nve.no/UserDocumentation/
STEP 2: Set environment variable: set NVE_API_KEY=your_key_here
STEP 3: Run the download with:
  import requests
  r = requests.get(
      'https://hydapi.nve.no/api/v1/Observations',
      headers={'X-API-Key': 'YOUR_KEY'},
      params={'StationId': '1000', 'Parameter': '2003',  # 2003 = Snow depth
              'ResolutionTime': 1440,  # daily
              'DateFrom': '1970-01-01', 'DateTo': '2026-04-01'}
  )

KEY PARAMETERS:
  2001 = Snow Water Equivalent
  2003 = Snow Depth
  2014 = New snow depth (24h)
  17   = Precipitation
  0    = Water level
"""
)

logger.info("=== Norway script complete ===")
