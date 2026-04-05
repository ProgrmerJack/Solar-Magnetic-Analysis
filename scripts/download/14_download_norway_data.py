"""
Script 14 — Norway NGI Avalanche Registry Contact + Public NVE Data
The Norwegian Geotechnical Institute (NGI) maintains avalanche occurrence data.
The Norwegian Water Resources and Energy Directorate (NVE) also has snow data.

Public data from NVE (Norges vassdrags- og energidirektorat) API is accessible
without registration and includes snow water equivalent, snow depth, and
avalanche-related hydrology data.
"""
import sys
import json
from pathlib import Path
import requests
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger, download_file, write_instructions

logger = get_logger("14_norway_data")
OUT_NGI   = DATA_DIR / "cryosphere" / "ngi_norway"
OUT_NVE   = DATA_DIR / "cryosphere" / "ngi_norway" / "nve_public"
OUT_NGI.mkdir(parents=True, exist_ok=True)
OUT_NVE.mkdir(exist_ok=True)

session = requests.Session()
session.headers.update({"User-Agent": "Solar-Magnetic-Analysis/1.0"})


# --------------------------------------------------------------------------- #
# 1. NVE (Norwegian Water Resources) Public API — free, no registration       #
# --------------------------------------------------------------------------- #
logger.info("=== NVE Norwegian Snow & Hydrology API ===")
NVE_API = "https://hydapi.nve.no/api/v1"

# NVE HydAPI — stations and observations
nve_endpoints = {
    "nve_stations.json": f"{NVE_API}/Stations?Active=true&ParameterName=Snødybde",
    "nve_snow_stations.json": f"{NVE_API}/Stations?Active=true&ParameterName=Snøvannekvivalent",
    "nve_parameters.json": f"{NVE_API}/Parameters",
}
# Note: NVE API requires API key for observations (free registration)
for fname, url in nve_endpoints.items():
    download_file(url, OUT_NVE / fname, desc=f"NVE {fname}", session=session)


# --------------------------------------------------------------------------- #
# 2. Norwegian Avalanche Warning Service (Varsom) — public API                #
# --------------------------------------------------------------------------- #
logger.info("=== Varsom Norwegian Avalanche Warning Service ===")
VARSOM_API = "https://api01.nve.no/hydrology/forecast/avalanche/v6.2.0/api"
varsom_out = OUT_NGI / "varsom"
varsom_out.mkdir(exist_ok=True)

from datetime import date, timedelta

# Get available regions
varsom_endpoints = {
    "varsom_regions.json": f"{VARSOM_API}/Region",
    "varsom_data_info.json": "https://api01.nve.no/doc/",
}
for fname, url in varsom_endpoints.items():
    download_file(url, varsom_out / fname, desc=f"Varsom {fname}", session=session)

# Download avalanche warnings by region (all Norwegian forecast regions)
# These are issued daily during winter season
logger.info("  Downloading Varsom avalanche bulletins...")
current_year = date.today().year
for year in range(2012, current_year + 1):  # Varsom started ~2012
    for region_id in range(3001, 3050):  # Norwegian region IDs
        url = (f"{VARSOM_API}/AvalancheWarningByRegion/Simple/"
               f"{region_id}/{year}-01-01/{year}-05-31/1")
        dest = varsom_out / str(year) / f"region_{region_id}_{year}_winter.json"
        if not dest.exists():
            try:
                r = session.get(url, timeout=30)
                if r.status_code == 200:
                    data = r.json()
                    if data:  # only save if data exists
                        dest.parent.mkdir(parents=True, exist_ok=True)
                        dest.write_text(json.dumps(data, indent=2), encoding="utf-8")
                        logger.info(f"  ✓  Varsom region {region_id} {year}: {len(data)} entries")
            except Exception as exc:
                logger.debug(f"  Varsom {region_id}/{year}: {exc}")


# --------------------------------------------------------------------------- #
# 3. Regobs Norwegian avalanche observations — public API                     #
# --------------------------------------------------------------------------- #
logger.info("=== Regobs Norwegian Avalanche Observations ===")
regobs_out = OUT_NGI / "regobs"
regobs_out.mkdir(exist_ok=True)

REGOBS_API = "https://api.regobs.no/app_v4"
regobs_endpoints = {
    "regobs_snow_profiles.json": f"{REGOBS_API}/Registration/Search",
    "regobs_api_info.json": "https://api.regobs.no/swagger/v1/swagger.json",
}
for fname, url in regobs_endpoints.items():
    download_file(url, regobs_out / fname, desc=f"Regobs {fname}", session=session)

# Search for avalanche observations (Norway, all years)
try:
    payload = {
        "GeoHazardTids": [10],  # 10 = Snow/avalanche
        "NumberOfRecords": 5000,
        "FromDate": "1970-01-01",
        "ToDate": str(date.today()),
    }
    r = session.post(
        f"{REGOBS_API}/Registration/Search",
        json=payload,
        timeout=120,
        headers={"Content-Type": "application/json"},
    )
    if r.status_code == 200:
        data = r.json()
        (regobs_out / "regobs_avalanche_obs_all.json").write_text(
            json.dumps(data, indent=2), encoding="utf-8"
        )
        count = len(data) if isinstance(data, list) else data.get("TotalMatches", "?")
        logger.info(f"  ✓  Regobs avalanche observations: {count} records")
    else:
        logger.warning(f"  Regobs search returned {r.status_code}")
except Exception as exc:
    logger.warning(f"  Regobs API: {exc}")


# --------------------------------------------------------------------------- #
# 4. NGI contact instructions                                                  #
# --------------------------------------------------------------------------- #
write_instructions(
    OUT_NGI,
    "Norway NGI Avalanche Registry — Access Instructions",
    """
Norwegian Geotechnical Institute (NGI) — www.ngi.no
Maintains the Norwegian avalanche registry (historical records).

DIRECT CONTACT:
  - General: info@ngi.no
  - Avalanche group: https://www.ngi.no/eng/Services/Avalanche
  - Dr. Karsten Müller or Dr. Regula Frauenfelder (avalanche researchers)

TEMPLATE EMAIL:
-----------------------------------------------------------------
Subject: Academic data request — Norwegian avalanche registry for solar forcing study

Dear NGI Avalanche Research Team,

I am conducting research on self-organized criticality in avalanche systems
and potential solar-atmospheric forcing of alpine snowpacks. I am investigating
whether geomagnetic storm activity and stratospheric sudden warming events
modulate avalanche frequency in high-latitude mountain ranges, through the
established atmospheric teleconnection chain (EPP → NOx → polar vortex → surface weather).

I would like to request access to the Norwegian avalanche registry with:
  - Event dates, location (region/fjord arm)
  - Size estimates (if available)
  - Trigger type (natural/artificial)
  Coverage: as far back as records exist

I am happy to sign a data agreement, acknowledge NGI, and offer co-authorship
for any publication using this dataset.

Kind regards,
[Your name and affiliation]
-----------------------------------------------------------------

PUBLIC ALTERNATIVES ALREADY DOWNLOADED:
  • Varsom bulletins (2012–present): data/cryosphere/ngi_norway/varsom/
  • Regobs observations: data/cryosphere/ngi_norway/regobs/
  • NVE snow data: data/cryosphere/ngi_norway/nve_public/
"""
)

logger.info("=== Script 14 complete ===")
