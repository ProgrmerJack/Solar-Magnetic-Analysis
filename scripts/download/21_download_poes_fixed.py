"""
NOAA POES / MEPED Energetic Particle Precipitation — CDAWeb archive
Downloads hemispheric power index + MEPED particle flux (NOAA-15 through -19).
CDAWeb is the correct access path; SWPC/NGDC direct URLs are stale.
"""
import sys, json, time, re
from pathlib import Path
from datetime import date
import requests
from bs4 import BeautifulSoup
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger, download_file

logger = get_logger("poes_meped")
OUT = DATA_DIR / "atmospheric" / "poes_meped"
OUT.mkdir(parents=True, exist_ok=True)

session = requests.Session()
session.headers["User-Agent"] = "Solar-Magnetic-Analysis/1.0"

CDAWEB = "https://cdaweb.gsfc.nasa.gov/pub/data"

# ── 1. NOAA POES Hemispheric Power Index (HPI) ───────────────────────────────
# NOAA SWPC publishes HPI as Ovation Prime model output
# Current JSON endpoint confirmed working:
logger.info("=== NOAA SWPC Ovation Aurora (HPI proxy) ===")
hpi_endpoints = {
    "ovation_aurora_latest": "https://services.swpc.noaa.gov/json/ovation_aurora_latest.json",
    "ovation_north_24h":     "https://services.swpc.noaa.gov/json/ovation_aurora_latest.json",
    "hemispheric_power":     "https://services.swpc.noaa.gov/json/hemispheric_power.json",
    "hemi_power_index":      "https://services.swpc.noaa.gov/json/geospace/geospace_hemi_power_index_1_hour.json",
    "hemi_power_index_7d":   "https://services.swpc.noaa.gov/json/geospace/geospace_hemi_power_index_7_day.json",
}
hpi_out = OUT / "swpc_hpi"
hpi_out.mkdir(exist_ok=True)
for name, url in hpi_endpoints.items():
    try:
        r = session.get(url, timeout=30)
        if r.status_code == 200:
            (hpi_out / f"{name}.json").write_bytes(r.content)
            logger.info(f"  ✓  {name}: {len(r.content)//1024} KB")
        else:
            logger.warning(f"  {r.status_code} {name}")
    except Exception as e:
        logger.warning(f"  FAIL {name}: {e}")

# ── 2. CDAWeb POES MEPED particle data ───────────────────────────────────────
logger.info("\n=== CDAWeb NOAA POES MEPED SEM-2 ===")
# Correct CDAWeb paths for POES:
poes_cdaweb_paths = {
    "noaa15": "noaa/noaa15/sem2_fluxes-2sec",
    "noaa16": "noaa/noaa16/sem2_fluxes-2sec",
    "noaa17": "noaa/noaa17/sem2_fluxes-2sec",
    "noaa18": "noaa/noaa18/sem2_fluxes-2sec",
    "noaa19": "noaa/noaa19/sem2_fluxes-2sec",
    "metop01": "noaa/metop01/sem2_fluxes-2sec",
    "metop02": "noaa/metop02/sem2_fluxes-2sec",
}
working = {}
for sat, path in poes_cdaweb_paths.items():
    url = f"{CDAWEB}/{path}/"
    try:
        r = session.get(url, timeout=15)
        if r.status_code == 200:
            soup = BeautifulSoup(r.text, "html.parser")
            years = sorted([a["href"].strip("/") for a in soup.find_all("a", href=True)
                           if re.match(r"^\d{4}/?$", a["href"])])
            logger.info(f"  ✓  {sat}: years {years[:5]}...{years[-3:] if len(years)>5 else ''}")
            working[sat] = (path, years)
        else:
            logger.info(f"  {r.status_code} {sat}")
    except Exception as e:
        logger.warning(f"  FAIL {sat}: {e}")
    time.sleep(0.3)

# ── 3. Download NOAA-15/18/19 recent years (2010–present) ────────────────────
logger.info("\n=== Downloading POES MEPED CDF files (2010–present) ===")
PRIORITY_SATS = ["noaa15", "noaa18", "noaa19"]
START_YEAR = 2010

for sat in PRIORITY_SATS:
    if sat not in working:
        logger.warning(f"  {sat} not found in CDAWeb, skipping")
        continue
    path, years = working[sat]
    sat_out = OUT / sat
    sat_out.mkdir(exist_ok=True)
    recent_years = [y for y in years if int(y) >= START_YEAR]
    logger.info(f"  {sat}: downloading {len(recent_years)} years")

    for yr in recent_years:
        yr_url = f"{CDAWEB}/{path}/{yr}/"
        yr_out = sat_out / yr
        done_flag = yr_out / "_complete.flag"
        if done_flag.exists():
            continue
        try:
            yr_r = session.get(yr_url, timeout=20)
            if yr_r.status_code != 200:
                continue
            yr_soup = BeautifulSoup(yr_r.text, "html.parser")
            cdf_files = [a["href"] for a in yr_soup.find_all("a", href=True)
                        if a["href"].endswith(".cdf") or a["href"].endswith(".CDF")]
            yr_out.mkdir(exist_ok=True)
            logger.info(f"    {sat}/{yr}: {len(cdf_files)} CDF files")
            success = True
            for cdf in cdf_files:
                dest = yr_out / cdf
                if not dest.exists():
                    ok = download_file(f"{CDAWEB}/{path}/{yr}/{cdf}", dest,
                                       desc=f"{sat} {yr}/{cdf[:20]}")
                    if not ok:
                        success = False
                    time.sleep(0.5)
            if success and cdf_files:
                done_flag.touch()
        except Exception as e:
            logger.warning(f"  {sat}/{yr}: {e}")

# ── 4. NOAA POES Hemispheric Power — NCEI archive ────────────────────────────
logger.info("\n=== NCEI POES Hemispheric Power Index (historical) ===")
# NCEI POES data portal
ncei_tests = [
    "https://www.ngdc.noaa.gov/stp/satellite/poes/dataaccess.html",
    "https://satdat.ngdc.noaa.gov/sem/poes/data/processed/ngdc/",
]
for url in ncei_tests:
    try:
        r = session.get(url, timeout=15)
        logger.info(f"  {r.status_code} [{len(r.content)//1024}KB] {url}")
    except Exception as e:
        logger.warning(f"  FAIL {url}: {e}")

# ── 5. Write instructions for full POES archive ───────────────────────────────
instructions = """
POES/MEPED Data Access — NOAA Energetic Particle Precipitation
==============================================================

1. CDAWeb (primary — auto-downloaded above):
   https://cdaweb.gsfc.nasa.gov/pub/data/noaa/noaa15/sem2_fluxes-2sec/
   NOAA-15, -16, -17, -18, -19 + MetOp-01, MetOp-02
   Format: CDF files, ~daily files per satellite

2. NCEI POES SEM-2 archive:
   https://satdat.ngdc.noaa.gov/sem/poes/data/processed/ngdc/
   30-second averages, txt format

3. Auroraweb hemispheric power index:
   https://www.swpc.noaa.gov/products/auroral-electrojet
   Or use OMNI combined dataset: HPI included in some products

4. For EPP-atmosphere coupling analysis (key for this research):
   Use Aura/MLS NOx data (downloaded via 11_setup_nasa_earthdata.py)
   together with POES MEPED to validate EPP flux → NOx pathway

5. NCEI Web Services for POES:
   https://www.ncei.noaa.gov/access/search/data-search/poes-l1b

Note: For the core SOC analysis, the OMNI HPI (Hemispheric Power Index) 
embedded in the OMNI2 dataset is sufficient for storm classification.
"""
(OUT / "INSTRUCTIONS.txt").write_text(instructions, encoding="utf-8")
logger.info("  Instructions written")
logger.info("=== POES/MEPED download complete ===")
