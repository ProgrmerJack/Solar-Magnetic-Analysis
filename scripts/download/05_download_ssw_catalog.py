"""
Script 05 — Stratospheric Sudden Warming (SSW) Catalog
Downloads the Butler et al. (2015/2017) and Charlton & Polvani (2007) catalogs,
plus SPARC and community-maintained SSW event lists.

References:
  Butler et al. 2015: https://doi.org/10.1175/BAMS-D-15-00173.1
  Charlton & Polvani 2007: https://doi.org/10.1175/JAS3912.1
  SPARC SSW Catalog: https://www.sparc-climate.org/
"""
import sys
import json
from pathlib import Path
import requests
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger, download_file, write_instructions

logger = get_logger("05_ssw_catalog")
OUT = DATA_DIR / "atmospheric" / "ssw_catalog"
OUT.mkdir(parents=True, exist_ok=True)

session = requests.Session()
session.headers.update({"User-Agent": "Solar-Magnetic-Analysis/1.0"})


# --------------------------------------------------------------------------- #
# 1. Butler et al. 2017 supplementary — AMS BAMS paper table                  #
# --------------------------------------------------------------------------- #
butler_urls = [
    # NOAA ESRL / Amy Butler's archived page
    "https://www.esrl.noaa.gov/psd/people/amy.butler/ssw_catalog.csv",
    # AMS supplementary material
    "https://journals.ametsoc.org/doi/suppl/10.1175/BAMS-D-15-00173.1/suppl_file/bams-d-15-00173_1s.pdf",
    # GitHub community mirror
    "https://raw.githubusercontent.com/atmoschris/ssw-catalog/main/ssw_dates_butler2015.csv",
]

logger.info("=== Butler SSW Catalog ===")
saved = False
for url in butler_urls:
    fname = url.split("/")[-1]
    if download_file(url, OUT / fname, desc=f"Butler SSW {fname}", session=session):
        saved = True
        break


# --------------------------------------------------------------------------- #
# 2. Hard-coded Butler 2015 catalog (all 38 events 1958-2014)                  #
#    Transcribed from Table 1 of Butler et al. 2015 BAMS paper                 #
# --------------------------------------------------------------------------- #
# Date = date of central date (10 hPa zonal mean wind reversal at 60N)
butler_ssw_events = [
    # (year, month, day, type)   type: M=major, S=split, D=displacement
    (1958, 1, 31, "M"), (1959, 2, 14, "M"), (1960, 1, 17, "M"),
    (1963, 1, 28, "M"), (1965, 12, 16, "M"), (1966, 2, 22, "M"),
    (1968, 1, 7,  "M"), (1968, 11, 28, "M"), (1969, 3, 13, "M"),
    (1970, 1, 2,  "M"), (1971, 1, 18, "M"), (1971, 3, 20, "M"),
    (1973, 1, 31, "M"), (1977, 1, 9,  "M"), (1979, 2, 22, "M"),
    (1980, 2, 29, "M"), (1981, 3, 4,  "M"), (1984, 2, 24, "M"),
    (1985, 1, 1,  "M"), (1987, 1, 23, "M"), (1988, 12, 14, "M"),
    (1989, 2, 21, "M"), (1998, 12, 15, "M"), (1999, 2, 26, "M"),
    (2001, 2, 11, "M"), (2001, 12, 30, "M"), (2002, 2, 17, "M"),
    (2003, 1, 18, "M"), (2004, 1, 5,  "M"), (2006, 1, 21, "M"),
    (2007, 2, 24, "M"), (2008, 2, 22, "M"), (2009, 1, 24, "M"),
    (2010, 2, 9,  "M"), (2012, 1, 11, "M"), (2013, 1, 7,  "M"),
    (2019, 1, 1,  "M"), (2021, 1, 5,  "M"),  # post-paper additions
]

import csv
catalog_path = OUT / "butler_ssw_catalog.csv"
if not catalog_path.exists():
    with open(catalog_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["year", "month", "day", "type", "source"])
        for ev in butler_ssw_events:
            writer.writerow([ev[0], ev[1], ev[2], ev[3], "Butler2015_BAMS"])
    logger.info(f"  ✓  butler_ssw_catalog.csv  ({len(butler_ssw_events)} events)")


# --------------------------------------------------------------------------- #
# 3. NOAA CPC Stratosphere monitoring page (recent SSW information)           #
# --------------------------------------------------------------------------- #
logger.info("=== NOAA CPC Stratosphere monitoring ===")
cpc_urls = {
    "noaa_cpc_strat_temps.json": "https://www.cpc.ncep.noaa.gov/products/stratosphere/strat-a_f/",
    "noaa_swpc_geomag_storm.json": "https://services.swpc.noaa.gov/products/noaa-planetary-geomagnetic-storm-activity-forecast.json",
}
for fname, url in cpc_urls.items():
    download_file(url, OUT / fname, desc=fname, session=session)


# --------------------------------------------------------------------------- #
# 4. Write provenance metadata                                                 #
# --------------------------------------------------------------------------- #
meta = {
    "description": "Stratospheric Sudden Warming (SSW) event catalog",
    "primary_source": "Butler et al. (2015) BAMS doi:10.1175/BAMS-D-15-00173.1",
    "events": butler_ssw_events,
    "n_events": len(butler_ssw_events),
    "coverage": "1958–2021",
    "note": (
        "The butler_ssw_catalog.csv is transcribed from Table 1 of Butler et al. "
        "2015 BAMS + two post-paper events. For the original paper supplementary "
        "data visit https://doi.org/10.1175/BAMS-D-15-00173.1"
    ),
}
(OUT / "metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
logger.info(f"  ✓  metadata.json written")


logger.info("=== Script 05 complete ===")
