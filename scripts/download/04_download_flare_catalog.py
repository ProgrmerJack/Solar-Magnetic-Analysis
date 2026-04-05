"""
Script 04 — NOAA SWPC Solar Flare / Geophysical Event Catalog
Downloads the complete NOAA/SWPC Solar & Geophysical Event Reports (1975–present).

SWPC publishes yearly event text files covering:
  FL = X-ray flares (GOES class A/B/C/M/X)
  RB = Radio Bursts, XRA = X-ray events, RSP = Radio Spectral events
  EPL = Eruptive Prominences, FLA = H-alpha flares, etc.

Sources tried in order:
  1. NOAA SWPC FTP mirror over HTTPS  (services.swpc.noaa.gov)
  2. NOAA NCEI GOES events archive
  3. HEK (Heliophysics Event Knowledgebase) JSON API at LMSAL
"""
import sys
import json
import re
from pathlib import Path
from datetime import date
import requests
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger, download_file

logger = get_logger("04_flare_catalog")
OUT = DATA_DIR / "solar" / "flare_catalog"
OUT.mkdir(parents=True, exist_ok=True)

session = requests.Session()
session.headers.update({"User-Agent": "Solar-Magnetic-Analysis/1.0"})


# --------------------------------------------------------------------------- #
# 1. NOAA SWPC yearly event report text files (1996–present)                  #
# --------------------------------------------------------------------------- #
def download_swpc_events():
    """Download SWPC yearly event report files."""
    logger.info("=== NOAA SWPC Yearly Event Reports ===")
    SWPC_EVENTS = "https://www.swpc.noaa.gov/ftpdir/indices/events"
    current_year = date.today().year
    for year in range(1996, current_year + 1):
        fname = f"{year}events.txt"
        url = f"{SWPC_EVENTS}/{fname}"
        dest = OUT / "swpc_yearly" / fname
        download_file(url, dest, desc=f"SWPC events {year}", session=session)


# --------------------------------------------------------------------------- #
# 2. NGDC/NCEI GOES X-ray event lists (1975–present)                          #
# --------------------------------------------------------------------------- #
def download_ngdc_xray_events():
    """Download NGDC GOES solar X-ray event catalog."""
    logger.info("=== NGDC GOES Solar X-ray Event Lists ===")
    NGDC_BASE = "https://www.ngdc.noaa.gov/stp/satellite/goes"
    current_year = date.today().year
    for year in range(1975, current_year + 1):
        url = f"{NGDC_BASE}/xray/GOES_XRS_EventLists/{year}_Xray_Events.txt"
        dest = OUT / "ngdc_xray" / f"{year}_Xray_Events.txt"
        download_file(url, dest, desc=f"NGDC X-ray events {year}", session=session)


# --------------------------------------------------------------------------- #
# 3. HEK (Heliophysics Event Knowledgebase) flare catalog via LMSAL API       #
# --------------------------------------------------------------------------- #
def download_hek_flares():
    """Download flare catalog from HEK API in chunks by year."""
    logger.info("=== HEK Flare Catalog (LMSAL) ===")
    HEK = "https://www.lmsal.com/hek/her"
    out_hek = OUT / "hek"
    out_hek.mkdir(exist_ok=True)

    current_year = date.today().year
    for year in range(2002, current_year + 1):  # RHESSI/GOES coverage from 2002
        out_file = out_hek / f"hek_flares_{year}.json"
        if out_file.exists():
            logger.info(f"  SKIP  hek_flares_{year}.json  (already exists)")
            continue

        params = {
            "cosec":        "2",
            "cmd":          "search",
            "type":         "column",
            "event_type":   "FL",
            "event_starttime": f"{year}-01-01T00:00:00",
            "event_endtime":   f"{year}-12-31T23:59:59",
            "result_limit": "100000",
            "page":         "1",
            "return":       ("kb_archivdate,ar_noaanum,frm_name,obs_observatory,"
                             "obs_instrument,fl_goescls,event_peaktime,"
                             "event_starttime,event_endtime,event_coord1,"
                             "event_coord2,fl_peakflux,hgs_x,hgs_y"),
        }
        try:
            r = session.get(HEK, params=params, timeout=120)
            r.raise_for_status()
            data = r.json()
            events = data.get("result", [])
            out_file.write_text(json.dumps(events, indent=2), encoding="utf-8")
            logger.info(f"  ✓  hek_flares_{year}.json  ({len(events)} events)")
        except Exception as exc:
            logger.warning(f"  HEK {year}: {exc}")


# --------------------------------------------------------------------------- #
# 4. NOAA SWPC JSON recent flare summary                                      #
# --------------------------------------------------------------------------- #
def download_swpc_recent():
    logger.info("=== SWPC Recent Flare Summary (JSON) ===")
    recent = {
        "flares_24h.json":  "https://services.swpc.noaa.gov/json/goes/primary/xray-flares-24-hours.json",
        "flares_latest.json": "https://services.swpc.noaa.gov/json/goes/primary/xray-flares-latest.json",
    }
    for fname, url in recent.items():
        download_file(url, OUT / fname, desc=f"SWPC {fname}", session=session)


# --------------------------------------------------------------------------- #
# 5. GOES event lists via NOAA NCEI Thredds catalog (alternate path)          #
# --------------------------------------------------------------------------- #
def download_goes_catalog_summary():
    logger.info("=== NOAA GOES XRS 1-day summary JSONs ===")
    SWPC = "https://services.swpc.noaa.gov"
    for sat in ("primary", "secondary"):
        for window in ("1-day", "3-day", "7-day"):
            fname = f"goes_{sat}_xrays_{window}.json"
            url = f"{SWPC}/json/goes/{sat}/xrays-{window}.json"
            download_file(url, OUT / fname, desc=fname, session=session)


if __name__ == "__main__":
    download_swpc_events()
    download_ngdc_xray_events()
    download_hek_flares()
    download_swpc_recent()
    download_goes_catalog_summary()
    logger.info("=== Script 04 complete ===")
