"""
GOES XRS Science-Quality Data — NCEI (correct paths)
Downloads annual 1-min avg NetCDF + flare summary NetCDF for GOES-8 through GOES-15.
Also checks for GOES-R (16/17/18) data under alternate NCEI paths.
"""
import sys, re, time
from pathlib import Path
import requests
from bs4 import BeautifulSoup
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger, download_file

logger = get_logger("goes_xrs_v2")
NCEI_XRS = "https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/science/xrs"
OUT = DATA_DIR / "solar" / "goes_xrs"
OUT.mkdir(parents=True, exist_ok=True)

session = requests.Session()
session.headers["User-Agent"] = "Solar-Magnetic-Analysis/1.0"

def ls(url):
    """Return list of hrefs in an NCEI directory listing."""
    try:
        r = session.get(url, timeout=20)
        if r.status_code != 200:
            return []
        soup = BeautifulSoup(r.text, "html.parser")
        return [a["href"] for a in soup.find_all("a", href=True)
                if not a["href"].startswith(("?", "http", "mailto", "/css", "/js", "/img", "/ngdc",
                                              "/instruments"))]
    except Exception as e:
        logger.warning(f"  ls({url}): {e}")
        return []

# ── GOES-8 through GOES-15 ────────────────────────────────────────────────────
SATELLITES = ["goes08", "goes09", "goes10", "goes11", "goes12", "goes13", "goes14", "goes15"]
PRODUCTS   = ["xrsf-l2-avg1m_science", "xrsf-l2-flsum_science"]   # 1-min avg + flare summary

logger.info("=== GOES XRS NCEI science archive (GOES-8 to GOES-15) ===")
for sat in SATELLITES:
    sat_out = OUT / sat
    sat_out.mkdir(exist_ok=True)
    products = ls(f"{NCEI_XRS}/{sat}/")
    if not products:
        logger.info(f"  {sat}: not found or empty")
        continue
    logger.info(f"  {sat}: products → {products}")

    for prod in PRODUCTS:
        prod_clean = prod.rstrip("/")
        if not any(prod_clean in p for p in products):
            continue
        prod_url = f"{NCEI_XRS}/{sat}/{prod_clean}/"
        items = ls(prod_url)
        # Prefer merged/full-mission file if present; otherwise yearly files
        nc_files = [f for f in items if f.endswith(".nc")]
        if not nc_files:
            logger.info(f"    {sat}/{prod_clean}: no NetCDF at top level, checking years …")
            years = sorted([d.strip("/") for d in items if re.match(r"^\d{4}/?$", d)])
            for yr in years:
                yr_ncs = [f for f in ls(f"{prod_url}{yr}/") if f.endswith(".nc")]
                nc_files += [f"{yr}/{f}" for f in yr_ncs]

        prod_out = sat_out / prod_clean
        prod_out.mkdir(exist_ok=True)
        for nc in nc_files:
            dest = prod_out / Path(nc).name
            if dest.exists() and dest.stat().st_size > 10_000:
                logger.info(f"    ✓  {sat}/{prod_clean}/{Path(nc).name} (exists)")
                continue
            url = f"{prod_url}{nc}"
            download_file(url, dest, desc=f"{sat} {prod_clean[-12:]} {Path(nc).name[:30]}")
            time.sleep(0.5)

# ── GOES-R series (16/17/18) — check NCEI alternate path ────────────────────
logger.info("\n=== GOES-R (16/17/18) XRS — searching NCEI ===")
NCEI_ALT = "https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access"
for sat_id, label in [("g16","goes16"), ("g17","goes17"), ("g18","goes18")]:
    for path_try in [f"{NCEI_ALT}/{sat_id}/", f"{NCEI_ALT}/{label}/",
                     f"{NCEI_ALT}/science/xrs/{label}/"]:
        try:
            r = session.get(path_try, timeout=10)
            if r.status_code == 200:
                items = ls(path_try)
                xrs_items = [i for i in items if "xrs" in i.lower()]
                logger.info(f"  ✓  {path_try} → {items[:8]}")
                if xrs_items:
                    logger.info(f"       XRS items: {xrs_items}")
                break
        except Exception:
            pass
    else:
        logger.info(f"  {label}: not found in NCEI (use DONKI for GOES-R era flares)")
    time.sleep(0.3)

logger.info("=== GOES XRS download complete ===")
