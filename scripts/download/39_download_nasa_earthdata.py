"""
Aura/MLS + MERRA-2 + MODIS Snow download using NASA Earthdata credentials.

Priority datasets:
  1. Aura/MLS L3 Daily Zonal HNO3 (ML3DZHNO3) -- EPP->NOx proxy, 2004-present
  2. Aura/MLS L3 Daily Zonal Temperature (ML3DZT) -- stratospheric temp
  3. Aura/MLS L3 Daily Zonal N2O (ML3DZN2O) -- dynamics tracer
  4. Aura/MLS L3 Monthly Binned T/O3/HNO3/N2O -- long-term trends
  5. MERRA-2 Monthly Mean Upper Atmosphere (M2IMNPASM) -- 1980-present
  6. MODIS MOD10A1 Daily Snow Cover -- Alps (h18v04), Colorado (h09v04), Norway (h16v02)

Auth: earthaccess with environment credentials (EARTHDATA_USERNAME / EARTHDATA_PASSWORD)
"""
import os
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parents[1]))
from utils import get_logger

LOG = get_logger("39_nasa_earthdata")

os.environ.setdefault("EARTHDATA_USERNAME", "andersonmark578")
os.environ.setdefault("EARTHDATA_PASSWORD", "Nasa.20080408@#")

import earthaccess

auth = earthaccess.login(strategy="environment", persist=False)
LOG.info("earthaccess auth: %s", auth.authenticated)

BASE = Path(__file__).parents[2] / "data"
MLS_OUT = BASE / "atmospheric" / "aura_mls"
MERRA_OUT = BASE / "atmospheric" / "merra2"
MODIS_OUT = BASE / "cryosphere" / "modis_snow"
for d in [MLS_OUT, MERRA_OUT, MODIS_OUT]:
    d.mkdir(parents=True, exist_ok=True)

TEMPORAL = ("2004-08-13", "2026-04-03")
WINTER_MONTHS = ["11", "12", "01", "02", "03", "04"]


def download_collection(short_name, out_dir, label, temporal=TEMPORAL,
                        version=None, batch=50, monthly_filter=None):
    """Search CMR and download all granules not already on disk.
    
    monthly_filter: if set, only download granules whose filename month is in this list
                    (e.g. ['10','11','12','01','02','03','04'] for winter months)
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    kw = dict(short_name=short_name, temporal=temporal, count=-1)
    if version:
        kw["version"] = version

    granules = earthaccess.search_data(**kw)
    LOG.info("  Found %d %s granules", len(granules), label)

    pending = []
    for g in granules:
        links = g.data_links(access="onprem")
        if not links:
            links = g.data_links()
        if not links:
            continue
        fname = Path(links[0]).name
        # Apply month filter for MODIS to limit to winter months only
        if monthly_filter:
            # MODIS filename: MOD10A1.AYYYYDDD.h18v04.006.*.hdf -> extract day-of-year
            # Filter by checking temporal metadata instead
            try:
                begin = g["umm"]["TemporalExtent"]["RangeDateTime"]["BeginningDateTime"]
                month = begin[5:7]  # "YYYY-MM-DD..." -> month
                if month not in monthly_filter:
                    continue
            except (KeyError, IndexError, TypeError):
                pass
        file_path = out_dir / fname
        if not file_path.exists():
            pending.append(g)

    LOG.info("  Pending: %d (skipping %d already downloaded)", len(pending), len(granules) - len(pending))
    ok = 0
    for i in range(0, len(pending), batch):
        chunk = pending[i:i + batch]
        try:
            files = earthaccess.download(chunk, local_path=str(out_dir), threads=4)
            ok += len([f for f in files if f])
            LOG.info("    %d/%d downloaded", ok, len(pending))
        except earthaccess.exceptions.EulaNotAccepted as e:
            LOG.error("  EULA not accepted for %s", short_name)
            LOG.error("  Approve app at https://urs.earthdata.nasa.gov/ then re-run")
            LOG.error("  %s", e)
            break
        except Exception as exc:
            LOG.warning("  batch %d/%d error: %s", i // batch + 1, (len(pending) + batch - 1) // batch, exc)
    LOG.info("  OK %s: %d/%d -> %s", label, ok, len(pending), out_dir)
    return ok


# 1. Aura/MLS L3 Daily Zonal HNO3 (ML3DZHNO3) -- KEY EPP-NOx proxy
#    HNO3 increases in polar upper stratosphere after energetic particle events.
#    Daily resolution enables storm-epoch superposed epoch analysis.
LOG.info("=== Aura/MLS L3 Daily Zonal HNO3 (EPP-NOx proxy) ===")
download_collection("ML3DZHNO3", MLS_OUT / "hno3_daily_zonal", "ML3DZHNO3 daily HNO3", version="005")

# 2. Aura/MLS L3 Daily Zonal Temperature (ML3DZT)
LOG.info("=== Aura/MLS L3 Daily Zonal Temperature ===")
download_collection("ML3DZT", MLS_OUT / "temperature_daily_zonal", "ML3DZT daily T", version="005")

# 3. Aura/MLS L3 Daily Zonal N2O (ML3DZN2O) -- atmospheric dynamics tracer
LOG.info("=== Aura/MLS L3 Daily Zonal N2O ===")
download_collection("ML3DZN2O", MLS_OUT / "n2o_daily_zonal", "ML3DZN2O daily N2O", version="005")

# 4. Aura/MLS L3 Monthly Binned products (long-term trend analysis)
LOG.info("=== Aura/MLS L3 Monthly Binned products ===")
for sn, label, subdir in [
    ("ML3MBT",    "monthly T",    "temperature_monthly"),
    ("ML3MBO3",   "monthly O3",   "ozone_monthly"),
    ("ML3MBHNO3", "monthly HNO3", "hno3_monthly"),
    ("ML3MBN2O",  "monthly N2O",  "n2o_monthly"),
]:
    download_collection(sn, MLS_OUT / subdir, sn, version="005")

# 5. MERRA-2 Monthly Mean Upper Atmosphere (M2IMNPASM) -- SKIPPED
#    Files are ~1.1 GB each (global 3D fields) × 554 months ≈ 600 GB total.
#    Instead we use:
#      - NCEP/NCAR stratosphere (already downloaded, 1948-present monthly)
#      - ERA5 polar stratosphere monthly means (script 47, polar subset, tiny)
#      - Aura/MLS daily for 2004-present (downloaded above)
LOG.info("=== MERRA-2: SKIPPED (too large; using NCEP + ERA5 polar strat instead) ===")

# 6. MODIS MOD10A1 Daily Snow Cover -- Swiss Alps tile, winter months only
#    Full 3-tile 26-year download = ~300 GB (too large).
#    Alps tile h18v04 (0-30E, 40-50N) + winter months Nov-Apr: ~22 yrs × 6 mo × ~30 days × 5 MB ≈ 20 GB
#    This validates IMS 4km snow cover at 500m resolution for avalanche-prone terrain.
LOG.info("=== MODIS MOD10A1 Daily Snow Cover (Alps winter only) ===")
# Use bounding_box to restrict to Swiss Alps tile h18v04 region
alps_kw = dict(
    short_name="MOD10A1",
    temporal=("2004-08-13", "2026-04-03"),
    count=-1,
    bounding_box=(5, 44, 11, 48),   # (W, S, E, N) Swiss Alps core → returns h18v04 only
    # Note: no version= filter here; earthaccess version param uses "61" not "061" and returns 0
)
granules_alps = earthaccess.search_data(**alps_kw)
LOG.info("  Found %d MOD10A1 Alps granules total", len(granules_alps))
alps_out = MODIS_OUT / "alps_h18v04"
alps_out.mkdir(parents=True, exist_ok=True)
pending_alps = []
for g in granules_alps:
    links = g.data_links(access="onprem")
    if not links:
        links = g.data_links()
    if not links:
        continue
    fname = Path(links[0]).name
    # Only keep Collection 6.1 granules
    if ".061." not in fname:
        continue
    # Only keep winter months (Nov-Apr) to avoid ~300 GB full archive
    try:
        begin = g["umm"]["TemporalExtent"]["RangeDateTime"]["BeginningDateTime"]
        month = begin[5:7]
        if month not in WINTER_MONTHS:
            continue
    except (KeyError, IndexError, TypeError):
        pass
    if not (alps_out / fname).exists():
        pending_alps.append(g)

LOG.info("  Pending (winter only): %d / %d", len(pending_alps), len(granules_alps))
ok_modis = 0
for i in range(0, len(pending_alps), 50):
    chunk = pending_alps[i:i+50]
    try:
        files = earthaccess.download(chunk, local_path=str(alps_out), threads=4)
        ok_modis += len([f for f in files if f])
        LOG.info("    %d/%d downloaded", ok_modis, len(pending_alps))
    except earthaccess.exceptions.EulaNotAccepted as e:
        LOG.error("  EULA not accepted for MOD10A1 - approve LP DAAC app at https://urs.earthdata.nasa.gov/")
        LOG.error("  %s", e)
        break
    except Exception as exc:
        LOG.warning("  MODIS batch error: %s", exc)
LOG.info("  OK MOD10A1 Alps winter: %d/%d -> %s", ok_modis, len(pending_alps), alps_out)

LOG.info("=== NASA Earthdata downloads complete ===")
