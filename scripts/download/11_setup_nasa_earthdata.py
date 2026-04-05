"""
Script 11 — NASA Earthdata Setup
Sets up credentials and download scripts for NASA Earthdata-hosted datasets:
  • MERRA-2 reanalysis (NASA GMAO)
  • Aura/MLS stratospheric composition (NOx, O3, temperature)
  • TIMED/SABER mesosphere/stratosphere temperature
  • MODIS Terra/Aqua MOD10A1 daily snow cover (500 m)
  • NSIDC Arctic/Antarctic snow and ice

REGISTRATION: https://urs.earthdata.nasa.gov/users/new (free)
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger, write_instructions

logger = get_logger("11_nasa_earthdata")

# --------------------------------------------------------------------------- #
# .netrc setup for NASA Earthdata (required for wget/curl/Python downloads)   #
# --------------------------------------------------------------------------- #
netrc_path = Path.home() / ".netrc"
if not netrc_path.exists():
    netrc_path.write_text(
        "machine urs.earthdata.nasa.gov\n"
        "    login YOUR_USERNAME\n"
        "    password YOUR_PASSWORD\n",
        encoding="utf-8"
    )
    logger.info("  Created template ~/.netrc — FILL IN YOUR NASA Earthdata credentials")
else:
    if "urs.earthdata.nasa.gov" not in netrc_path.read_text():
        with open(netrc_path, "a") as f:
            f.write("\nmachine urs.earthdata.nasa.gov\n    login YOUR_USERNAME\n    password YOUR_PASSWORD\n")
        logger.info("  Appended NASA Earthdata entry to ~/.netrc")
    else:
        logger.info("  ~/.netrc already has urs.earthdata.nasa.gov entry")


# --------------------------------------------------------------------------- #
# MERRA-2 download script                                                      #
# --------------------------------------------------------------------------- #
MERRA2_SCRIPT = '''\
"""
MERRA-2 Download Script
Downloads MERRA-2 reanalysis data via NASA GES DISC.
Requires NASA Earthdata credentials in ~/.netrc
"""
import requests
from pathlib import Path
from datetime import date, timedelta

OUT = Path(__file__).parents[2] / "data" / "atmospheric" / "merra2"
OUT.mkdir(parents=True, exist_ok=True)

session = requests.Session()
# Earthdata auth is handled by .netrc automatically

GES_DISC = "https://goldsmr4.gesdisc.eosdis.nasa.gov/data"

# MERRA-2 Monthly Mean datasets (smaller than daily, good for trend analysis)
DATASETS = {
    # Stratospheric temperature and wind (key for SSW analysis)
    "M2TMNXSLV": {
        "path": "MERRA2_MONTHLY/M2TMNXSLV.5.12.4",
        "desc": "Single level surface/atmosphere monthly mean",
    },
    "M2IMNPASM": {
        "path": "MERRA2_MONTHLY/M2IMNPASM.5.12.4",
        "desc": "Upper atmosphere monthly mean (strat/meso)",
    },
    # Daily data (for time series analysis)
    "M2T1NXSLV": {
        "path": "MERRA2/M2T1NXSLV.5.12.4",
        "desc": "Single level 1-hourly time-averaged",
    },
}

print("MERRA-2 structure ready. Visit GES DISC to browse:")
print("  https://disc.gsfc.nasa.gov/datasets?keywords=MERRA-2")
print("  Or use: pip install earthaccess && python -c \\"import earthaccess; earthaccess.login()\\"")
'''

merra2_script = DATA_DIR / "atmospheric" / "merra2" / "download_merra2.py"
merra2_script.parent.mkdir(parents=True, exist_ok=True)
if not merra2_script.exists():
    merra2_script.write_text(MERRA2_SCRIPT, encoding="utf-8")


# --------------------------------------------------------------------------- #
# Aura/MLS download script                                                     #
# --------------------------------------------------------------------------- #
MLS_SCRIPT = '''\
"""
Aura/MLS Download Script
Aura Microwave Limb Sounder — stratospheric NOx, O3, temperature
Coverage: 2004–present

MLS datasets used for EPP-NOx analysis:
  MLS/Aura_L2GP-NO2_v04: NO2 profiles (proxy for NOx)
  MLS/Aura_L2GP-N2O_v04: N2O (tracer)
  MLS/Aura_L2GP-Temperature_v04: Temperature profiles
  MLS/Aura_L2GP-O3_v04: Ozone profiles

Data access: https://disc.gsfc.nasa.gov/datasets?keywords=MLS
Or: https://mls.jpl.nasa.gov/data/
"""
import earthaccess  # pip install earthaccess
from pathlib import Path
from datetime import datetime

OUT = Path(__file__).parents[2] / "data" / "atmospheric" / "aura_mls"
OUT.mkdir(parents=True, exist_ok=True)

earthaccess.login()

results = earthaccess.search_data(
    short_name="ML2NO2",        # MLS NO2 Level 2
    version="004",
    temporal=("2004-08-13", "2026-04-03"),
)
print(f"Found {len(results)} MLS NO2 granules")
# earthaccess.download(results, str(OUT / "no2"))

results_t = earthaccess.search_data(
    short_name="ML2T",          # MLS Temperature Level 2
    version="004",
    temporal=("2004-08-13", "2026-04-03"),
)
print(f"Found {len(results_t)} MLS Temperature granules")
# earthaccess.download(results_t, str(OUT / "temperature"))
'''

mls_script = DATA_DIR / "atmospheric" / "aura_mls" / "download_aura_mls.py"
mls_script.parent.mkdir(parents=True, exist_ok=True)
if not mls_script.exists():
    mls_script.write_text(MLS_SCRIPT, encoding="utf-8")


# --------------------------------------------------------------------------- #
# MODIS Snow Cover download script                                             #
# --------------------------------------------------------------------------- #
MODIS_SCRIPT = '''\
"""
MODIS MOD10A1 Daily Snow Cover — 500 m resolution
Terra + Aqua (MOD10A1.061 / MYD10A1.061)
Coverage: 2000-present

Key tiles for analysis:
  h18v04 = Alps (Switzerland)
  h09v04, h10v04 = Colorado Rockies
  h16v02 = Scandinavia (Norway)
"""
import earthaccess
from pathlib import Path

OUT = Path(__file__).parents[2] / "data" / "cryosphere" / "modis_snow"
OUT.mkdir(parents=True, exist_ok=True)

earthaccess.login()

# Swiss Alps tile (h18v04)
results = earthaccess.search_data(
    short_name="MOD10A1",
    version="061",
    temporal=("2000-02-24", "2026-04-03"),
    granule_name="*h18v04*",
)
print(f"Found {len(results)} MODIS Terra snow granules for Swiss Alps tile")

# Colorado tile
results_co = earthaccess.search_data(
    short_name="MOD10A1",
    version="061",
    temporal=("2000-02-24", "2026-04-03"),
    granule_name="*h09v04*",
)
print(f"Found {len(results_co)} MODIS Terra snow granules for Colorado tile")

# To download (comment out for inspection only):
# earthaccess.download(results, str(OUT / "alps"))
# earthaccess.download(results_co, str(OUT / "colorado"))
'''

modis_script = DATA_DIR / "cryosphere" / "modis_snow" / "download_modis_snow.py"
modis_script.parent.mkdir(parents=True, exist_ok=True)
if not modis_script.exists():
    modis_script.write_text(MODIS_SCRIPT, encoding="utf-8")


# --------------------------------------------------------------------------- #
# TIMED/SABER download                                                         #
# --------------------------------------------------------------------------- #
SABER_SCRIPT = '''\
"""
TIMED/SABER Download
TIMED Sounding of the Atmosphere using Broadband Emission Radiometry
Coverage: 2002–present
Key data: Mesospheric/stratospheric temperature profiles (15-120 km)

Data access: https://saber.gats-inc.com/data.php
Direct FTP: ftp://saber.gats-inc.com/
"""
import ftplib
from pathlib import Path

OUT = Path(__file__).parents[2] / "data" / "atmospheric" / "timed_saber"
OUT.mkdir(parents=True, exist_ok=True)

FTP_HOST = "saber.gats-inc.com"
FTP_USER = "anonymous"
FTP_PASS = "your@email.com"  # replace with your email

try:
    with ftplib.FTP(FTP_HOST) as ftp:
        ftp.login(FTP_USER, FTP_PASS)
        ftp.cwd("/data/l2b/v2.0/")
        dirs = ftp.nlst()
        print(f"Available directories: {dirs[:10]}")
except Exception as e:
    print(f"FTP access: {e}")
    print("Alternative: https://saber.gats-inc.com/data.php")
    print("Or use NASA Earthdata: short_name=SABER_L2B_V2.0")
'''

saber_script = DATA_DIR / "atmospheric" / "timed_saber" / "download_saber.py"
saber_script.parent.mkdir(parents=True, exist_ok=True)
if not saber_script.exists():
    saber_script.write_text(SABER_SCRIPT, encoding="utf-8")


# --------------------------------------------------------------------------- #
# Write instructions                                                           #
# --------------------------------------------------------------------------- #
write_instructions(
    DATA_DIR / "atmospheric" / "aura_mls",
    "Aura/MLS Stratospheric Composition Data",
    """
REGISTRATION: https://urs.earthdata.nasa.gov/users/new (free)

STEP 1: pip install earthaccess
STEP 2: earthaccess.login()  (will prompt for credentials once, then save)
STEP 3: python download_aura_mls.py

KEY DATASETS:
  ML2NO2 v004  — NO2 profiles 2004–present (proxy for NOx from EPP)
  ML2T   v004  — Temperature profiles 2004–present (stratospheric temp)
  ML2O3  v004  — Ozone profiles 2004–present

ALTERNATIVE: NASA GES DISC web portal
  https://disc.gsfc.nasa.gov/datasets?keywords=MLS+Aura

ALTERNATIVE: MLS data portal at JPL
  https://mls.jpl.nasa.gov/data/
"""
)

write_instructions(
    DATA_DIR / "atmospheric" / "merra2",
    "MERRA-2 Reanalysis — NASA GMAO",
    """
REGISTRATION: https://urs.earthdata.nasa.gov/users/new (free)

STEP 1: pip install earthaccess
STEP 2: python download_merra2.py

KEY DATASETS (GES DISC):
  M2I3NPASM — Instantaneous 3-hourly assimilated upper atmosphere
              (1980–present, global, 0.5°x0.625°, 72 model levels)
  M2T1NXSLV — Single-level hourly time-averaged (surface fields)
  M2TMNXSLV — Monthly mean single-level

NASA GMAO direct access: https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/
GES DISC browser: https://disc.gsfc.nasa.gov/datasets?keywords=MERRA-2
"""
)

write_instructions(
    DATA_DIR / "cryosphere" / "modis_snow",
    "MODIS MOD10A1 Daily Snow Cover Product",
    """
REGISTRATION: https://urs.earthdata.nasa.gov/users/new (free)

STEP 1: pip install earthaccess
STEP 2: python download_modis_snow.py

KEY TILES:
  h18v04 — Swiss Alps / Central Europe
  h09v04 — Colorado Rockies (western)
  h10v04 — Colorado Rockies (eastern)
  h16v02 — Norway / Scandinavia
  h17v02 — Norway / Sweden

PRODUCTS:
  MOD10A1.061 — Terra daily snow cover, 500m, 2000–present
  MYD10A1.061 — Aqua daily snow cover, 500m, 2002–present

NSIDC viewer: https://nsidc.org/data/MOD10A1/versions/61
"""
)

write_instructions(
    DATA_DIR / "cryosphere" / "nsidc",
    "NSIDC Arctic & Antarctic Snow/Ice Datasets",
    """
NSIDC hosts multiple datasets relevant to this research:

1. NSIDC-0051: Sea Ice Concentrations (1979–present)
   https://nsidc.org/data/nsidc-0051
   
2. NSIDC-0192: MODIS/Terra Snow Cover Daily L3 Global 500m
   (same as MOD10A1, see modis_snow/)
   
3. GLiMS Glacier Database: https://www.glims.org/
4. IceBridge Airborne Snow: https://nsidc.org/data/icebridge

ACCESS: https://urs.earthdata.nasa.gov/users/new (free NASA Earthdata account)
Or: https://nsidc.org/data/user-resources/help-center/how-access-nsidc-data

pip install earthaccess
import earthaccess; earthaccess.login()
"""
)

logger.info("=== Script 11 complete (setup only — credentials required) ===")
logger.info("  Next steps: register at https://urs.earthdata.nasa.gov/ and run download scripts")
