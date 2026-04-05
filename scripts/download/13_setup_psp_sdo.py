"""
Script 13 — Parker Solar Probe & SDO/HMI+AIA Setup
Sets up access to the flagship solar physics datasets.

Parker Solar Probe (PSP):
  - Launch: August 2018
  - Instruments: FIELDS (magnetic field), SWEAP (solar wind), IS⊙IS (energetic particles)
  - Data: NASA SPDF (no strict authentication)
  - Best for: In-situ SOC analysis of solar wind magnetic fluctuations

SDO/HMI + AIA:
  - Launch: February 2010
  - HMI: Helioseismic and Magnetic Imager — photospheric magnetic flux maps
  - AIA: Atmospheric Imaging Assembly — EUV coronal imaging
  - Data: JSOC (Joint Science Operations Center, Stanford)
  - REGISTRATION: jsoc.stanford.edu (free)
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger, download_file, write_instructions

logger = get_logger("13_psp_sdo")
session = __import__("requests").Session()
session.headers.update({"User-Agent": "Solar-Magnetic-Analysis/1.0"})


# --------------------------------------------------------------------------- #
# 1. Parker Solar Probe — SPDF archive (no auth for most data)                #
# --------------------------------------------------------------------------- #
logger.info("=== Parker Solar Probe SPDF Data ===")
PSP_OUT = DATA_DIR / "solar" / "parker_solar_probe"
PSP_OUT.mkdir(parents=True, exist_ok=True)

SPDF_PSP = "https://spdf.gsfc.nasa.gov/pub/data/psp"

# PSP FIELDS instrument — magnetic field
# Level 2 fluxgate (1 min averages) for SOC analysis
PSP_FIELDS_URL = f"{SPDF_PSP}/fields/l2/mag_rtn_1min/"

# Try to get directory listing and download data
import re
try:
    r = session.get(PSP_FIELDS_URL, timeout=30)
    r.raise_for_status()
    years = re.findall(r'href="(\d{4})/"', r.text)
    logger.info(f"  PSP FIELDS available years: {years}")
    for yr in years:
        yr_url = f"{PSP_FIELDS_URL}{yr}/"
        yr_out = PSP_OUT / "fields_mag_rtn_1min" / yr
        yr_out.mkdir(parents=True, exist_ok=True)
        try:
            r2 = session.get(yr_url, timeout=30)
            r2.raise_for_status()
            files = re.findall(r'href="(psp_fld[^"]+\.cdf)"', r2.text)
            logger.info(f"  PSP FIELDS {yr}: {len(files)} files")
            for fname in files:
                download_file(f"{yr_url}{fname}", yr_out / fname,
                              desc=f"PSP FIELDS {yr} {fname}", session=session)
        except Exception as exc:
            logger.debug(f"  PSP {yr}: {exc}")
except Exception as exc:
    logger.warning(f"  PSP FIELDS listing failed: {exc}")


# PSP SWEAP — solar wind electrons, alphas, protons
PSP_SWEAP_URL = f"{SPDF_PSP}/sweap/spi/l3/spi_sf00_pad_l3/"
try:
    r = session.get(PSP_SWEAP_URL, timeout=30)
    r.raise_for_status()
    years = re.findall(r'href="(\d{4})/"', r.text)
    logger.info(f"  PSP SWEAP available years: {years}")
except Exception as exc:
    logger.debug(f"  PSP SWEAP: {exc}")


# --------------------------------------------------------------------------- #
# 2. PSP data via CDAWeb (more reliable SPDF mirror)                          #
# --------------------------------------------------------------------------- #
logger.info("=== PSP CDAWeb Access ===")
psp_cdaweb = {
    "psp_fld_l2_mag_rtn_1min": "https://cdaweb.gsfc.nasa.gov/pub/data/psp/fields/l2/mag_rtn_1min/",
    "psp_spi_sf00_l3":         "https://cdaweb.gsfc.nasa.gov/pub/data/psp/sweap/spi/l3/spi_sf00_pad_l3/",
}
for dataset, url in psp_cdaweb.items():
    try:
        r = session.get(url, timeout=20)
        r.raise_for_status()
        years = re.findall(r'href="(\d{4})/"', r.text)
        logger.info(f"  {dataset}: years available = {years}")
    except Exception as exc:
        logger.debug(f"  {dataset}: {exc}")


# --------------------------------------------------------------------------- #
# 3. SDO/HMI JSOC setup                                                       #
# --------------------------------------------------------------------------- #
write_instructions(
    DATA_DIR / "solar" / "sdo_hmi_aia",
    "SDO HMI + AIA Data — JSOC Stanford",
    """
SDO data is managed by JSOC (Joint Science Operations Center) at Stanford.

REGISTRATION: http://jsoc.stanford.edu/ajax/register_email.php (free)

ACCESSING HMI MAGNETIC FLUX MAPS:

  Method 1: JSOC Export Request (web UI)
    http://jsoc.stanford.edu/ajax/exportdata.html
    Series: hmi.M_720s   — Magnetogram, 720-second cadence
            hmi.B_720s   — Vector magnetic field
            hmi.sharp_720s — Solar Active Region Patches (SHARPs)
            aia.lev1_euv_12s — AIA EUV 12-second cadence

  Method 2: drms Python package (recommended)
    pip install drms
    
    import drms
    c = drms.Client(email='YOUR_EMAIL@domain.com')
    
    # SHARP patches (active regions) — for flare energy analysis
    ds = 'hmi.sharp_720s[2010.05.01_TAI/10y][!NOAA_ARS != "0"]{USFLUX,TOTUSJH,ABSNJZH}'
    keys = c.query(ds, key=['NOAA_ARS','USFLUX','TOTUSJH','T_REC'])
    
    # Full disk magnetogram for a day:
    ds = 'hmi.M_720s[2012.03.07_06:00:00_TAI/1d@6h]'
    r = c.export(ds, protocol='fits')
    r.wait()
    r.download('data/solar/sdo_hmi_aia/magnetograms/')

  Method 3: SunPy (high-level interface)
    pip install sunpy
    from sunpy.net import Fido, attrs as a
    result = Fido.search(
        a.Time('2012-01-01', '2012-12-31'),
        a.Instrument.hmi,
        a.Physobs.los_magnetic_field,
        a.Sample(720 * u.s)
    )
    files = Fido.fetch(result)

KEY DATASETS FOR THIS PROJECT:
  hmi.sharp_720s      — SHARPs: active region parameters (flux, helicity, etc.)
                         Use for flare energy proxy in power-law analysis
  hmi.M_720s          — Full-disk LOS magnetogram (720s cadence, 2010–present)
  aia.lev1_euv_12s    — AIA 12-second EUV images (multiple wavelengths)

JSOC CUTOUT SERVICE (for specific active regions):
  http://www.lmsal.com/hek/hcr?cmd=view-voevent&;ivorn=...
"""
)


# --------------------------------------------------------------------------- #
# 4. HMI SHARP parameter summary (key flare prediction metrics)               #
# --------------------------------------------------------------------------- #
logger.info("=== Attempting JSOC SHARP catalog (no auth needed for metadata) ===")
sharp_out = DATA_DIR / "solar" / "sdo_hmi_aia" / "sharp_catalog"
sharp_out.mkdir(parents=True, exist_ok=True)

# JSOC export of SHARP summary statistics (no download, just metadata)
jsoc_urls = {
    "hmi_sharp_info.html": "http://jsoc.stanford.edu/doc/data/hmi/sharp/sharp.htm",
    "jsoc_series_list.html": "http://jsoc.stanford.edu/ajax/lookdata.html",
}
for fname, url in jsoc_urls.items():
    download_file(url, sharp_out / fname, desc=f"JSOC {fname}", session=session)


# --------------------------------------------------------------------------- #
# 5. SunPy-based download script                                               #
# --------------------------------------------------------------------------- #
SUNPY_SCRIPT = '''\
"""
SunPy-based solar data download
Downloads GOES XRS, HEK flare catalog, and HMI data via SunPy/JSOC.
pip install sunpy drms
"""
import sunpy.net.attrs as a
from sunpy.net import Fido
import astropy.units as u
from pathlib import Path
from datetime import date

OUT = Path(__file__).parents[2] / "data" / "solar"

# ── GOES X-ray flux (SunPy handles authentication automatically) ──────────
result = Fido.search(
    a.Time("2017-01-01", str(date.today())),
    a.Instrument.xrs,
    a.goes.SatelliteNumber(16),
)
print(f"GOES-16 XRS results: {result}")
# files = Fido.fetch(result, path=str(OUT / "goes_xrs" / "{file}"))

# ── HEK Flare catalog ─────────────────────────────────────────────────────
result_fl = Fido.search(
    a.Time("2002-01-01", str(date.today())),
    a.hek.EventType("FL"),
    a.hek.OBS.Observatory == "GOES",
)
print(f"HEK flare events: {len(result_fl['hek'])}")
# files = Fido.fetch(result_fl, path=str(OUT / "flare_catalog" / "sunpy" / "{file}"))

# ── HMI SHARP (requires JSOC registration) ────────────────────────────────
import drms
c = drms.Client(email="YOUR_EMAIL@domain.com")
keys, segs = c.query(
    "hmi.sharp_720s[2010.05.01_TAI/15y@1d]{USFLUX,TOTBSQ,TOTUSJH}",
    key=["NOAA_ARS", "USFLUX", "TOTBSQ", "TOTUSJH", "T_REC"],
    seg=None,
)
print(f"SHARP records: {len(keys)}")
keys.to_csv(str(OUT / "sdo_hmi_aia" / "sharp_catalog" / "hmi_sharp_daily_2010_2025.csv"))
'''

script_path = DATA_DIR / "solar" / "sdo_hmi_aia" / "download_sdo_sunpy.py"
if not script_path.exists():
    script_path.write_text(SUNPY_SCRIPT, encoding="utf-8")
    logger.info(f"  ✓  SunPy download script written → {script_path}")


logger.info("=== Script 13 complete (setup + PSP direct downloads attempted) ===")
