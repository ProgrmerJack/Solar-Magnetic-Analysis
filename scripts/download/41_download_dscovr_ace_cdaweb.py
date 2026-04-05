"""
DSCOVR/ACE Historical Solar Wind Data via CDAWeb HAPI.
No authentication required — publicly accessible.

Products (via https://cdaweb.gsfc.nasa.gov/hapi/):
  - DSCOVR_H0_MAG: B-field GSE/RTN, 1-min avg (2016-present)
  - DSCOVR_H1_FC: Solar wind plasma speed/density/temp (2016-present)
  - AC_H0_MFI: ACE B-field GSE, 16-sec (1998-present)
  - AC_H0_SWE: ACE solar wind plasma (1998-present)
  - AC_H1_EPM: ACE energetic particles (1998-present)

Year-by-year downloads with resume support.
"""
import requests, time
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parents[1]))
from utils import get_logger

LOG = get_logger("41_dscovr_ace")
OUT = Path(__file__).parents[2] / "data" / "solar" / "dscovr_ace"
OUT.mkdir(parents=True, exist_ok=True)

HAPI = "https://cdaweb.gsfc.nasa.gov/hapi/data"
SESS = requests.Session()
SESS.headers["User-Agent"] = "Solar-Magnetic-Research/1.0"

def hapi_year(dataset, params, year, out_path):
    if out_path.exists() and out_path.stat().st_size > 500:
        LOG.info("    skip: %s", out_path.name)
        return True
    url = (f"{HAPI}?id={dataset}"
           f"&time.min={year}-01-01T00:00:00"
           f"&time.max={year}-12-31T23:59:59"
           f"&parameters={','.join(params)}"
           f"&format=csv")
    LOG.info("    %d …", year)
    for attempt in range(3):
        try:
            r = SESS.get(url, timeout=120, stream=True)
            if r.status_code == 200:
                data = r.content
                if len(data) > 200:
                    out_path.write_bytes(data)
                    LOG.info("      ✓ %d bytes", len(data))
                    return True
                else:
                    LOG.warning("      empty response (%d bytes)", len(data))
                    return False
            elif r.status_code == 404:
                LOG.warning("      404 (no data for %d)", year)
                return False
            else:
                LOG.warning("      HTTP %d attempt %d", r.status_code, attempt+1)
        except Exception as e:
            LOG.warning("      error attempt %d: %s", attempt+1, e)
        time.sleep(2)
    return False

def download_years(dataset, params, start, end, subdir):
    d = OUT / subdir
    d.mkdir(exist_ok=True)
    ok = 0
    for yr in range(start, end + 1):
        f = d / f"{dataset.lower()}_{yr}.csv"
        if hapi_year(dataset, params, yr, f):
            ok += 1
    LOG.info("  %s: %d/%d years", dataset, ok, end - start + 1)

# ─── 1. DSCOVR Magnetometer (2016-2025) ──────────────────────────────────────
LOG.info("=== DSCOVR Magnetometer (2016-2025) ===")
download_years("DSCOVR_H0_MAG", ["B1F1", "B1GSE"], 2016, 2025, "dscovr_mag")

# ─── 2. DSCOVR Solar Wind Plasma (2016-2025) ─────────────────────────────────
LOG.info("=== DSCOVR Plasma FC (2016-2025) ===")
download_years("DSCOVR_H1_FC", ["V_GSE", "Np", "THERMAL_TEMP"], 2016, 2025, "dscovr_plasma")

# ─── 3. ACE Magnetometer L2 (1998-2025) ──────────────────────────────────────
LOG.info("=== ACE MAG H0 (1998-2025) ===")
download_years("AC_H0_MFI", ["Magnitude", "BGSEc"], 1998, 2025, "ace_mag")

# ─── 4. ACE Solar Wind Plasma L2 (1998-2025) ─────────────────────────────────
LOG.info("=== ACE SWEPAM H0 (1998-2025) ===")
download_years("AC_H0_SWE", ["Np", "Vp", "Tpr"], 1998, 2025, "ace_swepam")

# ─── 5. ACE Energetic Particles H1 (1998-2025) ───────────────────────────────
LOG.info("=== ACE EPAM H1 (1998-2025) ===")
download_years("AC_H1_EPM", ["P1p", "P3p", "P5p", "P7p"], 1998, 2025, "ace_epam")

LOG.info("=== DSCOVR/ACE download complete ===")
LOG.info("Output: %s", OUT)
