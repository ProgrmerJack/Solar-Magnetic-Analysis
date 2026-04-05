"""
NOAA POES MEPED Historical — CDAWeb HAPI (1978-2004)
Downloads older NOAA satellites (NOAA-05 through NOAA-14) 1-minute MEPED
electron and proton flux via CDAWeb HAPI.

Together with 43_download_goes_particle.py (2006-2020) and
21_download_poes_fixed.py (2013-2025 pub server), this gives full EPP coverage:
  1978-2004  NOAA-05/06/07/08/10/12/14  (1-min MEPED via HAPI)
  2006-2020  GOES-13/14/15 MAGED        (1-min electron flux via HAPI)
  2013-2025  NOAA-15/18/19              (2-sec MEPED via CDAWeb pub)

Key parameters:
  Jperp_e  — perpendicular electron flux (3 energy bands) [primary EPP proxy]
  Jomni_p  — omnidirectional proton flux (5 energy bands)
  Jomni_e  — omnidirectional electron flux (3 energy bands)
"""
import sys, time, json, csv, gzip
from pathlib import Path
from datetime import date, timedelta
import requests

sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger

logger = get_logger("44_poes_hapi")
HAPI = "https://cdaweb.gsfc.nasa.gov/hapi"
OUT  = DATA_DIR / "atmospheric" / "poes_meped" / "hapi_historical"
OUT.mkdir(parents=True, exist_ok=True)

sess = requests.Session()
sess.headers["User-Agent"] = "Solar-Magnetic-Analysis/1.0"

# Satellites with 1-min MEPED data, confirmed via HAPI catalog
SATELLITES = [
    ("NOAA05_MEPED1MIN_SEM",  date(1978, 11,  2), date(1979,  8, 11)),
    ("NOAA06_MEPED1MIN_SEM",  date(1979,  6, 28), date(1986, 11, 18)),
    ("NOAA07_MEPED1MIN_SEM",  date(1981,  7, 11), date(1985,  2, 11)),
    ("NOAA08_MEPED1MIN_SEM",  date(1983,  5,  9), date(1985, 10, 14)),
    ("NOAA10_MEPED1MIN_SEM",  date(1986, 10, 11), date(1991,  8, 31)),
    ("NOAA12_MEPED1MIN_SEM",  date(1991,  6,  1), date(1999,  9, 30)),
    ("NOAA14_MEPED1MIN_SEM",  date(1995,  1, 11), date(2004, 12, 31)),
]

# Parameters to request (must match HAPI info order exactly)
PARAMS = "Jperp_e,Jomni_p,Jomni_e"


def fetch_month(dataset: str, year: int, month: int) -> bytes | None:
    """Fetch one month of MEPED data from HAPI, return gzip-compressed CSV bytes."""
    t0 = f"{year}-{month:02d}-01T00:00:00Z"
    # End = first day of next month
    if month == 12:
        t1 = f"{year+1}-01-01T00:00:00Z"
    else:
        t1 = f"{year}-{month+1:02d}-01T00:00:00Z"

    url = (f"{HAPI}/data?id={dataset}&parameters={PARAMS}"
           f"&time.min={t0}&time.max={t1}&format=csv")
    for attempt in range(3):
        try:
            r = sess.get(url, timeout=120, stream=True)
            if r.status_code != 200:
                logger.warning(f"    HTTP {r.status_code}: {dataset} {year}-{month:02d}")
                return None
            content = r.content
            if len(content) < 200:
                # HAPI error embedded in tiny response
                try:
                    err = json.loads(content)
                    logger.warning(f"    HAPI error {dataset} {year}-{month:02d}: {err}")
                except Exception:
                    pass
                return None
            return gzip.compress(content)
        except Exception as e:
            if attempt < 2:
                time.sleep(5)
            else:
                logger.warning(f"    Exception {dataset} {year}-{month:02d}: {e}")
    return None


for dataset, start, stop in SATELLITES:
    sat_out = OUT / dataset
    sat_out.mkdir(exist_ok=True)

    # Iterate month-by-month
    current = start.replace(day=1)
    end_m   = stop.replace(day=1)

    total_months = 0
    skipped = 0

    logger.info(f"=== {dataset}: {start} → {stop} ===")

    while current <= end_m:
        yr, mo = current.year, current.month
        dest = sat_out / f"{yr}" / f"{yr}{mo:02d}.csv.gz"
        dest.parent.mkdir(exist_ok=True)

        if dest.exists() and dest.stat().st_size > 200:
            skipped += 1
        else:
            data = fetch_month(dataset, yr, mo)
            if data:
                dest.write_bytes(data)
                total_months += 1
            time.sleep(0.3)

        # Advance to next month
        if mo == 12:
            current = current.replace(year=yr+1, month=1)
        else:
            current = current.replace(month=mo+1)

    logger.info(f"  Done: {total_months} downloaded, {skipped} skipped → {sat_out}")

logger.info("=== POES MEPED HAPI historical download complete ===")
