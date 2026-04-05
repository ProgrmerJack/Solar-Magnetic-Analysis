"""
GOES-R Series (GOES-16/17/18) SEISS MPSH Energetic Particle Data
(Magnetospheric Particle Sensor - High Energy), L2 1-minute science data.

Provides electron and proton flux continuity with GOES-13/14/15 coverage (which ends 2020-03-04).

Data source:
  https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes{N}/l2/data/mpsh-l2-avg1m_science/
  File pattern: sci_mpsh-l2-avg1m_g{N}_d{YYYYMMDD}_v2-0-3.nc
  File size: ~1 MB/day (1440 1-minute records × channels × variables)

Coverage:
  GOES-16 (g16): 2017-present (primary operational satellite)
  GOES-17 (g17): 2018-2023 (decommissioned early 2023)
  GOES-18 (g18): 2022-present (replaced GOES-17)

Size estimate: ~6 GB total (GOES-16: 3.3 GB, GOES-17: 1.5 GB, GOES-18: 1.1 GB)

Output: data/solar/goes_r_particle/goes{N}/mpsh_avg1m/{YYYY}/{MM}/
"""
import requests
import time
from pathlib import Path
from datetime import date, timedelta
import sys

sys.path.insert(0, str(Path(__file__).parents[1]))
from utils import get_logger

LOG = get_logger("48_goes_r_particle")

BASE_NCEI = "https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes"
OUT = Path(__file__).parents[2] / "data" / "solar" / "goes_r_particle"

SATELLITES = [
    ("goes16", "g16", date(2017, 2, 7),  date(2026, 1, 1)),
    ("goes17", "g17", date(2018, 4, 1),  date(2023, 1, 1)),
    ("goes18", "g18", date(2022, 8, 1),  date(2026, 1, 1)),
]

SESSION = requests.Session()
SESSION.headers["User-Agent"] = "SolarMagnetic-Research/1.0"


def download_day(sat_dir: str, tag: str, day: date, out_dir: Path) -> bool:
    yyyymmdd = day.strftime("%Y%m%d")
    fname = f"sci_mpsh-l2-avg1m_{tag}_d{yyyymmdd}_v2-0-3.nc"
    out_path = out_dir / fname
    if out_path.exists() and out_path.stat().st_size > 10_000:
        return True  # already downloaded

    url = f"{BASE_NCEI}/{sat_dir}/l2/data/mpsh-l2-avg1m_science/{day.year}/{day.month:02d}/{fname}"
    for attempt in range(3):
        try:
            r = SESSION.get(url, timeout=(10, 30), stream=True)  # connect=10s, read=30s
            if r.status_code == 404:
                return False  # data not available for this day
            r.raise_for_status()
            out_path.parent.mkdir(parents=True, exist_ok=True)
            with open(out_path, "wb") as fh:
                for chunk in r.iter_content(65536):
                    fh.write(chunk)
            return True
        except Exception as exc:
            if attempt < 2:
                LOG.warning("  Attempt %d/3 failed %s: %s", attempt + 1, fname, exc)
                time.sleep(5)
            else:
                LOG.error("  FAILED %s: %s", fname, exc)
                return False
    return False


def main():
    for sat_dir, tag, start, end in SATELLITES:
        sat_out = OUT / sat_dir / "mpsh_avg1m"
        LOG.info("=== %s MPSH 1-min (%s → %s) ===", sat_dir.upper(), start, end)

        total = ok = skipped = 0
        current = start
        while current < end:
            year_dir = sat_out / str(current.year) / f"{current.month:02d}"
            result = download_day(sat_dir, tag, current, year_dir)
            total += 1
            if result:
                ok += 1
            else:
                skipped += 1
            if total % 365 == 0:
                LOG.info("  %s: %d days done (%d ok, %d skipped/missing)", sat_dir, total, ok, skipped)
            current += timedelta(days=1)

        LOG.info("  DONE %s: %d ok / %d total -> %s", sat_dir, ok, total, sat_out)

    LOG.info("=== GOES-R MPSH download complete ===")


if __name__ == "__main__":
    main()
