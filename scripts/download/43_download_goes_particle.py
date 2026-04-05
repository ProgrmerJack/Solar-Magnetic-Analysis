"""
GOES Energetic Particle Data via CDAWeb HAPI (open access, no auth needed).

Downloads GOES MAGED (Magnetospheric Electron Detector) and EPEAD electron flux
as proxy for energetic particle precipitation (EPP) into the polar atmosphere.

Datasets:
  GOES13_EPS-MAGED_1MIN  : 2006-06-01 – 2017-12-14  (>0.1–4 MeV electrons)
  GOES14_EPS-MAGED_1MIN  : 2009-07-01 – 2020-03-04
  GOES15_EPS-MAGED_1MIN  : 2010-03-01 – 2020-03-04
  GOES15_EPEAD-SCIENCE-ELECTRONS-E13EW_1MIN : 2010-03-26 – 2020-03-04  (>0.8, >2 MeV)
  GOES14_EPEAD-SCIENCE-ELECTRONS-E13EW_1MIN : 2012-10-01 – 2020-03-04

Output: data/solar/goes_particle/<dataset>/<YYYY>/<dataset>_<YYYY>-<MM>.csv.gz
"""
import requests
import gzip
import io
import time
import calendar
from datetime import date, timedelta
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[1]))
from utils import get_logger

LOG = get_logger("43_goes_particle")

HAPI = "https://cdaweb.gsfc.nasa.gov/hapi/data"
OUT  = Path(__file__).parents[2] / "data" / "solar" / "goes_particle"
OUT.mkdir(parents=True, exist_ok=True)

# Dataset specs: (id, start, end, key_params)
DATASETS = [
    (
        # Total MAGED electron flux: stack9 (90° PA), stack1 (10° PA), stack5 (50° PA), total
        # HAPI requires params in info order: stack9 < stack1 < stack5 < dtc_cor_eflux
        "GOES15_EPS-MAGED_1MIN",
        date(2010, 3, 1), date(2020, 3, 4),
        "dtc_cor_eflux_stack9,dtc_cor_eflux_stack1,dtc_cor_eflux_stack5,dtc_cor_eflux",
    ),
    (
        "GOES14_EPS-MAGED_1MIN",
        date(2009, 7, 1), date(2020, 3, 4),
        "dtc_cor_eflux_stack9,dtc_cor_eflux_stack1,dtc_cor_eflux_stack5,dtc_cor_eflux",
    ),
    (
        "GOES13_EPS-MAGED_1MIN",
        date(2006, 6, 1), date(2017, 12, 14),
        "dtc_cor_eflux_stack9,dtc_cor_eflux_stack1,dtc_cor_eflux_stack5,dtc_cor_eflux",
    ),
    (
        "GOES15_EPEAD-SCIENCE-ELECTRONS-E13EW_1MIN",
        date(2010, 3, 26), date(2020, 3, 4),
        "E1W_COR_FLUX,E2W_COR_FLUX",        # >0.8 MeV, >2 MeV integral electrons
    ),
    (
        "GOES14_EPEAD-SCIENCE-ELECTRONS-E13EW_1MIN",
        date(2012, 10, 1), date(2020, 3, 4),
        "E1W_COR_FLUX,E2W_COR_FLUX",
    ),
]


def months_in_range(start: date, end: date):
    """Yield (year, month) tuples from start to end inclusive."""
    y, m = start.year, start.month
    while date(y, m, 1) <= end:
        yield y, m
        m += 1
        if m > 12:
            m = 1
            y += 1


def download_month(dataset_id: str, params: str, year: int, month: int,
                   out_dir: Path, retries: int = 3) -> bool:
    last_day = calendar.monthrange(year, month)[1]
    t_min = f"{year}-{month:02d}-01T00:00:00Z"
    t_max = f"{year}-{month:02d}-{last_day:02d}T23:59:59Z"
    fname = out_dir / f"{dataset_id}_{year}-{month:02d}.csv.gz"
    if fname.exists() and fname.stat().st_size > 500:
        return True

    for attempt in range(retries):
        try:
            r = requests.get(
                HAPI,
                params={
                    "id": dataset_id,
                    "parameters": params,
                    "time.min": t_min,
                    "time.max": t_max,
                    "format": "csv",
                },
                stream=True,
                timeout=120,
            )
            if r.status_code == 200:
                out_dir.mkdir(parents=True, exist_ok=True)
                raw = r.content
                with gzip.open(fname, "wb") as f:
                    f.write(raw)
                if fname.stat().st_size < 200:
                    fname.unlink(missing_ok=True)
                    return False  # empty response
                return True
            elif r.status_code == 404:
                return False   # no data for this period
            else:
                LOG.warning("  %d for %s %d-%02d (attempt %d)", r.status_code,
                            dataset_id, year, month, attempt + 1)
                time.sleep(5 * (attempt + 1))
        except Exception as exc:
            LOG.warning("  ERR %s %d-%02d: %s", dataset_id, year, month, exc)
            time.sleep(5 * (attempt + 1))
    return False


def download_dataset(dataset_id: str, params: str, start: date, end: date):
    LOG.info("=== %s  (%s – %s) ===", dataset_id, start, end)
    ds_dir = OUT / dataset_id
    ds_dir.mkdir(parents=True, exist_ok=True)

    total = ok = 0
    for y, m in months_in_range(start, end):
        yr_dir = ds_dir / str(y)
        yr_dir.mkdir(parents=True, exist_ok=True)
        total += 1
        if download_month(dataset_id, params, y, m, yr_dir):
            ok += 1
            if ok % 12 == 0:
                LOG.info("  %d months done (year %d)", ok, y)
        time.sleep(0.5)          # be polite to HAPI server

    LOG.info("  ✓ %s: %d/%d months downloaded → %s", dataset_id, ok, total, ds_dir)


# Run all datasets
for ds_id, t_start, t_end, params in DATASETS:
    download_dataset(ds_id, params, t_start, t_end)

LOG.info("=== GOES particle download complete ===")
