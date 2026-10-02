#!/usr/bin/env python3
"""
acquire_davis2022.py
====================
CESM2(WACCM6) forecasts of the January 2021 SSW with and without the SSW in the
initial stratosphere (Davis et al. 2022, Nat. Commun. 13, 1136), reduced for
07_physical_decomposition/davis2021_intervention_test.py (registered 738a422,
before this script first ran).

SOURCE  Zenodo record 5639805 (CC BY 4.0), anonymous: per member, Z3 (46 daily
        steps x 22 pressure levels, 192 x 288 grid; one compressed chunk, so each
        ~116 MB file is fetched whole) and tas_2m (46 x 181 x 360, 1 degree).
ARMS    ssfcst_04jan, 9to12kmRamp (04jan) | ssfcst_01feb, 9to12kmRamp_updated (01feb)
        | ssfcst_08feb, 9to12kmRamp_updated (08feb); members m00..m20.
REDUCTION
  zcap60   cos-latitude-weighted mean of 1000 hPa Z3 over the rows >= 60N (first row
           60.785N on this grid), all longitudes (m)
  regions  2 m temperature cos-lat means over acquire_s2s_reforecasts.T2M_REGIONS (K)
  Time stamps are used as given (daily, 00 UTC, the first equal to the initial
  time). Values are daily means (cell_methods 'time: mean'); the files do not say
  whether a stamp opens or closes its averaging day, a one-day ambiguity that an
  18-day window mean barely feels.
Raw files are kept under raw/davis2022 (read-only once written).
Output: davis2022_reduced.parquet  arm, init, member, time, zcap60, NEURASIA, ...
"""
import socket
import threading
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from acquire_s2s_reforecasts import regions_mean     # noqa: E402

RAW = HERE / "raw" / "davis2022"
OUT = HERE / "davis2022_reduced.parquet"
URL = "https://zenodo.org/api/records/5639805/files/{}/content"
ARMS = {  # arm label: (file prefix, date tag, init)
    "ssfcst_04jan": ("ssfcst_04jan", "04jan2021", "2021-01-04"),
    "9to12kmRamp_04jan": ("9to12kmRamp", "04jan2021", "2021-01-04"),
    "ssfcst_01feb": ("ssfcst_01feb", "01feb2021", "2021-02-01"),
    "9to12kmRamp_updated_01feb": ("9to12kmRamp_updated", "01feb2021", "2021-02-01"),
    "ssfcst_08feb": ("ssfcst_08feb", "08feb2021", "2021-02-08"),
    "9to12kmRamp_updated_08feb": ("9to12kmRamp_updated", "08feb2021", "2021-02-08"),
}
WORKERS = 6


def names(prefix, tag, m):
    return (f"{prefix}_Z3_70Lwaccm6_{tag}00z_d01_d46_m{m:02d}.nc",
            f"{prefix}_tas_2m_70Lwaccm6_{tag}_00z_d01_d46_m{m:02d}.nc")


socket.setdefaulttimeout(120)        # a stalled Zenodo connection hung for 5 h without one
# HDF5 is not thread-safe: concurrent opens from the download threads segfaulted the
# process (exit 139). Downloads stay parallel; every file open and reduction is serial.
H5 = threading.Lock()


def fetch(key):
    f = RAW / key
    if f.exists() and f.stat().st_size > 0:
        return f
    for attempt in range(8):
        try:
            tmp = f.with_suffix(".part")
            urllib.request.urlretrieve(URL.format(key), tmp)
            with H5:
                xr.open_dataset(tmp).close()                 # must open
            tmp.rename(f)
            return f
        except Exception as e:                               # noqa: BLE001 -- retried, then raised
            print(f"  {key} attempt {attempt + 1}: {type(e).__name__}: {str(e)[:120]}", flush=True)
            time.sleep(20 * (attempt + 1))
    raise RuntimeError(key)


def reduce(arm, m):
    prefix, tag, init = ARMS[arm]
    kz, kt = names(prefix, tag, m)
    fz, ft = fetch(kz), fetch(kt)
    with H5:
        return _reduce(arm, m, init, fz, ft)


def _reduce(arm, m, init, fz, ft):
    z = xr.open_dataset(fz)["Z3"].sel(lev_p=1000.0)
    z = z.sel(lat=z.lat[z.lat >= 60.0 - 1e-6])
    # the 192-row grid has no 60N row: the cap starts at 60.785N (row spacing 0.942)
    assert 60.0 <= float(z.lat.min()) < 60.0 + 0.95, float(z.lat.min())
    cap = z.weighted(np.cos(np.deg2rad(z.lat))).mean(("lat", "lon")).values
    t = xr.open_dataset(ft)["tas_2m"].rename({"lat": "latitude", "lon": "longitude"})
    reg = regions_mean(t).to_dataframe()
    times = [pd.Timestamp(str(x)[:10]) for x in z.time.values]
    df = pd.DataFrame({"arm": arm, "init": init, "member": f"m{m:02d}", "time": times, "zcap60": cap})
    for r in ("NEURASIA", "HI_EUROPE", "MID_EASIA", "MID_NAMER"):
        df[r] = reg[r].values
    assert len(df) == 46 and np.isfinite(df.zcap60).all(), (arm, m)
    print(f"  {arm} m{m:02d}: zcap60 mean {cap.mean():.1f} m", flush=True)
    return df


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    arms = sys.argv[sys.argv.index("--arms") + 1].split(",") if "--arms" in sys.argv else list(ARMS)
    jobs = [(a, m) for a in arms for m in range(21)]
    with ThreadPoolExecutor(WORKERS) as ex:
        parts = list(ex.map(lambda j: reduce(*j), jobs))
    d = pd.concat(parts, ignore_index=True)
    if OUT.exists():
        old = pd.read_parquet(OUT)
        d = pd.concat([old[~old.arm.isin(d.arm.unique())], d], ignore_index=True)
    d.to_parquet(OUT)
    print(d.groupby("arm").member.nunique().to_string())
    print(f"Saved -> {OUT.name} ({len(d):,} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
