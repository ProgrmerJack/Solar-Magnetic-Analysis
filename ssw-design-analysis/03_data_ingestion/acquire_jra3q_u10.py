#!/usr/bin/env python3
"""
acquire_jra3q_u10.py
====================
Input of 07_physical_decomposition/early_transfer_retest.py (registered 2026-10-09, committed
before any JRA-3Q value was read): daily-mean zonal-mean zonal wind at 10 hPa on the 60.0 N
row of JRA-3Q (Kosaka et al. 2024), November-April of the winters 1947/48 to 1979/80.
Source: NSF NCAR GDEX d640000 via OSDF (anonymous),
  https://osdf-director.osg-htc.org/ncar/gdex/d640000/anl_p125/YYYYMM/
  jra3q.anl_p125.0_2_2.ugrd-pres-an-ll125.YYYYMM0100_YYYYMMDD18.nc  (6-hourly, 1.25 deg)
The files are chunked one time step at a time, so each monthly file (~370 MB) is downloaded
to _jra3q_tmp/, reduced (zonal mean of the 60.0 N row at 10 hPa, mean of the four analysis
times of each UTC day; a day needs all four) and deleted. Per-month results are cached in
_jra3q_u10/ so an interrupted run resumes. Writes jra3q_u10_60N_daily.parquet (date, u10_60N).
"""
import re
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
TMP = HERE / "_jra3q_tmp"; CACHE = HERE / "_jra3q_u10"
OUT = HERE / "jra3q_u10_60N_daily.parquet"
BASE = "https://osdf-director.osg-htc.org/ncar/gdex/d640000/anl_p125/"
MONTHS = [(y, m) for w in range(1947, 1980) for (y, m) in [(w, 11), (w, 12), (w + 1, 1), (w + 1, 2), (w + 1, 3), (w + 1, 4)]]


def listing(ym):
    with urllib.request.urlopen(BASE + ym + "/", timeout=120) as r:
        t = r.read().decode()
    names = sorted(set(re.findall(r"jra3q\.anl_p125\.0_2_2\.ugrd-pres-an-ll125\.\d{10}_\d{10}\.nc", t)))
    if len(names) != 1:
        raise RuntimeError(f"{ym}: {len(names)} ugrd files")
    return names[0]


def one(ymt):
    y, m = ymt; ym = f"{y}{m:02d}"; c = CACHE / f"{ym}.csv"
    if c.exists():
        return ym, "cached"
    name = listing(ym); dest = TMP / name
    for a in range(5):
        try:
            with urllib.request.urlopen(BASE + ym + "/" + name, timeout=600) as r, open(dest, "wb") as f:
                while chunk := r.read(1 << 22):
                    f.write(chunk)
            with h5py.File(dest, "r") as h:
                lev = h["pressure_level"][:]; lat = h["lat"][:]
                k = int(np.flatnonzero(np.isclose(lev, 10.0))[0]); j = int(np.flatnonzero(np.isclose(lat, 60.0))[0])
                u = h["ugrd-pres-an-ll125"][:, k, j, :].astype(float).mean(axis=1)
                t = pd.to_datetime([str(x) for x in h["utc_date_int"][:]], format="%Y%m%d%H")
            break
        except Exception as e:                      # noqa: BLE001
            if dest.exists():
                dest.unlink()
            if a == 4:
                raise
            time.sleep(30 * (a + 1))
    dest.unlink()
    s = pd.Series(u, index=t)
    g = s.groupby(s.index.floor("D"))
    d = g.mean()[g.count() == 4]
    pd.DataFrame({"date": d.index.strftime("%Y-%m-%d"), "u10_60N": d.values, "n_times": 4}).to_csv(c, index=False, lineterminator="\n")
    return ym, f"{len(d)} days"


def main():
    TMP.mkdir(exist_ok=True); CACHE.mkdir(exist_ok=True)
    with ThreadPoolExecutor(3) as ex:
        for ym, msg in ex.map(one, MONTHS):
            print(ym, msg, flush=True)
    d = pd.concat([pd.read_csv(CACHE / f"{y}{m:02d}.csv") for y, m in MONTHS]).sort_values("date")
    d.to_parquet(OUT, index=False)
    print(f"{len(d)} days {d.date.min()} .. {d.date.max()} -> {OUT.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
