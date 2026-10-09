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
# The director redirects each request to a cache; on 2026-10-09 one (INFN) gave 6.8 MB/s and another
# (SURFsara) stalled at ~300 B/s. Attempts alternate between the INFN cache and the director, and a
# transfer slower than MIN_RATE after WATCH seconds is abandoned and retried.
FAST = "https://xsod14.cr.cnaf.infn.it:8443/ncar/gdex/d640000/anl_p125/"
MIN_RATE, WATCH = 500e3, 20.0


class Slow(Exception):
    pass
MONTHS = [(y, m) for w in range(1947, 1980) for (y, m) in [(w, 11), (w, 12), (w + 1, 1), (w + 1, 2), (w + 1, 3), (w + 1, 4)]]


def fname(y, m):
    """file names follow one pattern (first and last analysis time of the month), so no listing is needed"""
    last = pd.Timestamp(y, m, 1) + pd.offsets.MonthEnd(0)
    return f"jra3q.anl_p125.0_2_2.ugrd-pres-an-ll125.{y}{m:02d}0100_{y}{m:02d}{last.day:02d}18.nc"


def one(ymt):
    y, m = ymt; ym = f"{y}{m:02d}"; c = CACHE / f"{ym}.csv"
    if c.exists():
        return ym, "cached"
    name = fname(y, m); dest = TMP / name
    for a in range(8):
        try:
            src = (FAST if a % 2 == 0 else BASE) + ym + "/" + name
            t0, got = time.time(), 0
            with urllib.request.urlopen(src, timeout=120) as r, open(dest, "wb") as f:
                while chunk := r.read(1 << 20):
                    f.write(chunk); got += len(chunk)
                    el = time.time() - t0
                    if el > WATCH and got / el < MIN_RATE:
                        raise Slow(f"{got / el / 1e3:.0f} kB/s")
            with h5py.File(dest, "r") as h:
                lev = h["pressure_level"][:]; lat = h["lat"][:]
                k = int(np.flatnonzero(np.isclose(lev, 10.0))[0]); j = int(np.flatnonzero(np.isclose(lat, 60.0))[0])
                u = h["ugrd-pres-an-ll125"][:, k, j, :].astype(float).mean(axis=1)
                t = pd.to_datetime([str(x) for x in h["utc_date_int"][:]], format="%Y%m%d%H")
            break
        except Exception as e:                      # noqa: BLE001
            if dest.exists():
                dest.unlink()
            if a == 7:
                raise
            print(f"{ym}: retry {a + 1} after {type(e).__name__} {e}"[:120], flush=True)
            time.sleep(5 if isinstance(e, Slow) else min(600, 30 * 2 ** a))
    dest.unlink()
    s = pd.Series(u, index=t)
    g = s.groupby(s.index.floor("D"))
    d = g.mean()[g.count() == 4]
    pd.DataFrame({"date": d.index.strftime("%Y-%m-%d"), "u10_60N": d.values, "n_times": 4}).to_csv(c, index=False, lineterminator="\n")
    return ym, f"{len(d)} days"


def main():
    TMP.mkdir(exist_ok=True); CACHE.mkdir(exist_ok=True)
    def safe(ymt):
        try:
            return one(ymt)
        except Exception as e:                      # noqa: BLE001  a network outage must not end the run
            return f"{ymt[0]}{ymt[1]:02d}", f"FAILED this pass: {type(e).__name__}"
    for npass in range(6):                          # repeat passes until every month is cached
        todo = [ym for ym in MONTHS if not (CACHE / f"{ym[0]}{ym[1]:02d}.csv").exists()]
        if not todo:
            break
        print(f"pass {npass + 1}: {len(todo)} months", flush=True)
        with ThreadPoolExecutor(4) as ex:
            for ym, msg in ex.map(safe, todo):
                print(ym, msg, flush=True)
    missing = [ym for ym in MONTHS if not (CACHE / f"{ym[0]}{ym[1]:02d}.csv").exists()]
    if missing:
        sys.exit(f"{len(missing)} months still missing: {missing}")
    d = pd.concat([pd.read_csv(CACHE / f"{y}{m:02d}.csv") for y, m in MONTHS]).sort_values("date")
    d.to_parquet(OUT, index=False)
    print(f"{len(d)} days {d.date.min()} .. {d.date.max()} -> {OUT.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
