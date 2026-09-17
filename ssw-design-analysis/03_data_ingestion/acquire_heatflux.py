#!/usr/bin/env python3
"""
acquire_heatflux.py
===================
Acquires the gridded fields needed to compute the canonical SSW precursor
diagnostic — the 100 hPa eddy heat flux — plus polar-cap height, for Gate 8
(physical decomposition: what predicts design sensitivity?).

WHY THIS AND NOT SNAPSI
  SNAPSI is the preferred ground truth but sits behind CEDA authentication,
  which is unavailable today. Gate 8 is the one substantive gate that needs
  neither SNAPSI nor a second coder, and it needs per-event dynamical
  diagnostics that this repository does not have.

  `ncep_stratosphere.parquet` already holds zonal/area-mean air, height and
  zonal wind at six levels. It cannot give eddy heat flux, because v'T' is a
  covariance across LONGITUDE and is destroyed by zonal averaging. The full
  longitude grid is therefore required.

WHAT IS DOWNLOADED
  NCEP/NCAR Reanalysis 1 daily, via NOAA PSL OPeNDAP with server-side
  subsetting — open, no authentication.

    vwnd, air   level index 11 (100 hPa), lat indices 6..18 (75N..45N),
                all 144 longitudes, daily, 1958-2024

  Subsetting is what makes this tractable: the full pressure-level files are
  ~110-160 MB per variable per year (≈17 GB for both variables over the
  record). Restricting to one level and 13 latitudes gives ~2.7 MB per variable
  per year, ≈360 MB in total.

  Eddy heat flux is then  [v*T*]  averaged 45-75N, where * is the departure from
  the zonal mean — the standard wave-driving diagnostic used to characterise the
  forcing that precedes a warming.

Outputs (this directory):
  raw/ncep_heatflux/<var>.<year>.nc     subset files
  heatflux_daily.parquet                daily [v*T*] at 100 hPa, 45-75N
  heatflux_checksums.txt                sha256 per file
"""
import hashlib
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw" / "ncep_heatflux"
OUT = HERE / "heatflux_daily.parquet"
SUMS = HERE / "heatflux_checksums.txt"

BASE = ("https://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis/"
        "Dailies/pressure/{var}.{year}.nc.nc")
LEV = 11          # 100 hPa
LAT0, LAT1 = 6, 18   # 75.0N .. 45.0N inclusive
NLON = 144
YEARS = range(1958, 2025)
VARS = ("vwnd", "air")
TIMEOUT = 300


DODS = ("https://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis/"
        "Dailies/pressure/{var}.{year}.nc")


def year_heatflux(year):
    """Remote-slice one year at 100 hPa, 45-75N, and return daily [v*T*].

    netCDF4 opens the OPeNDAP endpoint and transfers only the requested slice,
    so the ~270 MB of full-grid files for a year become ~5 MB on the wire.
    A THREDDS `.nc` response is not offered (400); only `.dods` is, and letting
    netCDF4 speak DAP avoids parsing that format by hand.
    """
    import netCDF4
    out = {}
    for var in VARS:
        d = netCDF4.Dataset(DODS.format(var=var, year=year))
        try:
            a = np.asarray(d.variables[var][:, LEV, LAT0:LAT1 + 1, :])
            if var == VARS[0]:
                tv = netCDF4.num2date(d.variables["time"][:],
                                      d.variables["time"].units,
                                      only_use_cftime_datetimes=False)
                out["time"] = pd.to_datetime([pd.Timestamp(x) for x in tv])
                out["lat"] = np.asarray(d.variables["lat"][LAT0:LAT1 + 1])
        finally:
            d.close()
        out[var] = a
    v, T = out["vwnd"], out["air"]
    n = min(v.shape[0], T.shape[0], len(out["time"]))
    v, T = v[:n], T[:n]
    vp = v - v.mean(axis=2, keepdims=True)      # departure from zonal mean
    Tp = T - T.mean(axis=2, keepdims=True)
    vt = (vp * Tp).mean(axis=2)                 # zonal-mean covariance (time, lat)
    w = np.cos(np.deg2rad(out["lat"]))
    hf = (vt * w[None, :]).sum(axis=1) / w.sum()
    return pd.DataFrame({"date": out["time"][:n], "vT_100hPa_45_75N": hf})


def main():
    CACHE = RAW
    CACHE.mkdir(parents=True, exist_ok=True)
    frames, failed, cached = [], [], 0
    for year in YEARS:
        f = CACHE / f"hf_{year}.parquet"
        if f.exists():
            frames.append(pd.read_parquet(f)); cached += 1; continue
        try:
            df = year_heatflux(year)
            df.to_parquet(f, index=False)
            frames.append(df)
        except Exception as e:
            failed.append((year, f"{type(e).__name__}: {str(e)[:80]}"))
            continue
        if year % 10 == 0:
            print(f"  {year}: {len(df)} days, mean v'T' = "
                  f"{df['vT_100hPa_45_75N'].mean():+.1f} K m/s", flush=True)

    if not frames:
        print("no data acquired"); return 1
    df = pd.concat(frames, ignore_index=True).sort_values("date")
    df = df.drop_duplicates(subset="date")
    df.to_parquet(OUT, index=False)

    sums = [(f.name, hashlib.sha256(f.read_bytes()).hexdigest())
            for f in sorted(CACHE.glob("hf_*.parquet"))]
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    SUMS.write_text(
        "\n".join(f"{h}  {n}" for n, h in sums)
        + f"\n# NCEP/NCAR R1 100hPa v'T' 45-75N via NOAA PSL OPeNDAP, {stamp}\n",
        encoding="utf8")

    print(f"\n{len(df):,} daily values "
          f"{df['date'].min().date()}..{df['date'].max().date()}"
          f"  ({cached} years from cache)")
    w = df[df["date"].dt.month.isin((11, 12, 1, 2, 3, 4))]
    print(f"winter days {len(w):,}; v'T' mean {w['vT_100hPa_45_75N'].mean():+.1f} "
          f"sd {w['vT_100hPa_45_75N'].std():.1f} K m/s")
    if failed:
        print(f"failed {len(failed)}: {failed[:4]}")
    print(f"-> {OUT.name}, {SUMS.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
