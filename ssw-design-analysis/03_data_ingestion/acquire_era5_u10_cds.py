#!/usr/bin/env python3
"""
acquire_era5_u10_cds.py
=======================
ERA5 zonal-mean zonal wind at 10 hPa, 60 N, daily mean of 00/06/12/18 UTC,
1 January 1940 - 30 April 2026: the running variable of the ERA5 continuity test
and the series for detecting SSWs before 1958 (vortex_threshold_continuity.py,
shift_rule_forecast.py; designs in their docstrings).

WHY A SEPARATE SCRIPT
  Searched first (grep "u_component_of_wind", "reanalysis-era5-pressure-levels"):
  nothing in the repo retrieves ERA5 winds. ARCO-ERA5 stores each hour of a
  pressure-level variable as one 37-level global chunk (~150 MB), so an 86-year
  10 hPa series from it would read >1 TB; the Copernicus CDS returns one latitude
  row at one level.

SOURCE
  Copernicus Climate Data Store, reanalysis-era5-pressure-levels (0.25 deg),
  u_component_of_wind, 10 hPa, area [60, -180, 60, 180]; the CDS personal token
  is read from ~/.cdsapirc (same ECMWF account as the S2S retrievals) and never
  printed. One request per five calendar years.

REDUCTION
  Zonal mean over the 1440 distinct longitudes of the 60 N row (the duplicate
  180/-180 column dropped), then the mean of the four synoptic hours of each UTC
  day (days with fewer than four are dropped).

Output: era5_u10_60N_daily_cds.parquet   date, u10_60N (m/s)
"""
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw" / "era5_u10_cds"
OUT = HERE / "era5_u10_60N_daily_cds.parquet"
BLOCKS = [(y, y + 4, 12) for y in range(1940, 2025, 5)] + [(2025, 2025, 12), (2026, 2026, 8)]
WORKERS = 4


def client():
    from ecmwf.datastores import Client
    key = Path("~/.cdsapirc").expanduser().read_text().split("key:")[1].split()[0]
    return Client(url="https://cds.climate.copernicus.eu/api", key=key)


def fetch(block):
    a, b, last_month = block
    f = RAW / f"u10_60N_{a}_{b}.nc"
    if f.exists() and f.stat().st_size > 0:
        return f
    c = client()
    req = {"product_type": ["reanalysis"], "variable": ["u_component_of_wind"],
           "pressure_level": ["10"], "year": [str(y) for y in range(a, b + 1)],
           "month": [f"{m:02d}" for m in range(1, last_month + 1)], "day": [f"{d:02d}" for d in range(1, 32)],
           "time": ["00:00", "06:00", "12:00", "18:00"], "area": [60, -180, 60, 180],
           "data_format": "netcdf", "download_format": "unarchived"}
    for attempt in range(20):
        try:
            r = c.submit("reanalysis-era5-pressure-levels", req)
            while r.status not in ("successful", "failed", "rejected"):
                time.sleep(30)
            if r.status != "successful":
                raise RuntimeError(f"request {a}-{b} {r.status}")
            tmp = f.with_suffix(".part")
            r.download(str(tmp))
            tmp.rename(f)
            print(f"  {a}-{b}: {f.stat().st_size / 1e6:.1f} MB", flush=True)
            return f
        except Exception as e:                       # noqa: BLE001 -- retried, then raised
            print(f"  {a}-{b} attempt {attempt + 1}: {type(e).__name__}: {str(e)[:200]}", flush=True)
            time.sleep(60)
    raise RuntimeError(f"block {a}-{b} failed")


def reduce(f):
    ds = xr.open_dataset(f)
    u = ds["u"]
    tdim = "valid_time" if "valid_time" in u.dims else "time"
    lon = u["longitude"].values
    keep = ~np.isclose(lon, 180.0) if np.isclose(lon, -180.0).any() else np.ones(len(lon), bool)
    zm = u.isel(longitude=np.flatnonzero(keep)).squeeze(drop=True)
    if "pressure_level" in zm.dims:
        zm = zm.squeeze("pressure_level", drop=True)
    zm = zm.mean("longitude")
    s = pd.Series(zm.values, index=pd.to_datetime(zm[tdim].values))
    g = s.groupby(s.index.floor("D"))
    d = g.mean()[g.count() == 4]
    if int(np.sum(keep)) != 1440:
        raise ValueError(f"{f.name}: {int(np.sum(keep))} longitudes, expected 1440")
    return d


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(WORKERS) as ex:
        files = list(ex.map(fetch, BLOCKS))
    s = pd.concat([reduce(f) for f in files]).sort_index()
    s = s[~s.index.duplicated()]
    s = s[(s.index >= "1940-01-01") & (s.index <= "2026-04-30")]
    if not s.between(-60, 100).all():
        raise ValueError("implausible u(10 hPa, 60N)")
    gaps = pd.date_range(s.index.min(), s.index.max(), freq="D").difference(s.index)
    print(f"{len(s)} days {s.index.min().date()} .. {s.index.max().date()}; missing days: {len(gaps)}")
    pd.DataFrame({"date": s.index, "u10_60N": s.values}).to_parquet(OUT)
    print(f"-> {OUT.name}")


if __name__ == "__main__":
    sys.exit(main())
