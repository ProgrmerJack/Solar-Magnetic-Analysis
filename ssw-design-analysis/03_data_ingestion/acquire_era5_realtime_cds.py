#!/usr/bin/env python3
"""
acquire_era5_realtime_cds.py
============================
ERA5 / ERA5T (about five days behind real time) from the Copernicus CDS, reduced
to the series the prospective test needs (08_literature_audit/prospective_next_ssw.py):
  u10_60N      zonal-mean u at 10 hPa, 60N, mean of 00/06/12/18 UTC (m/s) -- onset
  z1000_m      polar-cap (>= 65N) cos-lat mean geopotential height at 1000 hPa, 00 UTC
  z150_m       the same at 150 hPa
  NEURASIA, HI_EUROPE, MID_EASIA, MID_NAMER   2 m temperature, cos-lat regional means,
               mean of 00/06/12/18 UTC (K)

WHY A SEPARATE SCRIPT
  Searched first (grep "reanalysis-era5-single-levels", "2m_temperature", "z1000_m"):
  the training series come from WeatherBench 2 (acquire_era5_nam.py,
  acquire_era5_t2m_regions.py), which ends on 10 January 2023, and from ARCO-ERA5,
  which is not updated within days; acquire_era5_u10_cds.py fetches fixed historical
  blocks of u only. A real-time test needs the same reductions from ERA5T.

SAME DEFINITIONS AS THE TRAINING DATA
  1.5-degree grid requested from the CDS (rows 90, 88.5, ... 66 are those of the
  WeatherBench 2 grid; the duplicate 180 column is dropped); geopotential / 9.80665;
  regions by acquire_s2s_reforecasts.regions_mean; u reduced by
  acquire_era5_u10_cds.reduce (0.25 degrees, 1440 longitudes). The CDS interpolates
  where WeatherBench 2 regridded conservatively, so --check measures the
  difference on a winter covered by both before the series is used (gate fixed in
  prospective_next_ssw.py: correlation >= 0.99 and RMS difference <= 0.10 of the
  day-of-year s.d., for each NAM level and northern-Eurasian temperature).

USAGE
  acquire_era5_realtime_cds.py --start 2026-11-01 --end 2027-04-30   (appends)
  acquire_era5_realtime_cds.py --check        (2021-12-01..2022-03-31, written to
                                                era5_realtime_cds_check.parquet)
The CDS key is read as in acquire_era5_u10_cds.client() and never printed.
Output: era5_realtime_cds.parquet (date index)
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import acquire_era5_u10_cds as U                    # noqa: E402  (client, reduce)
from acquire_s2s_reforecasts import regions_mean    # noqa: E402

RAW = HERE / "raw" / "era5_realtime_cds"
OUT = HERE / "era5_realtime_cds.parquet"
CHECK_OUT = HERE / "era5_realtime_cds_check.parquet"
CHECK_SPAN = ("2021-12-01", "2022-03-31")
G0 = 9.80665
CAP_LAT = 65.0
SYN = ["00:00", "06:00", "12:00", "18:00"]


def request(dataset, req, f):
    if f.exists() and f.stat().st_size > 0:
        return f
    c = U.client()
    for attempt in range(20):
        try:
            r = c.submit(dataset, {**req, "data_format": "netcdf", "download_format": "unarchived"})
            while r.status not in ("successful", "failed", "rejected"):
                time.sleep(20)
            if r.status != "successful":
                raise RuntimeError(r.status)
            tmp = f.with_suffix(".part")
            r.download(str(tmp))
            tmp.rename(f)
            print(f"  {f.name}: {f.stat().st_size / 1e6:.1f} MB", flush=True)
            return f
        except Exception as e:                       # noqa: BLE001 -- retried, then raised
            print(f"  {f.name} attempt {attempt + 1}: {type(e).__name__}: {str(e)[:200]}", flush=True)
            time.sleep(60)
    raise RuntimeError(f"{f.name} failed")


def month_days(y, m, a, b):
    days = pd.date_range(max(a, pd.Timestamp(y, m, 1)), min(b, pd.Timestamp(y, m, 1) + pd.offsets.MonthEnd(0)))
    return [f"{d.day:02d}" for d in days]


def fetch_month(y, m, a, b):
    days = month_days(y, m, a, b)
    tag = f"{y}{m:02d}_{days[0]}_{days[-1]}"
    base = {"product_type": ["reanalysis"], "year": [str(y)], "month": [f"{m:02d}"], "day": days}
    fz = request("reanalysis-era5-pressure-levels",
                 {**base, "variable": ["geopotential"], "pressure_level": ["150", "1000"], "time": ["00:00"],
                  "area": [90, -180, 64.5, 180], "grid": [1.5, 1.5]}, RAW / f"z_{tag}.nc")
    fu = request("reanalysis-era5-pressure-levels",
                 {**base, "variable": ["u_component_of_wind"], "pressure_level": ["10"], "time": SYN,
                  "area": [60, -180, 60, 180]}, RAW / f"u10_{tag}.nc")
    ft = request("reanalysis-era5-single-levels",
                 {**base, "variable": ["2m_temperature"], "time": SYN,
                  "area": [71, -180, 34, 180], "grid": [1.5, 1.5]}, RAW / f"t2m_{tag}.nc")
    return fz, fu, ft


def drop_dup_lon(da):
    lon = da["longitude"].values
    if np.isclose(lon, -180).any() and np.isclose(lon, 180).any():
        da = da.isel(longitude=np.flatnonzero(~np.isclose(lon, 180)))
    return da


def reduce_month(fz, fu, ft):
    z = xr.open_dataset(fz)["z"]
    tdim = "valid_time" if "valid_time" in z.dims else "time"
    z = drop_dup_lon(z)
    z = z.sel(latitude=z.latitude[z.latitude >= CAP_LAT])
    rows = z.latitude.values
    assert np.isclose(rows.min(), 66.0) and np.isclose(rows.max(), 90.0) and len(rows) == 17, rows
    assert z.sizes["longitude"] == 240, z.sizes["longitude"]
    cap = z.weighted(np.cos(np.deg2rad(z.latitude))).mean(("latitude", "longitude")) / G0
    plev = "pressure_level" if "pressure_level" in cap.dims else "level"
    t = pd.to_datetime(cap[tdim].values).floor("D")
    df = pd.DataFrame({f"z{int(p)}_m": cap.sel({plev: p}).values for p in (1000, 150)}, index=t)
    u = U.reduce(fu).rename("u10_60N")
    tt = xr.open_dataset(ft)["t2m"]
    tdt = "valid_time" if "valid_time" in tt.dims else "time"
    tt = drop_dup_lon(tt)
    r = regions_mean(tt).to_dataframe()
    r.index = pd.to_datetime(r.index if tdt not in r.columns else r[tdt])
    r = r[[c for c in ("NEURASIA", "HI_EUROPE", "MID_EASIA", "MID_NAMER")]]
    g = r.groupby(r.index.floor("D"))
    r = g.mean()[g.size() == 4]
    return df.join(u, how="outer").join(r, how="outer")


def acquire(a, b):
    RAW.mkdir(parents=True, exist_ok=True)
    a, b = pd.Timestamp(a), pd.Timestamp(b)
    parts = []
    for p in pd.period_range(a, b, freq="M"):
        parts.append(reduce_month(*fetch_month(p.year, p.month, a, b)))
    d = pd.concat(parts).sort_index()
    d.index.name = "date"
    return d[~d.index.duplicated()]


def main():
    argv = sys.argv[1:]
    if "--check" in argv:
        d = acquire(*CHECK_SPAN)
        d.to_parquet(CHECK_OUT)
        print(d.describe().round(2).to_string())
        print(f"Saved -> {CHECK_OUT.name}")
        return 0
    a, b = argv[argv.index("--start") + 1], argv[argv.index("--end") + 1]
    d = acquire(a, b)
    if OUT.exists():
        old = pd.read_parquet(OUT)
        # values are frozen as first retrieved (prospective_next_ssw registration):
        # a day already in the file keeps its earlier value; only new days are added
        d = pd.concat([old, d[~d.index.isin(old.index)]]).sort_index()
    if not d["u10_60N"].dropna().between(-60, 100).all():
        raise ValueError("implausible u(10 hPa, 60N)")
    d.to_parquet(OUT)
    print(d.tail().round(2).to_string())
    print(f"Saved -> {OUT.name} ({d.index.min().date()} .. {d.index.max().date()})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
