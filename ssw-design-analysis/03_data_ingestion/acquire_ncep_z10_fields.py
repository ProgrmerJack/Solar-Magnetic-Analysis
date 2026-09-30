#!/usr/bin/env python3
"""
acquire_ncep_z10_fields.py
==========================
Daily NCEP-NCAR reanalysis geopotential height at 10 hPa, 20-90N, every year
1958-2024, for the vortex-geometry (split/displacement) classification of
vortex_geometry_test.py (plan approved 2026-09-29).

Why a separate script: the other NCEP scripts reduce to polar-cap means on the
fly; the Seviour et al. (2013) moment method needs the 2-D field. Same transport
as extend_ncep_presatellite.py (THREDDS NCSS, one level, a latitude band) with
the whole-file fallback that script validated for years NCSS will not serve.

Output cache: _ncep_z10/z10_<year>.npz  z (time, lat, lon) float32 m; time (ns);
lat (descending as served); lon.
"""
import io
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import extend_ncep_presatellite as X                  # noqa: E402  (NCSS url, fallback)

CACHE = HERE / "_ncep_z10"
YEARS = range(1958, 2025)
SOUTH = 20.0


def fetch(year):
    out = CACHE / f"z10_{year}.npz"
    if out.exists():
        return year, "cached"
    url = (X.NCSS.format(var="hgt", year=year) +
           f"?var=hgt&north=90&south={SOUTH}&west=0&east=357.5&vertCoord=10"
           f"&time_start={year}-01-01T00:00:00Z&time_end={year}-12-31T23:59:59Z&accept=netcdf")
    ds, last = None, None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=600) as r:
                blob = r.read()
            if len(blob) > 1_000_000:
                ds = xr.open_dataset(io.BytesIO(blob))
                break
            last = f"short ({len(blob)} B)"
        except Exception as exc:
            last = exc
        time.sleep(3 * (attempt + 1))
    if ds is None:                                   # whole-file fallback (e.g. 1975)
        local = X._whole_file("hgt", year)
        ds = xr.open_dataset(local).sel(level=10.0, lat=slice(90, SOUTH))
        last = "whole-file"
    z = ds["hgt"].squeeze()
    if "level" in z.dims:
        z = z.isel(level=0)
    z = z.sel(lat=slice(90, SOUTH)) if z.lat[0] > z.lat[-1] else z.sel(lat=slice(SOUTH, 90))
    arr = z.values.astype(np.float32)
    if not (25000 < np.nanmean(arr) < 33000):
        raise ValueError(f"{year}: 10 hPa height mean {np.nanmean(arr):.0f} m is implausible")
    np.savez_compressed(out, z=arr, time=z.time.values.astype("datetime64[ns]").astype(np.int64),
                        lat=z.lat.values, lon=z.lon.values)
    ds.close()
    return year, f"fetched ({last or 'ncss'})"


def main():
    CACHE.mkdir(exist_ok=True)
    with ThreadPoolExecutor(4) as ex:
        for y, how in ex.map(fetch, YEARS):
            print(y, how, flush=True)


if __name__ == "__main__":
    main()
