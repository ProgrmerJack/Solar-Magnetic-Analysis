#!/usr/bin/env python3
"""
acquire_era5_nam.py
===================
Builds the NAM indices in the SOURCE VARIABLES the published criterion is defined
on, so the body count can be computed in the literature's own terms rather than in
CPC index proxies.

WHY
  `recompute_published_criterion.py` applied Karpechko et al. (2017) conditions 1-3
  using the CPC AO as a stand-in for the 1000 hPa NAM and the CPC NAO as the
  outcome. That establishes the DECOMPOSITION but cannot speak to any published
  number, because ACP 26, 3723 (2026) uses ERA5 850 hPa NAM for classification and
  a 1000 hPa-geopotential NAO for the outcome, reporting -0.620 (all DW) vs +0.088
  (NDW).

SOURCE
  WeatherBench2's ERA5 copy on Google Cloud Storage, public and anonymous:
  gs://weatherbench2/datasets/era5/1959-2023_01_10-6h-240x121_equiangular_with_poles_conservative.zarr
  1959-01-01..2023-01-10, 6-hourly, 1.5 deg, geopotential on 13 pressure levels
  including 1000, 850 and 150 hPa -- the exact three the criterion needs.

  No CDS registration, no API key. Resolution is 1.5 deg rather than 0.25 deg;
  for a polar-cap-mean index that is immaterial, and it is recorded here rather
  than glossed.

NAM DEFINITION
  Baldwin & Thompson (2009) polar-cap form, which is what the SSW literature means
  by "NAM index at level p": the cos-lat-weighted mean geopotential height anomaly
  over 65-90 N, negated and standardised by day of year. No EOF, no rotation, no
  analyst degrees of freedom.

  Geopotential (m^2/s^2) is divided by g0 = 9.80665 to give geopotential HEIGHT.

Output: era5_nam_daily.parquet with nam_1000, nam_850, nam_150 (standardised),
        plus the raw polar-cap heights for audit.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
OUT = HERE / "era5_nam_daily.parquet"

STORE = ("gs://weatherbench2/datasets/era5/"
         "1959-2023_01_10-6h-240x121_equiangular_with_poles_conservative.zarr")
LEVELS = [1000, 850, 150]
CAP_LAT = 65.0            # Baldwin & Thompson polar cap
G0 = 9.80665


def main():
    print(f"opening {STORE.split('/')[-1]} ...")
    ds = xr.open_zarr(STORE, storage_options={"token": "anon"}, chunks={"time": 500})
    z = ds["geopotential"].sel(level=LEVELS)
    z = z.sel(latitude=z.latitude[z.latitude >= CAP_LAT])
    # one timestep per day (00 UTC) -- a daily-mean polar-cap index is not
    # sensitive to the diurnal cycle at these levels, and this cuts the read 4x
    z = z.isel(time=slice(0, None, 4))
    print(f"  reading {z.sizes['time']:,} daily steps x {z.sizes['level']} levels "
          f"x {z.sizes['latitude']} lat x {z.sizes['longitude']} lon")

    w = np.cos(np.deg2rad(z.latitude))
    cap = (z.weighted(w).mean(dim=["latitude", "longitude"]) / G0).compute()
    print("  polar-cap means computed")

    t = pd.to_datetime(cap.time.values)
    df = pd.DataFrame(index=t)
    for lev in LEVELS:
        df[f"z{lev}_m"] = cap.sel(level=lev).values

    doy = df.index.dayofyear
    for lev in LEVELS:
        h = df[f"z{lev}_m"]
        # NAM = negated, day-of-year standardised polar-cap height anomaly
        mu = h.groupby(doy).transform("mean")
        sd = h.groupby(doy).transform("std")
        df[f"nam_{lev}"] = -((h - mu) / sd)

    df.index.name = "date"
    df.to_parquet(OUT)
    print(f"\n{df.index.min().date()} .. {df.index.max().date()}  {len(df):,} days")
    print(df[[f"nam_{l}" for l in LEVELS]].describe().round(3).to_string())
    print(f"\nSaved -> {OUT.name}")

    # sanity: NAM at 1000 hPa must track the CPC AO closely if it is right
    try:
        sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
        import multi_index_event_study as M
        ao = M.load("ao")["y"]
        j = pd.DataFrame({"nam": df["nam_1000"], "ao": ao.reindex(df.index)}).dropna()
        r = float(j["nam"].corr(j["ao"]))
        print(f"\nVALIDATION: corr(ERA5 NAM 1000 hPa, CPC AO) = {r:+.3f} "
              f"on {len(j):,} shared days")
        print("  (these are different constructions of the same mode; anything"
              " below ~0.7 means something is wrong)")
    except Exception as exc:
        print(f"validation skipped: {exc}")


if __name__ == "__main__":
    main()
