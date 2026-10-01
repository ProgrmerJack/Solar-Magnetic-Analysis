#!/usr/bin/env python3
"""
acquire_era5_t2m_regions.py
===========================
ERA5 daily-mean 2 m temperature over the four regions of the regional cold-risk
plan (approved 2026-09-26), the observation for the S2S regional test and the
observed-event regional analysis.

SOURCE AND REDUCTION
  WeatherBench 2 ERA5, the same public store as acquire_era5_psl_cap.py
  (6-hourly, 1.5 degree, 1959-01-01 .. 2023-01-10, anonymous), variable
  `2m_temperature`. Every November-April 1959-2023 in the store. Regions are
  taken from acquire_s2s_reforecasts.T2M_REGIONS (one definition for forecasts
  and observations) and reduced with the same cos-lat mean (regions_mean).
  Daily mean = mean of 00, 06, 12, 18 UTC of the UTC day; days without all four
  steps are dropped. Winters 2023/24 onward are not in this store.

Output: era5_t2m_regions_daily.parquet  date, NEURASIA, HI_EUROPE, MID_EASIA,
        MID_NAMER (K)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import acquire_era5_psl_cap as E                     # noqa: E402  (store)
import acquire_s2s_reforecasts as A                  # noqa: E402  (regions, reduction)

OUT = HERE / "era5_t2m_regions_daily.parquet"
PERIODS = [(f"{y}-11-01", f"{y + 1}-04-30T18") for y in range(1958, 2023)]


def main():
    ds = xr.open_zarr(E.STORE, storage_options={"token": "anon"})
    v = ds["2m_temperature"]
    if v.attrs.get("units") not in ("K", None):
        raise ValueError(f"unexpected units {v.attrs.get('units')!r}")
    lat = v.latitude
    v = v.sel(latitude=lat[(lat >= 34) & (lat <= 71)])
    frames = []
    for a, b in PERIODS:
        s = v.sel(time=slice(a, b))
        if s.time.size == 0:
            continue
        r = A.regions_mean(s.load()).to_dataframe()
        r["date"] = r.index.floor("D")
        g = r.groupby("date")
        d = g[list(A.T2M_REGIONS)].mean()[g.size() == 4]
        frames.append(d)
        print(f"  {a[:4]}/{b[:4]}: {len(d)} days, NEURASIA mean {d['NEURASIA'].mean():.1f} K", flush=True)
    out = pd.concat(frames).reset_index()
    if out["date"].duplicated().any():
        raise ValueError("duplicate dates")
    vals = out[list(A.T2M_REGIONS)]
    if vals.isna().any().any() or not vals.stack().between(200, 310).all():
        raise ValueError("implausible regional 2 m temperature")
    out.to_parquet(OUT)
    print(f"{len(out):,} days ({out.date.min().date()} .. {out.date.max().date()}) -> {OUT.name}")


if __name__ == "__main__":
    # --full: November-May (adds May) 1959-2022, for the revision-3 ERA5 continuity
    # test and weeks 3-6 (commit add7ba9); a separate file
    if "--full" in sys.argv:
        PERIODS = [(f"{y}-11-01", f"{y + 1}-05-31T18") for y in range(1958, 2023)]
        OUT = HERE / "era5_t2m_regions_daily_full.parquet"
    main()
