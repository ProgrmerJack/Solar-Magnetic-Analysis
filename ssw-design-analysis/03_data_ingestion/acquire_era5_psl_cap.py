#!/usr/bin/env python3
"""
acquire_era5_psl_cap.py
=======================
Observed polar-cap mean sea-level pressure, reduced EXACTLY as the SNAPSI members
are, so the observed outcome of each SNAPSI event can be placed inside its own
ensemble.

WHY A SEPARATE SCRIPT
  Searched first: nothing in the repo extracts ERA5 sea-level pressure
  ("mean_sea_level_pressure", "msl"). `acquire_era5_nam.py` reads the same store
  but builds a DIFFERENT product -- daily 00 UTC geopotential height over
  65-90N, standardised by day of year. The SNAPSI comparison needs raw 6-hourly
  psl over 60-90N, cos-lat weighted, in Pa, which is what
  `acquire_snapsi_surface.py` computes per member. Adding it to the NAM script
  would restamp that script's hash and mark its downstream results stale.

SOURCE
  WeatherBench2 ERA5 on public GCS, anonymous, no CDS account:
  gs://weatherbench2/datasets/era5/1959-2023_01_10-6h-240x121_equiangular_with_poles_conservative.zarr
  6-hourly, 1.5 deg. The cap mean over 60-90N at 1.5 deg is recorded as a
  resolution difference from the model grids, not glossed.

PERIODS
  Both NH SNAPSI winters with margin (2017-12-01 .. 2018-04-30 and
  2018-12-01 .. 2019-04-30) and the SH case (2019-08-01 .. 2019-11-30). Both caps
  are computed for every time, as in the SNAPSI reduction, so the question can
  never be asked of the wrong hemisphere.

Output: era5_psl_cap_6h.parquet  time, psl_cap_N, psl_cap_S  (Pa)
"""
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
OUT = HERE / "era5_psl_cap_6h.parquet"
STORE = ("gs://weatherbench2/datasets/era5/"
         "1959-2023_01_10-6h-240x121_equiangular_with_poles_conservative.zarr")
CAP_LAT = 60.0
PERIODS = [("2017-12-01", "2018-04-30T18"), ("2018-12-01", "2019-04-30T18"),
           ("2019-08-01", "2019-11-30T18")]


# WeatherBench2 stores the 60N row as 59.999999999999986, so `lat >= 60` silently
# dropped it and built a 61.5-90N cap (the SH cap kept -60 exactly). Seven of
# the nine SNAPSI grids carry an exact 60.0 row and keep it, KMA and SNU on this
# same 1.5-degree grid among them, so the row belongs in. Tolerance, not rounding.
LAT_TOL = 1e-6


def cap(da, north):
    lat = da.latitude
    sub = da.sel(latitude=lat[lat >= CAP_LAT - LAT_TOL] if north
                 else lat[lat <= -CAP_LAT + LAT_TOL])
    if abs(float(abs(sub.latitude).min()) - CAP_LAT) > 1e-3:
        raise ValueError(f"cap edge at {float(abs(sub.latitude).min())}, not {CAP_LAT}")
    w = np.cos(np.deg2rad(sub.latitude))
    return sub.weighted(w).mean(dim=["latitude", "longitude"])


def main():
    ds = xr.open_zarr(STORE, storage_options={"token": "anon"})
    v = ds["mean_sea_level_pressure"]
    if v.attrs.get("units") != "Pa":
        raise ValueError(f"unexpected units {v.attrs.get('units')!r}")
    frames = []
    for a, b in PERIODS:
        s = v.sel(time=slice(a, b)).load()
        n, so = cap(s, True).values, cap(s, False).values
        frames.append(pd.DataFrame({"time": s.time.values,
                                    "psl_cap_N": n, "psl_cap_S": so}))
        print(f"  {a} .. {b}: {len(s.time)} steps, "
              f"N mean {np.nanmean(n):.0f} Pa, S mean {np.nanmean(so):.0f} Pa", flush=True)
    out = pd.concat(frames, ignore_index=True)
    if out[["psl_cap_N", "psl_cap_S"]].isna().any().any():
        raise ValueError("NaN in the ERA5 polar cap")
    for c in ("psl_cap_N", "psl_cap_S"):
        m = float(out[c].mean())
        if not 95000.0 < m < 105000.0:
            raise ValueError(f"{c} mean {m:.0f} Pa is not sea-level pressure")
    out.to_parquet(OUT)
    print(f"{len(out):,} rows -> {OUT.name}")


if __name__ == "__main__":
    main()
