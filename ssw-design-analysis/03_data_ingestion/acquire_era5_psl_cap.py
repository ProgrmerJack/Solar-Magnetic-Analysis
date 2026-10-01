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
  Every NH winter Nov-Apr 1998/99 .. 2021/22 (SNAPSI's two winters and the S2S
  reforecast span) and the SH case (2019-08-01 .. 2019-11-30). Both caps are
  computed for every time, as in the SNAPSI reduction, so the question can never
  be asked of the wrong hemisphere.

Output: era5_psl_cap_6h.parquet  time, psl_cap_N, psl_cap_S  (Pa)

--z100 (plan approved 2026-09-30, operational lead-lag diagnosis): geopotential
height at 100 hPa (geopotential / 9.80665, gpm), 00 UTC, the same NH cap and the
same winters -> era5_z100_cap_00utc.parquet  time, z100_cap_N (gpm)
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
# Every NH winter (Nov-Apr) 1998/99-2021/22 -- the span of the ECMWF S2S
# reforecasts used by the forecast-reliability test, with margin for its
# leave-one-year-out climatology -- plus the SH SNAPSI case. The first version
# held only the two SNAPSI winters; the reduction is unchanged, so their values
# are identical.
PERIODS = ([(f"{y}-11-01", f"{y + 1}-04-30T18") for y in range(1998, 2022)]
           + [("2019-08-01", "2019-11-30T18")])


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
    # The SH period (Aug-Nov 2019) overlaps the NH winter 2019/20 from 1 Nov, so
    # 2019-11 was written twice until 2026-09-25 (identical values). Keep one
    # copy, and refuse if overlapping periods ever disagree.
    dup = out[out["time"].duplicated(keep=False)]
    if len(dup) and (dup.groupby("time")[["psl_cap_N", "psl_cap_S"]].nunique() > 1).any().any():
        raise ValueError("overlapping periods give different values for the same time")
    out = out.drop_duplicates("time").sort_values("time", ignore_index=True)
    if out[["psl_cap_N", "psl_cap_S"]].isna().any().any():
        raise ValueError("NaN in the ERA5 polar cap")
    for c in ("psl_cap_N", "psl_cap_S"):
        m = float(out[c].mean())
        if not 95000.0 < m < 105000.0:
            raise ValueError(f"{c} mean {m:.0f} Pa is not sea-level pressure")
    out.to_parquet(OUT)
    print(f"{len(out):,} rows -> {OUT.name}")


G0 = 9.80665
OUT_Z100 = HERE / "era5_z100_cap_00utc.parquet"


def main_z100():
    """Chunk-aligned parallel reads (the store's chunks hold 8 time steps x 13
    levels; a single sequential read stalled for over an hour, 2026-09-30), each
    with a timeout and retries."""
    import concurrent.futures as cf
    import aiohttp
    # the same public store over plain HTTPS with socket timeouts: through gcsfs a
    # stalled request hung the whole read twice (2026-09-30)
    url = STORE.replace("gs://", "https://storage.googleapis.com/")
    ds = xr.open_zarr(url, chunks=None, storage_options={
        "client_kwargs": {"timeout": aiohttp.ClientTimeout(total=180, sock_read=90)}})
    v = ds["geopotential"]
    if v.attrs.get("units") != "m**2 s**-2":
        raise ValueError(f"unexpected units {v.attrs.get('units')!r}")
    li = int(np.flatnonzero(ds.level.values == 100)[0])
    tt = pd.DatetimeIndex(ds.time.values)
    want = np.zeros(len(tt), bool)
    for a, b in PERIODS[:-1]:                       # NH winters only
        want |= (tt >= pd.Timestamp(a)) & (tt <= pd.Timestamp(b)) & (tt.hour == 0)
    blocks = sorted(set(np.flatnonzero(want) // 8))

    def read(k):
        sub = v.isel(time=slice(8 * k, 8 * k + 8), level=li)
        keep = want[8 * k: 8 * k + 8]
        x = sub.isel(time=np.flatnonzero(keep)).load() / G0
        return pd.DataFrame({"time": x.time.values, "z100_cap_N": cap(x, True).values})

    frames, todo = [], list(blocks)
    for attempt in range(5):
        failed = []
        with cf.ThreadPoolExecutor(16) as ex:
            futs = {ex.submit(read, k): k for k in todo}
            for f in cf.as_completed(futs):
                try:
                    frames.append(f.result(timeout=300))
                except Exception:
                    failed.append(futs[f])
        print(f"  pass {attempt + 1}: {len(todo) - len(failed)} of {len(todo)} blocks", flush=True)
        if not failed:
            break
        todo = failed
    if failed:
        raise IOError(f"{len(failed)} blocks failed")
    out = pd.concat(frames, ignore_index=True).sort_values("time", ignore_index=True)
    if out["time"].duplicated().any() or out["z100_cap_N"].isna().any():
        raise ValueError("duplicate or missing 100 hPa cap values")
    if not out["z100_cap_N"].between(14000, 18000).all():
        raise ValueError("100 hPa cap height outside 14-18 km")
    if len(out) != int(want.sum()):
        raise ValueError(f"{len(out)} rows for {int(want.sum())} wanted days")
    out.to_parquet(OUT_Z100)
    print(f"{len(out):,} rows -> {OUT_Z100.name}")


if __name__ == "__main__":
    import sys
    # --full: every November-May 1958/59-2022/23 (the WB2 span), for the revision-3
    # ERA5 continuity test (vortex_threshold_continuity.py, commit add7ba9); a
    # separate file, so the registered inputs above are untouched
    if "--full" in sys.argv:
        PERIODS = [(f"{y}-11-01", f"{y + 1}-05-31T18") for y in range(1958, 2022)] + [("2022-11-01", "2023-01-10T18")]
        PERIODS = [("1959-01-01", "1959-05-31T18")] + PERIODS[1:]
        OUT = HERE / "era5_psl_cap_6h_full.parquet"
    main_z100() if "--z100" in sys.argv else main()
