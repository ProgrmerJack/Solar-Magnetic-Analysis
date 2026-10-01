#!/usr/bin/env python3
"""
acquire_era5_arco_extension.py
==============================
ERA5 observations for the held-out 2023-2024 SSWs (plan approved 2026-09-30;
design in 07_physical_decomposition/s2s_heldout_test.py), which the WeatherBench 2
copy used everywhere else does not reach (it ends 2023-01-10; the polar-cap
extraction ends 2022-04-30).

WHY A SEPARATE SCRIPT
  Searched first: nothing in the repo reads ARCO-ERA5 (grep "arco",
  "gcp-public-data"). The WeatherBench 2 scripts read a different store, on a
  different grid; folding this into them would restamp their hashes.

SOURCE
  ARCO-ERA5, gs://gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3,
  anonymous, hourly, 0.25 deg (one chunk per time step and variable).

REDUCTION (made to match the WeatherBench 2 series, then checked on an overlap)
  Each field is first averaged onto the WeatherBench 2 grid (1.5 deg, poles
  included, longitudes 0..358.5): a target cell centred at (phi, lam) averages the
  0.25 deg points within +-0.75 deg in each direction, points on a cell edge
  weighted 1/2, all weighted by cos(latitude) -- a conservative average for these
  nested grids. Then:
    mean_sea_level_pressure, 00 UTC -> NH polar cap exactly as
      acquire_era5_psl_cap.cap (cos-lat weighted 60-90N, the 60N row included)
    2m_temperature, 00/06/12/18 UTC -> the four regions of acquire_s2s_reforecasts
      (regions_mean), daily mean of the four steps of the UTC day
  Winters (Nov-Apr) 2020/21 .. 2023/24. 2020/21-2021/22 (both variables) and
  2022-11-01 .. 2023-01-10 (temperature) overlap the WeatherBench 2 series; the
  held-out test removes the mean ARCO-minus-WB2 difference over that overlap
  (one constant per variable/region, fixed in the design) and reports it with the
  daily RMS difference and correlation.

Output: era5_psl_cap_00utc_arco.parquet   time, psl_cap_N (Pa)
        era5_t2m_regions_daily_arco.parquet  date, NEURASIA, HI_EUROPE, MID_EASIA,
                                            MID_NAMER (K)
"""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import acquire_era5_psl_cap as E                     # noqa: E402  (cap reduction)
import acquire_s2s_reforecasts as A                  # noqa: E402  (regions, reduction)

STORE = "gs://gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
OUT_PSL = HERE / "era5_psl_cap_00utc_arco.parquet"
OUT_T2M = HERE / "era5_t2m_regions_daily_arco.parquet"
WINTERS = [(f"{y}-11-01", f"{y + 1}-04-30") for y in range(2020, 2024)]
T_LAT = np.linspace(-90, 90, 121)
T_LON = np.arange(240) * 1.5
WORKERS = int(__import__("os").environ.get("ARCO_WORKERS", 16))
PARTS = HERE / "raw" / "arco_parts"


def weights_1d(src, tgt, cyclic):
    """Rows: target cells; columns: source points; trapezoid weights within
    +-0.75 deg (edge points 1/2)."""
    W = np.zeros((len(tgt), len(src)))
    for i, c in enumerate(tgt):
        d = np.abs(src - c)
        if cyclic:
            d = np.minimum(d, 360.0 - d)
        W[i, d < 0.75 - 1e-9] = 1.0
        W[i, np.abs(d - 0.75) <= 1e-9] = 0.5
    return W


def regridder(lat, lon):
    wl = weights_1d(lat, T_LAT, False) * np.cos(np.deg2rad(lat))[None, :]
    wl = wl / wl.sum(1, keepdims=True)
    wo = weights_1d(lon, T_LON, True)
    wo = wo / wo.sum(1, keepdims=True)

    def f(x):                                         # x: (lat, lon)
        return wl @ x @ wo.T
    return f


def to_da(a, t):
    return xr.DataArray(a, dims=("latitude", "longitude"),
                        coords={"latitude": T_LAT, "longitude": T_LON}).expand_dims(time=[t])


def main():
    global WINTERS, OUT_PSL, OUT_T2M
    ds = xr.open_zarr(STORE, chunks=None, storage_options={"token": "anon"})
    lat, lon = ds.latitude.values, ds.longitude.values
    rg = regridder(lat, lon)
    msl, t2m = ds["mean_sea_level_pressure"], ds["2m_temperature"]

    def psl_at(t):
        x = rg(msl.sel(time=t).values)
        return t, float(E.cap(to_da(x, t), north=True).values[0])

    def t2m_at(t):
        x = rg(t2m.sel(time=t).values)
        r = A.regions_mean(to_da(x, t))
        return t, {k: float(r[k].values[0]) for k in A.T2M_REGIONS}

    ps, ts = [], []
    PARTS.mkdir(parents=True, exist_ok=True)
    for a, b in WINTERS:
        fp, ft = PARTS / f"psl_{a}_{b}.parquet", PARTS / f"t2m_{a}_{b}.parquet"
        if fp.exists() and ft.exists():                # resume: a winter already reduced
            ps += [tuple(r) for r in pd.read_parquet(fp)[["time", "psl_cap_N"]].itertuples(index=False)]
            ts.append(pd.read_parquet(ft).set_index("date"))
            print(f"  {a[:4]}/{b[:4]}: cached", flush=True)
            continue
        days = pd.date_range(a, b, freq="D")
        with ThreadPoolExecutor(WORKERS) as ex:
            p = list(ex.map(psl_at, days))
            steps = [d + pd.Timedelta(hours=h) for d in days for h in (0, 6, 12, 18)]
            t = list(ex.map(t2m_at, steps))
        ps += p
        tt = pd.DataFrame([{"time": s, **v} for s, v in t])
        tt["date"] = tt["time"].dt.floor("D")
        g = tt.groupby("date")
        tw = g[list(A.T2M_REGIONS)].mean()[g.size() == 4]
        ts.append(tw)
        pd.DataFrame(p, columns=["time", "psl_cap_N"]).to_parquet(fp)
        tw.reset_index().to_parquet(ft)
        print(f"  {a[:4]}/{b[:4]}: {len(p)} days", flush=True)
    psl = pd.DataFrame(ps, columns=["time", "psl_cap_N"])
    t2 = pd.concat(ts).reset_index()
    if not psl["psl_cap_N"].between(95000, 106000).all():
        raise ValueError("ARCO polar cap outside 950-1060 hPa")
    if t2[list(A.T2M_REGIONS)].isna().any().any() or not t2[list(A.T2M_REGIONS)].stack().between(200, 310).all():
        raise ValueError("implausible ARCO regional temperature")
    psl.to_parquet(OUT_PSL)
    t2.to_parquet(OUT_T2M)
    print(f"{len(psl)} days -> {OUT_PSL.name}; {len(t2)} days -> {OUT_T2M.name}")


if __name__ == "__main__":
    # --late: winters 2024/25 and 2025/26 for the 2026 held-out event
    # (s2s_heldout2026_test.py), written to separate files so the inputs of the
    # registered 2023-24 test are untouched
    # --early: 1940-1958 (November-May, and January-May 1940) for the revision-3
    # out-of-sample SSWs and the ERA5 continuity test (designs registered in
    # shift_rule_forecast.py and vortex_threshold_continuity.py, commit add7ba9)
    if "--early" in sys.argv:
        WINTERS = ([("1940-01-01", "1940-05-31")] + [(f"{y}-11-01", f"{y + 1}-05-31") for y in range(1940, 1958)]
                   + [("1958-11-01", "1958-12-31")])
        OUT_PSL = HERE / "era5_psl_cap_00utc_arco_early.parquet"
        OUT_T2M = HERE / "era5_t2m_regions_daily_arco_early.parquet"
    if "--late" in sys.argv:
        WINTERS = [("2024-11-01", "2025-04-30"), ("2025-11-01", "2026-04-30")]
        OUT_PSL = HERE / "era5_psl_cap_00utc_arco_late.parquet"
        OUT_T2M = HERE / "era5_t2m_regions_daily_arco_late.parquet"
    main()
