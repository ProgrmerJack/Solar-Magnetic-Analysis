#!/usr/bin/env python3
"""
acquire_enso_regime_fields.py
=============================
Monthly surface temperature over the tropical band (20S-20N, all longitudes), box-averaged onto a
common 2.5-degree grid, for the registered test 07_physical_decomposition/enso_regime_relative.py
(committed with it, before this script first ran). Over ocean cells ts is SST on the atmospheric grid.

SOURCES
  CESM1-LE   AWS s3://ncar-cesm-lens/atm/monthly/cesmLE-{20C,RCP85}-TS.zarr (anonymous), 40 members
  GFDL-CM3   MMLEA v2 (NCAR GDEX d651039) ts_Amon_GFDL-CM3_historical_rcp85_r*i1p1_g025_192001-210012.nc
  CMIP6      Pangeo store (as acquire_cmip6_nino3.py): MIROC6, CanESM5, ACCESS-ESM1-5, ts, historical + ssp585,
             members present in both
  MASK       NOAA OISST v2 0.25-degree land-sea mask (lsmask.oisst.nc, PSL), ocean fraction per 2.5-degree box
WINDOWS  1950-12..1980-11, 2020-12..2050-11, 2069-12..2099-11 (360 months each; RCP8.5 for the CMIP5-era
  ensembles, SSP5-8.5 for CMIP6)
GRID  centres at multiples of 2.5 degrees (lat -20..20, lon 0..357.5); a native cell belongs to the box
  that contains its centre (box edges +-1.25); cos-latitude weights. Plausibility: 270-315 K over the
  ocean cells of the band.
Output: _enso_fields/<model>__<member>.npz (ts[1080, 17, 144], year, month), _enso_fields/ocean_frac.npy
"""
import socket
import sys
import time
import urllib.request
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
OUTD = HERE / "_enso_fields"
RAW = HERE / "raw" / "gfdl_cm3_ts"
WINDOWS = [("1950-12-01", "1980-11-30"), ("2020-12-01", "2050-11-30"), ("2069-12-01", "2099-11-30")]
LAT = np.arange(-20.0, 20.01, 2.5); LON = np.arange(0.0, 360.0, 2.5)
CMIP6 = ("MIROC6", "CanESM5", "ACCESS-ESM1-5")
MASK_URL = "https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/lsmask.oisst.nc"
GFDL_URL = "https://osdf-director.osg-htc.org/ncar/gdex/d651039/gfdl_cm3_lens/Amon/ts/ts_Amon_GFDL-CM3_historical_rcp85_r{r}i1p1_g025_192001-210012.nc"
socket.setdefaulttimeout(600)


def to_grid(da):
    """box-average a (time, lat, lon) DataArray with 0..360 longitudes onto LAT x LON."""
    lat = da.lat.values; lon = np.mod(da.lon.values, 360.0)
    il = np.round(lat / 2.5).astype(int); jl = np.mod(np.round(lon / 2.5).astype(int), 144)
    v = da.values.astype(np.float64); w = np.cos(np.deg2rad(lat))
    out = np.zeros((v.shape[0], LAT.size, LON.size)); wsum = np.zeros((LAT.size, LON.size))
    for a, i in enumerate(il):
        if not (-8 <= i <= 8):
            continue
        for b, j in enumerate(jl):
            out[:, i + 8, j] += w[a] * v[:, a, b]; wsum[i + 8, j] += w[a]
    return (out / np.where(wsum > 0, wsum, np.nan)).astype(np.float32)


def window_months(t):
    yrs = np.array([x.year for x in t]); mon = np.array([x.month for x in t])
    keep = np.zeros(len(t), bool)
    for a, b in WINDOWS:
        ya, ma = int(a[:4]), int(a[5:7]); yb, mb = int(b[:4]), int(b[5:7])
        k = yrs * 12 + mon; keep |= (k >= ya * 12 + ma) & (k <= yb * 12 + mb)
    return keep, yrs, mon


def save(model, member, da):
    sub = da.sel(lat=slice(-22, 22))
    keep, yrs, mon = window_months(sub.time.values)
    sub = sub.isel(time=np.where(keep)[0]).load()
    g = to_grid(sub)
    if g.shape[0] != 1080:
        raise ValueError(f"{model} {member}: {g.shape[0]} months, expected 1080")
    of = np.load(OUTD / "ocean_frac.npy")
    oc = g[:, of >= 0.9]
    if not (np.nanmin(oc) > 270 and np.nanmax(oc) < 315):
        raise ValueError(f"{model} {member}: implausible ocean ts {np.nanmin(oc):.1f}-{np.nanmax(oc):.1f}")
    np.savez_compressed(OUTD / f"{model}__{member}.npz", ts=g, year=yrs[keep], month=mon[keep])


def mask():
    import xarray as xr
    f = OUTD / "lsmask.oisst.nc"
    if not f.exists():
        urllib.request.urlretrieve(MASK_URL, f)
    m = xr.open_dataset(f).lsmask.squeeze()
    lat = m.lat.values; lon = m.lon.values; v = m.values.astype(float)
    il = np.round(lat / 2.5).astype(int); jl = np.mod(np.round(lon / 2.5).astype(int), 144)
    s = np.zeros((LAT.size, LON.size)); n = np.zeros((LAT.size, LON.size))
    for a, i in enumerate(il):
        if -8 <= i <= 8:
            np.add.at(s[i + 8], jl, v[a]); np.add.at(n[i + 8], jl, 1)
    np.save(OUTD / "ocean_frac.npy", s / n)


def cesm(member):
    import s3fs
    import xarray as xr
    if (OUTD / f"CESM1-LE__{member}.npz").exists():
        return "ok"
    fs = s3fs.S3FileSystem(anon=True)
    parts = []
    for e in ("20C", "RCP85"):
        ds = xr.open_zarr(fs.get_mapper(f"ncar-cesm-lens/atm/monthly/cesmLE-{e}-TS.zarr"), consolidated=True)
        parts.append(ds.TS.sel(member_id=member))
    da = xr.concat(parts, "time")
    save("CESM1-LE", str(member), da)
    return "ok"


def gfdl(r):
    import xarray as xr
    if (OUTD / f"GFDL-CM3__r{r}i1p1.npz").exists():
        return "ok"
    RAW.mkdir(parents=True, exist_ok=True)
    f = RAW / f"r{r}.nc"
    if not f.exists():
        urllib.request.urlretrieve(GFDL_URL.format(r=r), f.with_suffix(".part")); f.with_suffix(".part").rename(f)
    da = xr.open_dataset(f, use_cftime=True).ts
    save("GFDL-CM3", f"r{r}i1p1", da)
    return "ok"


def cmip6(job):
    import gcsfs
    import xarray as xr
    model, member, urls = job
    if (OUTD / f"{model}__{member}.npz").exists():
        return "ok"
    fs = gcsfs.GCSFileSystem(token="anon")
    parts = [xr.open_zarr(fs.get_mapper(u), consolidated=True, use_cftime=True).ts for u in urls]
    da = xr.concat(parts, "time")
    if da.lat.values[0] > da.lat.values[-1]:
        da = da.sortby("lat")
    save(model, member, da)
    return "ok"


def cmip6_jobs():
    from acquire_cmip6_monthly_boxes import CATALOGUE
    c = pd.read_csv(CATALOGUE)
    c = c[(c.table_id == "Amon") & (c.variable_id == "ts") & c.experiment_id.isin(["historical", "ssp585"]) & c.source_id.isin(CMIP6)]
    out = []
    for m in CMIP6:
        cm = c[c.source_id == m]
        common = sorted(set(cm[cm.experiment_id == "historical"].member_id) & set(cm[cm.experiment_id == "ssp585"].member_id))
        for mem in common:
            urls = [cm[(cm.experiment_id == e) & (cm.member_id == mem)].sort_values("version").zstore.iloc[-1]
                    for e in ("historical", "ssp585")]
            out.append((m, mem, urls))
    return out


def run(fn, arg):
    for attempt in range(4):
        try:
            return fn(arg)
        except Exception as e:                       # noqa: BLE001 -- retried, then reported
            last = e; time.sleep(20 * (attempt + 1))
    return f"FAILED {fn.__name__} {arg if not isinstance(arg, tuple) else arg[:2]}: {type(last).__name__}: {str(last)[:200]}"


def _job(a):
    return run(*a)


def main():
    OUTD.mkdir(exist_ok=True)
    if not (OUTD / "ocean_frac.npy").exists():
        mask()
    import s3fs
    import xarray as xr
    ds = xr.open_zarr(s3fs.S3FileSystem(anon=True).get_mapper("ncar-cesm-lens/atm/monthly/cesmLE-RCP85-TS.zarr"), consolidated=True)
    jobs = [(cesm, int(m)) for m in ds.member_id.values] + [(gfdl, r) for r in range(1, 21)] + [(cmip6, j) for j in cmip6_jobs()]
    print(f"{len(jobs)} members", flush=True)
    if "--list" in sys.argv:
        return 0
    import multiprocessing as mp
    with ProcessPoolExecutor(4, mp_context=mp.get_context("spawn"), max_tasks_per_child=5) as ex:
        for i, r in enumerate(ex.map(_job, jobs)):
            if r != "ok":
                print(r, flush=True)
            if (i + 1) % 10 == 0:
                print(f"  {i + 1}/{len(jobs)}", flush=True)
    print("files:", len(list(OUTD.glob("*__*.npz"))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
