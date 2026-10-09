#!/usr/bin/env python3
"""
acquire_cmip6_nino3.py
======================
Monthly Nino3 (5S-5N, 210-270E) area means of precipitation (pr, mm/day) and surface
temperature (ts, K; all-ocean box, so SST on the atmospheric grid) for the registered test
07_physical_decomposition/enso_extreme_class.py (committed with it, before this script
first ran). Pangeo CMIP6 store, latest version per member, members present in both
historical and ssp585.
MODELS  pr: CanESM5, MIROC6, ACCESS-ESM1-5, EC-Earth3; ts: CanESM5, MIROC6, ACCESS-ESM1-5
PERIODS historical 1950-12 .. 1980-02; ssp585 2020-12 .. 2099-02
Same reduction as acquire_cmip6_general.py (cos-latitude weights of cell centres inside
the box); plausibility checks pr in [0, 60] mm/day, ts in [285, 315] K.
Output: cmip6_nino3.parquet  model, experiment, var, member, year, month, NINO3
        per-member cache in _cmip6_nino3/
"""
import socket
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from acquire_cmip6_monthly_boxes import CATALOGUE   # noqa: E402

CACHE = HERE / "_cmip6_nino3"
OUT = HERE / "cmip6_nino3.parquet"
BOX = (-5.0, 5.0, 210.0, 270.0)
PERIODS = {"historical": ("1950-12-01", "1980-02-28"), "ssp585": ("2020-12-01", "2099-02-28")}
MODELS = {"pr": ("CanESM5", "MIROC6", "ACCESS-ESM1-5", "EC-Earth3"), "ts": ("CanESM5", "MIROC6", "ACCESS-ESM1-5")}
socket.setdefaulttimeout(300)


def jobs():
    c = pd.read_csv(CATALOGUE)
    c = c[(c.table_id == "Amon") & c.variable_id.isin(["pr", "ts"]) & c.experiment_id.isin(PERIODS)]
    out = []
    for var, models in MODELS.items():
        for m in models:
            cm = c[(c.variable_id == var) & (c.source_id == m)]
            common = set(cm[cm.experiment_id == "historical"].member_id) & set(cm[cm.experiment_id == "ssp585"].member_id)
            for e in PERIODS:
                df = cm[(cm.experiment_id == e) & cm.member_id.isin(common)].sort_values("version").groupby("member_id").tail(1)
                out += [(m, e, var, r.member_id, r.zstore) for r in df.itertuples()]
    return out


def reduce(job):
    model, exp, var, member, url = job
    import dask
    import gcsfs
    import xarray as xr
    dask.config.set(scheduler="synchronous")
    CACHE.mkdir(exist_ok=True)
    cf = CACHE / f"{model}_{exp}_{var}_{member}.parquet"
    if cf.exists():
        return str(cf)
    last = None
    for attempt in range(5):
        try:
            ds = xr.open_zarr(gcsfs.GCSFileSystem(token="anon").get_mapper(url), consolidated=True, use_cftime=True)
            t = ds[var].sel(time=slice(*PERIODS[exp]))
            if t.lat.values[0] > t.lat.values[-1]:
                t = t.sortby("lat")
            if float(t.lon.min()) < 0:
                t = t.assign_coords(lon=np.mod(t.lon, 360)).sortby("lon")
            la0, la1, lo0, lo1 = BOX
            sub = t.sel(lat=slice(la0, la1), lon=slice(lo0, lo1)).load()
            if sub.sizes["lat"] < 1 or sub.sizes["lon"] < 1:
                raise ValueError(f"box empty: {dict(sub.sizes)}")
            v = sub.weighted(np.cos(np.deg2rad(sub.lat))).mean(("lat", "lon")).values
            v = v * 86400.0 if var == "pr" else v
            tt = sub["time"].values
            df = pd.DataFrame({"model": model, "experiment": exp, "var": var, "member": member,
                               "year": [x.year for x in tt], "month": [x.month for x in tt], "NINO3": v})
            lo, hi = (0.0, 60.0) if var == "pr" else (285.0, 315.0)
            if len(df) < 300 or not df.NINO3.between(lo, hi).all():
                raise ValueError(f"short or implausible: {len(df)} months, range {df.NINO3.min():.2f}-{df.NINO3.max():.2f}")
            df.to_parquet(cf)
            return str(cf)
        except Exception as e:                       # noqa: BLE001 -- retried, then reported
            last = e
            time.sleep(10 * (attempt + 1))
    return f"FAILED {model} {exp} {var} {member}: {type(last).__name__}: {str(last)[:160]}"


def main():
    js = jobs()
    import collections
    for k, v in sorted(collections.Counter((j[0], j[1], j[2]) for j in js).items()):
        print(k, v, flush=True)
    if "--list" in sys.argv:
        return 0
    import multiprocessing as mp
    done = 0
    with ProcessPoolExecutor(6, mp_context=mp.get_context("spawn"), max_tasks_per_child=20) as ex:
        for r in ex.map(reduce, js):
            done += 1
            if r.startswith("FAILED"):
                print(r, flush=True)
            if done % 50 == 0:
                print(f"  {done}/{len(js)}", flush=True)
    d = pd.concat([pd.read_parquet(p) for p in sorted(CACHE.glob("*.parquet"))], ignore_index=True)
    d.to_parquet(OUT)
    print(d.groupby(["model", "experiment", "var"]).member.nunique().to_string())
    print(f"Saved -> {OUT.name} ({len(d):,} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
