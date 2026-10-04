#!/usr/bin/env python3
"""
acquire_cmip6_general.py
========================
Monthly box means of 2 m temperature (tas, K) and precipitation (pr, mm/day) for
07_physical_decomposition/class_law_general.py (registered before this script first
ran). Reuses the box definitions and reduction of acquire_cmip6_monthly_boxes.py;
writes a SEPARATE file so the registered D-MM dataset is untouched.

SETS (Pangeo; latest version per member)
  CanESM5 (p1, p2), MIROC6 : hist-nat 1984-12..2014-12, historical 1950-12..2014-12 and
                             ssp585 2070-12..2100-12, members with a hist-nat variant
  EC-Earth3, ACCESS-ESM1-5 : historical 1950-12..2014-12 and ssp585 2070-12..2100-12,
                             members present in both experiments
  HadGEM3-GC31-LL          : pr, historical and hist-nat 1984-12..2014-12 (ESGF CEDA);
                             its tas is in cmip6_monthly_boxes.parquet
Output: cmip6_general_boxes.parquet  model, experiment, var, member, year, month, <boxes>
        per-member cache in _cmip6_general/
"""
import json
import socket
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from acquire_cmip6_monthly_boxes import BOXES, CATALOGUE, SEARCH, variant   # noqa: E402

CACHE = HERE / "_cmip6_general"
OUT = HERE / "cmip6_general_boxes.parquet"
PERIODS = {"hist-nat": ("1984-12-01", "2014-12-30"), "historical": ("1950-12-01", "2014-12-30"),
           "ssp585": ("2070-12-01", "2100-12-30")}
socket.setdefaulttimeout(300)


def pangeo_jobs():
    c = pd.read_csv(CATALOGUE)
    c = c[(c.table_id == "Amon") & c.variable_id.isin(["tas", "pr"])]
    jobs = []
    def latest(df):
        return df.sort_values("version").groupby("member_id").tail(1)
    for var in ("tas", "pr"):
        cv = c[c.variable_id == var]
        for m in ("CanESM5", "MIROC6"):
            cm = cv[cv.source_id == m]
            nat = cm[cm.experiment_id == "hist-nat"]; vs = set(nat.member_id.map(variant))
            for e in ("hist-nat", "historical", "ssp585"):
                df = latest(cm[(cm.experiment_id == e) & cm.member_id.map(variant).isin(vs)])
                jobs += [("pangeo", m, e, var, r.member_id, r.zstore) for r in df.itertuples()]
        for m in ("EC-Earth3", "ACCESS-ESM1-5"):
            cm = cv[cv.source_id == m]
            common = set(cm[cm.experiment_id == "historical"].member_id) & set(cm[cm.experiment_id == "ssp585"].member_id)
            for e in ("historical", "ssp585"):
                df = latest(cm[(cm.experiment_id == e) & cm.member_id.isin(common)])
                jobs += [("pangeo", m, e, var, r.member_id, r.zstore) for r in df.itertuples()]
    return jobs


def esgf_jobs():
    jobs = []
    for exp in ("historical", "hist-nat"):
        q = {"type": "File", "project": "CMIP6", "source_id": "HadGEM3-GC31-LL", "experiment_id": exp,
             "variable_id": "pr", "table_id": "Amon", "format": "application/solr+json", "limit": 2000,
             "fields": "title,url,member_id", "distrib": "false"}
        docs = json.load(urllib.request.urlopen(SEARCH + "?" + urllib.parse.urlencode(q), timeout=120))["response"]["docs"]
        best = {}
        for doc in docs:
            m = doc["member_id"][0] if isinstance(doc["member_id"], list) else doc["member_id"]
            if not m.endswith("i1p1f3"):
                continue
            a, b = doc["title"].rsplit("_", 1)[-1].replace(".nc", "").split("-")
            if a[:4] <= "1984" and b[:4] >= "2014":
                for u in doc["url"]:
                    if "OPENDAP" in u and "esgf.ceda.ac.uk" in u:
                        best[m] = u.split("|")[0].removesuffix(".html")
        jobs += [("esgf", "HadGEM3-GC31-LL", exp, "pr", m, u) for m, u in sorted(best.items())]
    return jobs


def reduce(job):
    src, model, exp, var, member, url = job
    import dask
    import xarray as xr
    dask.config.set(scheduler="synchronous")           # one chunk in memory at a time per worker
    CACHE.mkdir(exist_ok=True)
    cf = CACHE / f"{model}_{exp}_{var}_{member}.parquet"
    if cf.exists():
        return str(cf)
    t0, t1 = ("1984-12-01", "2014-12-30") if src == "esgf" else PERIODS[exp]
    last = None
    for attempt in range(5):
        try:
            if src == "pangeo":
                import gcsfs
                ds = xr.open_zarr(gcsfs.GCSFileSystem(token="anon").get_mapper(url), consolidated=True, use_cftime=True)
            else:
                ds = xr.open_dataset(url, use_cftime=True)
            t = ds[var].sel(time=slice(t0, t1))
            if t.lat.values[0] > t.lat.values[-1]:
                t = t.sortby("lat")
            t = t.sel(lat=slice(30, 80)).load()           # one read per member; every box lies in 35-75N
            if float(t.lon.min()) < 0:
                t = t.assign_coords(lon=np.mod(t.lon, 360)).sortby("lon")
            cols = {}
            for b, (la0, la1, lo0, lo1) in BOXES.items():
                sub = t.sel(lat=slice(la0, la1), lon=slice(lo0, lo1))
                if sub.sizes["lat"] < 1 or sub.sizes["lon"] < 1:
                    raise ValueError(f"box {b}: {dict(sub.sizes)}")
                v = sub.weighted(np.cos(np.deg2rad(sub.lat))).mean(("lat", "lon")).values
                cols[b] = v * 86400.0 if var == "pr" else v
            tt = t["time"].values
            df = pd.DataFrame({"model": model, "experiment": exp, "var": var, "member": member,
                               "year": [x.year for x in tt], "month": [x.month for x in tt], **cols})
            ok = df[list(BOXES)].stack()
            if len(df) < 360 or not (ok.between(180, 340).all() if var == "tas" else ok.between(0, 60).all()):
                raise ValueError(f"short or implausible: {len(df)} months")
            df.to_parquet(cf)
            return str(cf)
        except Exception as e:                       # noqa: BLE001 -- retried, then reported
            last = e
            time.sleep(10 * (attempt + 1))
    return f"FAILED {model} {exp} {var} {member}: {type(last).__name__}: {str(last)[:160]}"


def main():
    jobs = esgf_jobs() if "--only-esgf" in sys.argv else (pangeo_jobs() if "--only-pangeo" in sys.argv
                                                           else pangeo_jobs() + esgf_jobs())
    import collections
    for k, v in sorted(collections.Counter((j[1], j[2], j[3]) for j in jobs).items()):
        print(k, v, flush=True)
    import multiprocessing as mp
    done = 0
    # spawn: gcsfs (FAILURES 2026-10-04); 6 recycled workers: with 10 the workers vanished on EC-Earth3 (2026-10-05)
    with ProcessPoolExecutor(6, mp_context=mp.get_context("spawn"), max_tasks_per_child=20) as ex:
        for r in ex.map(reduce, jobs):
            done += 1
            if r.startswith("FAILED"):
                print(r, flush=True)
            if done % 50 == 0:
                print(f"  {done}/{len(jobs)}", flush=True)
    d = pd.concat([pd.read_parquet(p) for p in sorted(CACHE.glob("*.parquet"))], ignore_index=True)
    d.to_parquet(OUT)
    print(d.groupby(["model", "experiment", "var"]).member.nunique().to_string())
    print(f"Saved -> {OUT.name} ({len(d):,} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
