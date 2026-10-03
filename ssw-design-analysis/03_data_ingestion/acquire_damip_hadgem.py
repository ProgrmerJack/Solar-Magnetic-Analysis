#!/usr/bin/env python3
"""
acquire_damip_hadgem.py
=======================
Daily 2 m temperature box means from CMIP6 HadGEM3-GC31-LL `historical` (55
members) and DAMIP `hist-nat` (60 members), 1984-12-01 .. 2014-12-30, for
07_physical_decomposition/damip_class_law.py (registered 60f3d8b before this script
first ran). Recipe and access test: research note cross_domain_interventions.md.

SOURCE  ESGF search API at esgf.ceda.ac.uk (type=File, frequency=day, variable_id=tas),
  OPeNDAP URLs on esgf.ceda.ac.uk only, files 19500101-2014/2020; anonymous,
  CC BY 4.0 (CMIP6 terms). Only the four boxes are requested from the server.
BOXES (cos-latitude-weighted means over grid cells whose centres lie inside)
  CEUR 45-55N, 5-15E; WRUS 50-60N, 30-60E; CNAM 35-50N, 260-280E; CEAS 45-60N, 60-120E
CALENDAR 360-day; year, month and day are taken from the model's own time axis.
PARALLELISM processes, not threads (NetCDF/HDF5 is not thread-safe: FAILURES
  2026-10-03); a socket timeout; the Python exit status is recorded by the caller.
Output: damip_hadgem_boxes.parquet  arm, member, year, month, day, CEUR, WRUS, CNAM, CEAS
        (per-member cache in _damip_hadgem/)
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
CACHE = HERE / "_damip_hadgem"
OUT = HERE / "damip_hadgem_boxes.parquet"
SEARCH = "https://esgf.ceda.ac.uk/esg-search/search"
BOXES = {"CEUR": (45, 55, 5, 15), "WRUS": (50, 60, 30, 60), "CNAM": (35, 50, 260, 280), "CEAS": (45, 60, 60, 120)}
ARMS = {"historical": "CMIP", "hist-nat": "DAMIP"}
T0, T1 = "1984-12-01", "2014-12-30"
socket.setdefaulttimeout(300)


def files(exp):
    q = {"type": "File", "project": "CMIP6", "source_id": "HadGEM3-GC31-LL", "experiment_id": exp,
         "variable_id": "tas", "frequency": "day", "format": "application/solr+json", "limit": 1000,
         "fields": "title,url,member_id"}
    d = json.load(urllib.request.urlopen(SEARCH + "?" + urllib.parse.urlencode(q), timeout=120))["response"]["docs"]
    out = {}
    for doc in d:
        m = doc["member_id"][0] if isinstance(doc["member_id"], list) else doc["member_id"]
        if not m.endswith("i1p1f3") or "_19500101-" not in doc["title"]:
            continue
        for u in doc["url"]:
            if "OPENDAP" in u and "esgf.ceda.ac.uk" in u:
                out[m] = u.split("|")[0].removesuffix(".html")
    return out


def reduce_member(job):
    exp, member, url = job
    import xarray as xr
    CACHE.mkdir(exist_ok=True)
    cf = CACHE / f"{exp}_{member}.parquet"
    if cf.exists():
        return str(cf)
    last = None
    for attempt in range(5):
        try:
            with xr.open_dataset(url, decode_times=True, use_cftime=True) as ds:
                t = ds["tas"].sel(time=slice(T0, T1))
                cols = {}
                for b, (la0, la1, lo0, lo1) in BOXES.items():
                    sub = t.sel(lat=slice(la0, la1), lon=slice(lo0, lo1))
                    if sub.sizes["lat"] < 2 or sub.sizes["lon"] < 2:
                        raise ValueError(f"box {b} has {sub.sizes} cells")
                    w = np.cos(np.deg2rad(sub.lat))
                    cols[b] = sub.weighted(w).mean(("lat", "lon")).values
                tt = t["time"].values
            df = pd.DataFrame({"arm": exp, "member": member, "year": [x.year for x in tt],
                               "month": [x.month for x in tt], "day": [x.day for x in tt], **cols})
            if not (df[list(BOXES)].stack().between(180, 340).all()) or len(df) < 10500:
                raise ValueError(f"implausible or short: {len(df)} days")
            df.to_parquet(cf)
            return str(cf)
        except Exception as e:                       # noqa: BLE001 -- retried, then reported
            last = e
            time.sleep(15 * (attempt + 1))
    return f"FAILED {exp} {member}: {type(last).__name__}: {str(last)[:160]}"


def main():
    CACHE.mkdir(exist_ok=True)
    jobs = []
    for exp in ARMS:
        f = files(exp)
        print(f"{exp}: {len(f)} members with an i1p1f3 1950- file on CEDA", flush=True)
        jobs += [(exp, m, u) for m, u in sorted(f.items())]
    done = 0
    import multiprocessing as mp           # 'fork': the parent never initialises NetCDF; forkserver
    with ProcessPoolExecutor(6, mp_context=mp.get_context("fork")) as ex:   # failed when detached
        for r in ex.map(reduce_member, jobs):
            done += 1
            if r.startswith("FAILED"):
                print(r, flush=True)
            if done % 10 == 0:
                print(f"  {done}/{len(jobs)}", flush=True)
    parts = [pd.read_parquet(p) for p in sorted(CACHE.glob("*.parquet"))]
    d = pd.concat(parts, ignore_index=True)
    d.to_parquet(OUT)
    print(d.groupby("arm").member.nunique().to_string())
    print(f"Saved -> {OUT.name} ({len(d):,} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
