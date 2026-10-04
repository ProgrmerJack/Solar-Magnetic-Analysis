#!/usr/bin/env python3
"""
acquire_cmip6_monthly_boxes.py
==============================
Monthly 2 m temperature box means, 1984-12 .. 2014-12, CMIP6 historical vs DAMIP
hist-nat, for 07_physical_decomposition/damip_class_law_multimodel.py (registered
a280a10 before this script first ran).

SOURCES  Pangeo CMIP6 cloud store (anonymous; catalogue
  https://storage.googleapis.com/cmip6/pangeo-cmip6.csv) for CanESM5, MIROC6 and
  GISS-E2-1-G; ESGF CEDA OPeNDAP (anonymous) for HadGEM3-GC31-LL, Amon tas.
MEMBERS  every hist-nat member; historical members with the same physics and forcing
  variant (i*p*f* part of the member id) as the hist-nat members.
BOXES (cos-latitude means over cell centres inside; longitudes 0-360)
  CEUR 45-55N 5-15E; WRUS 50-60N 30-60E; CNAM 35-50N 260-280E; CEAS 45-60N 60-120E;
  NSIB 60-75N 60-140E; NWNA 55-70N 200-260E
Output: cmip6_monthly_boxes.parquet  model, arm, member, year, month, <boxes>
        per-member cache in _cmip6_monthly_boxes/
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
CACHE = HERE / "_cmip6_monthly_boxes"
OUT = HERE / "cmip6_monthly_boxes.parquet"
CATALOGUE = "https://storage.googleapis.com/cmip6/pangeo-cmip6.csv"
SEARCH = "https://esgf.ceda.ac.uk/esg-search/search"
BOXES = {"CEUR": (45, 55, 5, 15), "WRUS": (50, 60, 30, 60), "CNAM": (35, 50, 260, 280),
         "CEAS": (45, 60, 60, 120), "NSIB": (60, 75, 60, 140), "NWNA": (55, 70, 200, 260)}
PANGEO_MODELS = ["CanESM5", "MIROC6", "GISS-E2-1-G"]
T0, T1 = "1984-12-01", "2014-12-31"
socket.setdefaulttimeout(300)


def variant(m):
    return m[m.index("i"):]


def pangeo_jobs():
    c = pd.read_csv(CATALOGUE)
    c = c[(c.variable_id == "tas") & (c.table_id == "Amon") & c.experiment_id.isin(["historical", "hist-nat"])
          & c.source_id.isin(PANGEO_MODELS)]
    jobs = []
    for m in PANGEO_MODELS:
        cm = c[c.source_id == m]
        nat = cm[cm.experiment_id == "hist-nat"]
        vs = set(nat.member_id.map(variant))
        his = cm[(cm.experiment_id == "historical") & cm.member_id.map(variant).isin(vs)]
        for arm, df in (("hist-nat", nat), ("historical", his)):
            df = df.sort_values("version").groupby("member_id").tail(1)       # latest version per member
            jobs += [("pangeo", m, arm, r.member_id, r.zstore) for r in df.itertuples()]
    return jobs


def esgf_jobs():
    jobs = []
    for exp in ("historical", "hist-nat"):
        q = {"type": "File", "project": "CMIP6", "source_id": "HadGEM3-GC31-LL", "experiment_id": exp,
             "variable_id": "tas", "table_id": "Amon", "format": "application/solr+json", "limit": 2000,
             "fields": "title,url,member_id"}
        docs = json.load(urllib.request.urlopen(SEARCH + "?" + urllib.parse.urlencode(q), timeout=120))["response"]["docs"]
        best = {}
        for doc in docs:
            m = doc["member_id"][0] if isinstance(doc["member_id"], list) else doc["member_id"]
            if not m.endswith("i1p1f3"):
                continue
            span = doc["title"].rsplit("_", 1)[-1].replace(".nc", "")
            a, b = span.split("-")
            if not (a[:4] <= "1984" and b[:4] >= "2014"):
                continue
            for u in doc["url"]:
                if "OPENDAP" in u and "esgf.ceda.ac.uk" in u:
                    best[m] = u.split("|")[0].removesuffix(".html")
        jobs += [("esgf", "HadGEM3-GC31-LL", exp, m, u) for m, u in sorted(best.items())]
    return jobs


def reduce(job):
    src, model, arm, member, url = job
    import xarray as xr
    CACHE.mkdir(exist_ok=True)
    cf = CACHE / f"{model}_{arm}_{member}.parquet"
    if cf.exists():
        return str(cf)
    last = None
    for attempt in range(5):
        try:
            if src == "pangeo":
                import gcsfs
                ds = xr.open_zarr(gcsfs.GCSFileSystem(token="anon").get_mapper(url), consolidated=True, use_cftime=True)
            else:
                ds = xr.open_dataset(url, use_cftime=True)
            t = ds["tas"].sel(time=slice(T0, T1))
            if float(t.lon.min()) < 0:
                t = t.assign_coords(lon=np.mod(t.lon, 360)).sortby("lon")
            cols = {}
            for b, (la0, la1, lo0, lo1) in BOXES.items():
                sub = t.sel(lat=slice(la0, la1), lon=slice(lo0, lo1))
                if sub.sizes["lat"] < 1 or sub.sizes["lon"] < 1:
                    raise ValueError(f"box {b}: {dict(sub.sizes)}")
                cols[b] = sub.weighted(np.cos(np.deg2rad(sub.lat))).mean(("lat", "lon")).values
            tt = t["time"].values
            df = pd.DataFrame({"model": model, "arm": arm, "member": member, "year": [x.year for x in tt],
                               "month": [x.month for x in tt], **cols})
            if len(df) < 360 or not df[list(BOXES)].stack().between(180, 340).all():
                raise ValueError(f"short or implausible: {len(df)} months")
            df.to_parquet(cf)
            return str(cf)
        except Exception as e:                       # noqa: BLE001 -- retried, then reported
            last = e
            time.sleep(10 * (attempt + 1))
    return f"FAILED {model} {arm} {member}: {type(last).__name__}: {str(last)[:160]}"


def main():
    jobs = pangeo_jobs() + esgf_jobs()
    for m in PANGEO_MODELS + ["HadGEM3-GC31-LL"]:
        print(m, {a: sum(1 for j in jobs if j[1] == m and j[2] == a) for a in ("historical", "hist-nat")}, flush=True)
    import multiprocessing as mp
    done = 0
    with ProcessPoolExecutor(8, mp_context=mp.get_context("fork")) as ex:
        for r in ex.map(reduce, jobs):
            done += 1
            if r.startswith("FAILED"):
                print(r, flush=True)
            if done % 25 == 0:
                print(f"  {done}/{len(jobs)}", flush=True)
    d = pd.concat([pd.read_parquet(p) for p in sorted(CACHE.glob("*.parquet"))], ignore_index=True)
    d.to_parquet(OUT)
    print(d.groupby(["model", "arm"]).member.nunique().to_string())
    print(f"Saved -> {OUT.name} ({len(d):,} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
