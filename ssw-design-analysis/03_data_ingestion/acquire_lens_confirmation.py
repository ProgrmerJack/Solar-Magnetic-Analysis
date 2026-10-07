#!/usr/bin/env python3
"""
acquire_lens_confirmation.py
============================
Monthly box means of 2 m temperature (K) and precipitation (mm/day) from two large
ensembles never analysed in this project for surface outcomes, for the registered
confirmation 07_physical_decomposition/lens_confirmation.py (committed before this
script first fetched any value).

SOURCES (anonymous)
  CESM2-LE  AWS s3://ncar-cesm2-lens (catalogue github.com/NCAR/cesm2-le-aws), monthly
            TREFHT and PRECC + PRECL; historical and SSP3-7.0; forcing variants "cmip6"
            (50 members, f1) and "smbb" (50 members, f2); a member keeps its id across
            historical and SSP3-7.0
  MPI-GE    DKRZ ESGF (project MPI-GE), Amon tas and pr; historical r###i1850p3 and
            RCP8.5 r###i2005p3 (run ### continues), 100 members; OPeNDAP subsets
PERIODS   historical 1950-12 .. 1980-12; future 2039-12 .. 2099-12
BOXES     as acquire_cmip6_monthly_boxes.py (cos-latitude means of cell centres inside)
GATE (--gate, run before the MPI-GE download, as FAILURES/ssw-data-access-traps require):
  the OPeNDAP subset of member 001 tas (historical and RCP8.5) must equal the values
  of the full HTTP file (max |diff| < 1e-4 K) or nothing is downloaded.
Output: lens_confirmation_boxes.parquet  model, experiment, var, member, year, month, <boxes>
        per-member cache in _lens_confirmation/
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
from acquire_cmip6_monthly_boxes import BOXES   # noqa: E402

CACHE = HERE / "_lens_confirmation"
OUT = HERE / "lens_confirmation_boxes.parquet"
PERIODS = {"historical": ("1950-12-01", "1980-12-31"), "future": ("2039-12-01", "2099-12-31")}
S3 = "s3://ncar-cesm2-lens/atm/monthly/cesm2LE-{exp}-{variant}-{var}.zarr"
SEARCH = "https://esgf-data.dkrz.de/esg-search/search"
socket.setdefaulttimeout(300)


def boxes(t, scale, kind):
    if t.lat.values[0] > t.lat.values[-1]:
        t = t.sortby("lat")
    t = t.sel(lat=slice(30, 80)).load()
    if float(t.lon.min()) < 0:
        t = t.assign_coords(lon=np.mod(t.lon, 360)).sortby("lon")
    cols = {}
    for b, (la0, la1, lo0, lo1) in BOXES.items():
        sub = t.sel(lat=slice(la0, la1), lon=slice(lo0, lo1))
        if sub.sizes["lat"] < 1 or sub.sizes["lon"] < 1:
            raise ValueError(f"box {b}")
        cols[b] = sub.weighted(np.cos(np.deg2rad(sub.lat))).mean(("lat", "lon")).values * scale
    v = pd.DataFrame(cols).stack()
    ok = v.between(180, 340).all() if kind == "tas" else v.between(0, 60).all()
    if not ok or not np.isfinite(v).all():
        raise ValueError(f"implausible {kind}: {float(v.min())}..{float(v.max())}")
    return cols


def cesm_job(job):
    variant, member = job
    import fsspec
    import xarray as xr
    CACHE.mkdir(exist_ok=True)
    cf = CACHE / f"CESM2-{variant}_{member}.parquet"
    if cf.exists():
        return str(cf)
    last = None
    for attempt in range(4):
        try:
            frames = []
            for exp, lab in (("historical", "historical"), ("ssp370", "future")):
                t0, t1 = PERIODS[lab]
                op = lambda v: xr.open_zarr(fsspec.get_mapper(S3.format(exp=exp, variant=variant, var=v), anon=True),
                                            consolidated=True)[v].sel(member_id=member, time=slice(t0, t1))
                for var, parts, scale in (("tas", ["TREFHT"], 1.0), ("pr", ["PRECC", "PRECL"], 86400.0 * 1000.0)):
                    x = sum(op(p) for p in parts)
                    cols = boxes(x, scale, var)
                    yrs = x["time"].dt.year.values; mos = x["time"].dt.month.values      # mid-month stamps (catalogue)
                    frames.append(pd.DataFrame({"model": f"CESM2-{variant}", "experiment": lab, "var": var, "member": member,
                                                "year": yrs, "month": mos, **cols}))
            df = pd.concat(frames, ignore_index=True)
            df.to_parquet(cf)
            return str(cf)
        except Exception as e:                       # noqa: BLE001 -- retried, then reported
            last = e
            time.sleep(15 * (attempt + 1))
    return f"FAILED CESM2-{variant} {member}: {type(last).__name__}: {str(last)[:160]}"


def mpi_urls():
    out = {}
    for var in ("tas", "pr"):
        for exp in ("historical", "rcp85"):
            q = {"project": "MPI-GE", "type": "File", "variable": var, "experiment": exp, "time_frequency": "mon",
                 "cmor_table": "Amon", "distrib": "false", "format": "application/solr+json", "limit": 500,
                 "fields": "title,url,ensemble"}
            docs = json.load(urllib.request.urlopen(SEARCH + "?" + urllib.parse.urlencode(q), timeout=300))["response"]["docs"]
            for d in docs:
                if "_Amon_" not in d["title"]:
                    continue
                ens = d["ensemble"][0] if isinstance(d["ensemble"], list) else d["ensemble"]
                run = ens.split("i")[0]
                urls = {u.split("|")[2]: u.split("|")[0] for u in d["url"]}
                out[(var, exp, run)] = (urls["OPENDAP"].removesuffix(".html"), urls["HTTPServer"])
    return out


def mpi_job(job):
    run, urls = job
    import xarray as xr
    CACHE.mkdir(exist_ok=True)
    cf = CACHE / f"MPI-GE_{run}.parquet"
    if cf.exists():
        return str(cf)
    last = None
    for attempt in range(4):
        try:
            frames = []
            for exp, lab in (("historical", "historical"), ("rcp85", "future")):
                t0, t1 = PERIODS[lab]
                for var, scale in (("tas", 1.0), ("pr", 86400.0)):
                    ds = xr.open_dataset(urls[(var, exp)][0], use_cftime=True)
                    x = ds[var].sel(time=slice(t0, t1))
                    cols = boxes(x, scale, var)
                    frames.append(pd.DataFrame({"model": "MPI-GE", "experiment": lab, "var": var, "member": run,
                                                "year": x["time"].dt.year.values, "month": x["time"].dt.month.values, **cols}))
            pd.concat(frames, ignore_index=True).to_parquet(cf)
            return str(cf)
        except Exception as e:                       # noqa: BLE001
            last = e
            time.sleep(15 * (attempt + 1))
    return f"FAILED MPI-GE {run}: {type(last).__name__}: {str(last)[:160]}"


def gate(urls):
    import xarray as xr
    res = {}
    for exp in ("historical", "rcp85"):
        dap, http = urls[("tas", exp, "r001")]
        local = CACHE / f"_gate_{exp}.nc"
        CACHE.mkdir(exist_ok=True)
        if not local.exists():
            urllib.request.urlretrieve(http, local)
        t0, t1 = PERIODS["historical" if exp == "historical" else "future"]
        a = xr.open_dataset(dap, use_cftime=True)["tas"].sel(time=slice(t0, t1)).sortby("lat").sel(lat=slice(30, 80)).load()
        b = xr.open_dataset(local, use_cftime=True)["tas"].sel(time=slice(t0, t1)).sortby("lat").sel(lat=slice(30, 80)).load()
        d = float(np.abs(a.values - b.values).max())
        res[exp] = {"max_abs_diff": d, "n": int(a.size), "mean_dap": float(a.mean()), "mean_http": float(b.mean()),
                    "zeros_dap": int((a.values == 0).sum())}
        local.unlink()
    res["passed"] = all(v["max_abs_diff"] < 1e-4 and v["zeros_dap"] == 0 for v in res.values() if isinstance(v, dict))
    print(json.dumps(res, indent=1))
    return res


def main():
    import multiprocessing as mp
    urls_raw = mpi_urls() if ("--gate" in sys.argv or "--only-cesm" not in sys.argv) else {}
    if "--gate" in sys.argv:
        r = gate(urls_raw)
        (CACHE / "gate.json").write_text(json.dumps(r, indent=1))
        return 0 if r["passed"] else 1
    jobs = []
    if "--only-mpi" not in sys.argv:
        jobs += [("cesm", (v, m)) for v in ("cmip6", "smbb") for m in _cesm_members(v)]
    if "--only-cesm" not in sys.argv:
        g = CACHE / "gate.json"
        if not g.exists() or not json.loads(g.read_text())["passed"]:
            sys.exit("MPI-GE: run --gate first (OPeNDAP must match the full file)")
        runs = sorted({k[2] for k in urls_raw})
        for run in runs:
            u = {(v, e): urls_raw[(v, e, run)] for v in ("tas", "pr") for e in ("historical", "rcp85") if (v, e, run) in urls_raw}
            if len(u) == 4:
                jobs.append(("mpi", (run, u)))
    print(f"{len(jobs)} jobs", flush=True)
    done = 0
    with ProcessPoolExecutor(6, mp_context=mp.get_context("spawn"), max_tasks_per_child=10) as ex:
        futs = [ex.submit(cesm_job if k == "cesm" else mpi_job, a) for k, a in jobs]
        for f in futs:
            r = f.result(); done += 1
            if r.startswith("FAILED"):
                print(r, flush=True)
            if done % 20 == 0:
                print(f"  {done}/{len(jobs)}", flush=True)
    d = pd.concat([pd.read_parquet(p) for p in sorted(CACHE.glob("*.parquet"))], ignore_index=True)
    d.to_parquet(OUT)
    print(d.groupby(["model", "experiment", "var"]).member.nunique().to_string())
    print(f"Saved -> {OUT.name} ({len(d):,} rows)")
    return 0


def _cesm_members(variant):
    import fsspec
    import xarray as xr
    ds = xr.open_zarr(fsspec.get_mapper(S3.format(exp="historical", variant=variant, var="TREFHT"), anon=True), consolidated=True)
    fut = xr.open_zarr(fsspec.get_mapper(S3.format(exp="ssp370", variant=variant, var="TREFHT"), anon=True), consolidated=True)
    return sorted(set(map(str, ds.member_id.values)) & set(map(str, fut.member_id.values)))


if __name__ == "__main__":
    sys.exit(main())
