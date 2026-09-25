#!/usr/bin/env python3
"""
acquire_s2s_reforecasts.py
==========================
ECMWF S2S reforecasts of polar-cap mean sea-level pressure, reduced EXACTLY as
the ERA5 observation is (`acquire_era5_psl_cap.cap`: cos-lat weighted 60-90N,
the 60N row included), for the forecast-reliability test
(07_physical_decomposition/s2s_forecast_test.py).

WHY A SEPARATE SCRIPT
  Searched first: nothing in the repo retrieves S2S data. The ERA5 cap script
  was extended to 1998-2022 for this test but reads a different store.

SOURCE
  ECMWF Data Store (ECDS), dataset `s2s-reforecasts`, via `cdsapi` with
  ~/.cdsapirc (url https://ecds.ecmwf.int/api). The old Web-API route for S2S
  was decommissioned on 27 May 2026. Data licence CC BY-NC 4.0; acknowledgement
  wording and Vitart et al. (2017, BAMS 98, 163) are required in publications.

WHAT IS RETRIEVED (fixed before any download; analysis plan approved 2026-09-25)
  origin ecmwf, model year 2022 (one model version; hindcast years 2002-2021,
  checked with the ECDS constraint service), every Monday/Thursday model date
  in December-March, control + 10 perturbed members, 00 UTC, lead days 10-42
  (the +8..+25 day window after onset for starts 2-17 days before onset),
  area 60-90N. One request per model date: 20 hindcast years x 11 members x 33
  leads = 7,260 fields (ECDS size limit 1,000,000).

  Each GRIB is reduced in memory to the cap mean and deleted; only the reduced
  table is kept. A model date already in the cache is not requested again.

Output: s2s_ecmf_psl_cap.parquet
        model_date, init (hindcast start, 00 UTC), member (0 = control),
        lead_day, valid, psl_cap_N (Pa)
"""
import hashlib
import json
import os
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import acquire_era5_psl_cap as E                     # noqa: E402

OUT = HERE / "s2s_ecmf_psl_cap.parquet"
CACHE = HERE / "_s2s_reduced"
MANIFEST = HERE / "s2s_manifest.json"
DATASET = "s2s-reforecasts"
MODEL_YEAR = 2022
MONTHS = ("12", "01", "02", "03")
HYEARS = [str(y) for y in range(2002, 2022)]
LEAD_DAYS = list(range(10, 43))
AREA = [90, -180, 60, 180]
N_MEMBERS = 11
WORKERS = 3                  # model dates downloaded at once
RANGES = 16                  # parallel byte-ranges per file
CHUNK = 4 * 1024 * 1024


def model_dates():
    """Monday/Thursday model dates of MODEL_YEAR in MONTHS, from ECDS itself."""
    from ecmwf.datastores import Client
    key = Path("~/.cdsapirc").expanduser().read_text().split("key:")[1].split()[0]
    p = Client(url="https://ecds.ecmwf.int/api", key=key).get_process(DATASET)
    out = []
    for m in MONTHS:
        r = p.apply_constraints({"origin": ["ecmwf"], "year": [str(MODEL_YEAR)],
                                 "month": [m]})
        hy = set(r.get("hyear", []))
        if not set(HYEARS) <= hy:
            raise ValueError(f"model month {m}: hindcast years {sorted(hy)} "
                             f"do not cover {HYEARS[0]}-{HYEARS[-1]}")
        out += [pd.Timestamp(f"{MODEL_YEAR}-{m}-{d}") for d in r["day"]]
    return sorted(out)


def request(md):
    return {"origin": ["ecmwf"], "year": [f"{md.year}"], "month": [f"{md.month:02d}"],
            "day": [f"{md.day:02d}"], "time": ["00:00"], "hyear": HYEARS,
            "hmonth": [f"{md.month:02d}"], "hday": [f"{md.day:02d}"],
            "level_type": "single_level", "variable": ["mean_sea_level_pressure"],
            "forecast_type": ["control_forecast", "perturbed_forecast"],
            "leadtime_hour": [str(24 * d) for d in LEAD_DAYS],
            "data_format": "grib", "area": AREA}


def reduce(grib, md):
    frames = []
    for dtype, member_of in (("cf", None), ("pf", "number")):
        ds = xr.open_dataset(grib, engine="cfgrib",
                             backend_kwargs={"indexpath": "",
                                             "filter_by_keys": {"dataType": dtype}})
        v = ds[list(ds.data_vars)[0]]
        c = E.cap(v, north=True)                      # identical to the ERA5 reduction
        df = c.to_dataframe(name="psl_cap_N").reset_index()
        df["member"] = 0 if member_of is None else df["number"].astype(int)
        frames.append(df)
        ds.close()
    df = pd.concat(frames, ignore_index=True)
    df["init"] = pd.to_datetime(df["time"])
    df["lead_day"] = (pd.to_timedelta(df["step"]) / pd.Timedelta(days=1)).round().astype(int)
    df["valid"] = df["init"] + pd.to_timedelta(df["lead_day"], unit="D")
    df["model_date"] = md
    df = df[["model_date", "init", "member", "lead_day", "valid", "psl_cap_N"]]
    want = len(HYEARS) * N_MEMBERS * len(LEAD_DAYS)
    if len(df) != want or df["psl_cap_N"].isna().any():
        raise ValueError(f"{md.date()}: {len(df)} rows, {df['psl_cap_N'].isna().sum()} NaN; "
                         f"expected {want} complete")
    if not df["psl_cap_N"].between(95000, 106000).all():
        raise ValueError(f"{md.date()}: cap mean outside 950-1060 hPa")
    return df


def client():
    from ecmwf.datastores import Client
    key = Path("~/.cdsapirc").expanduser().read_text().split("key:")[1].split()[0]
    return Client(url="https://ecds.ecmwf.int/api", key=key)


def ranged_download(url, size, path):
    """ECDS serves each connection at ~60 KB/s (measured 2026-09-25: one
    connection 60 KB/s, 16 parallel ranges 915 KB/s), so the file is fetched as
    CHUNK-sized byte ranges on RANGES connections and assembled in place.
    Each range must come back as 206 with the exact Content-Range and length
    asked for, and the bytes written must add up to the server's size (the file
    is pre-allocated, so its size alone proves nothing -- review 2026-09-25)."""
    import requests
    with open(path, "wb") as f:
        f.truncate(size)
    spans = [(a, min(a + CHUNK, size) - 1) for a in range(0, size, CHUNK)]

    def one(span):
        a, b = span
        for attempt in range(8):
            try:
                r = requests.get(url, headers={"Range": f"bytes={a}-{b}"}, timeout=120)
                if (r.status_code == 206 and len(r.content) == b - a + 1
                        and r.headers.get("Content-Range") == f"bytes {a}-{b}/{size}"):
                    with open(path, "r+b") as f:
                        f.seek(a)
                        f.write(r.content)
                    return len(r.content)
            except requests.RequestException:
                pass
            time.sleep(5 * (attempt + 1))
        raise IOError(f"range {a}-{b} failed")
    with ThreadPoolExecutor(RANGES) as ex:
        written = sum(ex.map(one, spans))
    if written != size:
        raise IOError(f"wrote {written} of {size} bytes")


def fetch(md):
    out = CACHE / f"{md.date()}.parquet"
    if out.exists():
        return md, "cached", None
    with tempfile.TemporaryDirectory() as td:
        g = Path(td) / "req.grib"
        remote = client().submit(DATASET, request(md))
        while remote.status not in ("successful", "failed", "rejected", "dismissed", "deleted"):
            time.sleep(10)
        if remote.status != "successful":
            raise IOError(f"ECDS request {remote.request_id} ended {remote.status}")
        res = remote.get_results()
        ranged_download(res.location, int(res.content_length), g)
        size = g.stat().st_size
        sha = hashlib.sha256(g.read_bytes()).hexdigest()
        df = reduce(str(g), md)
    df.to_parquet(out)
    return md, "fetched", {"bytes": size, "sha256": sha}


def main():
    CACHE.mkdir(exist_ok=True)
    mds = model_dates()
    print(f"{len(mds)} model dates in {MODEL_YEAR} ({MONTHS}); "
          f"{len(HYEARS)} hindcast years; {len(LEAD_DAYS)} leads", flush=True)
    man = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {}
    failed = []
    with ThreadPoolExecutor(WORKERS) as ex:
        futs = {ex.submit(fetch, md): md for md in mds}
        for f in as_completed(futs):
            md = futs[f]
            try:
                _, how, meta = f.result()
                if meta:
                    man[str(md.date())] = meta
                    MANIFEST.write_text(json.dumps(man, indent=2, sort_keys=True),
                                        newline="\n")
                print(f"  {md.date()} {how}", flush=True)
            except Exception as e:
                failed.append(str(md.date()))
                print(f"  {md.date()} FAILED: {str(e)[:200]}", flush=True)
    if failed:
        print(f"\n{len(failed)} model dates failed: {failed}; re-run to resume")
        return 1
    df = pd.concat([pd.read_parquet(CACHE / f"{md.date()}.parquet") for md in mds],
                   ignore_index=True).sort_values(["init", "member", "lead_day"])
    df.to_parquet(OUT)
    print(f"\n{len(df):,} rows, {df['init'].nunique()} starts "
          f"({df['init'].min().date()} .. {df['init'].max().date()}), "
          f"members {sorted(df['member'].unique())}")
    print(f"Saved -> {OUT.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
