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

SECOND VARIABLE: `--var t2m` (plan approved 2026-09-26, regional cold risk)
  Daily-mean 2 m temperature ("2_m_temperature", requested as 24 h lead windows
  "24k_24(k+1)", k = 10..33 -- post-onset days +8..+24 for starts 2-9 d before
  onset; JMA's 12 h-offset windows are not requested), area 35-70N, reduced to the
  four regions of snapsi_regional_test.py (cos-lat means, land and sea):
  NEURASIA 50-65N 10-130E, HI_EUROPE 55-70N 0-60E, MID_EASIA 35-55N 90-150E,
  MID_NAMER 35-55N 120-60W. lead_day is the START day of the 24 h window; valid is
  init + lead_day. Centres define the daily average differently (true daily mean
  or the mean of 00/06/12/18 UTC); anomalies are taken within each system.
  Retrieved by the packed path only, ECMWF included (hindcasts on the fly).
  Output: s2s_<origin>_t2m_regions.parquet (ECMWF: s2s_ecmf_t2m_regions.parquet)
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
VAR = "msl"                  # "t2m" with --var t2m
T2M_CACHE = HERE / "_s2s_reduced_t2m"
T2M_AREA = [70, -180, 35, 180]
T2M_LEADS = list(range(10, 34))
T2M_REGIONS = {"NEURASIA": (50, 65, 10, 130), "HI_EUROPE": (55, 70, 0, 60),
               "MID_EASIA": (35, 55, 90, 150), "MID_NAMER": (35, 55, 240, 300)}
WORKERS = int(os.environ.get("S2S_WORKERS", 6))   # requests in flight (queue-limit rejections are retried)
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


def request(md, origin="ecmwf", hm=None, hd=None, hyears=None, leads=None):
    hm, hd = hm or md.month, hd or md.day
    if VAR == "t2m":
        return {"origin": [origin], "year": [f"{md.year}"], "month": [f"{md.month:02d}"],
                "day": [f"{md.day:02d}"], "time": ["00:00"], "hyear": hyears or HYEARS,
                "hmonth": [f"{hm:02d}"], "hday": [f"{hd:02d}"],
                "level_type": "single_level", "variable": ["2_m_temperature"],
                "forecast_type": ["control_forecast", "perturbed_forecast"],
                "leadtime_hour": [f"{24 * d}_{24 * (d + 1)}" for d in T2M_LEADS],
                "data_format": "grib", "area": T2M_AREA}
    return {"origin": [origin], "year": [f"{md.year}"], "month": [f"{md.month:02d}"],
            "day": [f"{md.day:02d}"], "time": ["00:00"], "hyear": hyears or HYEARS,
            "hmonth": [f"{hm:02d}"], "hday": [f"{hd:02d}"],
            "level_type": "single_level", "variable": ["mean_sea_level_pressure"],
            "forecast_type": ["control_forecast", "perturbed_forecast"],
            "leadtime_hour": [str(24 * d) for d in (leads or LEAD_DAYS)],
            "data_format": "grib", "area": AREA}


def regions_mean(da):
    """cos-lat means over T2M_REGIONS (bounds inclusive; longitudes 0-360)."""
    da = da.assign_coords(longitude=np.mod(da["longitude"], 360.0))
    out = {}
    for r, (la0, la1, lo0, lo1) in T2M_REGIONS.items():
        m = ((da.latitude >= la0) & (da.latitude <= la1)
             & (da.longitude >= lo0) & (da.longitude <= lo1))
        sub = da.where(m, drop=True)
        out[r] = sub.weighted(np.cos(np.deg2rad(sub.latitude))).mean(("latitude", "longitude"))
    return xr.Dataset(out)


def reduce(grib, md, origin="ecmwf", n_hyears=None, leads=None):
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
    # model_date: the calendar start key the anomaly climatology groups on. For
    # ECMWF (and every on-the-fly system) it is the model date; for a fixed
    # hindcast set it is the hindcast month-day (the model version date is the
    # same for every start and would lump them together).
    df["model_date"] = md if origin == "ecmwf" else pd.Timestamp(
        year=2000, month=int(df["init"].dt.month.iloc[0]), day=int(df["init"].dt.day.iloc[0]))
    df = df[["model_date", "init", "member", "lead_day", "valid", "psl_cap_N"]]
    if origin == "ecmwf":
        want = len(HYEARS) * N_MEMBERS * len(LEAD_DAYS)
        if len(df) != want or df["psl_cap_N"].isna().any():
            raise ValueError(f"{md.date()}: {len(df)} rows, {df['psl_cap_N'].isna().sum()} NaN; "
                             f"expected {want} complete")
    else:
        # complete and rectangular: every start has the same members at every lead
        per = df.groupby(["init", "member"])["lead_day"].nunique()
        mem = df.groupby("init")["member"].nunique()
        if (df["psl_cap_N"].isna().any() or (per != len(leads)).any()
                or mem.nunique() != 1 or df["init"].nunique() != n_hyears):
            raise ValueError(f"{origin} {md.date()}: incomplete -- starts {df['init'].nunique()}/"
                             f"{n_hyears}, members per start {sorted(mem.unique())}, "
                             f"leads {sorted(per.unique())}/{len(leads)}, NaN {df['psl_cap_N'].isna().sum()}")
    if not df["psl_cap_N"].between(95000, 106000).all():
        raise ValueError(f"{md.date()}: cap mean outside 950-1060 hPa")
    return df


# Other S2S systems (multi-model test, plan approved 2026-09-25). Short-lead
# window only (starts 2-9 d before onset -> leads 10-34 d); one model version per
# centre, the one that best covers the 2003-2021 events; Dec-Mar hindcast starts.
# "otf": hindcasts follow the model date (same day and month in each hindcast
# year); "fixed": a fixed hindcast set addressed by hmonth/hday. Daily-start
# fixed sets are thinned to every THIN-th day. BoM (2014, poorly resolved
# stratosphere, very large hindcast ensemble) and UKMO / IAP-CAS (no msl) are not
# used. ECMWF keeps its own configuration above, unchanged.
CENTRES = {
    "eccc":     {"year": "2025", "kind": "otf",   "thin": 1},
    "cma":      {"year": "2022", "kind": "otf",   "thin": 1},
    "hmcr":     {"year": "2025", "kind": "otf",   "thin": 1},
    "kma":      {"year": "2026", "kind": "otf",   "thin": 1},
    "cnrm":     {"year": "2025", "kind": "fixed", "thin": 1},
    # JMA lists two 2022 model versions (31 Mar, 30 Sep); MARS holds no msl for
    # the 31 Mar one ("MARS returned no data", 2026-09-26), and one version per
    # centre is the design anyway: the 30 Sep version only.
    "jma":      {"year": "2022", "kind": "fixed", "thin": 1, "model_md": ("09", "30")},
    "cnr_isac": {"year": "2023", "kind": "fixed", "thin": 1},
    "ncep":     {"year": "2011", "kind": "fixed", "thin": 3},
    "cptec":    {"year": "2023", "kind": "fixed", "thin": 3},
}
MULTI_LEAD_DAYS = list(range(10, 35))


def jobs(origin):
    """One job per (model date, hindcast month-day), all hindcast years, from the
    ECDS constraint service. For ECMWF this reproduces the original model dates."""
    p = client().get_process(DATASET)
    if origin == "ecmwf":
        root = CACHE if VAR == "msl" else T2M_CACHE / "ecmwf"
        return [{"origin": "ecmwf", "md": md, "hm": md.month, "hd": md.day,
                 "hyears": HYEARS, "key": str(md.date()),
                 "cache": root / f"{md.date()}.parquet"} for md in model_dates()]
    cfg = CENTRES[origin]
    q = {"origin": [origin], "year": [cfg["year"]],
         "variable": ["mean_sea_level_pressure" if VAR == "msl" else "2_m_temperature"]}
    # constraint queries in parallel: one sequential query per (date, hindcast
    # date) took most of an hour when ECDS was returning 502s
    def days_of(m):
        return [(m, d) for d in p.apply_constraints({**q, "month": [m]}).get("day", [])]

    def combos_of(md):
        m, d = md
        q2 = {**q, "month": [m], "day": [d]}
        out_ = []
        for hm in p.apply_constraints(q2).get("hmonth", []):
            if hm not in MONTHS or (cfg["kind"] == "otf" and hm != m):
                continue
            for hd in p.apply_constraints({**q2, "hmonth": [hm]}).get("hday", []):
                if (cfg["kind"] == "otf" and hd != d) or (int(hd) - 1) % cfg["thin"]:
                    continue
                out_.append((m, d, hm, hd))
        return out_

    def job_of(c):
        m, d, hm, hd = c
        hy = p.apply_constraints({**q, "month": [m], "day": [d], "hmonth": [hm],
                                  "hday": [hd]}).get("hyear", [])
        key = f"{cfg['year']}-{m}-{d}_h{hm}-{hd}"
        return {"origin": origin, "md": pd.Timestamp(f"{cfg['year']}-{m}-{d}"),
                "hm": int(hm), "hd": int(hd), "hyears": sorted(hy), "key": key,
                "cache": (CACHE if VAR == "msl" else T2M_CACHE) / origin / f"{key}.parquet"}

    months = p.apply_constraints(q).get("month", [])
    with ThreadPoolExecutor(8) as ex:
        mds = [x for lst in ex.map(days_of, months) for x in lst]
        if "model_md" in cfg:
            mds = [md for md in mds if md == cfg["model_md"]]
        combos = [x for lst in ex.map(combos_of, mds) for x in lst]
        out = sorted(ex.map(job_of, combos), key=lambda j: j["key"])
    return out


def client():
    # ECDS runs ONE request at a time per account (measured 2026-09-26: each job
    # started as the previous finished, 2-11 min of server time each), so client
    # threads cannot add throughput; a second account can. ECDS_RC names another
    # rc file (same format as ~/.cdsapirc) for a second process on other origins.
    from ecmwf.datastores import Client
    rc = Path(os.environ.get("ECDS_RC", "~/.cdsapirc")).expanduser()
    key = rc.read_text().split("key:")[1].split()[0]
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


def fetch(job):
    out = job["cache"]
    if out.exists():
        return job, "cached", None
    out.parent.mkdir(parents=True, exist_ok=True)
    leads = LEAD_DAYS if job["origin"] == "ecmwf" else MULTI_LEAD_DAYS
    with tempfile.TemporaryDirectory() as td:
        g = Path(td) / "req.grib"
        # ECDS limits queued requests per user and dataset ("Number queued
        # requests for this dataset is temporarily limited"): a rejection is
        # retried with back-off, never treated as data failure.
        for attempt in range(40):
            remote = client().submit(DATASET, request(job["md"], job["origin"], job["hm"], job["hd"],
                                                      job["hyears"], leads))
            # read the status ONCE per poll: every access re-queries ECDS, and
            # checking it three times raced ("FAILED ... ended successful")
            while True:
                st = remote.status
                if st in ("successful", "failed", "rejected", "dismissed", "deleted"):
                    break
                time.sleep(10)
            if st != "rejected":
                break
            time.sleep(60 + 30 * min(attempt, 8))
        if st != "successful":
            raise IOError(f"ECDS request {remote.request_id} ended {st}")
        res = remote.get_results()
        ranged_download(res.location, int(res.content_length), g)
        size = g.stat().st_size
        sha = hashlib.sha256(g.read_bytes()).hexdigest()
        df = reduce(str(g), job["md"], job["origin"], len(job["hyears"]), leads)
    df.to_parquet(out)
    return job, "fetched", {"bytes": size, "sha256": sha}


def main(origin="ecmwf"):
    CACHE.mkdir(exist_ok=True)
    js = jobs(origin)
    if VAR == "t2m":
        missing = [j["key"] for j in js if not j["cache"].exists()]
        if missing:
            print(f"{origin} t2m: {len(missing)} starts not cached; run --packed first")
            return 1
        tag = "ecmf" if origin == "ecmwf" else origin
        df = pd.concat([pd.read_parquet(j["cache"]) for j in js],
                       ignore_index=True).sort_values(["init", "member", "lead_day"])
        out_path = HERE / f"s2s_{tag}_t2m_regions.parquet"
        df.to_parquet(out_path)
        print(f"{len(df):,} rows, {df['init'].nunique()} starts -> {out_path.name}")
        return 0
    out_path = OUT if origin == "ecmwf" else HERE / f"s2s_{origin}_psl_cap.parquet"
    man_path = MANIFEST if origin == "ecmwf" else HERE / f"s2s_manifest_{origin}.json"
    print(f"{origin}: {len(js)} requests; hindcast years "
          f"{min(min(j['hyears']) for j in js)}-{max(max(j['hyears']) for j in js)}", flush=True)
    man = json.loads(man_path.read_text()) if man_path.exists() else {}
    failed = []
    with ThreadPoolExecutor(WORKERS) as ex:
        futs = {ex.submit(fetch, j): j for j in js}
        for fu in as_completed(futs):
            j = futs[fu]
            try:
                _, how, meta = fu.result()
                if meta:
                    man[j["key"]] = meta
                    man_path.write_text(json.dumps(man, indent=2, sort_keys=True), newline="\n")
                print(f"  {j['key']} {how}", flush=True)
            except Exception as ex_:
                failed.append(j["key"])
                print(f"  {j['key']} FAILED: {str(ex_)[:200]}", flush=True)
    if failed:
        print(f"\n{len(failed)} requests failed: {failed}; re-run to resume")
        return 1
    df = pd.concat([pd.read_parquet(j["cache"]) for j in js],
                   ignore_index=True).sort_values(["init", "member", "lead_day"])
    df.to_parquet(out_path)
    print(f"\n{len(df):,} rows, {df['init'].nunique()} starts "
          f"({df['init'].min().date()} .. {df['init'].max().date()}), "
          f"members {sorted(int(m) for m in df['member'].unique())}")
    print(f"Saved -> {out_path.name}")
    return 0


# ---------------------------------------------------------------------------
# PACKED retrieval. ECDS runs about four requests per user and dataset at a
# time and rejects the rest (measured 2026-09-25: 74 of 80 rejected when 27 were
# in flight), so throughput is set by the NUMBER of requests. A request may list
# several dates: ECDS returns only the valid (date, hindcast date) combinations
# (checked on ECCC and CNRM). One request per centre, model date and hindcast
# month; the result is split back into the per-start cache files above, and each
# start is checked for completeness on its own.
QUEUE_SLOTS = 4


def packs(origin):
    out = {}
    for j in jobs(origin):
        if j["cache"].exists():
            continue
        cfg = CENTRES.get(origin, {"kind": "otf"})         # ECMWF: on the fly
        k = (j["md"].month, j["hm"]) if cfg["kind"] == "otf" else (j["md"], j["hm"])
        out.setdefault(k, []).append(j)
    return [(origin, v) for v in out.values()]


def fetch_pack(item):
    origin, js = item
    leads = MULTI_LEAD_DAYS if VAR == "msl" else T2M_LEADS
    col = "psl_cap_N" if VAR == "msl" else list(T2M_REGIONS)
    req = request(js[0]["md"], origin, js[0]["hm"], js[0]["hd"],
                  sorted({h for j in js for h in j["hyears"]}), leads)
    req["day"] = sorted({f"{j['md'].day:02d}" for j in js})
    req["hday"] = sorted({f"{j['hd']:02d}" for j in js})
    with tempfile.TemporaryDirectory() as td:
        g = Path(td) / "req.grib"
        for attempt in range(60):
            remote = client().submit(DATASET, req)
            while True:
                st = remote.status
                if st in ("successful", "failed", "rejected", "dismissed", "deleted"):
                    break
                time.sleep(10)
            if st != "rejected":
                break
            time.sleep(60)
        if st != "successful":
            raise IOError(f"ECDS request {remote.request_id} ended {st}")
        res = remote.get_results()
        ranged_download(res.location, int(res.content_length), g)
        sha = hashlib.sha256(g.read_bytes()).hexdigest()
        size = g.stat().st_size
        frames = []
        for dtype, member_of in (("cf", None), ("pf", "number")):
            ds = xr.open_dataset(str(g), engine="cfgrib",
                                 backend_kwargs={"indexpath": "", "filter_by_keys": {"dataType": dtype}})
            if VAR == "msl":
                c = E.cap(ds[list(ds.data_vars)[0]], north=True)
                df = c.to_dataframe(name="psl_cap_N").reset_index()
            else:
                df = regions_mean(ds[list(ds.data_vars)[0]]).to_dataframe().reset_index()
            df["member"] = 0 if member_of is None else df["number"].astype(int)
            frames.append(df); ds.close()
    df = pd.concat(frames, ignore_index=True)
    df["init"] = pd.to_datetime(df["time"])
    df["lead_day"] = (pd.to_timedelta(df["step"]) / pd.Timedelta(days=1)).round().astype(int)
    if VAR == "t2m":
        # an averaged field carries the window END as its step: "240_264" -> 11 d.
        # lead_day is the window START; the requested set must come back exactly
        df["lead_day"] -= 1
        if set(df["lead_day"]) != set(T2M_LEADS):
            raise ValueError(f"{origin}: lead days {sorted(set(df['lead_day']))[:4]}... "
                             f"!= requested windows {T2M_LEADS[0]}..{T2M_LEADS[-1]}")
        if not df[col].stack().between(200, 310).all():
            raise ValueError(f"{origin}: regional 2 m temperature outside 200-310 K")
    df["valid"] = df["init"] + pd.to_timedelta(df["lead_day"], unit="D")
    if VAR == "msl" and not df["psl_cap_N"].between(95000, 106000).all():
        raise ValueError(f"{origin}: cap mean outside 950-1060 hPa")
    done = []
    for j in js:
        sub = df[(df["init"].dt.month == j["hm"]) & (df["init"].dt.day == j["hd"])
                 & df["init"].dt.year.astype(str).isin(j["hyears"])].copy()
        sub["model_date"] = pd.Timestamp(year=2000, month=j["hm"], day=j["hd"])
        cols = [col] if isinstance(col, str) else col
        sub = sub[["model_date", "init", "member", "lead_day", "valid", *cols]]
        per = sub.groupby(["init", "member"])["lead_day"].nunique()
        mem = sub.groupby("init")["member"].nunique()
        if (sub.empty or sub[cols].isna().any().any() or (per != len(leads)).any()
                or mem.nunique() != 1 or sub["init"].nunique() != len(j["hyears"])):
            raise ValueError(f"{origin} {j['key']}: incomplete in packed result")
        j["cache"].parent.mkdir(parents=True, exist_ok=True)
        sub.to_parquet(j["cache"])
        done.append(j["key"])
    return origin, done, {"bytes": size, "sha256": sha, "starts": done}


def main_packed(origins):
    with ThreadPoolExecutor(len(origins)) as ex:
        allp = [p for lst in ex.map(packs, origins) for p in lst]
    print(f"{len(allp)} packed requests for {origins}", flush=True)
    failed = []
    with ThreadPoolExecutor(QUEUE_SLOTS) as ex:
        futs = {ex.submit(fetch_pack, p): p for p in allp}
        for fu in as_completed(futs):
            origin, js = futs[fu]
            try:
                _, done, meta = fu.result()
                mp_ = HERE / (f"s2s_manifest_{origin}.json" if VAR == "msl"
                              else f"s2s_manifest_t2m_{origin}.json")
                man = json.loads(mp_.read_text()) if mp_.exists() else {}
                man["packed:" + "+".join(done)] = meta
                mp_.write_text(json.dumps(man, indent=2, sort_keys=True), newline="\n")
                print(f"  {origin} {len(done)} starts: {done[0]} .. {done[-1]} fetched", flush=True)
            except Exception as ex_:
                failed.append((origin, [j["key"] for j in js]))
                print(f"  {origin} pack {js[0]['key']} FAILED: {str(ex_)[:200]}", flush=True)
    print(f"packed done; {len(failed)} packs failed", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    if "--var" in sys.argv:
        i = sys.argv.index("--var")
        VAR = sys.argv[i + 1]
        del sys.argv[i:i + 2]
        assert VAR in ("msl", "t2m"), VAR
    if len(sys.argv) > 2 and sys.argv[1] == "--packed":
        sys.exit(main_packed(sys.argv[2:]))
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "ecmwf"))
