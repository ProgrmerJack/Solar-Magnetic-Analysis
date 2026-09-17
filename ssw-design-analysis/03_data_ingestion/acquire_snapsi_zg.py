#!/usr/bin/env python3
"""
acquire_snapsi_zg.py
====================
SNAPSI geopotential height at 100 hPa, reduced to a polar-cap mean, fetched by
OPeNDAP subsetting rather than by downloading whole files.

WHY A SECOND ACQUISITION SCRIPT
  `acquire_snapsi_surface.py` downloads whole `psl` files and reduces them in
  memory. That cannot work for `zg`: it is 3-D on 34 pressure levels and ~34x
  the size of psl per file, putting `nudged`+`control` for all centres at
  roughly **5 TB** against 398 GB of disk. The transport has to differ, so the
  script does. Searched before writing this: nothing in the repo speaks CEDA
  OPeNDAP -- the existing dodsC users (`acquire_heatflux.py`,
  `extend_ncep_presatellite.py`, and the superseded `scripts/acquisition/*`) all
  target NOAA PSL, a different server with no authentication.

WHAT IT IS FOR
  Loeffel et al. (2026, WCD 7, 895-913) report r = 0.85 between the week-2
  100 hPa polar-cap GPH anomaly and the surface response over weeks 3-7, across
  18 events, and conclude that SSWs differ in their capacity to couple downward.
  That correlation is computed ACROSS EVENTS on ensemble means. In a SNAPSI
  `nudged` ensemble the event is identical in every member by construction, so
  the same relation can be measured BETWEEN MEMBERS, where no event-to-event
  difference exists to explain it.
  This script only acquires the field. The test is a separate analysis.

WHY 100 hPa IS NOT NUDGED, WHICH IS WHAT MAKES THE TEST POSSIBLE
  The SNAPSI protocol (Hitchcock et al., GMD 15, 5073, 2022) nudges only the
  ZONAL-MEAN temperature and zonal wind, at full strength above 50 hPa, tapering
  to no nudging below 90 hPa. Eddies are free at every level. So 100 hPa
  polar-cap zg is not directly constrained and members genuinely differ there --
  but that has to be MEASURED, not assumed, and the analysis gates on it.

TRANSPORT, AND THE TRAP IT IS GATED AGAINST
  `https://dap.ceda.ac.uk/thredds/dodsC/<archive path>` with the same bearer JWT
  the file download uses, opened through pydap with a requests.Session carrying
  the header. Subsetting to one level and one 30-degree cap is **152x smaller**:
  204.8 MB -> 1.35 MB in 6 s.

  NCEP OPeNDAP once returned exactly 0.0 for every day of a year through this
  same class of interface -- no error, no NaN (see the project's data-access
  notes). So every centre is gated on first use: its subset is compared against
  a full-file download of the same field, and the centre is refused unless they
  agree exactly. CanESM5 passing does not license the rest.

  Two dead ends, recorded so they are not retried: hand-parsing the `.dods`
  binary returns the coordinate MAPS rather than the array, because DAP2 Grid
  payloads append them after the data; and `.nc` responses return HTTP 400 here.

Output: snapsi_polarcap_zg100.parquet
        centre, model, experiment, init, member, lead_days,
        zg100_cap_N, zg100_cap_S, hemisphere, zg100_cap
"""
from __future__ import annotations

import io
import json
import os
import pathlib
import re
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import requests
import xarray as xr

HERE = pathlib.Path(__file__).resolve().parent
MANIFEST = HERE / "snapsi_manifest.csv"
OUT = HERE / "snapsi_polarcap_zg100.parquet"
CACHE = HERE / "_snapsi_reduced"          # shared; keys carry the variable
GATE = HERE / "snapsi_zg_opendap_gate.json"
ENV = pathlib.Path(os.environ.get("CEDA_ENV", HERE.parents[1] / ".env"))

DAP_ROOT = "https://dap.ceda.ac.uk/thredds/dodsC"
FILE_ROOT = "https://dap.ceda.ac.uk"

NH_INITS = ["s20180125", "s20180208", "s20181213", "s20190108"]
SH_INITS = ["s20190829", "s20191001"]
ALL_INITS = NH_INITS + SH_INITS
EXPERIMENTS = ["nudged", "control"]
ALL_CENTRES = ["CCCma", "CNR-ISAC", "ECCC", "ECMWF", "KMA", "Meteo-France",
               "NCAR", "NRL", "SNU", "UKMO"]

TARGET_PA = 10000.0      # 100 hPa, the Loeffel level
CAP_LAT = 60.0
WORKERS = 8              # gentler than the psl run: each request is a server-side subset


# ------------------------------------------------------------------- auth ---

def token() -> str:
    txt = ENV.read_text(encoding="utf8", errors="ignore")
    m = re.search(r"^\s*CEDA\s*=\s*['\"]?([A-Za-z0-9_\-\.]+)", txt, re.M)
    if not m:
        sys.exit(f"no CEDA token in {ENV}")
    tok = m.group(1)
    import base64
    pay = json.loads(base64.urlsafe_b64decode(
        tok.split(".")[1] + "=" * (-len(tok.split(".")[1]) % 4)))
    left = pay.get("exp", 0) - time.time()
    if left <= 0:
        sys.exit("CEDA token EXPIRED -- refresh it before running")
    print(f"CEDA token valid, {left / 3600:.1f} h remaining "
          f"(user {pay.get('preferred_username', '?')})")
    return tok


def session(tok: str) -> requests.Session:
    s = requests.Session()
    s.headers.update({"Authorization": f"Bearer {tok}"})
    return s


def zg_path(row) -> str:
    """The zg twin of a psl archive path.

    Layout is .../<freq>/<var>/<grid>/<version>/<var>_<freq>_... so only the
    variable token changes. Verified to resolve for all 10 centres by HEAD.
    """
    return (row.archive_path
            .replace("/psl/", "/zg/")
            .replace("psl_6hrPt", "zg_6hrPt"))


# ------------------------------------------------------------- reduction ---

def cap_means(da, lat):
    north = da.sel(lat=lat[lat >= CAP_LAT].values)
    south = da.sel(lat=lat[lat <= -CAP_LAT].values)
    n = north.weighted(np.cos(np.deg2rad(north.lat))).mean(dim=["lat", "lon"]).values
    s = south.weighted(np.cos(np.deg2rad(south.lat))).mean(dim=["lat", "lon"]).values
    return np.asarray(n, dtype=float), np.asarray(s, dtype=float)


_UNIT_TO_DAYS = {"day": 1.0, "days": 1.0,
                 "hour": 1 / 24, "hours": 1 / 24,
                 "minute": 1 / 1440, "minutes": 1 / 1440,
                 "second": 1 / 86400, "seconds": 1 / 86400}


def lead_days(tv, units=None):
    """Forecast lead in days, without reconciling the calendar.

    Absolute dates are never needed -- lead is a difference within one file --
    so no 365-vs-366 drift from the model calendars can creep in. CCCma runs a
    365_day calendar, which is exactly why decoding is avoided.

    Opened with decode_times=False, the axis arrives as raw numbers in the
    file's own units ("days since 1850-01-01" for CanESM5, 6-hourly, so values
    step by 0.25). The first version assumed cftime objects and called
    .total_seconds() on a float, failing every member. The units attribute is
    honoured rather than assumed, and an unrecognised one RAISES: silently
    treating hours as days would scale every lead by 24 and quietly misplace
    the analysis windows.
    """
    arr = np.asarray(tv)
    t0 = arr[0]
    if np.issubdtype(arr.dtype, np.datetime64):
        return (pd.to_datetime(arr) - pd.to_datetime(t0)).total_seconds() / 86400.0
    if np.issubdtype(arr.dtype, np.number):
        if not units:
            raise ValueError("numeric time axis with no units attribute")
        unit = str(units).strip().split()[0].lower()
        if unit not in _UNIT_TO_DAYS:
            raise ValueError(f"unhandled time unit {unit!r} in {units!r}")
        return (arr.astype(float) - float(t0)) * _UNIT_TO_DAYS[unit]
    # cftime objects
    return np.array([(x - t0).total_seconds() / 86400.0 for x in arr], dtype=float)


def plausible(z) -> bool:
    """100 hPa geopotential height sits near 15-17 km. Anything else is not it."""
    if z is None or len(z) == 0 or not np.isfinite(z).any():
        return False
    return 13000.0 < float(np.nanmean(z)) < 18000.0


def cache_ok(df) -> bool:
    if df is None or len(df) == 0:
        return False
    need = {"centre", "experiment", "init", "member", "lead_days",
            "zg100_cap_N", "zg100_cap_S", "hemisphere", "zg100_cap"}
    if not need <= set(df.columns):
        return False
    return plausible(df["zg100_cap_N"].values) and plausible(df["zg100_cap_S"].values)


# ------------------------------------------------- the per-centre gate ------

def verify_centre(row, tok, sess) -> dict:
    """Compare an OPeNDAP subset against the full file for the SAME field.

    Refuses the centre unless they agree exactly. This is the guard against the
    silent-zeros failure mode, in which a subsetting server returns a perfectly
    shaped array of zeros with no error and no NaN.
    """
    path = zg_path(row)
    dap = f"{DAP_ROOT}/{path.lstrip('/')}"

    req = urllib.request.Request(f"{FILE_ROOT}{path}?download=1",
                                 headers={"Authorization": f"Bearer {tok}"})
    t0 = time.time()
    with urllib.request.urlopen(req, timeout=1800) as r:
        blob = r.read()
    dl = time.time() - t0

    with xr.open_dataset(io.BytesIO(blob)) as ds:
        lat = ds["lat"]
        truth = ds["zg"].sel(plev=TARGET_PA).sel(lat=lat[lat >= CAP_LAT].values).values

    t0 = time.time()
    with xr.open_dataset(dap, engine="pydap", session=sess, decode_times=False) as rem:
        rlat = rem["lat"]
        got = (rem["zg"].sel(plev=TARGET_PA)
               .sel(lat=rlat[rlat >= CAP_LAT].values).load().values)
    sub = time.time() - t0

    a = np.asarray(truth, float)
    b = np.asarray(got, float).reshape(a.shape)
    res = {
        "centre": row.centre,
        "file": pathlib.Path(path).name,
        "full_MB": round(len(blob) / 1e6, 2),
        "full_seconds": round(dl, 1),
        "subset_MB": round(b.size * 4 / 1e6, 3),
        "subset_seconds": round(sub, 1),
        "saving_x": round(len(blob) / max(b.size * 4, 1), 1),
        "all_zero": bool(np.all(b == 0)),
        "any_nan": bool(np.isnan(b).any()),
        "max_abs_diff": float(np.max(np.abs(b - a))),
        "correlation": float(np.corrcoef(b.ravel(), a.ravel())[0, 1]),
        "truth_mean": float(a.mean()),
        "subset_mean": float(b.mean()),
        "plausible": plausible(b),
    }
    res["passes"] = bool(
        res["max_abs_diff"] == 0.0 and not res["all_zero"]
        and not res["any_nan"] and res["plausible"])
    return res


# ------------------------------------------------------------- the fetch ---

def fetch(row, sess):
    key = f"zg100_{row.centre}_{row.experiment}_{row.start_date}_{row.member}"
    cf = CACHE / f"{key}.parquet"
    if cf.exists():
        try:
            cached = pd.read_parquet(cf)
        except Exception:
            cached = None
        if cache_ok(cached):
            return cached
        print(f"  cache REJECTED, refetching: {cf.name}", flush=True)
        cf.unlink(missing_ok=True)

    dap = f"{DAP_ROOT}/{zg_path(row).lstrip('/')}"
    last = None
    for attempt in range(4):
        try:
            with xr.open_dataset(dap, engine="pydap", session=sess,
                                 decode_times=False) as rem:
                lat = rem["lat"]
                da = rem["zg"].sel(plev=TARGET_PA)
                cap_n, cap_s = cap_means(da.load(), lat)
                lead = lead_days(rem["time"].values,
                                 rem["time"].attrs.get("units"))
            if not (plausible(cap_n) and plausible(cap_s)):
                raise ValueError(
                    f"implausible 100 hPa GPH: N={np.nanmean(cap_n):.0f} "
                    f"S={np.nanmean(cap_s):.0f} m")
            hemi = "S" if str(row.start_date) in SH_INITS else "N"
            out = pd.DataFrame({
                "centre": row.centre, "model": row.model,
                "experiment": row.experiment, "init": str(row.start_date),
                "member": row.member, "lead_days": lead,
                "zg100_cap_N": cap_n, "zg100_cap_S": cap_s,
                "hemisphere": hemi,
                "zg100_cap": cap_s if hemi == "S" else cap_n})
            tmp = cf.with_suffix(".parquet.tmp")
            out.to_parquet(tmp)
            tmp.replace(cf)
            return out
        except (ImportError, ModuleNotFoundError) as exc:
            raise SystemExit(
                f"missing dependency, not a transfer error: {exc}\n"
                f"pydap is required; check with environment/check_environment.py") from exc
        except Exception as exc:
            last = exc
            if attempt == 0:
                print(f"  retrying {row.centre}/{row.start_date}/{row.member}: "
                      f"{type(exc).__name__}: {str(exc)[:120]}", flush=True)
            time.sleep(2 * (attempt + 1))
    raise IOError(f"{row.centre}/{row.experiment}/{row.start_date}/{row.member}: {last}")


def main() -> int:
    argv = sys.argv[1:]

    def opt(name, default=None):
        return argv[argv.index(name) + 1] if name in argv else default

    centres = opt("--centres")
    centres = [c.strip() for c in centres.split(",")] if centres else ALL_CENTRES
    inits_arg = (opt("--inits") or "all").lower()
    inits = {"all": ALL_INITS, "nh": NH_INITS, "sh": SH_INITS}.get(
        inits_arg, [i.strip() for i in inits_arg.split(",")])

    CACHE.mkdir(exist_ok=True)
    d = pd.read_csv(MANIFEST)
    q = d[(d.variable == "psl")                      # psl rows carry the paths
          & (d.experiment.isin(EXPERIMENTS))
          & (d.start_date.astype(str).isin(inits))
          & (d.centre.isin(centres))].reset_index(drop=True)

    def cached(r):
        return (CACHE / f"zg100_{r.centre}_{r.experiment}_"
                        f"{r.start_date}_{r.member}.parquet").exists()

    have = q.apply(cached, axis=1) if len(q) else pd.Series(dtype=bool)
    print(f"centres    {centres}")
    print(f"inits      {inits}")
    print(f"matched    {len(q):,} members")
    print(f"cached     {int(have.sum()):,}")
    print(f"to fetch   {len(q) - int(have.sum()):,}", flush=True)
    if "--dry-run" in argv:
        return 0

    tok = token()
    sess = session(tok)

    # ---- gate every centre before fetching any of its members ----
    gate = json.loads(GATE.read_text()) if GATE.exists() else {}
    todo_centres = [c for c in q.centre.unique()
                    if not gate.get(c, {}).get("passes")]
    if "--skip-gate" in argv and todo_centres:
        print(f"WARNING --skip-gate: {todo_centres} unverified", flush=True)
        todo_centres = []
    for c in todo_centres:
        row = q[q.centre == c].iloc[0]
        print(f"\ngating {c} -- full file vs OPeNDAP subset ...", flush=True)
        try:
            res = verify_centre(row, tok, sess)
        except Exception as exc:
            res = {"centre": c, "passes": False, "error": f"{type(exc).__name__}: {exc}"}
        gate[c] = res
        GATE.write_text(json.dumps(gate, indent=1))
        if res.get("passes"):
            print(f"  PASS  max|diff|={res['max_abs_diff']:g} "
                  f"corr={res['correlation']:.10f} "
                  f"mean={res['subset_mean']:.1f} m  "
                  f"{res['saving_x']}x smaller", flush=True)
        else:
            print(f"  FAIL  {res}", flush=True)

    ok_centres = [c for c in q.centre.unique() if gate.get(c, {}).get("passes")]
    refused = sorted(set(q.centre.unique()) - set(ok_centres))
    if refused:
        print(f"\nREFUSED (failed the subset gate): {refused}", flush=True)
    q = q[q.centre.isin(ok_centres)].reset_index(drop=True)
    if q.empty:
        sys.exit("no centre passed the gate")

    rows, failed = [], []
    t0 = time.time()
    with ThreadPoolExecutor(WORKERS) as ex:
        futs = {ex.submit(fetch, r, sess): i for i, r in enumerate(q.itertuples())}
        for n, f in enumerate(as_completed(futs), 1):
            try:
                rows.append(f.result())
            except Exception as exc:
                failed.append(exc)
                print(f"  FAILED {exc}", flush=True)
            if n % 100 == 0:
                el = time.time() - t0
                print(f"  {n}/{len(q)}  {el/60:.1f} min  "
                      f"({n/el*60:.0f} members/min)", flush=True)

    files = sorted(CACHE.glob("zg100_*.parquet"))
    for s in CACHE.glob("zg100_*.parquet.tmp"):
        s.unlink(missing_ok=True)
    frames, dropped = [], 0
    for f in files:
        try:
            df = pd.read_parquet(f)
        except Exception:
            dropped += 1
            continue
        if cache_ok(df):
            frames.append(df)
        else:
            dropped += 1
    if dropped:
        print(f"excluded {dropped} member(s) that failed validation")
    if not frames:
        sys.exit("nothing usable")
    out = pd.concat(frames, ignore_index=True)
    out.to_parquet(OUT)
    print(f"\n{len(out):,} rows from {len(frames):,} members -> {OUT.name} "
          f"({len(failed)} failed this run)")
    print(out.groupby(["centre", "experiment"]).init.nunique()
          .unstack(fill_value=0).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
