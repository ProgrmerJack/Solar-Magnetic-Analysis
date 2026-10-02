#!/usr/bin/env python3
"""
acquire_snapsi_surface.py
=========================
Downloads the SNAPSI surface field and reduces it to a polar-cap index on the fly,
so the CAUSAL stratosphere-to-surface effect can be measured by experimental
design rather than inferred from composites.

WHY SNAPSI SETTLES SOMETHING NOTHING ELSE CAN
  Every number in this project so far -- the beta x dS projection law, the forced
  variance ceiling -- is inferred from observed or free-running data, where the
  stratospheric anomaly is never manipulated. SNAPSI manipulates it. The
  stratosphere is nudged, so the difference between two nudged ensembles is a
  controlled experiment, not a composite.

THE CONTRAST, AND THE TRAP THAT WAS AVOIDED
  Verified from the protocol paper (Hitchcock et al., GMD 15, 5073, 2022):
    nudged  -- "the zonally symmetric stratospheric state is nudged globally to
                the observed time evolution of the stratospheric event of
                interest", full strength above 50 hPa tapering to none
                below 90 hPa, 6-hour relaxation timescale
    control -- "nudged globally to a time-evolving climatological state"
               (1979-2019 climatology), SAME nudging machinery
    free    -- "the atmosphere evolves freely after initialization", NO nudging

  The causal contrast is therefore **nudged minus control**, NOT nudged minus
  free. Both nudged arms carry identical nudging artefacts, which cancel in the
  difference, leaving only the effect of the stratospheric ANOMALY. Differencing
  against `free` would confound the anomaly with the act of nudging itself, which
  suppresses internal stratospheric variability. An earlier plan for this script
  used free+nudged and would have measured the wrong thing.

  Note also that the six SNAPSI dates are TWO INITIALISATIONS PER EVENT, not six
  events: Feb 2018 (init 25 Jan, 8 Feb), Jan 2019 (init 13 Dec, 8 Jan), and the
  Sep 2019 SOUTHERN hemisphere case (init 29 Aug, 1 Oct). So the NH sample is
  2 events at 2 lead times, never 4 or 6 independent events, and any spread
  between them is reported as such.

WHAT IS DOWNLOADED
  psl (sea level pressure), 6hrPt, 45-60 day forecasts, for all six
  initialisations (four NH, two SH), experiments `nudged` and `control`, across
  every centre that carries both. Each file is reduced IN MEMORY to cos-weighted
  60-90N and 60-90S polar-cap means and discarded; only the reduced series are written, so the
  ~12 GB of transfer leaves a few MB on disk.

AUTH
  CEDA issues an RS256 JWT (accounts.ceda.ac.uk realm). It is read from the .env
  supplied by the user and sent as `Authorization: Bearer <token>`. Tokens are
  short-lived -- the one this was built against had 72 h of life -- so expiry is
  checked up front and the run refuses to start rather than failing halfway.

Output: snapsi_polarcap_psl.parquet
        columns centre, model, experiment, init, member, time, lead_days, psl_cap
"""
import base64
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
import xarray as xr

HERE = pathlib.Path(__file__).resolve().parent
MANIFEST = HERE / "snapsi_manifest.csv"
OUT = HERE / "snapsi_polarcap_psl.parquet"
CACHE = HERE / "_snapsi_reduced"
# Repo-root .env. This was a hardcoded absolute Windows path into an unrelated
# project (global_supply_chain_inflation_analysis), so the script could not run on
# this machine at all. CEDA_ENV overrides it if the token lives elsewhere.
ENV = pathlib.Path(os.environ.get("CEDA_ENV", HERE.parents[1] / ".env"))

# SNAPSI has six initialisations: four NH (two events x two lead times) and two
# SH (the Sep 2019 Antarctic warming). The SH pair is what tests whether the
# result is an NH-specific artefact.
NH_INITS = ["s20180125", "s20180208", "s20181213", "s20190108"]
SH_INITS = ["s20190829", "s20191001"]
ALL_INITS = NH_INITS + SH_INITS

EXPERIMENTS = ["nudged", "control"]

# Every centre in the archive carrying psl for both experiments. NRL is retained
# in the manifest but its submission is defective -- nudged and control members
# differ by at most 0.023 Pa at three of four initialisations -- and the analysis scripts
# exclude it by their own guard, not by this list.
ALL_CENTRES = ["CCCma", "CNR-ISAC", "ECCC", "ECMWF", "KMA", "Meteo-France",
               "NCAR", "NRL", "SNU", "UKMO"]
# The original four: the cheapest grids, NH only. Kept as a named subset so the
# earlier 3-usable-centre result can be reproduced exactly.
PILOT_CENTRES = ["NRL", "CCCma", "SNU", "KMA"]

CENTRES = ALL_CENTRES
INITS = ALL_INITS
CAP_LAT = 60.0
WORKERS = 12


def token():
    txt = ENV.read_text(encoding="utf8", errors="ignore")
    m = re.search(r"^\s*CEDA\s*=\s*['\"]?([A-Za-z0-9_\-\.]+)", txt, re.M)
    if not m:
        sys.exit("no CEDA token in .env")
    tok = m.group(1)
    pay = json.loads(base64.urlsafe_b64decode(
        tok.split(".")[1] + "=" * (-len(tok.split(".")[1]) % 4)))
    left = pay.get("exp", 0) - time.time()
    if left <= 0:
        sys.exit("CEDA token EXPIRED -- refresh it before running")
    print(f"CEDA token valid, {left / 3600:.1f} h remaining "
          f"(user {pay.get('preferred_username', '?')})")
    return tok


def lead_days(tv):
    """Forecast lead in days, WITHOUT converting the calendar.

    CCCma runs a 365-day NoLeap calendar, so its time axis decodes to
    cftime.DatetimeNoLeap objects that pandas refuses to convert ("is not
    convertible to datetime"). Every CCCma file failed on that. Nothing here needs
    absolute wall-clock dates though: lead is a difference within one file, and
    post-onset day is recovered downstream as lead + (init - onset), an offset
    computed once from the two known real dates. So the calendar never has to be
    reconciled at all, which also removes any silent 365-vs-366 day drift.
    """
    t0 = tv[0]
    if np.issubdtype(np.asarray(tv).dtype, np.datetime64):
        return (pd.to_datetime(tv) - pd.to_datetime(t0)).total_seconds() / 86400.0
    return np.array([(x - t0).total_seconds() / 86400.0 for x in tv], dtype=float)


# ------------------------------------------------------------ time origin ---
# `lead_days()` measures lead from the FILE'S FIRST TIME STEP. Every consumer then
# computes post-onset day as lead - (onset - init date), i.e. it assumes lead 0 is
# 00 UTC on the initialisation date. UKMO files start at 06 UTC (their names say
# 201801250600-...), so every UKMO window in K, L, M, N and O sat 6 h early --
# found in review 2026-09-24. The origin is MEASURED per ensemble from one
# member's time axis (`--measure-time-origin`, recorded in TIME_ORIGIN) and every
# cached member is rebased so that lead_days = days since 00 UTC on the init
# date. Rebased frames carry lead_origin == LEAD_ORIGIN; the smoke test checks.
TIME_ORIGIN = HERE / "snapsi_time_origin.json"
LEAD_ORIGIN = "init_00UTC"


def measure_time_origin(tok, sample=3):
    """First time step of each (centre, experiment, init), in hours after 00 UTC
    of the init date, read from `sample` members' files. Refuses a centre whose
    sampled members disagree -- the rebase applies one offset per ensemble."""
    import cftime
    import h5py
    d = pd.read_csv(MANIFEST)
    q = d[(d.variable == "psl") & d.experiment.isin(EXPERIMENTS)]
    rows = q.groupby(["centre", "experiment", "start_date"]).head(sample)

    import fsspec
    fs = fsspec.filesystem("http", client_kwargs={
        "headers": {"Authorization": f"Bearer {tok}"}})

    def first(r):
        # Byte-range read of the time variable only (4 KB blocks): the first
        # version downloaded 240 whole files to read one number from each and
        # starved the concurrent zg fetches of bandwidth for 40 minutes.
        path = r.archive_path if bool(r.verified) else resolve_archive_path(r.archive_path, tok)
        with fs.open(f"{FILE_ROOT}{path}", "rb", block_size=2 ** 12) as fh, \
                h5py.File(fh, "r") as f:
            t = f["time"]
            u = t.attrs["units"]; u = u.decode() if isinstance(u, bytes) else str(u)
            cal = t.attrs.get("calendar", b"standard")
            cal = cal.decode() if isinstance(cal, bytes) else str(cal)
            t0 = cftime.num2date(t[0], u, cal)
        init = pd.Timestamp(str(r.start_date)[1:])
        hrs = (pd.Timestamp(year=t0.year, month=t0.month, day=t0.day, hour=t0.hour,
                            minute=t0.minute) - init).total_seconds() / 3600
        return (r.centre, r.experiment, str(r.start_date)), hrs, str(t0)

    def first_retry(r):
        for attempt in range(6):
            try:
                return first(r)
            except Exception as exc:
                if attempt == 5:
                    raise
                print(f"  retry {attempt + 1} {r.centre}/{r.start_date}/{r.member}: "
                      f"{type(exc).__name__}", flush=True)
                time.sleep(3 * (attempt + 1))

    with ThreadPoolExecutor(8) as ex:
        got = list(ex.map(first_retry, list(rows.itertuples())))
    # Merge into the existing table: measuring one set of experiments (e.g. --full)
    # must not drop the origins already measured for the others.
    table = json.loads(TIME_ORIGIN.read_text()) if TIME_ORIGIN.exists() else {}
    fresh = set()
    for key, hrs, t0 in got:
        k = "|".join(key)
        if k in fresh and table[k]["offset_hours"] != hrs:
            raise ValueError(f"{k}: members start at different times "
                             f"({table[k]['first_time']} vs {t0})")
        table[k] = {"offset_hours": hrs, "first_time": t0}
        fresh.add(k)
    TIME_ORIGIN.write_text(json.dumps(table, indent=1, sort_keys=True),
                           encoding="utf8", newline="\n")
    return table


class MissingTimeOrigin(Exception):
    """An ensemble with no measured lead-0 time. Not a transfer error: fetch
    loops re-raise it immediately instead of retrying the download."""


def origin_days(centre, experiment, init):
    """Measured offset of lead 0 from 00 UTC on the init date, in days."""
    table = json.loads(TIME_ORIGIN.read_text())
    k = f"{centre}|{experiment}|{init}"
    if k not in table:
        raise MissingTimeOrigin(f"no measured time origin for {k}; "
                                f"run --measure-time-origin")
    return table[k]["offset_hours"] / 24.0


def rebase_cache(cache_dir, pattern="*.parquet"):
    """Rebase every cached member not yet marked. Idempotent; atomic per file."""
    n = 0
    for f in sorted(pathlib.Path(cache_dir).glob(pattern)):
        df = pd.read_parquet(f)
        if "lead_origin" in df.columns and (df["lead_origin"] == LEAD_ORIGIN).all():
            continue
        r = df.iloc[0]
        df["lead_days"] = df["lead_days"] + origin_days(r.centre, r.experiment, r.init)
        df["lead_origin"] = LEAD_ORIGIN
        tmp = f.with_suffix(f".parquet.{os.getpid()}.rebase")
        df.to_parquet(tmp)
        tmp.replace(f)
        n += 1
    return n


def cache_ok(df):
    """Is a cached member usable? Cheap structural + physical plausibility test."""
    if df is None or len(df) == 0:
        return False
    need = {"centre", "experiment", "init", "member", "lead_days", "psl_cap"}
    if not need <= set(df.columns):
        return False
    # Entries written before 2026-09-17 have no `hemisphere` column and hold a
    # NORTHERN cap regardless of the event. That is correct for the four NH
    # initialisations and wrong for the two SH ones, so reject those and refetch.
    if "hemisphere" not in df.columns:
        if str(df["init"].iloc[0]) in SH_INITS:
            return False
    v = df["psl_cap"]
    if v.isna().all():
        return False
    # polar-cap sea-level pressure in Pa; anything outside this is not psl
    return 90000.0 < float(v.mean()) < 110000.0


FILE_ROOT = "https://dap.ceda.ac.uk"


def resolve_archive_path(archive_path: str, tok: str) -> str:
    """The path that actually exists, for a manifest row marked verified=False.

    `acquire_snapsi.py` builds the manifest by discovering one member's layout
    and substituting the other member ids, then samples the result. When the
    sample fails it flags the whole node `verified=False` -- and until
    2026-09-23 no consumer read the flag. CNR-ISAC nudged s20190108 is such a
    node: 10 of its 50 members live under version v20230307, not the v20230110
    the manifest constructed, so they 404ed and were silently absent.

    Lists the member's <var>/<grid>/ directory and requires EXACTLY one version
    holding a file of the manifest's name. Anything else raises: guessing
    between versions would mix data vintages without saying so.
    """
    parts = archive_path.rstrip("/").split("/")
    fname, grid_dir = parts[-1], "/".join(parts[:-2]) + "/"

    def listing(path):
        req = urllib.request.Request(f"{FILE_ROOT}{path}",
                                     headers={"Authorization": f"Bearer {tok}"})
        with urllib.request.urlopen(req, timeout=120) as r:
            return r.read().decode("utf8", "replace")

    versions = sorted(set(re.findall(r'href="(v\d{8})/"', listing(grid_dir))))
    hits = [v for v in versions
            if f'href="{fname}"' in listing(f"{grid_dir}{v}/")]
    if len(hits) != 1:
        raise FileNotFoundError(
            f"{fname}: {len(hits)} versions hold it (listed {versions}); refusing to guess")
    return f"{grid_dir}{hits[0]}/{fname}"


def fetch_reduce(row, tok):
    """Download one file, reduce to a polar-cap mean, discard the raw bytes."""
    # CENTRE-FIRST, and it must stay that way. snapsi_selection_test.py and
    # snapsi_distribution_test.py read this directory by filename --
    # `p.name.split("_")[0]` for the centre and `{centre}_{exp}_{init}_*` to
    # load -- so the cache name is a PUBLIC CONTRACT, not an internal detail.
    # Prefixing it with the variable on 2026-09-17 broke both silently: they
    # parsed "psl" as a centre name and quietly analysed only the 2,026
    # legacy-named members while ignoring 4,135 new ones. Other variables get
    # their own DIRECTORY (see acquire_snapsi_zg.py), never a prefix here.
    key = f"{row.centre}_{row.experiment}_{row.start_date}_{row.member}"
    cf = CACHE / f"{key}.parquet"

    if cf.exists():
        # A cached file is NOT trusted on existence alone. A process killed
        # mid-write leaves a truncated parquet that exists() happily accepts and
        # that would then be believed forever. Validate, and re-fetch if bad.
        try:
            cached = pd.read_parquet(cf)
        except Exception:
            cached = None
        if cache_ok(cached):
            return cached
        print(f"  cache REJECTED, refetching: {cf.name}", flush=True)
        cf.unlink(missing_ok=True)
        cf = CACHE / f"{key}.parquet"

    url = row.download_url
    if not bool(row.verified):
        url = FILE_ROOT + resolve_archive_path(row.archive_path, tok)
    if not url.endswith("?download=1"):
        url = url.rstrip("/") + "?download=1"
    last = None
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"Authorization": f"Bearer {tok}"})
            with urllib.request.urlopen(req, timeout=300) as r:
                blob = r.read()
            if len(blob) < 50_000:
                last = f"short file ({len(blob)} B)"
                continue
            with xr.open_dataset(io.BytesIO(blob)) as ds:
                v = ds["psl"]
                lat = ds["lat"]
                # BOTH caps, always. The original code took lat >= 60 for every
                # initialisation, which is the wrong hemisphere for the two SH
                # cases (s20190829, s20191001 -- the Sep 2019 Antarctic warming)
                # and silently gave them a NORTHERN cap. Computing both costs
                # nothing here and means the raw file never has to be fetched
                # again to answer a question about the other hemisphere.
                north = v.sel(lat=lat[lat >= CAP_LAT])
                south = v.sel(lat=lat[lat <= -CAP_LAT])
                cap_n = north.weighted(np.cos(np.deg2rad(north.lat))).mean(
                    dim=["lat", "lon"]).values
                cap_s = south.weighted(np.cos(np.deg2rad(south.lat))).mean(
                    dim=["lat", "lon"]).values
                lead = lead_days(ds.time.values)
            if not (np.isfinite(cap_n).any() and np.isfinite(cap_s).any()):
                raise ValueError("all-NaN polar cap")
            hemi = "S" if str(row.start_date) in SH_INITS else "N"
            out = pd.DataFrame({
                "centre": row.centre, "model": row.model,
                "experiment": row.experiment, "init": str(row.start_date),
                "member": row.member,
                "lead_days": lead + origin_days(row.centre, row.experiment,
                                                str(row.start_date)),
                "lead_origin": LEAD_ORIGIN,
                "psl_cap_N": cap_n, "psl_cap_S": cap_s,
                # the cap of the hemisphere the event is IN -- what an analysis
                # of that event should use
                "hemisphere": hemi,
                "psl_cap": cap_s if hemi == "S" else cap_n})
            # Atomic: write beside the target, then rename. A kill during
            # to_parquet() would otherwise leave a truncated file that the next
            # run treats as a completed download.
            tmp = cf.with_suffix(".parquet.tmp")
            out.to_parquet(tmp)
            tmp.replace(cf)
            return out
        except MissingTimeOrigin:
            raise
        except (ImportError, ModuleNotFoundError) as exc:
            # Not transient. Retrying a missing backend just re-downloads the
            # file four times and reports it as a network problem: this cost a
            # 5 GB pilot on 2026-09-17 (h5py absent -> every reduce failed).
            raise SystemExit(
                f"missing dependency, not a transfer error: {exc}\n"
                f"install it into the project venv and re-run; "
                f"check with environment/check_environment.py") from exc
        except Exception as exc:
            last = exc
            if attempt == 0:
                print(f"  retrying {row.centre}/{row.start_date}/{row.member}: "
                      f"{type(exc).__name__}: {str(exc)[:120]}", flush=True)
            time.sleep(2 * (attempt + 1))
    raise IOError(f"{row.centre}/{row.experiment}/{row.start_date}/{row.member}: {last}")


def main():
    """
    Usage:
      acquire_snapsi_surface.py                      every centre, every init
      acquire_snapsi_surface.py --centres UKMO,NCAR  only those centres
      acquire_snapsi_surface.py --inits sh           only the two SH inits
      acquire_snapsi_surface.py --dry-run            print the budget, transfer nothing
      acquire_snapsi_surface.py --measure-time-origin  record each ensemble's lead-0 time
      acquire_snapsi_surface.py --rebase             rebase cached leads, then continue
      acquire_snapsi_surface.py --full --centres ECMWF,UKMO,Meteo-France --inits nh
                                                     the nudged-full / control-full arms, into
                                                     their OWN cache (_snapsi_reduced_full) and
                                                     output (snapsi_polarcap_psl_full.parquet),
                                                     so no consumer of the zonal-nudging cache
                                                     can pick them up (FAILURES 2026-09-17)

    Already-reduced members are served from _snapsi_reduced/ and cost no
    transfer, so re-running to widen the scope only fetches what is new.
    """
    argv = sys.argv[1:]
    if "--full" in argv:
        global EXPERIMENTS, CACHE, OUT
        EXPERIMENTS = ["nudged-full", "control-full"]
        CACHE = HERE / "_snapsi_reduced_full"
        OUT = HERE / "snapsi_polarcap_psl_full.parquet"

    def opt(name, default=None):
        if name in argv:
            return argv[argv.index(name) + 1]
        return default

    centres = opt("--centres")
    centres = ([c.strip() for c in centres.split(",")] if centres
               else (PILOT_CENTRES if "--pilot" in argv else CENTRES))
    inits_arg = (opt("--inits") or "all").lower()
    inits = {"all": ALL_INITS, "nh": NH_INITS, "sh": SH_INITS}.get(
        inits_arg, [i.strip() for i in inits_arg.split(",")])

    CACHE.mkdir(exist_ok=True)
    if "--measure-time-origin" in argv:
        t = measure_time_origin(token())
        print(json.dumps({k: v["offset_hours"] for k, v in t.items()}, indent=1))
        return 0
    if "--rebase" in argv:
        print(f"rebased {rebase_cache(CACHE)} psl cache file(s)")
    d = pd.read_csv(MANIFEST)
    q = d[(d.variable == "psl")
          & (d.experiment.isin(EXPERIMENTS))
          & (d.start_date.astype(str).isin(inits))
          & (d.centre.isin(centres))].reset_index(drop=True)

    # Separate what must be transferred from what is already reduced on disk,
    # so the printed budget is the real one rather than the notional total.
    def cached(r):
        return (CACHE / f"{r.centre}_{r.experiment}_{r.start_date}_{r.member}.parquet").exists()

    have = q.apply(cached, axis=1) if len(q) else pd.Series(dtype=bool)
    todo = q[~have] if len(q) else q
    print(f"centres    {centres}")
    print(f"inits      {inits}")
    print(f"matched    {len(q):,} files, {q.size_bytes_est.sum()/1e9:.1f} GB total")
    print(f"cached     {int(have.sum()):,} files (no transfer)")
    print(f"to fetch   {len(todo):,} files, "
          f"{todo.size_bytes_est.sum()/1e9:.1f} GB", flush=True)
    if "--dry-run" in argv:
        by = todo.groupby("centre").size_bytes_est.agg(["count", "sum"])
        by["GB"] = (by["sum"] / 1e9).round(1)
        print(by[["count", "GB"]].to_string())
        return 0
    if todo.empty:
        print("nothing to fetch; rebuilding the parquet from cache")
    tok = token()
    q = q.reset_index(drop=True)

    rows, failed = [], []
    t0 = time.time()
    with ThreadPoolExecutor(WORKERS) as ex:
        futs = {ex.submit(fetch_reduce, r, tok): i for i, r in enumerate(q.itertuples())}
        for n, f in enumerate(as_completed(futs), 1):
            try:
                rows.append(f.result())
            except Exception as exc:
                failed.append(exc)
                print(f"  FAILED {exc}", flush=True)
            if n % 100 == 0:
                el = time.time() - t0
                print(f"  {n}/{len(q)}  {el/60:.1f} min  "
                      f"({n/el*60:.0f} files/min)", flush=True)

    # Rebuild the parquet from the WHOLE cache, never from this run's selection.
    # Writing only the current subset made the output scope-dependent: running
    # `--centres KMA,SNU --inits sh` overwrote the file with SH data alone and
    # silently dropped every NH member the downstream analyses read. The cache is
    # the source of truth; this file is a materialised view of it.
    # psl members only. The cache is now keyed by variable, so once zg or ta
    # members land beside these a bare *.parquet glob would concatenate
    # different fields into one column. Legacy entries have no variable prefix
    # and are psl by construction -- they predate the key change.
    # This directory holds psl and nothing else; other variables live elsewhere.
    rb = rebase_cache(CACHE)
    if rb:
        print(f"rebased {rb} cached member(s) to lead 0 = 00 UTC on the init date")
    cached_files = sorted(CACHE.glob("*.parquet"))
    stray = list(CACHE.glob("*.parquet.tmp"))
    for s in stray:
        s.unlink(missing_ok=True)
    if stray:
        print(f"removed {len(stray)} interrupted partial write(s)")
    if not cached_files:
        sys.exit("nothing downloaded and cache is empty")

    # Validate on the way in. A cached file that fetch_reduce() would reject
    # must not reach the output either -- otherwise an interrupted run leaves
    # wrong-hemisphere SH members concatenated into the parquet the analyses
    # read, which is exactly how the NH-cap-for-SH-events bug propagated.
    frames, dropped = [], 0
    for f in cached_files:
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
        print(f"excluded {dropped} cached member(s) that failed validation "
              f"(re-run to refetch them)")
    if not frames:
        sys.exit("every cached member failed validation")
    out = pd.concat(frames, ignore_index=True)
    out.to_parquet(OUT)
    print(f"\n{len(out):,} rows from {len(cached_files):,} cached members "
          f"-> {OUT.name}  ({len(failed)} files failed this run)")
    print(out.groupby(["centre", "experiment", "init"]).member.nunique()
          .unstack(fill_value=0).to_string())


if __name__ == "__main__":
    main()
