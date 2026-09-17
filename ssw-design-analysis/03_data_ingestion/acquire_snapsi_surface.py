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
                interest", 50-90 hPa, 6-hourly relaxation
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
  psl (sea level pressure), 6hrPt, 45-day forecasts, global 1 deg, for the four
  NH initialisations, experiments `nudged` and `control`, across the cheapest
  centres that carry both. Each file is reduced IN MEMORY to a cos-weighted
  60-90N polar-cap mean and discarded; only the reduced series are written, so the
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
# in the manifest but its submission is corrupt -- nudged and control are
# byte-identical at three of four initialisations -- and the analysis scripts
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


def fetch_reduce(row, tok):
    """Download one file, reduce to a polar-cap mean, discard the raw bytes."""
    key = f"{row.centre}_{row.experiment}_{row.start_date}_{row.member}"
    cf = CACHE / f"{key}.parquet"
    if cf.exists():
        return pd.read_parquet(cf)
    url = row.download_url
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
                v = v.sel(lat=lat[lat >= CAP_LAT])
                w = np.cos(np.deg2rad(v.lat))
                cap = v.weighted(w).mean(dim=["lat", "lon"]).values
                lead = lead_days(ds.time.values)
            if not np.isfinite(cap).any():
                raise ValueError("all-NaN polar cap")
            out = pd.DataFrame({
                "centre": row.centre, "model": row.model,
                "experiment": row.experiment, "init": str(row.start_date),
                "member": row.member, "lead_days": lead, "psl_cap": cap})
            out.to_parquet(cf)          # so a mid-run failure costs nothing
            return out
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

    Already-reduced members are served from _snapsi_reduced/ and cost no
    transfer, so re-running to widen the scope only fetches what is new.
    """
    argv = sys.argv[1:]

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
    cached_files = sorted(CACHE.glob("*.parquet"))
    if not cached_files:
        sys.exit("nothing downloaded and cache is empty")
    out = pd.concat([pd.read_parquet(f) for f in cached_files], ignore_index=True)
    out.to_parquet(OUT)
    print(f"\n{len(out):,} rows from {len(cached_files):,} cached members "
          f"-> {OUT.name}  ({len(failed)} files failed this run)")
    print(out.groupby(["centre", "experiment", "init"]).member.nunique()
          .unstack(fill_value=0).to_string())


if __name__ == "__main__":
    main()
