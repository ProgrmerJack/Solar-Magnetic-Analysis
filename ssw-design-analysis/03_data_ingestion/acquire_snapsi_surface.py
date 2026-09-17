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
ENV = pathlib.Path(r"C:\Users\Jack0\GitHub\global_supply_chain_inflation_analysis\.env")

NH_INITS = ["s20180125", "s20180208", "s20181213", "s20190108"]
EXPERIMENTS = ["nudged", "control"]
CENTRES = ["NRL", "CCCma", "SNU", "KMA"]      # cheapest carrying both experiments
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
        except Exception as exc:
            last = exc
            time.sleep(2 * (attempt + 1))
    raise IOError(f"{row.centre}/{row.experiment}/{row.start_date}/{row.member}: {last}")


def main():
    tok = token()
    CACHE.mkdir(exist_ok=True)
    d = pd.read_csv(MANIFEST)
    q = d[(d.variable == "psl")
          & (d.experiment.isin(EXPERIMENTS))
          & (d.start_date.astype(str).isin(NH_INITS))
          & (d.centre.isin(CENTRES))].reset_index(drop=True)
    gb = q.size_bytes_est.sum() / 1e9
    print(f"{len(q):,} files, {gb:.2f} GB to transfer, reduced on the fly")
    print(f"centres {CENTRES}, experiments {EXPERIMENTS}, inits {NH_INITS}")

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

    if not rows:
        sys.exit("nothing downloaded")
    out = pd.concat(rows, ignore_index=True)
    out.to_parquet(OUT)
    print(f"\n{len(out):,} rows -> {OUT.name}  ({len(failed)} files failed)")
    print(out.groupby(["centre", "experiment", "init"]).member.nunique()
          .unstack(fill_value=0).to_string())


if __name__ == "__main__":
    main()
