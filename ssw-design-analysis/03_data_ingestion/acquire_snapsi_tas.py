#!/usr/bin/env python3
"""
acquire_snapsi_tas.py
=====================
SNAPSI near-surface air temperature (tas, 6hrPt), reduced on the fly to DAILY
MEANS on a common 2.5-degree grid north of 20N, for the regional cold-risk
analysis (plan approved 2026-09-26).

WHY A GENERIC FIELD, NOT REGIONAL MEANS
  The regions and cold thresholds are fixed in the analysis plan from the
  literature, BEFORE any regional result is seen. Reducing to a coarse field
  rather than to boxes keeps the download independent of that choice, so the
  plan can be written after the (token-limited) transfer without the data having
  steered it. 29 lat x 144 lon x ~60 days x 4 B is ~1 MB per member.

WHAT IS DOWNLOADED
  tas for the four NH initialisations, experiments nudged and control, every
  centre in acquire_snapsi_surface.ALL_CENTRES except NRL (its nudged and
  control members are identical to < 0.03 Pa in psl; excluded by every analysis).
  Whole files are fetched and discarded after reduction. On 2026-09-26 CEDA served
  ~0.06 MB/s per connection (psl had run at ~15 MB/s in total); CEDA honours
  byte ranges (HTTP 206), and 32 parallel ranges on one file gave 1.5 MB/s, so
  each file is fetched as RANGES parallel ranges, each checked for status 206
  and exact length, the whole against the Content-Range total.

DUPLICATE STEPS
  UKMO s20181213 files repeat 2019-01-01 00 UTC. Repeated steps are dropped (the
  first copy kept). The 100 UKMO s20181213 members cached before this fix carry
  the repeat in lead day 19 (5 steps), outside every analysis window; consumers
  require exactly 4 steps on each window day.

LEAD TIME
  Measured from each file's OWN time stamps against 00 UTC on the initialisation
  date (start_date sYYYYMMDD), in the file's own calendar -- not from the psl
  time-origin table: CCCma tas files start at 03 UTC, so the two variables need
  not share an origin. Day d holds the mean of the steps with lead in [d, d+1);
  the count per day is kept (4 for 6-hourly data) and consumers require it.

REGRIDDING
  Linear interpolation of each 6-hourly field to the 2.5-degree grid (lat 20..90,
  lon 0..357.5), after sorting latitude and wrapping longitude. Adequate for means
  over regions of several hundred km; not for grid-point extremes.

Output cache: _snapsi_tas/<centre>_<experiment>_<init>_<member>.npz
  tas (day, lat, lon) float32 K; n_steps (day); lead_day (day); lat; lon;
  first_lead_hours (lead of the first step, hours after 00 UTC of the init date)
"""
import io
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import acquire_snapsi_surface as S                  # noqa: E402  (token, manifest, paths)

CACHE = HERE / "_snapsi_tas"
CENTRES = [c for c in S.ALL_CENTRES if c != "NRL"]
INITS = S.NH_INITS
EXPERIMENTS = ["nudged", "control"]
LAT = np.arange(20.0, 90.0 + 1e-9, 2.5)
LON = np.arange(0.0, 360.0, 2.5)
WORKERS = int(os.environ.get("SNAPSI_WORKERS", 8))
RANGES = int(os.environ.get("SNAPSI_RANGES", 24))


def ranged(url, tok):
    import urllib.request
    h = {"Authorization": f"Bearer {tok}"}
    with urllib.request.urlopen(urllib.request.Request(url, headers={**h, "Range": "bytes=0-0"}),
                                timeout=120) as x:
        if x.status != 206:
            raise IOError(f"range probe HTTP {x.status}")
        size = int(x.headers["Content-Range"].split("/")[1])
    step = -(-size // RANGES)

    def get(i):
        a, b = i * step, min(size - 1, (i + 1) * step - 1)
        for attempt in range(6):
            try:
                with urllib.request.urlopen(urllib.request.Request(
                        url, headers={**h, "Range": f"bytes={a}-{b}"}), timeout=900) as x:
                    data = x.read()
                    if x.status != 206 or len(data) != b - a + 1:
                        raise IOError(f"HTTP {x.status}, {len(data)} of {b - a + 1} B")
                    return data
            except Exception:
                if attempt == 5:
                    raise
                time.sleep(3 * (attempt + 1))
    with ThreadPoolExecutor(RANGES) as ex:
        blob = b"".join(ex.map(get, range(-(-size // step))))
    if len(blob) != size:
        raise IOError(f"assembled {len(blob)} B != {size}")
    return blob


def lead_hours(times, start_date):
    """Hours after 00 UTC of the init date, in the file's own calendar."""
    y, m, d = int(start_date[1:5]), int(start_date[5:7]), int(start_date[7:9])
    t = np.asarray(times)
    if np.issubdtype(t.dtype, np.datetime64):
        return (t - np.datetime64(f"{y:04d}-{m:02d}-{d:02d}T00:00")) / np.timedelta64(1, "h")
    import cftime
    t0 = cftime.datetime(y, m, d, calendar=t[0].calendar)
    return np.array([(x - t0).total_seconds() / 3600.0 for x in t])


def reduce(ds, start_date):
    v = ds["tas"]
    lat, lon = v["lat"].values, v["lon"].values
    v = v.assign_coords(lon=np.mod(lon, 360.0)).sortby("lat").sortby("lon")
    v = v.sel(lat=slice(17.0, 90.0))
    # wrap one column each side so interpolation across 0/360 is defined
    v = xr.concat([v.isel(lon=[-1]).assign_coords(lon=v.lon[-1:].values - 360.0), v,
                   v.isel(lon=[0]).assign_coords(lon=v.lon[:1].values + 360.0)], dim="lon")
    g = v.interp(lat=LAT, lon=LON).values.astype(np.float32)        # (time, lat, lon)
    if np.isnan(g[:, -1, :]).all():                                  # grid without a pole row
        g[:, -1, :] = np.nanmean(g[:, -2, :], axis=1)[:, None]
    h = lead_hours(ds["time"].values, start_date)
    # UKMO s20181213 repeats 2019-01-01 00 UTC (year boundary; found 2026-09-26):
    # keep the first copy of any repeated step
    keep = np.r_[True, np.diff(h) > 0]
    g, h = g[keep], h[keep]
    day = np.floor(h / 24.0 + 1e-9).astype(int)
    days = np.unique(day[day >= 0])
    out = np.stack([np.nanmean(g[day == k], axis=0) for k in days]).astype(np.float32)
    n = np.array([(day == k).sum() for k in days])
    return out, n, days, float(h[0])


def plausible(a):
    m = float(np.nanmean(a))
    return 200.0 < m < 300.0 and np.isfinite(a).mean() > 0.99


def fetch(row, tok):
    cf = CACHE / f"{row.centre}_{row.experiment}_{row.start_date}_{row.member}.npz"
    if cf.exists():
        try:
            with np.load(cf) as z:
                if plausible(z["tas"]):
                    return "cached"
        except Exception:
            pass
        cf.unlink(missing_ok=True)
    url = row.download_url
    if not bool(row.verified):
        url = S.FILE_ROOT + S.resolve_archive_path(row.archive_path, tok)
    if not url.endswith("?download=1"):
        url = url.rstrip("/") + "?download=1"
    last = None
    for attempt in range(4):
        try:
            blob = ranged(url, tok)
            if len(blob) < 50_000:
                last = f"short file ({len(blob)} B)"
                continue
            with xr.open_dataset(io.BytesIO(blob)) as ds:
                tas, n, days, h0 = reduce(ds, str(row.start_date))
            if not plausible(tas):
                raise ValueError(f"implausible tas mean {np.nanmean(tas):.1f} K")
            tmp = cf.with_name(cf.name + ".tmp.npz")
            np.savez_compressed(tmp, tas=tas, n_steps=n, lead_day=days, lat=LAT, lon=LON,
                                first_lead_hours=h0)
            tmp.replace(cf)
            return "fetched"
        except (ImportError, ModuleNotFoundError) as exc:
            raise SystemExit(f"missing dependency, not a transfer error: {exc}") from exc
        except Exception as exc:
            last = exc
            time.sleep(2 * (attempt + 1))
    raise IOError(f"{row.centre}/{row.experiment}/{row.start_date}/{row.member}: {last}")


def main():
    argv = sys.argv[1:]
    d = pd.read_csv(S.MANIFEST)
    q = d[(d.variable == "tas") & d.experiment.isin(EXPERIMENTS)
          & d.start_date.astype(str).isin(INITS) & d.centre.isin(CENTRES)].reset_index(drop=True)
    if "--centres" in argv:
        q = q[q.centre.isin(argv[argv.index("--centres") + 1].split(","))]
    if "--limit" in argv:
        q = q.groupby("centre").head(int(argv[argv.index("--limit") + 1]))
    CACHE.mkdir(exist_ok=True)
    have = q.apply(lambda r: (CACHE / f"{r.centre}_{r.experiment}_{r.start_date}_{r.member}.npz")
                   .exists(), axis=1)
    print(f"matched {len(q):,} files, {q.size_bytes_est.sum()/1e9:.1f} GB; cached {int(have.sum()):,}; "
          f"to fetch {int((~have).sum()):,} ({q[~have].size_bytes_est.sum()/1e9:.1f} GB)", flush=True)
    if "--dry-run" in argv:
        return 0
    tok = S.token()
    t0, done, failed = time.time(), 0, []
    with ThreadPoolExecutor(WORKERS) as ex:
        futs = {ex.submit(fetch, r, tok): r for r in q.itertuples()}
        for n, f in enumerate(as_completed(futs), 1):
            try:
                f.result(); done += 1
            except Exception as exc:
                failed.append(str(exc)); print(f"  FAILED {exc}", flush=True)
            if n % 100 == 0:
                el = time.time() - t0
                print(f"  {n}/{len(q)}  {el/60:.1f} min  ({n/el*60:.0f} files/min)", flush=True)
    print(f"done {done}, failed {len(failed)} of {len(q)}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
