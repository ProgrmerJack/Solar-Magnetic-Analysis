#!/usr/bin/env python3
"""
acquire_cmip6_ensemble.py
=========================
Acquires a large ensemble of daily stratosphere + surface data, to break the
binding constraint on every open question in this project.

THE CONSTRAINT BEING BROKEN
  Everything unresolved here is unresolved for one reason: 47 winters of
  observations is one realisation of history, and it does not contain enough
  independent information.
    - the pre-onset AO anomaly: -0.82 on 43 events, loses significance under
      every restriction, but every restricted interval still CONTAINS -0.82
    - stratospheric mediation: 45-54%, disattenuated CI spans -11% to +140%
  No estimator fixes that. Only more independent winters do.

WHY THIS SOURCE
  SNAPSI is the ideal experiment but sits behind CEDA authentication. CMIP6 daily
  output is on Google Cloud as consolidated zarr with NO authentication at all,
  and daily `ua` is archived on plev8 = 1000/850/700/500/250/100/50/**10 hPa** --
  10 hPa is present, which is exactly the level Charlton-Polvani SSW detection
  needs. `psl` gives the surface annular mode.

  CanESM5 has 35 members with daily ua at 64x128, the smallest arrays and so the
  fastest to pull, and is a SNAPSI participating model. Each member is 165 years,
  so ten members give ~1,650 winters against 47 observed.

WHAT IS EXTRACTED, AND WHY IT IS SMALL
  The zarr chunks span all levels/lats/lons, so the full array must be READ even
  though only a slice is wanted. But what is KEPT is tiny: zonal means.
      u_zm    (time, plev8, lat)   -> vortex at 10 hPa 60N, and the descent profile
      psl_zm  (time, lat)          -> annular mode index from the meridional dipole
  ~140 MB retained per member against ~18 GB read.

  Reads are threaded because a single stream gives 16 MB/s and 24 threads give
  55 MB/s, which is the difference between 17 minutes and 5 minutes per member.

Outputs (this directory):
  raw/cmip6/<model>_<member>_zm.nc     zonal-mean fields per member
  cmip6_ensemble_manifest.csv          what was fetched, with sizes and timings
"""
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw" / "cmip6"
MANIFEST = HERE / "cmip6_ensemble_manifest.csv"

CATALOG = "https://storage.googleapis.com/cmip6/cmip6-zarr-consolidated-stores.csv"
# usage: acquire_cmip6_ensemble.py [n_members] [model]
# A second, independent model is not optional for the headline claim: one model's
# stratosphere-troposphere coupling strength is a property of that model.
N_MEMBERS = int(sys.argv[1]) if len(sys.argv) > 1 else 10
MODEL = sys.argv[2] if len(sys.argv) > 2 else "CanESM5"
EXPERIMENT = "historical"
N_THREADS = 24
# Throughput to Google Cloud is capped PER CONNECTION, not in aggregate: one
# stream gives 0.8 MB/s while 4 parallel streams give 27 MB/s. Concurrency IS the
# bandwidth here. An earlier version cut threads to 9 on large-chunk stores to
# bound memory and threw the bandwidth away with it. Observed working set at 24
# threads on 157 MB chunks was 2.4 GB, which is fine, so the cap is set well above
# that and only guards genuinely pathological chunk sizes.
IN_FLIGHT_MB = 6000


def zonal_mean_series(store, var, fs):
    """Read a whole zarr array in chunk blocks, keeping only the zonal mean."""
    import zarr
    g = zarr.open_group(fs.get_mapper(store), mode="r")
    a = g[var]
    C = a.chunks[0]
    nblk = int(np.ceil(a.shape[0] / C))

    # CMIP6 arrays carry a _FillValue (typically 1e20) wherever data is missing --
    # for pressure-level fields, below-ground points and, in some models, whole
    # regions. A plain .mean() averages that sentinel in and destroys the result:
    # MIROC6 came back with a 10 hPa 60N "wind" of 2.6e18 m/s and 26% of its
    # (time, plev, lat) cells corrupted, across every level and latitude, with no
    # clean timestep anywhere. CanESM5 happened to have no fill and was unaffected,
    # which is exactly why this went unnoticed for ten members.
    # Mask before averaging, and use nanmean so partially-filled rows degrade
    # gracefully instead of exploding.
    # Kept deliberately cheap: this runs over every element of a ~20 GB array and
    # the first version was CPU-bound, not network-bound (1048 CPU-seconds in 17
    # minutes on one core, stalling a 15.8 GB model). Two np.isclose scans against
    # _FillValue/missing_value were the bulk of it and are redundant -- CMIP6
    # sentinels are ~1e20 and are already caught by the magnitude test. Take the
    # fast path when a chunk has no missing data, which is the common case.
    def grab(i):
        blk = np.asarray(a[i * C:(i + 1) * C], dtype="float32")
        bad = ~np.isfinite(blk)
        np.logical_or(bad, np.abs(blk) > 1e10, out=bad)
        if not bad.any():
            return blk.mean(axis=-1)                             # fast path
        blk = np.where(bad, np.nan, blk)
        return np.nanmean(blk, axis=-1)                          # mean over lon

    # Bound in-flight memory instead of using a fixed thread count. Chunk size
    # varies 3x across models (CanESM5 53.7 MB, BCC-ESM1 157.3 MB) and 24 threads
    # on the large-chunk stores puts ~3.8 GB in flight, which pages and shows up
    # as CPU: BCC-ESM1 sat at two cores pegged and made no progress in 17 minutes
    # while CanESM5 ran at 35-42 MB/s. Same dtype and compressor in both, so chunk
    # size was the only difference.
    chunk_mb = np.prod(a.chunks) * a.dtype.itemsize / 1e6
    nthreads = int(np.clip(IN_FLIGHT_MB / max(chunk_mb, 1), 4, N_THREADS))
    print(f"    chunks {chunk_mb:.0f} MB x {nblk} -> {nthreads} threads", flush=True)
    with ThreadPoolExecutor(nthreads) as ex:
        parts = list(ex.map(grab, range(nblk)))
    out = np.concatenate(parts, axis=0)
    coords = {k: np.asarray(g[k][:]) for k in g.array_keys()
              if k in ("time", "plev", "lat")}
    # Decode CF time. CanESM5 uses calendar="noleap", so the raw values are
    # integer day counts in a 365-day calendar and casting them to datetime64
    # silently yields nanoseconds-since-epoch -- which is what broke the first
    # two members (every date collapsed into 1970, giving 1 "winter" per member).
    ta = dict(g["time"].attrs)
    coords["time"] = decode_time(coords["time"], ta.get("units", ""),
                                 ta.get("calendar", "standard"))
    return out, coords, a.nbytes


def decode_time(vals, units, calendar):
    """CF time -> pandas Timestamps, via cftime so non-standard calendars work.

    A noleap date is always a valid Gregorian date (it simply omits 29 Feb), so
    mapping (year, month, day) across is lossless for day-of-year and month
    logic. Day differences within a winter are unaffected.
    """
    import cftime
    d = cftime.num2date(np.asarray(vals), units, calendar=calendar)
    return pd.to_datetime([f"{x.year:04d}-{x.month:02d}-{x.day:02d}" for x in d])


def main():
    import gcsfs
    RAW.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(CATALOG)
    fs = gcsfs.GCSFileSystem(token="anon")

    ua = df[(df.table_id == "day") & (df.variable_id == "ua")
            & (df.source_id == MODEL) & (df.experiment_id == EXPERIMENT)]
    psl = df[(df.table_id == "day") & (df.variable_id == "psl")
             & (df.source_id == MODEL) & (df.experiment_id == EXPERIMENT)]
    members = sorted(set(ua.member_id) & set(psl.member_id))[:N_MEMBERS]
    print(f"{MODEL} {EXPERIMENT}: {len(members)} members with both ua and psl daily")
    print(f"  {members}\n")

    rows = []
    for i, mem in enumerate(members, 1):
        dest = RAW / f"{MODEL}_{mem}_zm.nc"
        if dest.exists():
            print(f"[{i}/{len(members)}] {mem}: cached")
            continue
        t0 = time.time()
        zu = ua[ua.member_id == mem].iloc[0].zstore
        zp = psl[psl.member_id == mem].iloc[0].zstore
        try:
            u, cu, nb_u = zonal_mean_series(zu, "ua", fs)
            p, cp, nb_p = zonal_mean_series(zp, "psl", fs)
        except Exception as e:
            print(f"[{i}/{len(members)}] {mem}: FAILED {type(e).__name__} {str(e)[:70]}")
            rows.append({"member": mem, "status": f"FAILED {type(e).__name__}"})
            continue
        dt = time.time() - t0

        import xarray as xr
        n = min(u.shape[0], p.shape[0])
        ds = xr.Dataset(
            {"u_zm": (("time", "plev", "lat"), u[:n]),
             "psl_zm": (("time", "lat"), p[:n])},
            coords={"time": cu["time"][:n], "plev": cu["plev"], "lat": cu["lat"]})
        ds.attrs["source"] = f"{MODEL} {EXPERIMENT} {mem} (CMIP6 GC zarr, anon)"
        ds.to_netcdf(dest)
        gb = (nb_u + nb_p) / 1e9
        print(f"[{i}/{len(members)}] {mem}: {n:,} days, read {gb:.1f} GB in "
              f"{dt/60:.1f} min ({gb*1000/dt:.0f} MB/s) -> {dest.name} "
              f"({dest.stat().st_size/1e6:.0f} MB)", flush=True)
        rows.append({"member": mem, "status": "ok", "n_days": n,
                     "gb_read": round(gb, 2), "minutes": round(dt / 60, 2),
                     "file": dest.name})

    if rows:
        pd.DataFrame(rows).to_csv(MANIFEST, index=False, lineterminator="\n")
    have = sorted(RAW.glob(f"{MODEL}_*_zm.nc"))
    print(f"\n{len(have)} members on disk, "
          f"{sum(f.stat().st_size for f in have)/1e9:.2f} GB retained")
    print(f"  ~{len(have) * 165:,} model winters available "
          f"(observations provide 47)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
