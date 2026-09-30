#!/usr/bin/env python3
"""
extend_ncep_presatellite.py
===========================
Extends the NCEP/NCAR stratospheric polar-cap record back from 1979 to 1958 to
match the 100 hPa eddy heat flux record, taking the observational event count for
every stratospheric stratifier from 29 to ~42.

WHY THIS IS NOT FREE POWER, AND IS TESTED RATHER THAN ASSUMED
  `scripts/download/27_download_ncep_stratosphere.py` starts at 1979 and that was
  very likely deliberate: 1979 is the start of the satellite era. NCEP/NCAR R1
  stratospheric analyses before it rest on radiosondes with almost no Arctic
  coverage above 100 hPa, and 10 hPa fields in particular are widely regarded as
  unreliable in that period.

  This project has ALREADY been bitten on exactly this seam. The original event
  catalogue rule (">=4 of 6 reanalyses") demanded unanimity pre-1979 because only
  a few products cover it, producing era-dependent event selection at Fisher
  P=0.011. That is why the catalogue rule is now ">=2/3 of *covering* reanalyses".

  So the extra 13 events are acquired here, but the decomposition is then run
  SEPARATELY on 1958-1978 and 1979-2024. If the two eras disagree, the pooled
  number is not used and the observational arm stays at n=29 with that stated.

TRANSPORT: THREDDS NetcdfSubset (NCSS), not OPeNDAP and not whole files.
  OPeNDAP returns SILENT ZEROS here for large aggregate requests (see
  `_reduce` docstring). Whole-file downloads are correct but cost ~8 GB for
  21 years x 3 variables, because every file is global on all 17 levels.
  NCSS serves one level over the 60-90N cap as a real netCDF file: 2.6 MB and
  ~5 s per variable-year-level, and it is BIT-IDENTICAL to the whole-file route
  (validated on uwnd.1958 at 10 hPa: max|diff| = 0.0, corr = 1.00000000).

  Only the REDUCED series are cached, so a rerun is instant and costs no disk.

Output: ncep_presatellite_extension.parquet, holding the pre-1979 rows with the
        polar-cap means matching the target schema plus the 60N zonal means
        needed to detect pre-satellite SSW onsets. The merge into
        data/processed/atmospheric/ncep_stratosphere.parquet stays a separate,
        reversible step.
"""
import io
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TARGET = ROOT / "data" / "processed" / "atmospheric" / "ncep_stratosphere.parquet"
OUT = HERE / "ncep_presatellite_extension.parquet"

NCSS = ("https://psl.noaa.gov/thredds/ncss/grid/Datasets/"
        "ncep.reanalysis.dailyavgs/pressure/{var}.{year}.nc")
VARS = {"air": "air_K", "uwnd": "uwnd_ms", "hgt": "hgt_m"}
YEARS = range(1958, 1979)
LEVELS = [100, 70, 50, 30, 20, 10]
POLAR_MIN_LAT = 60.0
CACHE = HERE / "_ncep_reduced"
WORKERS = 4

# Plausible winter ranges for a 60-90N polar-cap mean, used to REFUSE bad data
# rather than average it into the record.
SANE = {"uwnd_ms": (-40.0, 80.0), "air_K": (180.0, 280.0), "hgt_m": (10.0, 32000.0)}


def _fetch(var, year, level):
    """One NCSS request -> raw netCDF bytes, with backoff.

    PSL's host intermittently DNS-fails under a fast sequential loop; the first
    attempt at this task lost every year from 1959 to 1978 to "getaddrinfo
    failed" while 1958 succeeded.
    """
    url = (NCSS.format(var=var, year=year) +
           f"?var={var}&north=90&south={POLAR_MIN_LAT}&west=0&east=357.5"
           f"&vertCoord={level}"
           f"&time_start={year}-01-01T00:00:00Z&time_end={year}-12-31T23:59:59Z"
           "&accept=netcdf")
    last = None
    for attempt in range(6):
        try:
            with urllib.request.urlopen(url, timeout=300) as r:
                blob = r.read()
            if len(blob) >= 100_000:
                return blob
            last = f"short response ({len(blob)} B)"
        except Exception as exc:
            last = exc
        time.sleep(2 * (attempt + 1))
    raise IOError(f"{var}.{year} {level}hPa: {last}")


def _whole_file(var, year):
    """Fallback transport: the whole global file, reduced locally.

    NCSS is not uniformly healthy across the archive -- every request touching
    1975 returns HTTP 500, at every level, for both variables, and for half-year
    sub-ranges too, while the plain file for the same year downloads fine. So one
    bad year is not allowed to cost the whole extension. This route is the one
    validated bit-identical to NCSS on uwnd.1958 at 10 hPa (max|diff| = 0.0), it
    is just 60x heavier, so it is only used when NCSS fails.
    """
    CACHE.mkdir(exist_ok=True)
    local = CACHE / f"{var}.{year}.whole.nc"
    if not local.exists() or local.stat().st_size < 1_000_000:
        last = None
        for attempt in range(4):
            try:
                urllib.request.urlretrieve(
                    f"https://downloads.psl.noaa.gov/Datasets/"
                    f"ncep.reanalysis.dailyavgs/pressure/{var}.{year}.nc", local)
                if local.stat().st_size >= 1_000_000:
                    break
                last = f"short file ({local.stat().st_size} B)"
            except Exception as exc:
                last = exc
            time.sleep(3 * (attempt + 1))
        else:
            raise IOError(f"{var}.{year} whole-file fallback: {last}")
    return local


def _reduce(blob, var, year, level):
    """Cosine-weighted 60-90N cap mean, plus the 60N zonal mean, computed LOCALLY.

    The OPeNDAP route this script originally used returned SILENT ZEROS for large
    aggregate requests: a single-timestep read of uwnd.1965 at 10 hPa gave a
    correct 30.593 m/s polar-cap value, while the mean over the same year's full
    array came back exactly 0.0 for every day, with no error and no NaN. Both the
    weighted and unweighted means were affected, so it was a transport failure,
    not an aggregation one. Averaging that in would have silently added 21 years
    of zeros -- the same silent-contamination class as the CMIP6 _FillValue trap
    already documented in this project.

    Every year-level is therefore checked against SANE before it is accepted, and
    a violation raises instead of warning.
    """
    with xr.open_dataset(io.BytesIO(blob)) as ds:
        da = ds[var].squeeze(drop=True).load()
    return _reduce_da(da, var, year, level)


def _reduce_whole(path, var, year, level):
    """Same reduction, from a whole global file: subset first, then reduce."""
    with xr.open_dataset(path, engine="netcdf4") as ds:
        da = ds[var].sel(level=level)
        da = da.sel(lat=da.lat[da.lat >= POLAR_MIN_LAT]).squeeze(drop=True).load()
    return _reduce_da(da, var, year, level)


def _reduce_da(da, var, year, level):
    lat = da.lat
    raw = da.values
    if not np.isfinite(raw).any() or np.nanstd(raw) == 0.0:
        raise ValueError(f"{var}.{year} {level}hPa: degenerate array")

    w = np.cos(np.deg2rad(lat)).clip(min=0.0)
    cap = da.weighted(w).mean(dim=["lat", "lon"])
    # WMO SSW definition uses the zonal-mean zonal wind at 60N exactly, not the
    # cap mean, so it is carried alongside rather than reconstructed later.
    at60 = da.sel(lat=POLAR_MIN_LAT, method="nearest").mean(dim="lon")

    t = pd.to_datetime(cap.time.values).normalize()
    lo, hi = SANE[VARS[var]]
    m = float(np.nanmean(cap.values))
    if not (lo <= m <= hi):
        raise ValueError(f"{var}.{year} {level}hPa: cap mean {m:.2f} outside "
                         f"plausible range [{lo}, {hi}] -- refusing")
    return pd.DataFrame({f"{VARS[var]}_{level}hPa": cap.values,
                         f"{VARS[var]}_60N_{level}hPa": at60.values}, index=t)


def series(var, year, level):
    CACHE.mkdir(exist_ok=True)
    p = CACHE / f"{var}.{year}.{level}.parquet"
    if p.exists():
        return pd.read_parquet(p)
    try:
        df = _reduce(_fetch(var, year, level), var, year, level)
    except IOError:
        df = _reduce_whole(_whole_file(var, year), var, year, level)
    df.to_parquet(p)
    return df


def main():
    ref = pd.read_parquet(TARGET)
    tz = ref.index.tz
    print(f"existing: {ref.index.min()} .. {ref.index.max()}  ({len(ref):,} days)")

    jobs = [(v, y, l) for y in YEARS for v in VARS for l in LEVELS]
    print(f"{len(jobs)} variable-year-level requests, {WORKERS} workers")

    done, failed = {}, []
    with ThreadPoolExecutor(WORKERS) as ex:
        futs = {ex.submit(series, *j): j for j in jobs}
        for n, (f, j) in enumerate(futs.items(), 1):
            try:
                done[j] = f.result()
            except Exception as exc:
                failed.append((j, exc))
                print(f"  FAILED {j}: {exc}", flush=True)
            if n % 25 == 0:
                print(f"  {n}/{len(jobs)}", flush=True)

    if failed:
        print(f"\n{len(failed)} requests failed -- refusing to write a partial "
              f"record. Rerun; cached successes are reused.")
        sys.exit(1)

    ext = pd.concat([pd.concat([done[(v, y, l)] for v in VARS for l in LEVELS],
                               axis=1) for y in YEARS]).sort_index()
    if tz is not None:
        ext.index = ext.index.tz_localize("UTC")
    ext.to_parquet(OUT)
    print(f"\nextension: {ext.index.min()} .. {ext.index.max()} "
          f"({len(ext):,} days, {len(ext.columns)} cols) -> {OUT.name}")

    missing = [c for c in ref.columns if c not in ext.columns]
    print(f"target columns not produced: {missing or 'none'}")

    # sanity: does the pre-satellite era look physically like the satellite era?
    print("\n=== PRE-SATELLITE PLAUSIBILITY CHECK (winter DJF) ===")
    print(f"{'column':24s} {'1958-1978':>11s} {'1979-2024':>11s} {'ratio':>8s}")
    print("-" * 58)
    for c in ["uwnd_ms_10hPa", "uwnd_ms_100hPa", "air_K_10hPa", "hgt_m_100hPa"]:
        if c not in ref.columns:
            continue
        a = ext[c][np.isin(ext.index.month, (12, 1, 2))].mean()
        b = ref[c][np.isin(ref.index.month, (12, 1, 2))].mean()
        print(f"{c:24s} {a:11.2f} {b:11.2f} {a / b if b else np.nan:8.3f}")
    print()
    for c in ["uwnd_ms_10hPa", "air_K_10hPa", "hgt_m_100hPa"]:
        if c not in ref.columns:
            continue
        a = ext[c][np.isin(ext.index.month, (12, 1, 2))].std()
        b = ref[c][np.isin(ref.index.month, (12, 1, 2))].std()
        print(f"{c:20s} SD {a:11.2f} {b:11.2f} {a / b if b else np.nan:8.3f}")
    print("\n  A variance ratio far below 1 at 10 hPa is the signature of a")
    print("  radiosonde-era analysis being nudged toward climatology, and would")
    print("  suppress event detection in the pre-satellite arm.")


if __name__ == "__main__":
    main()
