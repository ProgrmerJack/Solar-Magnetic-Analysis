#!/usr/bin/env python3
"""
check_environment.py
====================
Verifies that this environment can actually do the project's work -- not that
packages import, but that the four things every analysis here depends on
really function.

WHY THIS EXISTS
  On 2026-09-17 the system interpreter had xarray with only the `scipy`
  backend, so `xr.open_dataset()` failed on every .nc file in the repo with a
  message about missing dependencies rather than a missing file. Nothing in
  the repo detected that; the analyses simply could not run, and the reason
  was three lines deep in an exception. An import check would have passed.

WHAT IT CHECKS
  1. imports            every direct dependency, with its version
  2. netcdf             opens a real .nc file from the repo
  3. parquet            opens a real .parquet file from the repo
  4. remote             the two anonymous endpoints the ingestion uses

Exit status is 0 only if every required check passes. Remote checks are
advisory: they fail on a plane, which is not an environment defect.

    .venv/bin/python ssw-design-analysis/environment/check_environment.py
    .venv/bin/python ssw-design-analysis/environment/check_environment.py --no-remote
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

REQUIRED = [
    "numpy", "pandas", "scipy", "xarray", "netCDF4", "h5netcdf",
    "cftime", "pyarrow", "sklearn", "statsmodels", "matplotlib",
    "requests", "zarr", "gcsfs", "dask",
]

# Real files, committed or re-acquirable, that the analyses actually read.
NETCDF_PROBE = ROOT / "ssw-design-analysis/03_data_ingestion/_ncep_cache/uwnd.1958.nc"
PARQUET_PROBE = ROOT / "ssw-design-analysis/03_data_ingestion/era5_nam_daily.parquet"

# Anonymous endpoints. CEDA browsing needs no credentials; only downloads do.
REMOTES = [
    ("CEDA / SNAPSI",
     "https://dap.ceda.ac.uk/badc/snap/data/post-cmip6/SNAPSI/"),
    ("NCEP THREDDS",
     "https://psl.noaa.gov/thredds/catalog/Datasets/ncep.reanalysis/catalog.html"),
]

rows: list[tuple[str, str, str]] = []


def record(status: str, name: str, detail: str) -> None:
    rows.append((status, name, detail))


def check_imports() -> bool:
    ok = True
    for mod in REQUIRED:
        try:
            m = importlib.import_module(mod)
            record("PASS", f"import {mod}", getattr(m, "__version__", "(no __version__)"))
        except Exception as exc:  # noqa: BLE001 -- we want the class name reported
            record("FAIL", f"import {mod}", f"{type(exc).__name__}: {exc}")
            ok = False
    return ok


def check_netcdf() -> bool:
    try:
        import xarray as xr
    except Exception as exc:  # noqa: BLE001
        record("FAIL", "netcdf", f"xarray unavailable: {exc}")
        return False

    engines = sorted(xr.backends.list_engines().keys())
    record("INFO", "xarray backends", ", ".join(engines))
    if not ({"netcdf4", "h5netcdf"} & set(engines)):
        record("FAIL", "netcdf backend",
               "neither netcdf4 nor h5netcdf is installed -- no .nc in this repo can be opened")
        return False

    if not NETCDF_PROBE.exists():
        record("SKIP", "netcdf open", f"probe absent: {NETCDF_PROBE.relative_to(ROOT)}")
        return True
    try:
        with xr.open_dataset(NETCDF_PROBE) as ds:
            record("PASS", "netcdf open",
                   f"{NETCDF_PROBE.name}: vars={list(ds.data_vars)} dims={dict(ds.sizes)}")
        return True
    except Exception as exc:  # noqa: BLE001
        record("FAIL", "netcdf open", f"{type(exc).__name__}: {str(exc)[:200]}")
        return False


def check_parquet() -> bool:
    try:
        import pandas as pd
    except Exception as exc:  # noqa: BLE001
        record("FAIL", "parquet", f"pandas unavailable: {exc}")
        return False
    if not PARQUET_PROBE.exists():
        record("SKIP", "parquet open", f"probe absent: {PARQUET_PROBE.relative_to(ROOT)}")
        return True
    try:
        df = pd.read_parquet(PARQUET_PROBE)
        record("PASS", "parquet open", f"{PARQUET_PROBE.name}: shape={df.shape}")
        return True
    except Exception as exc:  # noqa: BLE001
        record("FAIL", "parquet open", f"{type(exc).__name__}: {str(exc)[:200]}")
        return False


def check_remote() -> None:
    """Advisory only -- never changes the exit status."""
    try:
        import requests
    except Exception as exc:  # noqa: BLE001
        record("WARN", "remote", f"requests unavailable: {exc}")
        return
    for name, url in REMOTES:
        try:
            r = requests.get(url, timeout=20)
            status = "PASS" if r.status_code == 200 else "WARN"
            record(status, f"remote {name}", f"HTTP {r.status_code}")
        except Exception as exc:  # noqa: BLE001
            record("WARN", f"remote {name}", f"{type(exc).__name__}: {str(exc)[:120]}")


def main() -> int:
    do_remote = "--no-remote" not in sys.argv

    print(f"python      {sys.version.split()[0]}  ({sys.executable})")
    print(f"repo root   {ROOT}")
    print()

    ok_imports = check_imports()
    ok_netcdf = check_netcdf()
    ok_parquet = check_parquet()
    if do_remote:
        check_remote()

    width = max(len(n) for _, n, _ in rows)
    for status, name, detail in rows:
        print(f"[{status:4s}] {name:<{width}}  {detail}")

    required_ok = ok_imports and ok_netcdf and ok_parquet
    n_fail = sum(1 for s, _, _ in rows if s == "FAIL")
    n_warn = sum(1 for s, _, _ in rows if s == "WARN")
    print()
    print(f"{'ENVIRONMENT OK' if required_ok else 'ENVIRONMENT BROKEN'}"
          f"  --  {n_fail} failure(s), {n_warn} advisory warning(s)")
    if not required_ok:
        print("Fix with: python3 -m venv .venv && "
              ".venv/bin/pip install -r ssw-design-analysis/environment/requirements.lock.txt")
    return 0 if required_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
