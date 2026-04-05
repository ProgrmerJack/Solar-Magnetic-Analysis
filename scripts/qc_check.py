"""
Quality Control check for all processed datasets.
Validates Parquet and NetCDF4 files against scientific community standards.
"""
import json
import sys
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import netCDF4 as nc4
from pathlib import Path

PROC = Path("data/processed")
results = []  # (status, name, message)


def ok(name, msg):
    results.append(("PASS", name, msg))


def warn(name, msg):
    results.append(("WARN", name, msg))


def fail(name, msg):
    results.append(("FAIL", name, msg))


# ── 1. PARQUET FILES ──────────────────────────────────────────────────────────
print("Checking parquet files...")
for fp in sorted(PROC.rglob("*.parquet")):
    name = fp.relative_to(PROC).as_posix()
    try:
        df = pd.read_parquet(fp)
        nrows, ncols = df.shape
        if nrows == 0:
            fail(name, "EMPTY — 0 rows")
            continue
        idx = df.index
        if not isinstance(idx, pd.DatetimeIndex):
            # Auxiliary datasets legitimately use RangeIndex (e.g. HARP-indexed SHARP data)
            import pyarrow.parquet as _pq
            _meta = (_pq.read_metadata(fp).schema.to_arrow_schema().metadata or {})
            _role = _meta.get(b"data_role", b"").decode()
            if "auxiliary" in _role:
                ok(name, f"{nrows:,}r x {ncols}c | auxiliary — RangeIndex expected (no timestamps)")
            else:
                fail(name, "Non-datetime index: " + type(idx).__name__)
            continue
        if idx.tz is None:
            fail(name, "DatetimeIndex has no timezone — must be UTC for reproducibility")
        elif str(idx.tz) not in ("UTC", "utc", "UTC+00:00"):
            warn(name, "Timezone=" + str(idx.tz) + " expected UTC")
        else:
            dups = int(idx.duplicated().sum())
            nan_r = float(df.isnull().values.mean())
            tmin = idx.min().date()
            tmax = idx.max().date()
            msg = f"{nrows:,}r x {ncols}c | {tmin}/{tmax} | NaN={nan_r:.0%}"
            if dups > 0:
                warn(name, msg + f" | {dups} duplicate timestamps")
            elif nan_r > 0.95:
                warn(name, msg + " | >95% NaN — verify fill-value handling")
            else:
                ok(name, msg)
    except Exception as e:
        fail(name, str(e)[:120])

print(f"  {len([r for r in results if r[0]=='PASS'])} pass, "
      f"{len([r for r in results if r[0]=='WARN'])} warn, "
      f"{len([r for r in results if r[0]=='FAIL'])} fail so far")


# ── 2. NETCDF4 / CF COMPLIANCE ────────────────────────────────────────────────
print("Checking NetCDF4 files...")
for fp in sorted(PROC.rglob("*.nc")):
    name = fp.relative_to(PROC).as_posix()
    try:
        ds = nc4.Dataset(str(fp))
        conv = getattr(ds, "Conventions", "")
        has_cf = "CF" in str(conv)

        # Required coordinate variables
        # lon/longitude is optional for zonal-mean products (MLS, MERRA-2 polar)
        varnames = list(ds.variables.keys())
        has_time = "time" in varnames
        has_lat  = any(v in varnames for v in ("lat", "latitude"))
        has_lon  = any(v in varnames for v in ("lon", "longitude"))
        is_zonal = (not has_lon) and any(v in str(getattr(ds, "title", "")).lower()
                                         for v in ("zonal", "mls", "polar", "mean"))
        coords_ok = has_time and has_lat and (has_lon or is_zonal)

        # Time must have units (CF requirement); xarray may store in encoding not attrs
        time_units = ""
        if has_time:
            time_units = getattr(ds.variables["time"], "units", "")

        # Data variables should have units
        coord_set = {"time", "lat", "lon", "lev", "level", "pressure",
                     "latitude", "longitude"}
        data_vars = [v for v in varnames if v not in coord_set]
        no_units  = [v for v in data_vars if not getattr(ds.variables[v], "units", "")]

        dims = {k: v.size for k, v in ds.dimensions.items()}
        sz_mb = fp.stat().st_size / 1e6

        issues = []
        if not has_cf:
            issues.append("Conventions attr missing/not CF-1.x")
        if not coords_ok:
            issues.append("Missing time/lat/lon coordinate variables")
        if not time_units:
            issues.append("time variable lacks units attr")
        if no_units:
            issues.append(f"{len(no_units)} data vars lack units: {no_units[:4]}")

        ds.close()

        msg = f"{sz_mb:.0f} MB | dims={dims}"
        if issues:
            warn(name, msg + " | " + "; ".join(issues))
        else:
            ok(name, msg + " | CF-compliant")
    except Exception as e:
        fail(name, str(e)[:120])

print(f"  {len([r for r in results if r[0]=='PASS'])} pass, "
      f"{len([r for r in results if r[0]=='WARN'])} warn, "
      f"{len([r for r in results if r[0]=='FAIL'])} fail so far")


# ── 3. MANIFEST COMPLETENESS ─────────────────────────────────────────────────
mf = PROC / "manifest.json"
if mf.exists():
    with open(mf) as f:
        manifest = json.load(f)
    # Manifest may be a list of records or a dict
    if isinstance(manifest, list):
        n_entries = len(manifest)
        missing = [e.get("key", "?") for e in manifest
                   if not Path(str(e.get("path", ""))).exists()]
    else:
        n_entries = len(manifest)
        missing = [k for k, v in manifest.items()
                   if isinstance(v, dict) and not Path(str(v.get("path", ""))).exists()]
    ok("manifest.json", f"{n_entries} registered datasets")
    if missing:
        warn("manifest.json", f"{len(missing)} registered paths not found: {missing[:5]}")
else:
    fail("manifest.json", "MISSING — required for data provenance")


# ── 4. PHYSICAL RANGE CHECKS ─────────────────────────────────────────────────
print("Checking physical value ranges...")
range_checks = [
    # (file_rel, column, lo, hi, units)
    ("solar/goes_xrs.parquet",           "xrsb_flux",     1e-10, 1e-2,  "W/m²"),
    ("geomagnetic/dst_index.parquet",    "dst",           -600,  100,   "nT"),
    ("geomagnetic/kp_index.parquet",     "kp",             0,    9,     "0-9"),
    ("solar/omni_1min.parquet",          "Vx",             -2500, -200, "km/s"),   # Vx is negative (sunward)
    ("solar/omni_1min.parquet",          "Bz",            -200,  200,   "nT"),
    ("solar/ace_dscovr.parquet",         "proton_speed",   200,  2500,  "km/s"),
    ("solar/ace_dscovr.parquet",         "bz_gsm",        -100,  100,   "nT"),
    ("solar/psp_mag.parquet",            "B_mag_mean",     0,    1000,  "nT"),
    ("cryosphere/snotel_daily.parquet",  "snwd_mm",        0,    15000, "mm"),
]

for rel, col, lo, hi, units in range_checks:
    fp = PROC / rel
    if not fp.exists():
        warn(rel, f"File missing — cannot check {col}")
        continue
    try:
        df = pd.read_parquet(fp, columns=[col])
        vals = df[col].dropna()
        if len(vals) == 0:
            warn(rel, f"{col}: all NaN")
            continue
        vmin, vmax = float(vals.min()), float(vals.max())
        outliers = int(((vals < lo) | (vals > hi)).sum())
        pct = outliers / len(vals) * 100
        msg = f"{col} [{vmin:.3g}, {vmax:.3g}] {units}"
        if pct > 1.0:
            warn(rel, msg + f" — {pct:.1f}% outside expected [{lo},{hi}]")
        else:
            ok(rel, msg + " ✓")
    except KeyError:
        # Column not found — list available
        df2 = pd.read_parquet(fp)
        warn(rel, f"Column {col!r} not found. Available: {list(df2.columns)[:8]}")
    except Exception as e:
        warn(rel, f"Range check error: {e!s:.80}")


# ── 5. KEY SCIENCE DATASETS PRESENT ──────────────────────────────────────────
required = {
    # SOC core datasets
    "solar/flares.parquet":                     "GOES flare catalog (SOC power-law fit)",
    "cryosphere/slf_snow_events.parquet":        "SLF individual avalanches (SOC power-law fit)",
    # Mechanistic chain
    "atmospheric/poes_noaa15_2013.parquet":      "POES EPP (EPP→NOx link)",
    "atmospheric/mls_hno3_polar.parquet":        "MLS HNO3 polar (NOx proxy)",
    "atmospheric/era5_polar_strat_gridded.nc":   "ERA5 polar strat (SSW detection)",
    "atmospheric/ssw_catalog.parquet":           "SSW event catalog (Butler 2017)",
    # Surface forcing
    "geomagnetic/dst_index.parquet":             "Dst index (storm classification)",
    "geomagnetic/kp_index.parquet":              "Kp index (storm classification)",
    "solar/goes_xrs.parquet":                    "GOES XRS (flare energy proxy)",
    "solar/omni_1min.parquet":                   "OMNI solar wind (L1 coupling)",
    # Avalanche records
    "cryosphere/snotel_daily.parquet":           "SNOTEL snowpack (meteorological control)",
    "cryosphere/norway_avalanche.parquet":        "Norway avalanche record",
}
print("Checking required science datasets...")
for rel, desc in required.items():
    fp = PROC / rel
    if fp.exists():
        sz = fp.stat().st_size / 1e6
        ok(f"REQUIRED: {fp.name}", f"{desc} — {sz:.1f} MB ✓")
    else:
        fail(f"REQUIRED: {fp.name}", f"MISSING — needed for: {desc}")


# ── PRINT SUMMARY ─────────────────────────────────────────────────────────────
passes = [r for r in results if r[0] == "PASS"]
warns  = [r for r in results if r[0] == "WARN"]
fails  = [r for r in results if r[0] == "FAIL"]

print("\n" + "=" * 70)
print(f"  QC RESULTS:  {len(passes)} PASS  |  {len(warns)} WARN  |  {len(fails)} FAIL")
print("=" * 70)

if fails:
    print("\n[FAILURES — must resolve before submission]")
    for _, n, m in fails:
        print(f"  FAIL  {n}")
        print(f"        {m}")

if warns:
    print("\n[WARNINGS — review before submission]")
    for _, n, m in warns:
        print(f"  WARN  {n}")
        print(f"        {m}")

print("\n[PASSES]")
for _, n, m in passes:
    print(f"  PASS  {n}: {m}")

print("\n" + "=" * 70)
print(f"VERDICT: ", end="")
if fails:
    print(f"{len(fails)} FAILURE(S) — fix before data submission")
elif warns:
    print(f"ACCEPTABLE with {len(warns)} warning(s) — review noted items")
else:
    print("ALL CLEAR — data meets scientific community standards")
print("=" * 70)
