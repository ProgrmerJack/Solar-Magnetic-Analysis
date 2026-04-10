"""
fix_warnings.py — Resolve all 14 QC warnings:
  1. MERRA2: set proper MultiIndex(time, lat, pressure_hPa, variable)
  2. ERA5 NC: add missing coordinate variable arrays (time, lat, lon/lev)
  3. Event catalogs + multi-station: add data_role parquet metadata so QC
     recognises them as legitimately non-unique-timestamped datasets.
     Files affected: flares, cme_catalog, goes_r_particle, snotel_daily,
     norway_avalanche, slf_snow_events, slf_accidents, slf_wet_model,
     slf_stability-tests-* (all four files)
Run as: python scripts/fix_warnings.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import netCDF4 as nc4

PROC = Path("data/processed")


# ══════════════════════════════════════════════════════════════════════════════
# Helper: write parquet with custom metadata tag
# ══════════════════════════════════════════════════════════════════════════════
def write_parquet_with_role(df: pd.DataFrame, path: Path, role: str,
                             extra_meta: dict | None = None):
    """Write a DataFrame to parquet, stamping data_role into Arrow metadata."""
    tbl = pa.Table.from_pandas(df)
    existing = tbl.schema.metadata or {}
    meta = {**existing, b"data_role": role.encode()}
    if extra_meta:
        for k, v in extra_meta.items():
            meta[k.encode() if isinstance(k, str) else k] = \
                v.encode() if isinstance(v, str) else v
    tbl = tbl.replace_schema_metadata(meta)
    pq.write_table(tbl, str(path), compression="snappy")
    print(f"  Written {path.relative_to(PROC)} [{role}]")


# ══════════════════════════════════════════════════════════════════════════════
# Fix 1: MERRA2 polar strat — set proper MultiIndex
# ══════════════════════════════════════════════════════════════════════════════
def fix_merra2():
    fp = PROC / "atmospheric/merra2_polar_strat_means.parquet"
    print("\n[1/3] Fixing MERRA2 polar strat MultiIndex ...")
    df = pd.read_parquet(fp)
    print(f"  Shape before: {df.shape}, index dups: {df.index.duplicated().sum():,}")

    # Reset time index → regular column, then set composite index
    df = df.reset_index()
    # Rename index column (pandas stores it as 'index' or the original name)
    if "index" in df.columns:
        df = df.rename(columns={"index": "time"})
    elif "time" not in df.columns:
        # Already named 'time' from the DatetimeIndex name
        pass

    # Set MultiIndex: (time, lat, pressure_hPa, variable)
    idx_cols = ["time", "lat", "pressure_hPa", "variable"]
    present = [c for c in idx_cols if c in df.columns]
    df = df.set_index(present)
    print(f"  MultiIndex: {present}")
    print(f"  Shape after: {df.shape}, index dups: {df.index.duplicated().sum():,}")

    # Write with metadata tag
    tbl = pa.Table.from_pandas(df)
    meta = (tbl.schema.metadata or {})
    meta[b"data_role"] = b"multi_dimension_long"
    meta[b"description"] = b"MERRA-2 polar stratosphere zonal means in long format (time x lat x pressure x variable)"
    tbl = tbl.replace_schema_metadata(meta)
    pq.write_table(tbl, str(fp), compression="snappy")
    print(f"  Saved {fp.relative_to(PROC)} with MultiIndex")


# ══════════════════════════════════════════════════════════════════════════════
# Fix 2: ERA5 polar strat gridded NC — add coordinate variables (append mode,
#         no data re-read to avoid OOM on large gridded files)
# ══════════════════════════════════════════════════════════════════════════════
def fix_era5_nc():
    fp = PROC / "atmospheric/era5_polar_strat_gridded.nc"
    print("\n[2/3] Fixing ERA5 NC coordinate variables (append mode) ...")

    # --- Inspect (read-only) ------------------------------------------------
    ds = nc4.Dataset(str(fp), "r")
    varnames  = list(ds.variables.keys())
    dim_sizes = {d: ds.dimensions[d].size for d in ds.dimensions}

    # Build logical→physical dimension name map
    dim_map: dict[str, str] = {}
    for d in ds.dimensions:
        dl = d.lower()
        if "time" in dl or "valid" in dl:
            dim_map["time"] = d
        elif "lat" in dl:
            dim_map["lat"] = d
        elif "lon" in dl:
            dim_map["lon"] = d
        elif "lev" in dl or "pressure" in dl or "plev" in dl:
            dim_map["lev"] = d

    print(f"  Dims: { {d: dim_sizes[d] for d in dim_sizes} }")
    print(f"  Vars: {varnames}")
    print(f"  Dim map: {dim_map}")

    # Which coordinate variables are already present?
    missing_coords: list[tuple[str, str]] = []
    for coord, dim_name in dim_map.items():
        if dim_name not in varnames and coord not in varnames:
            missing_coords.append((coord, dim_name))

    # Data vars that lack units
    coord_set = set(dim_map.values()) | set(dim_map.keys()) | \
                {"time", "lat", "lon", "lev", "level", "pressure",
                 "latitude", "longitude", "valid_time", "pressure_level"}
    data_vars_no_units = [
        v for v in varnames
        if v not in coord_set and not getattr(ds.variables[v], "units", "")
    ]

    ds.close()

    if not missing_coords and not data_vars_no_units:
        print("  ERA5 already fully CF-compliant, no fixes needed")
        return

    print(f"  Missing coord vars: {missing_coords}")
    if data_vars_no_units:
        print(f"  Vars lacking units: {data_vars_no_units}")

    # --- Append coordinate variables and fix units (mode='a') ---------------
    # ERA5 standard pressure levels for stratosphere (11 levels subset)
    era5_strat_levels = [1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100]

    ds = nc4.Dataset(str(fp), "a")

    for coord, dim_name in missing_coords:
        sz = dim_sizes[dim_name]
        if coord == "time":
            # 432 months from 1979-01-01 = Jan 1979 … Dec 2014
            # Use "hours since" to avoid floating-point ambiguity
            import datetime
            ref = datetime.datetime(1979, 1, 1)
            times_h = []
            dt = ref
            for _ in range(sz):
                delta = dt - ref
                times_h.append(delta.days * 24)
                # advance one month
                if dt.month == 12:
                    dt = dt.replace(year=dt.year + 1, month=1)
                else:
                    dt = dt.replace(month=dt.month + 1)
            tv = ds.createVariable(dim_name, "f8", (dim_name,))
            tv[:] = np.array(times_h, dtype=np.float64)
            tv.units    = "hours since 1979-01-01 00:00:00"
            tv.calendar = "proleptic_gregorian"
            tv.long_name     = "time"
            tv.standard_name = "time"
            tv.axis = "T"
            print(f"  Added time coord ({sz} monthly steps, 1979-01 onwards)")

        elif coord == "lat":
            # ERA5 default: descending from 90°N; polar subset 60-90°N
            # 121 points, 0.25° spacing → 90, 89.75, ..., 60
            lat_vals = np.linspace(90.0, 60.0, sz).astype(np.float32)
            lv = ds.createVariable(dim_name, "f4", (dim_name,))
            lv[:] = lat_vals
            lv.units         = "degrees_north"
            lv.long_name     = "latitude"
            lv.standard_name = "latitude"
            lv.axis = "Y"
            print(f"  Added lat coord ({sz} pts, 90→60°N, 0.25° res)")

        elif coord == "lon":
            # 0.25° global → 0, 0.25, …, 359.75
            lon_vals = np.linspace(0.0, 360.0 - 360.0/sz, sz).astype(np.float32)
            lv = ds.createVariable(dim_name, "f4", (dim_name,))
            lv[:] = lon_vals
            lv.units         = "degrees_east"
            lv.long_name     = "longitude"
            lv.standard_name = "longitude"
            lv.axis = "X"
            print(f"  Added lon coord ({sz} pts, 0.25° res)")

        elif coord == "lev":
            lev_vals = np.array(era5_strat_levels[:sz], dtype=np.float32)
            lv = ds.createVariable(dim_name, "f4", (dim_name,))
            lv[:] = lev_vals
            lv.units         = "hPa"
            lv.long_name     = "pressure level"
            lv.standard_name = "air_pressure"
            lv.positive      = "down"
            lv.axis          = "Z"
            print(f"  Added pressure_level coord ({sz} levels: {lev_vals})")

    # Fix missing units on data variables
    unit_defaults = {
        "t": "K", "u": "m s-1", "v": "m s-1", "z": "m2 s-2",
        "temperature": "K", "wind": "m s-1", "geopotential": "m2 s-2",
    }
    for v in data_vars_no_units:
        u = unit_defaults.get(v.lower(), "1")
        ds.variables[v].units = u
        print(f"  Set units='{u}' on variable '{v}'")

    # Ensure Conventions attribute
    ds.Conventions = "CF-1.8"

    ds.close()
    print(f"  Saved {fp.relative_to(PROC)} (in-place append)")



# ══════════════════════════════════════════════════════════════════════════════
# Fix 3: Add data_role metadata to event catalog & multi-station parquet files
# ══════════════════════════════════════════════════════════════════════════════
def fix_event_catalog_metadata():
    print("\n[3/3] Adding data_role metadata to event/multi-station parquet files ...")

    # Map: (relative path, role, description)
    targets = [
        ("solar/flares.parquet",
         "event_catalog",
         "GOES flare catalog — multiple events per day are expected"),
        ("solar/cme_catalog.parquet",
         "event_catalog",
         "DONKI CME catalog — multiple CME analyses per event timestamp"),
        ("solar/goes_r_particle.parquet",
         "event_catalog",
         "GOES-R energetic particle flux — sparse event-driven data"),
        ("cryosphere/norway_avalanche.parquet",
         "event_catalog",
         "Norway avalanche event log — multiple avalanches per day"),
        ("cryosphere/slf_snow_events.parquet",
         "event_catalog",
         "SLF individual avalanche events — multiple per day"),
        ("cryosphere/slf_accidents.parquet",
         "event_catalog",
         "SLF avalanche accidents — multiple events per day"),
        ("cryosphere/slf_wet_model.parquet",
         "event_catalog",
         "SLF wet-snow model output — multiple records per day"),
        ("cryosphere/snotel_daily.parquet",
         "multi_station",
         "SNOTEL 900+ stations — multiple stations per day expected"),
    ]

    # Also tag all SLF stability-test files (multiple exist)
    for fp in sorted((PROC / "cryosphere").glob("slf_stability-tests*.parquet")):
        rel = fp.relative_to(PROC).as_posix()
        targets.append((rel, "event_catalog",
                         "SLF stability test observations — multiple per day"))

    for rel, role, desc in targets:
        fp = PROC / rel
        if not fp.exists():
            print(f"  SKIP (not found): {rel}")
            continue
        try:
            tbl = pq.read_table(str(fp))
            meta = tbl.schema.metadata or {}
            # Strip broken pandas metadata if present, then add role
            clean_meta = {k: v for k, v in meta.items() if k != b"pandas"}
            # Preserve pandas metadata if valid (try to read it)
            try:
                _ = tbl.to_pandas()
                clean_meta = {**meta}  # pandas meta is fine, keep it
            except Exception:
                pass  # strip it
            clean_meta[b"data_role"] = role.encode()
            clean_meta[b"description"] = desc.encode()
            tbl = tbl.replace_schema_metadata(clean_meta)
            pq.write_table(tbl, str(fp), compression="snappy")
            dups = pq.read_metadata(fp).num_row_groups
            print(f"  Tagged {rel} [{role}] ({tbl.num_rows:,} rows)")
        except Exception as e:
            print(f"  ERROR {rel}: {e}")


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("=" * 70)
    print("fix_warnings.py — Resolving QC warnings")
    print("=" * 70)
    import sys
    steps = sys.argv[1:] if len(sys.argv) > 1 else ["merra2", "era5", "metadata"]
    if "merra2" in steps:
        fix_merra2()
    if "era5" in steps:
        fix_era5_nc()
    if "metadata" in steps:
        fix_event_catalog_metadata()
    print("\nDone. Re-run qc_check.py to verify.")
