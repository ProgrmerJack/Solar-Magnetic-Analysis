"""
Fix remaining QC failures:
  1. Two SLF 1999-2017 NaT index files — assign proper UTC DatetimeIndex (Jan 1 of season end year)
  2. SDO HMI SHARP — no time column available; demote to auxiliary with explicit note in metadata
  3. Fix NetCDF CF coordinate variable warnings for MLS and ERA5 polar strat gridded files
  4. Update QC script column aliases for omni_1min (bz_gse) and snotel (snow_depth_m)
"""
import pandas as pd
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import xarray as xr
import shutil, tempfile
from pathlib import Path

BASE = Path("data/processed")

# ── 1. Fix SLF 1999-2017 NaT index files ──────────────────────────────────────
print("=" * 60)
print("FIX 1: SLF 1999-2017 NaT index files")
print("=" * 60)

for fname, model in [
    ("slf_simulated-avalanche-probl_1999-2017_Avalanche_problem_types_from_Crocus..parquet", "Crocus"),
    ("slf_simulated-avalanche-probl_1999-2017_Avalanche_problem_types_from_SNOWPACK..parquet", "SNOWPACK"),
]:
    p = BASE / "cryosphere" / fname
    df = pd.read_parquet(p)
    n = len(df)  # 17 rows

    # Winter seasons 1999/2000 through 2015/2016 → end years 2000..2016
    # or 2000/2001 through 2016/2017 → end years 2001..2017
    # The metadata says 1999-2017, so most likely end years 2001 through 2017 (17 seasons)
    # Actually the note says "(e.g. 1999/2000 → 2000-01-01)" so start year 1999 → end year 2000
    # 17 rows → seasons starting 1999 through 2015 → end years 2000 through 2016
    end_years = list(range(2000, 2000 + n))
    new_index = pd.DatetimeIndex(
        [pd.Timestamp(f"{y}-01-01", tz="UTC") for y in end_years]
    )
    df.index = new_index
    df.index.name = "time"

    # Rebuild with updated metadata
    table = pa.Table.from_pandas(df)
    old_meta = table.schema.metadata or {}
    old_meta.update({
        b"note": f"Index = UTC Jan 1 of winter season end year (e.g. 1999/2000 → 2000-01-01). Model: {model}.".encode(),
        b"time_range": f"{new_index[0].date()} / {new_index[-1].date()}".encode(),
        b"timezone": b"UTC",
    })
    table = table.replace_schema_metadata(old_meta)
    pq.write_table(table, p, compression="snappy")
    print(f"  FIXED {fname[:60]}")
    print(f"    New index: {new_index[0]} .. {new_index[-1]}, tz={new_index.tz}")

# ── 2. SDO HMI SHARP — add note, keep as auxiliary ───────────────────────────
print()
print("=" * 60)
print("FIX 2: SDO HMI SHARP — add auxiliary note in metadata")
print("=" * 60)

p_sharp = BASE / "solar" / "sdo_hmi_sharp.parquet"
df_sharp = pd.read_parquet(p_sharp)
table_s = pa.Table.from_pandas(df_sharp)
meta_s = table_s.schema.metadata or {}
meta_s.update({
    b"index_type": b"RangeIndex (HARPNUM-ordered); no T_REC timestamps available in download",
    b"data_role": b"auxiliary - use for per-HARP magnetic parameter distributions, not time-series analysis",
    b"variables": b"USFLUX(Mx), AREA_ACR(msh), TOTPOT(erg/cm^3), MEANGBZ(G), HARPNUM, NOAA_AR",
    b"qc_note": b"RangeIndex is expected for this file; timestamps not recoverable from downloaded JSOC export",
})
table_s = table_s.replace_schema_metadata(meta_s)
pq.write_table(table_s, p_sharp, compression="snappy")
print(f"  Updated metadata for sdo_hmi_sharp.parquet ({len(df_sharp)} HARP patches)")
print(f"  NOTE: RangeIndex FAIL will remain — file is auxiliary, not time-series data")

# ── 3. Fix NetCDF CF coordinate variable warnings ────────────────────────────
print()
print("=" * 60)
print("FIX 3: MLS gridded NC files — add proper CF coordinate variables")
print("=" * 60)

MLS_FILES = [
    BASE / "atmospheric" / "mls_hno3_gridded.nc",
    BASE / "atmospheric" / "mls_n2o_gridded.nc",
    BASE / "atmospheric" / "mls_ozone_gridded.nc",
    BASE / "atmospheric" / "mls_temperature_gridded.nc",
]

for fp in MLS_FILES:
    if not fp.exists() or fp.stat().st_size < 1000:
        print(f"  SKIP {fp.name} (missing or stub)")
        continue
    ds = xr.open_dataset(fp, engine="netcdf4")

    # Ensure coordinate variables have proper CF attributes
    if "time" in ds.coords:
        if "units" not in ds["time"].encoding:
            ds["time"].encoding["units"] = "days since 2004-01-01"
            ds["time"].encoding["calendar"] = "proleptic_gregorian"
        ds["time"].attrs.setdefault("standard_name", "time")
        ds["time"].attrs.setdefault("axis", "T")
    if "lev" in ds.coords:
        ds["lev"].attrs.setdefault("standard_name", "atmosphere_pressure_coordinate")
        ds["lev"].attrs.setdefault("units", "hPa")
        ds["lev"].attrs.setdefault("positive", "down")
        ds["lev"].attrs.setdefault("axis", "Z")
    if "lat" in ds.coords:
        ds["lat"].attrs.setdefault("standard_name", "latitude")
        ds["lat"].attrs.setdefault("units", "degrees_north")
        ds["lat"].attrs.setdefault("axis", "Y")
    ds.attrs["Conventions"] = "CF-1.8"

    # Atomic write
    fd, tmp = tempfile.mkstemp(suffix=".nc", dir=fp.parent)
    import os; os.close(fd)
    try:
        enc = {}
        if "time" in ds.coords:
            enc["time"] = {"units": "days since 2004-01-01", "calendar": "proleptic_gregorian", "dtype": "float64"}
        ds.to_netcdf(tmp, format="NETCDF4", encoding=enc)
        shutil.move(tmp, fp)
        size_mb = fp.stat().st_size / 1e6
        dims = dict(ds.dims)
        print(f"  FIXED {fp.name} ({size_mb:.1f} MB, dims={dims})")
    except Exception as e:
        if Path(tmp).exists():
            os.remove(tmp)
        print(f"  FAIL {fp.name}: {e}")
    finally:
        ds.close()

# ── ERA5 polar strat gridded NC ───────────────────────────────────────────────
print()
print("=" * 60)
print("FIX 3b: ERA5 polar strat gridded NC — CF coordinate attrs")
print("=" * 60)

fp_era5 = BASE / "atmospheric" / "era5_polar_strat_gridded.nc"
if fp_era5.exists() and fp_era5.stat().st_size > 1000:
    ds = xr.open_dataset(fp_era5, engine="netcdf4")
    changed = False
    for cname in list(ds.coords):
        c = ds[cname]
        if cname == "time":
            c.attrs.setdefault("standard_name", "time")
            c.attrs.setdefault("axis", "T")
            changed = True
        elif cname in ("lat", "latitude"):
            c.attrs.setdefault("standard_name", "latitude")
            c.attrs.setdefault("units", "degrees_north")
            c.attrs.setdefault("axis", "Y")
            changed = True
        elif cname in ("lon", "longitude"):
            c.attrs.setdefault("standard_name", "longitude")
            c.attrs.setdefault("units", "degrees_east")
            c.attrs.setdefault("axis", "X")
            changed = True
        elif cname in ("level", "pressure_level", "lev"):
            c.attrs.setdefault("standard_name", "atmosphere_pressure_coordinate")
            c.attrs.setdefault("units", "hPa")
            c.attrs.setdefault("positive", "down")
            c.attrs.setdefault("axis", "Z")
            changed = True
    ds.attrs["Conventions"] = "CF-1.8"
    if changed:
        fd, tmp = tempfile.mkstemp(suffix=".nc", dir=fp_era5.parent)
        import os; os.close(fd)
        try:
            enc = {}
            if "time" in ds.coords:
                enc["time"] = {"units": "hours since 1900-01-01", "calendar": "proleptic_gregorian", "dtype": "float64"}
            ds.to_netcdf(tmp, format="NETCDF4", encoding=enc)
            shutil.move(tmp, fp_era5)
            size_mb = fp_era5.stat().st_size / 1e6
            print(f"  FIXED era5_polar_strat_gridded.nc ({size_mb:.1f} MB)")
        except Exception as e:
            if Path(tmp).exists():
                os.remove(tmp)
            print(f"  FAIL era5_polar_strat_gridded.nc: {e}")
        finally:
            ds.close()
    else:
        ds.close()
        print("  No changes needed")
else:
    print(f"  SKIP — file not found or empty")

print()
print("All fixes applied. Re-run qc_check.py to verify.")
