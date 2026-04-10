"""
Add missing CF units attributes to MLS gridded NC data variables
and fix ERA5 polar strat time units attr.
"""
import netCDF4 as nc4
import shutil, tempfile, os
import xarray as xr
from pathlib import Path

PROC = Path("data/processed/atmospheric")

# Units for each MLS product variable
MLS_UNITS = {
    "hno3":        "ppbv",
    "temperature": "K",
    "n2o":         "ppbv",
    "ozone":       "ppmv",
}

print("Adding units to MLS gridded NC data variables...")
for fname, units in [
    ("mls_hno3_gridded.nc",        "ppbv"),
    ("mls_n2o_gridded.nc",         "ppbv"),
    ("mls_ozone_gridded.nc",       "ppmv"),
    ("mls_temperature_gridded.nc", "K"),
]:
    fp = PROC / fname
    if not fp.exists() or fp.stat().st_size < 1000:
        print(f"  SKIP {fname}")
        continue

    # Load into memory, add units, atomic rewrite
    with xr.open_dataset(fp, engine="netcdf4") as ds_raw:
        ds = ds_raw.load()

    for var in ds.data_vars:
        if not ds[var].attrs.get("units"):
            ds[var].attrs["units"] = units
            print(f"  {fname}: set {var}.units = {units!r}")

    # Ensure time has units attr (copy from encoding if needed)
    if "time" in ds.coords:
        ds["time"].attrs.setdefault("standard_name", "time")
        ds["time"].attrs.setdefault("axis", "T")

    fd, tmp = tempfile.mkstemp(suffix=".nc", dir=fp.parent)
    os.close(fd)
    tmp_p = Path(tmp)
    try:
        enc = {}
        if "time" in ds.coords:
            enc["time"] = {
                "units": "days since 2004-01-01",
                "calendar": "proleptic_gregorian",
                "dtype": "float64",
            }
        for var in ds.data_vars:
            enc[var] = {"zlib": True, "complevel": 4}
        ds.to_netcdf(str(tmp_p), format="NETCDF4", encoding=enc)

        # Verify magic before committing
        magic = tmp_p.read_bytes()[:8]
        if magic[1:4] != b"HDF":
            raise RuntimeError(f"Bad magic {magic.hex()}")

        # Open and add time units attr directly via netCDF4 (encoding doesn't appear in attrs)
        with nc4.Dataset(str(tmp_p), "a") as nc:
            if "time" in nc.variables and not getattr(nc.variables["time"], "units", ""):
                nc.variables["time"].units = "days since 2004-01-01"
                nc.variables["time"].calendar = "proleptic_gregorian"

        shutil.move(str(tmp_p), str(fp))
        print(f"  FIXED {fname} ({fp.stat().st_size/1e6:.1f} MB)")
    except Exception as e:
        if tmp_p.exists():
            tmp_p.unlink()
        print(f"  FAIL {fname}: {e}")
    finally:
        ds.close()

# Fix ERA5 polar strat time units attr
print()
print("Checking ERA5 polar strat time units...")
fp_era5 = PROC / "era5_polar_strat_gridded.nc"
if fp_era5.exists():
    with nc4.Dataset(str(fp_era5), "r") as nc:
        t_units = getattr(nc.variables.get("time", None), "units", "") if "time" in nc.variables else ""
        print(f"  Current time units: {t_units!r}")
    if not t_units:
        # Add time units in-place (append mode)
        with nc4.Dataset(str(fp_era5), "a") as nc:
            if "time" in nc.variables:
                nc.variables["time"].units = "hours since 1900-01-01"
                nc.variables["time"].calendar = "proleptic_gregorian"
        print("  FIXED: added time units to era5_polar_strat_gridded.nc")
    else:
        print("  OK: time units already present")

print()
print("Done. Re-run qc_check.py to verify.")
