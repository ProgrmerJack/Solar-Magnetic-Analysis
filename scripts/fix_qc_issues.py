"""
Fix all QC failures and real warnings identified by qc_check.py.

Issues to fix:
1. SNOTEL: set 'date' column as UTC DatetimeIndex
2. SDO HMI: no date column — need T_REC from raw or use HARPNUM as index with note
3. SLF Crocus/SNOWPACK 1999-2017: 'season' col (e.g. '1999/2000') → parse to Jan 1 of second year
4. Dst: 98% NaN — values are there but tiny; the fill sentinel is 9999.0 → mask it
5. MLS polar: 100% NaN — investigate what's stored
6. OMNI fill values: mask known OMNI fill sentinels (9999, 99999, 999999)
7. NetCDF4 CF compliance: add Conventions, units, rename coords where needed
8. slf_wet_model: index is 1970 (epoch=0) — parse 'datum' column instead
"""
import warnings
warnings.filterwarnings("ignore")
import os
import pandas as pd
import numpy as np
import netCDF4 as nc4
import xarray as xr
from pathlib import Path
import sys

sys.path.insert(0, "scripts/process")
from _utils import save_parquet, save_netcdf4, PROCESSED_ROOT

PROC = PROCESSED_ROOT
fixed = []
skipped = []


def report(tag, name, msg):
    icon = {"FIX": "✓ FIX", "SKIP": "  SKIP", "FAIL": "✗ FAIL"}[tag]
    print(f"  {icon}  {name}: {msg}")
    if tag == "FIX":
        fixed.append(name)
    else:
        skipped.append(name)


# ── 1. SNOTEL: date column → DatetimeIndex ────────────────────────────────────
print("\n[1] SNOTEL daily")
fp = PROC / "cryosphere/snotel_daily.parquet"
df = pd.read_parquet(fp)
if not isinstance(df.index, pd.DatetimeIndex):
    df["date"] = pd.to_datetime(df["date"], utc=True)
    df = df.set_index("date").sort_index()
    meta = {
        "title": "SNOTEL Daily Snow/Weather Observations",
        "source": "USDA NRCS SNOTEL Network",
        "references": "https://www.nrcs.usda.gov/wps/portal/wcc/home/",
        "columns": "snwd_mm=snow depth; prec_mm=precipitation; tavg/tmax/tmin_c=temperature; wteq_mm=SWE",
        "time_range": f"{df.index.min()} / {df.index.max()}",
    }
    save_parquet(df, fp, meta)
    report("FIX", "snotel_daily.parquet", f"{len(df):,} rows, UTC DatetimeIndex set")
else:
    report("SKIP", "snotel_daily.parquet", "already has DatetimeIndex")


# ── 2. SDO HMI: add T_REC-based datetime index ────────────────────────────────
print("\n[2] SDO HMI SHARP")
fp = PROC / "solar/sdo_hmi_sharp.parquet"
df = pd.read_parquet(fp)
if not isinstance(df.index, pd.DatetimeIndex):
    # T_REC column holds timestamp strings like "2010.05.03_00:00:00_TAI"
    t_cols = [c for c in df.columns if "t_rec" in c.lower() or "time" in c.lower()]
    if t_cols:
        col = t_cols[0]
        # TAI timestamps: replace '_TAI' and '_' separators
        ts = df[col].astype(str).str.replace("_TAI", "").str.replace("_", " ")
        df.index = pd.to_datetime(ts, format="%Y.%m.%d %H:%M:%S", utc=True, errors="coerce")
        df.index.name = "time"
        df = df.drop(columns=[col]).sort_index()
        report_msg = f"T_REC → UTC DatetimeIndex, {df.index.isna().sum()} unparsed"
    else:
        # No time column — create synthetic index from HARPNUM+sequential, flag as auxiliary
        df.index = pd.RangeIndex(len(df))
        report("SKIP", "sdo_hmi_sharp.parquet", "No T_REC column — keeping RangeIndex (auxiliary data)")
        df = None
    if df is not None:
        meta = {
            "title": "SDO HMI Active Region SHARP Parameters",
            "source": "NASA SDO/HMI via JSOC",
            "references": "https://doi.org/10.1007/s11207-014-0529-3",
            "columns": "USFLUX=total unsigned flux Mx; AREA_ACR=area in CRS; TOTPOT=total potential energy; MEANGBZ=mean Bz",
            "time_range": f"{df.index.min()} / {df.index.max()}",
        }
        save_parquet(df, fp, meta)
        report("FIX", "sdo_hmi_sharp.parquet", f"{len(df):,} rows, UTC DatetimeIndex set")


# ── 3. SLF Simulated 1999–2017 (season string → datetime) ────────────────────
print("\n[3] SLF Simulated 1999-2017 Crocus & SNOWPACK")
for fname in [
    "slf_simulated-avalanche-probl_1999-2017_Avalanche_problem_types_from_Crocus..parquet",
    "slf_simulated-avalanche-probl_1999-2017_Avalanche_problem_types_from_SNOWPACK..parquet",
]:
    fp = PROC / "cryosphere" / fname
    df = pd.read_parquet(fp)
    if not isinstance(df.index, pd.DatetimeIndex) and "season" in df.columns:
        # 'season' like '1999/2000' → Jan 1 of the winter's end year as representative date
        def season_to_ts(s):
            try:
                year = int(str(s).split("/")[-1])
                return pd.Timestamp(f"{year}-01-01", tz="UTC")
            except Exception:
                return pd.NaT
        df.index = df["season"].apply(season_to_ts)
        df.index.name = "time"
        df = df.drop(columns=["season"]).sort_index()
        meta = {
            "title": "SLF Simulated Avalanche Problem Types 1999-2017",
            "source": "WSL/SLF Davos — Crocus/SNOWPACK model output",
            "references": "https://www.slf.ch/",
            "note": "Index = Jan 1 of winter season end year (e.g. 1999/2000 → 2000-01-01)",
            "time_range": f"{df.index.min()} / {df.index.max()}",
        }
        save_parquet(df, fp, meta)
        report("FIX", fname[:50], f"{len(df)} rows, season → UTC DatetimeIndex")
    else:
        report("SKIP", fname[:50], "already OK or no season column")


# ── 4. Dst: mask fill sentinels (9999, 99999) ─────────────────────────────────
print("\n[4] Dst index fill values")
fp = PROC / "geomagnetic/dst_index.parquet"
df = pd.read_parquet(fp)
col = "dst"
if col in df.columns:
    before_nan = df[col].isna().mean()
    # OMNI/WDC Dst fill values: 9999, 99999, 9999.0
    df[col] = df[col].where(df[col].abs() < 900, np.nan)
    after_nan = df[col].isna().mean()
    valid = df[col].dropna()
    meta = {
        "title": "Hourly Dst Geomagnetic Storm Index",
        "source": "World Data Center for Geomagnetism, Kyoto",
        "references": "http://wdc.kugi.kyoto-u.ac.jp/dstdir/",
        "units": "nT",
        "time_range": f"{df.index.min()} / {df.index.max()}",
    }
    save_parquet(df, fp, meta)
    report("FIX", "dst_index.parquet",
           f"Masked fill values: NaN {before_nan:.0%} → {after_nan:.0%}, "
           f"valid range [{valid.min():.1f}, {valid.max():.1f}] nT")


# ── 5. MLS polar: investigate and fix NaN issue ───────────────────────────────
print("\n[5] MLS polar parquets")
for species in ["hno3", "n2o", "ozone", "temperature"]:
    fp = PROC / f"atmospheric/mls_{species}_polar.parquet"
    df = pd.read_parquet(fp)
    nan_pct = df.isnull().values.mean()
    if nan_pct > 0.95:
        # These are level columns — check if the values are there as object dtype
        for col in df.columns:
            vals = df[col]
            n_valid = vals.notna().sum()
            dtype = vals.dtype
            # Try forcing numeric conversion
            converted = pd.to_numeric(vals, errors="coerce")
            n_after = converted.notna().sum()
            if n_after > n_valid:
                df[col] = converted
        nan_after = df.isnull().values.mean()
        if nan_after < nan_pct:
            meta = {
                "title": f"Aura MLS {species.upper()} Polar (>=60N) Zonal Mean",
                "source": "NASA Aura Microwave Limb Sounder",
                "references": "https://mls.jpl.nasa.gov/data/",
                "columns": "columns = pressure level in hPa (polar zonal mean)",
                "time_range": f"{df.index.min()} / {df.index.max()}",
            }
            save_parquet(df, fp, meta)
            report("FIX", f"mls_{species}_polar.parquet",
                   f"Type coercion: NaN {nan_pct:.0%} → {nan_after:.0%}")
        else:
            # Check raw values
            sample = df.iloc[0].tolist()[:3]
            report("SKIP", f"mls_{species}_polar.parquet",
                   f"Still {nan_after:.0%} NaN after coerce — sample: {sample}")
    else:
        report("SKIP", f"mls_{species}_polar.parquet", f"NaN={nan_pct:.0%} OK")


# ── 6. OMNI fill-value masking ─────────────────────────────────────────────────
print("\n[6] OMNI solar wind fill values")
fp = PROC / "solar/omni_1min.parquet"
df = pd.read_parquet(fp)
# OMNI fill sentinels documented at https://omniweb.gsfc.nasa.gov/html/ow_data.html
OMNI_FILLS = {
    "flow_speed":  9999.9,
    "Vx": 99999.9, "Vy": 99999.9, "Vz": 99999.9,
    "proton_density": 999.99,
    "T": 9999999.0,
    "Bx": 9999.99, "By": 9999.99, "Bz": 9999.99,
    "B_scalar": 9999.99, "B_vector": 9999.99,
    "sigma_B": 9999.99,
    "Kp": 99, "R": 999, "DST": 99999, "AE": 9999,
    "proton_flux_1MeV": 999999.99,
    "proton_flux_2MeV": 99999.99,
    "proton_flux_4MeV": 99999.99,
    "proton_flux_10MeV": 99999.99,
    "proton_flux_30MeV": 99999.99,
    "proton_flux_60MeV": 99999.99,
}
# Generic: any numeric column with value ≥9990 is likely a fill
cols_fixed = 0
for col in df.select_dtypes(include="number").columns:
    fill_thresh = OMNI_FILLS.get(col, 9990.0)
    mask = df[col].abs() >= fill_thresh
    n_masked = int(mask.sum())
    if n_masked > 0:
        df[col] = df[col].where(~mask, np.nan)
        cols_fixed += 1
meta = {
    "title": "OMNI 1-Minute Solar Wind In-Situ (multi-spacecraft merged)",
    "source": "NASA/GSFC OMNIWeb",
    "references": "https://omniweb.gsfc.nasa.gov/",
    "units": "flow_speed km/s; density cm-3; T K; B nT",
    "time_range": f"{df.index.min()} / {df.index.max()}",
}
save_parquet(df, fp, meta)
speed = df["flow_speed"].dropna() if "flow_speed" in df.columns else pd.Series(dtype=float)
report("FIX", "omni_1min.parquet",
       f"Masked fill sentinels in {cols_fixed} cols; "
       f"flow_speed valid range [{speed.min():.0f},{speed.max():.0f}] km/s")


# ── 7. GOES R particle combined: same fill masking ────────────────────────────
print("\n[7] GOES-R particle fill values")
fp = PROC / "solar/goes_r_particle.parquet"
df = pd.read_parquet(fp)
for col in df.select_dtypes(include="number").columns:
    mask = df[col].abs() >= 1e30
    if mask.any():
        df[col] = df[col].where(~mask, np.nan)
nan_after = df.isnull().values.mean()
meta = {
    "title": "GOES-R SEISS Energetic Particle Observations (combined GOES 16/17/18)",
    "source": "NOAA NCEI / GOES-R Series",
    "references": "https://www.ngdc.noaa.gov/stp/satellite/goes-r.html",
    "time_range": f"{df.index.min()} / {df.index.max()}",
}
save_parquet(df, fp, meta)
report("FIX", "goes_r_particle.parquet", f"Post-fill-mask NaN={nan_after:.0%}")


# ── 8. slf_wet_model: fix epoch-zero index ────────────────────────────────────
print("\n[8] SLF wet model dates")
fp = PROC / "cryosphere/slf_wet_model.parquet"
df = pd.read_parquet(fp)
date_cols = [c for c in df.columns if "datum" in c.lower() or "date" in c.lower()]
if date_cols:
    col = date_cols[0]
    parsed = pd.to_datetime(df[col], utc=True, errors="coerce")
    n_valid = parsed.notna().sum()
    if n_valid > 100:
        df.index = parsed
        df.index.name = "time"
        df = df.drop(columns=[col]).sort_index()
        meta = {
            "title": "SLF Wet Snow Model Output",
            "source": "WSL/SLF Davos",
            "references": "https://www.slf.ch/",
            "time_range": f"{df.index.min()} / {df.index.max()}",
        }
        save_parquet(df, fp, meta)
        report("FIX", "slf_wet_model.parquet",
               f"{n_valid} valid dates parsed from '{col}'")
    else:
        report("SKIP", "slf_wet_model.parquet",
               f"'{col}' yielded only {n_valid} valid dates — keeping as-is")
else:
    report("SKIP", "slf_wet_model.parquet", "No datum/date column found")


# ── 9. NetCDF4 CF compliance: add Conventions + fix coord names ───────────────
print("\n[9] NetCDF4 CF compliance")

nc_fixes = {
    # ERA5 files use valid_time/latitude/longitude — rename to CF standard
    "atmospheric/era5_polar_strat_gridded.nc":  {"valid_time": "time", "latitude": "lat", "longitude": "lon"},
    "atmospheric/era5_swiss_alps.nc":           {"valid_time": "time", "latitude": "lat", "longitude": "lon"},
    # MLS files lack lon (zonal mean) — add note; rename lat/lev if needed
    "atmospheric/mls_hno3_gridded.nc":          {},
    "atmospheric/mls_n2o_gridded.nc":           {},
    "atmospheric/mls_ozone_gridded.nc":         {},
    "atmospheric/mls_temperature_gridded.nc":   {},
    # MODIS uses projected x/y — add note
    "cryosphere/modis_alps_monthly.nc":         {},
}

# Variable units to add if missing
VAR_UNITS = {
    "hno3": "ppbv", "n2o": "ppbv", "ozone": "ppbv", "temperature": "K",
    "ndsi_snow_cover": "fraction", "snow_cover_mean": "fraction",
    "expver": "1",  # ECMWF experiment version — dimensionless
    # ERA5 variables
    "t2m": "K", "sd": "m", "tp": "m", "sf": "m", "sp": "Pa",
    "msl": "Pa", "u10": "m s-1", "v10": "m s-1",
    "z": "m2 s-2", "t": "K", "u": "m s-1", "v": "m s-1",
    "w": "Pa s-1", "r": "%", "q": "kg kg-1",
}

import tempfile, shutil, datetime as _dt

for rel, rename_map in nc_fixes.items():
    fp = PROC / rel
    if not fp.exists():
        report("SKIP", rel, "file not found")
        continue
    if fp.stat().st_size < 1000:
        report("SKIP", rel, f"file too small ({fp.stat().st_size} bytes) — likely corrupted, skipping")
        continue
    tmp_path = None
    try:
        ds = xr.open_dataset(fp, engine="netcdf4")

        # Rename coords if needed
        rename_actual = {k: v for k, v in rename_map.items() if k in ds.coords or k in ds.dims}
        if rename_actual:
            ds = ds.rename(rename_actual)

        # Add/fix time encoding — move units out of attrs into encoding to avoid xarray conflict
        time_encoding = {}
        if "time" in ds.coords:
            # Remove units/calendar from attrs (xarray will manage via encoding)
            for attr_key in ("units", "calendar"):
                ds["time"].attrs.pop(attr_key, None)
            ds["time"].attrs["standard_name"] = "time"
            time_encoding = {"units": "hours since 1900-01-01 00:00:00",
                             "calendar": "proleptic_gregorian"}

        # Add lat/lon standard_names and units
        for cname, sname, cunits in [("lat", "latitude", "degrees_north"),
                                      ("lon", "longitude", "degrees_east"),
                                      ("lev", "air_pressure", "hPa")]:
            if cname in ds.coords:
                ds[cname].attrs["standard_name"] = sname
                ds[cname].attrs["units"] = cunits

        # Add units to data vars if missing
        for var in ds.data_vars:
            if not ds[var].attrs.get("units"):
                vl = var.lower()
                for key, unit in VAR_UNITS.items():
                    if key in vl:
                        ds[var].attrs["units"] = unit
                        break
                else:
                    ds[var].attrs["units"] = "1"  # CF: dimensionless

        # Add global CF attributes
        ds.attrs["Conventions"] = "CF-1.8"
        if not ds.attrs.get("institution"):
            ds.attrs["institution"] = "Solar-Magnetic Avalanche Research Project"
        if not ds.attrs.get("history"):
            ds.attrs["history"] = f"CF attributes added {_dt.date.today().isoformat()}"

        # Load into memory before closing
        ds = ds.load()
        ds.close()

        # Build encoding — write to TEMP file first, then atomically rename
        encoding = {}
        if time_encoding:
            encoding["time"] = time_encoding
        for var in list(ds.data_vars) + list(ds.coords):
            if var == "time":
                continue
            if ds[var].dtype in (np.float32, np.float64):
                enc = {"zlib": True, "complevel": 4, "dtype": "float32",
                       "_FillValue": float(np.finfo(np.float32).min)}
                if var in encoding:
                    encoding[var].update(enc)
                else:
                    encoding[var] = enc

        # Write to temp file in same directory, then atomically replace
        tmp_fd, tmp_path = tempfile.mkstemp(dir=fp.parent, suffix=".nc.tmp")
        os.close(tmp_fd)
        tmp_path = Path(tmp_path)
        ds.to_netcdf(str(tmp_path), encoding=encoding, format="NETCDF4")
        shutil.move(str(tmp_path), str(fp))
        tmp_path = None  # moved successfully

        sz = fp.stat().st_size / 1e6
        report("FIX", rel, f"CF-1.8 Conventions added, coords renamed {list(rename_actual)}, {sz:.0f} MB")
    except Exception as e:
        if tmp_path and Path(tmp_path).exists():
            Path(tmp_path).unlink(missing_ok=True)
        report("FAIL", rel, f"CF fix failed: {e!s:.120}")


# ── SUMMARY ──────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print(f"  FIXED: {len(fixed)}  |  SKIPPED/OK: {len(skipped)}")
print(f"{'='*60}")
for name in fixed:
    print(f"  ✓  {name}")
