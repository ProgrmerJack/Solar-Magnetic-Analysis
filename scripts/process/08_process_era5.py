"""
08_process_era5.py
Process ERA5 NetCDF files → merged compressed NetCDF4 + polar Parquet.
Swiss Alps: concatenate all years → era5_swiss_alps.nc (delete year files).
Polar stratosphere: compute polar means → parquet + gridded NC (delete year files).
Output: data/processed/atmospheric/
"""
import logging
import shutil
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

import sys
sys.path.insert(0, str(Path(__file__).parent))
from _utils import (
    DATA_ROOT, PROCESSED_ROOT, LOG, setup_logging,
    save_parquet, save_netcdf4, register_output, safe_delete, disk_free_gb,
)

MANIFEST = PROCESSED_ROOT / "manifest.json"
ATM_OUT = PROCESSED_ROOT / "atmospheric"

ERA5_DIR = DATA_ROOT / "atmospheric" / "era5"


def _area_weights(lats: np.ndarray) -> np.ndarray:
    """Cosine-of-latitude area weights, normalised to sum=1."""
    w = np.cos(np.deg2rad(lats))
    w = np.where(np.isnan(w), 0.0, w)
    total = w.sum()
    return w / total if total > 0 else w


def _open_alps_zip(fp: Path) -> "xr.Dataset | None":
    """
    Open a Swiss Alps ZIP archive (CDS API 2024+ format).
    Each ZIP contains two inner NetCDF4 files:
      - data_stream-oper_stepType-instant.nc  (t2m, u10, v10, sp, msl, d2m, sd)
      - data_stream-oper_stepType-accum.nc    (tp, sf)
    Extract to temp dir, merge, load into memory, clean up.
    """
    tmpdir = None
    try:
        tmpdir = Path(tempfile.mkdtemp())
        inner_dsets = []
        with zipfile.ZipFile(fp) as zf:
            for inner_name in zf.namelist():
                tmp_path = tmpdir / inner_name
                tmp_path.write_bytes(zf.read(inner_name))
                inner_dsets.append(
                    xr.open_dataset(str(tmp_path), engine="netcdf4", use_cftime=False)
                )
        merged = xr.merge(inner_dsets).load()   # load into RAM before temp cleanup
        for ds in inner_dsets:
            ds.close()
        return merged
    except Exception as exc:
        LOG.warning("Error opening Alps zip %s: %s", fp.name, exc)
        return None
    finally:
        if tmpdir and tmpdir.exists():
            shutil.rmtree(tmpdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Swiss Alps surface data
# ---------------------------------------------------------------------------
def process_swiss_alps() -> None:
    out = ATM_OUT / "era5_swiss_alps.nc"
    if out.exists():
        LOG.info("SKIP era5_swiss_alps.nc")
        return

    src_dir = ERA5_DIR / "swiss_alps"
    if not src_dir.exists():
        LOG.warning("Missing directory: %s", src_dir)
        return

    year_files = sorted(src_dir.glob("era5_swiss_alps_*.nc"))
    if not year_files:
        LOG.warning("No era5_swiss_alps_YYYY.nc files found in %s", src_dir)
        return

    LOG.info("Opening %d Swiss Alps ERA5 ZIP files …", len(year_files))

    frames: list[xr.Dataset] = []
    for fp in year_files:
        ds = _open_alps_zip(fp)
        if ds is not None:
            frames.append(ds)

    if not frames:
        LOG.warning("No Swiss Alps data loaded — check file format")
        return

    # Concatenate along time dimension
    time_dim = "valid_time" if "valid_time" in frames[0].coords else "time"
    ds_all = xr.concat(frames, dim=time_dim)

    # Add CF units to known variables
    cf_units = {"t2m": "K", "u10": "m s-1", "v10": "m s-1",
                "tp": "m", "sf": "m", "d2m": "K", "sp": "Pa",
                "msl": "Pa", "sd": "m of water equivalent"}
    for var, unit in cf_units.items():
        if var in ds_all:
            ds_all[var].attrs.setdefault("units", unit)

    meta = {
        "title": "ERA5 Reanalysis — Swiss Alps Surface (1959–2007)",
        "source": "ECMWF ERA5",
        "references": "https://doi.org/10.1002/qj.3803",
        "Conventions": "CF-1.8",
    }
    save_netcdf4(ds_all, out, meta)
    ds_all.close()

    if out.exists() and out.stat().st_size > 0:
        safe_delete(year_files)
        register_output(MANIFEST, "era5_swiss_alps", out, True, meta)
    else:
        LOG.error("Swiss Alps NC output missing/empty — year files NOT deleted")


# ---------------------------------------------------------------------------
# Polar stratosphere
# ---------------------------------------------------------------------------
def _compute_polar_means(ds: xr.Dataset, year: int) -> pd.DataFrame:
    """
    Compute area-weighted (lat≥60°N) + zonal mean for each pressure level and variable.
    Returns long-form DataFrame: year, month, pressure_level, var, value.
    """
    # Identify latitude coordinate name
    lat_name = None
    for name in ("latitude", "lat", "XLAT"):
        if name in ds.coords or name in ds:
            lat_name = name
            break
    if lat_name is None:
        LOG.warning("No latitude coordinate found in year %d", year)
        return pd.DataFrame()

    lats = ds[lat_name].values
    polar_mask = lats >= 60.0

    if not polar_mask.any():
        LOG.warning("No polar latitudes in year %d (lat range %.1f–%.1f)",
                    year, lats.min(), lats.max())
        return pd.DataFrame()

    weights = _area_weights(lats[polar_mask])

    # Identify time and pressure level coordinates
    time_name = None
    for name in ("time", "valid_time", "forecast_reference_time"):
        if name in ds.coords:
            time_name = name
            break

    lev_name = None
    for name in ("level", "pressure_level", "plev", "isobaricInhPa", "lev"):
        if name in ds.coords or name in ds:
            lev_name = name
            break

    if time_name is None:
        LOG.warning("No time coordinate in year %d", year)
        return pd.DataFrame()

    time_vals = pd.DatetimeIndex(ds[time_name].values, tz="UTC")
    records = []

    target_vars = [v for v in ds.data_vars if v.lower() in ("t", "u", "v", "z", "o3",
                                                               "temperature", "wind_u", "wind_v")]
    if not target_vars:
        target_vars = list(ds.data_vars)[:5]   # fallback to first 5

    for t_idx, t in enumerate(time_vals):
        row: dict = {"year": t.year, "month": t.month, "datetime": t}

        for var in target_vars:
            da = ds[var]
            try:
                # Select time
                t_slice = da.isel({time_name: t_idx})

                # Average over longitude if present
                lon_name = next((n for n in ("longitude", "lon") if n in t_slice.dims), None)
                if lon_name:
                    t_slice = t_slice.mean(dim=lon_name)

                # Select polar latitudes
                lat_dim = next((n for n in (lat_name, "latitude", "lat") if n in t_slice.dims), None)
                if lat_dim is None:
                    continue
                polar_slice = t_slice.sel({lat_dim: lats[polar_mask]})
                polar_vals = polar_slice.values

                if lev_name and lev_name in t_slice.dims:
                    levs = ds[lev_name].values
                    # polar_vals shape: (npolar, nlev) or (nlev, npolar)
                    if polar_vals.ndim == 2:
                        if polar_vals.shape[0] == len(weights):
                            # (npolar, nlev)
                            col_mean = np.average(polar_vals, axis=0, weights=weights)
                        else:
                            col_mean = np.average(polar_vals, axis=1, weights=weights)
                        for li, lv in enumerate(levs):
                            row[f"{var}_{int(lv)}hPa"] = float(col_mean[li])
                    else:
                        mean_val = np.average(polar_vals, weights=weights)
                        row[var] = float(mean_val)
                else:
                    if polar_vals.ndim == 1:
                        mean_val = np.average(polar_vals, weights=weights)
                    else:
                        mean_val = float(np.nanmean(polar_vals))
                    row[var] = float(mean_val)

            except Exception as exc:
                LOG.debug("Polar mean error var=%s t=%s: %s", var, t, exc)

        records.append(row)

    return pd.DataFrame(records)


def process_polar_strat() -> None:
    out_parquet = ATM_OUT / "era5_polar_strat_means.parquet"
    out_grid = ATM_OUT / "era5_polar_strat_gridded.nc"

    src_dir = ERA5_DIR / "polar_strat_monthly"
    if not src_dir.exists():
        LOG.warning("Missing directory: %s", src_dir)
        return

    year_files = sorted(src_dir.glob("era5_polar_strat_*.nc"))
    if not year_files:
        LOG.warning("No era5_polar_strat_YYYY.nc files found in %s", src_dir)
        return

    # --- Polar means Parquet ---
    if not out_parquet.exists():
        all_means: list[pd.DataFrame] = []
        success_files: list[Path] = []

        for fp in year_files:
            try:
                year_str = fp.stem.split("_")[-1]
                year = int(year_str)
            except ValueError:
                year = 0

            try:
                ds = xr.open_dataset(str(fp), engine="netcdf4", use_cftime=False)
                ds.close()
                if not mean_df.empty:
                    all_means.append(mean_df)
                success_files.append(fp)
            except Exception as exc:
                LOG.warning("Could not process %s: %s", fp.name, exc)

        if all_means:
            combined = pd.concat(all_means, ignore_index=True)
            if "datetime" in combined.columns:
                combined = combined.set_index("datetime").sort_index()

            meta = {
                "title": "ERA5 Polar Stratosphere Zonal+Meridional Means (≥60°N)",
                "source": "ECMWF ERA5",
                "references": "https://doi.org/10.1002/qj.3803",
                "time_range": f"{combined.index.min()} / {combined.index.max()}",
                "units": "T K; U/V m s-1; Z m2 s-2; O3 kg kg-1",
            }
            save_parquet(combined, out_parquet, meta)
            register_output(MANIFEST, "era5_polar_strat_means", out_parquet, False, meta)
        else:
            LOG.warning("No polar means computed")

    # --- Gridded NetCDF (polar ≥60°N only, loaded one file at a time) ---
    if not out_grid.exists():
        LOG.info("Building gridded ERA5 polar stratosphere NC (≥60°N slice) …")
        polar_chunks = []
        time_name = None
        lat_name = None
        for fp in year_files:
            try:
                ds = xr.open_dataset(str(fp), engine="netcdf4", use_cftime=False)
                # Determine coordinate names on first file
                if time_name is None:
                    time_name = next(
                        (n for n in ("valid_time", "time") if n in ds.coords), "time"
                    )
                    lat_name = next(
                        (n for n in ("latitude", "lat") if n in ds.coords), "latitude"
                    )
                # Select polar latitudes only (≥60°N)
                polar = ds.sel({lat_name: ds[lat_name] >= 60}).load()
                ds.close()
                polar_chunks.append(polar)
            except Exception as exc:
                LOG.warning("Could not open %s: %s", fp.name, exc)

        if polar_chunks:
            try:
                ds_all = xr.concat(polar_chunks, dim=time_name)
                del polar_chunks
                # Convert to float32 to reduce file size
                for var in ds_all.data_vars:
                    if ds_all[var].dtype == np.float64:
                        ds_all[var] = ds_all[var].astype(np.float32)
                meta_nc = {
                    "title": "ERA5 Polar Stratosphere Monthly Gridded ≥60°N (1975–2023)",
                    "source": "ECMWF ERA5",
                    "references": "https://doi.org/10.1002/qj.3803",
                    "Conventions": "CF-1.8",
                }
                save_netcdf4(ds_all, out_grid, meta_nc)
                ds_all.close()
                register_output(MANIFEST, "era5_polar_strat_gridded", out_grid, False, meta_nc)
            except Exception as exc:
                LOG.warning("Could not build gridded NC: %s", exc)
            except Exception as exc:
                LOG.warning("Could not build gridded NC: %s", exc)

    # Delete year files only if both outputs exist and are non-empty
    if (out_parquet.exists() and out_parquet.stat().st_size > 0 and
            out_grid.exists() and out_grid.stat().st_size > 0):
        safe_delete(year_files)
    else:
        LOG.warning("Not deleting ERA5 polar strat year files — one or both outputs missing/empty")


def main() -> None:
    setup_logging()
    LOG.info("=== 08_process_era5.py | disk free=%.1f GB ===", disk_free_gb())
    ATM_OUT.mkdir(parents=True, exist_ok=True)

    process_swiss_alps()
    LOG.info("After Swiss Alps | disk free=%.1f GB", disk_free_gb())

    process_polar_strat()
    LOG.info("=== 08 complete | disk free=%.1f GB ===", disk_free_gb())


if __name__ == "__main__":
    main()
