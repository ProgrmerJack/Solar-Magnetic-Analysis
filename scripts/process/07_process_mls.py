"""
07_process_mls.py
Process Aura/MLS HDF5 files → polar-mean Parquet + gridded NetCDF4.
Does NOT delete raw MLS files (critical EPP proxy data).
Output: data/processed/atmospheric/
"""
import logging
from datetime import datetime, timedelta, timezone
from pathlib import Path

import netCDF4 as nc4
import numpy as np
import pandas as pd
import xarray as xr

import sys
sys.path.insert(0, str(Path(__file__).parent))
from _utils import (
    DATA_ROOT, PROCESSED_ROOT, LOG, setup_logging,
    save_parquet, save_netcdf4, register_output, disk_free_gb,
)

MANIFEST = PROCESSED_ROOT / "manifest.json"
ATM_OUT = PROCESSED_ROOT / "atmospheric"
MLS_DIR = DATA_ROOT / "atmospheric" / "aura_mls"

# MLS reference epoch: days since 1993-01-01
MLS_EPOCH = datetime(1993, 1, 1, tzinfo=timezone.utc)

# EPP-relevant pressure levels (hPa)
EPP_LEVELS_HPA = [1.0, 2.0, 3.2, 4.6, 6.8, 10.0]

# Product configuration: (subdir_hint, group_name, short_name)
MLS_PRODUCTS = [
    ("HNO3",       "HNO3 PressureZM",        "hno3"),
    ("Temperature","Temperature PressureZM",  "temperature"),
    ("N2O",        "N2O PressureZM",          "n2o"),
    ("O3",         "O3 PressureZM",           "ozone"),
]


def _find_mls_files(product: str) -> list[Path]:
    """Return sorted files for a given MLS product (NetCDF4 .nc or legacy HDF5)."""
    import re
    # Use word-boundary match so "O3" doesn't match inside "HNO3"
    pattern = re.compile(r'(?<![A-Za-z0-9])' + re.escape(product) + r'(?![A-Za-z0-9])',
                         re.IGNORECASE)
    candidates: list[Path] = []
    for subdir in sorted(MLS_DIR.iterdir()):
        if not subdir.is_dir():
            continue
        if pattern.search(subdir.name):
            for ext in ("*.nc", "*.he5", "*.h5", "*.hdf5"):
                candidates.extend(sorted(subdir.glob(ext)))
    # Also search MLS_DIR root directly
    for ext in ("*.nc", "*.he5", "*.h5", "*.hdf5"):
        for fp in MLS_DIR.glob(ext):
            if pattern.search(fp.name) and fp not in candidates:
                candidates.append(fp)
    return sorted(set(candidates))


def _mls_time_to_utc(time_arr: np.ndarray) -> pd.DatetimeIndex:
    """Convert MLS time (days since 1993-01-01) to UTC DatetimeIndex."""
    timestamps = [MLS_EPOCH + timedelta(days=float(t)) for t in time_arr]
    return pd.DatetimeIndex(timestamps, tz="UTC")


def _weighted_polar_mean(data: np.ndarray, lats: np.ndarray, min_lat: float = 60.0) -> np.ndarray:
    """
    Compute area-weighted mean over latitudes ≥ min_lat, ignoring NaN values.
    data shape: (nlat, ntime) or (ntime, nlat) or (nlat, 1) for single timestep
    Returns: (ntime,)
    """
    polar_mask = lats >= min_lat
    if not polar_mask.any():
        LOG.warning("No latitudes ≥ %.0f°N found (lat range: %.1f–%.1f)",
                    min_lat, lats.min(), lats.max())
        return np.full(data.shape[-1], np.nan)

    base_weights = np.cos(np.deg2rad(lats[polar_mask]))

    if data.ndim == 2:
        if data.shape[0] == len(lats):
            polar_data = data[polar_mask, :]          # (npolar, ntime)
        elif data.shape[1] == len(lats):
            polar_data = data[:, polar_mask].T        # (npolar, ntime)
        else:
            return np.full(1, np.nan)

        ntime = polar_data.shape[1]
        result = np.full(ntime, np.nan)
        for t in range(ntime):
            col = polar_data[:, t]
            valid = ~np.isnan(col)
            if valid.sum() >= 2:                      # at least 2 polar lats with data
                w = base_weights[valid]
                w = w / w.sum()
                result[t] = float(np.sum(w * col[valid]))
        return result

    return np.full(1, np.nan)


def _find_closest_levels(lev: np.ndarray, targets: list[float]) -> dict[float, int]:
    """Return {target_hPa: index_in_lev} for closest matches."""
    mapping = {}
    for t in targets:
        idx = int(np.argmin(np.abs(lev - t)))
        mapping[t] = idx
    return mapping


def _open_mls_group(fp: Path, group_name: str) -> tuple | None:
    """
    Open a NetCDF4 MLS file and return (lat, lev, time_index, value).
    value shape: (ntime, nlev, nlat) — confirmed from Aura/MLS L3DZ/L3MB files.
    time_index: pd.DatetimeIndex (UTC), decoded from CF units in file.
    """
    ds = None
    try:
        ds = nc4.Dataset(str(fp), "r")

        # Locate the PressureZM group (prefer exact match, then fuzzy)
        grp = None
        if group_name in ds.groups:
            grp = ds.groups[group_name]
        else:
            product_word = group_name.split()[0].lower()
            for gname, g in ds.groups.items():
                gname_norm = gname.lower().replace(" ", "")
                # Use startswith to avoid "o3" matching inside "hno3pressurezm"
                if gname_norm.startswith(product_word) and "pressurezm" in gname_norm:
                    grp = g
                    break
            if grp is None:
                for gname, g in ds.groups.items():
                    if gname.lower().startswith(product_word):
                        grp = g
                        break

        if grp is None:
            LOG.debug("Group '%s' not found in %s (available: %s)",
                      group_name, fp.name, list(ds.groups))
            return None

        lat  = np.asarray(grp.variables["lat"][:], dtype=float)
        lev  = np.asarray(grp.variables["lev"][:], dtype=float)

        # Decode CF time using units stored in file (e.g. "days since 1950-01-01")
        time_var  = grp.variables["time"]
        time_num  = np.asarray(time_var[:])
        time_units = time_var.units
        time_cal   = getattr(time_var, "calendar", "standard")
        dates_cf   = nc4.num2date(time_num, time_units, time_cal)
        timestamps = []
        for d in dates_cf:
            try:
                timestamps.append(pd.Timestamp(
                    year=d.year, month=d.month, day=d.day,
                    hour=getattr(d, "hour", 0), minute=getattr(d, "minute", 0),
                    tz="UTC",
                ))
            except Exception:
                timestamps.append(pd.NaT)
        time_index = pd.DatetimeIndex(timestamps)

        # value shape: (ntime, nlev, nlat)
        val_raw = grp.variables["value"][:]
        value   = np.ma.filled(np.ma.asarray(val_raw, dtype=float), fill_value=np.nan)
        # Mask MLS fill sentinels: any value ≤ -900 is below-detection-limit or fill
        value[value <= -900.0] = np.nan

        return (lat, lev, time_index, value)

    except Exception as exc:
        LOG.warning("Error opening %s: %s", fp.name, exc)
        return None
    finally:
        if ds is not None:
            try:
                ds.close()
            except Exception:
                pass


def process_mls_product(
    product_dir_hint: str,
    group_name: str,
    short_name: str,
) -> None:
    out_polar = ATM_OUT / f"mls_{short_name}_polar.parquet"
    out_grid = ATM_OUT / f"mls_{short_name}_gridded.nc"

    if out_polar.exists() and out_grid.exists():
        LOG.info("SKIP mls_%s (both outputs exist)", short_name)
        return

    files = _find_mls_files(product_dir_hint)
    if not files:
        LOG.warning("No MLS files found for product %s", product_dir_hint)
        return

    LOG.info("Processing MLS %s: %d files …", product_dir_hint, len(files))

    polar_records: list[dict] = []
    grid_lats: np.ndarray | None = None
    grid_levs: np.ndarray | None = None
    grid_times: list[np.datetime64] = []
    grid_values: list[np.ndarray] = []

    for fp in files:
        result = _open_mls_group(fp, group_name)
        if result is None:
            continue
        lat, lev, time_index, value = result   # time_index is pd.DatetimeIndex

        # Store coordinate arrays from first file
        if grid_lats is None:
            grid_lats = lat
            grid_levs = lev

        # value shape: (ntime, nlev, nlat)
        level_map = _find_closest_levels(lev, EPP_LEVELS_HPA)

        for t_idx, t in enumerate(time_index):
            if pd.isna(t):
                continue
            row: dict = {"time": t}
            for target_hPa, lev_idx in level_map.items():
                try:
                    # Extract lat profile at (t_idx, lev_idx)
                    lat_profile = value[t_idx, lev_idx, :]     # (nlat,)
                    polar_mean = _weighted_polar_mean(
                        lat_profile.reshape(-1, 1), lat
                    )
                    row[f"lev_{target_hPa:.1f}hpa".replace(".", "p")] = float(polar_mean[0])
                except Exception as exc:
                    LOG.debug("Level %s extraction error: %s", target_hPa, exc)
                    row[f"lev_{target_hPa:.1f}hpa".replace(".", "p")] = np.nan

            polar_records.append(row)

        # Accumulate gridded time slices: one (nlev, nlat) array per timestep
        for t_idx in range(len(time_index)):
            if not pd.isna(time_index[t_idx]):
                grid_times.append(time_index[t_idx].to_datetime64())
                grid_values.append(value[t_idx, :, :])   # (nlev, nlat)

    # --- Save polar Parquet ---
    if not out_polar.exists() and polar_records:
        df = pd.DataFrame(polar_records).set_index("time").sort_index()
        meta = {
            "title": f"Aura/MLS {product_dir_hint} Polar Cap Mean (≥60°N)",
            "source": "NASA Aura/MLS Level 3",
            "references": "https://disc.gsfc.nasa.gov/datasets?project=Aura",
            "time_range": f"{df.index.min()} / {df.index.max()}",
            "units": "ppbv (mixing ratio) or K (Temperature)",
        }
        save_parquet(df, out_polar, meta)
        register_output(MANIFEST, f"mls_{short_name}_polar", out_polar, False, meta)

    # --- Save gridded NetCDF4 ---
    if not out_grid.exists() and grid_values and grid_lats is not None and grid_levs is not None:
        try:
            data_arr = np.stack(grid_values, axis=0)    # (ntime, nlev, nlat)
            ds = xr.Dataset(
                {short_name: (["time", "lev", "lat"], data_arr)},
                coords={
                    "time": (["time"], np.array(grid_times, dtype="datetime64[ns]")),
                    "lev":  (["lev"],  grid_levs),
                    "lat":  (["lat"],  grid_lats),
                },
            )
            ds["time"].attrs["standard_name"] = "time"
            ds["lat"].attrs.update({"standard_name": "latitude", "units": "degrees_north"})
            ds["lev"].attrs.update({"standard_name": "air_pressure", "units": "hPa",
                                     "positive": "down"})
            ds[short_name].attrs["long_name"] = f"MLS {product_dir_hint}"

            meta_nc = {
                "title": f"Aura/MLS {product_dir_hint} Gridded Zonal Mean",
                "source": "NASA Aura/MLS Level 3",
                "references": "https://disc.gsfc.nasa.gov/datasets?project=Aura",
            }
            save_netcdf4(ds, out_grid, meta_nc)
            register_output(MANIFEST, f"mls_{short_name}_gridded", out_grid, False, meta_nc)
        except Exception as exc:
            LOG.warning("Could not save gridded NC for %s: %s", short_name, exc)


def main() -> None:
    setup_logging()
    LOG.info("=== 07_process_mls.py | disk free=%.1f GB ===", disk_free_gb())
    ATM_OUT.mkdir(parents=True, exist_ok=True)

    for prod_hint, group_name, short_name in MLS_PRODUCTS:
        process_mls_product(prod_hint, group_name, short_name)

    LOG.info("=== 07 complete | disk free=%.1f GB ===", disk_free_gb())


if __name__ == "__main__":
    main()
