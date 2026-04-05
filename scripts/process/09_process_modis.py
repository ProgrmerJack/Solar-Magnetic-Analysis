"""
09_process_modis.py
Process MODIS MOD10A1 HDF4 snow cover files → daily stats Parquet + monthly grid NetCDF4.
Deletes raw HDF4 files after processing.
Output: data/processed/cryosphere/
"""
import logging
import re
from datetime import datetime, timezone
from itertools import groupby
from pathlib import Path

import netCDF4 as nc4
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
CRYO_OUT = PROCESSED_ROOT / "cryosphere"

MODIS_DIR = DATA_ROOT / "cryosphere" / "modis_snow" / "alps_h18v04"

# NDSI special flag values
NDSI_MISSING = 200
NDSI_NO_DECISION = 201
NDSI_NIGHT = 211
NDSI_INLAND_WATER = 237
NDSI_OCEAN = 239
NDSI_CLOUD = 250
NDSI_SATURATED = 254
NDSI_FILL = 255
NDSI_SNOW_THRESHOLD = 10   # NDSI ≥ 10 = snow-covered

# Valid observation mask: integer value 0-100 AND QA ≤ 2
QA_MAX_VALID = 2


def _parse_date_from_filename(fname: str) -> datetime | None:
    """Parse YYYYDDD from MOD10A1.A{YYYYDDD}.h18v04.*.hdf filename."""
    m = re.search(r"\.A(\d{7})\.", fname)
    if not m:
        return None
    try:
        return datetime.strptime(m.group(1), "%Y%j").replace(tzinfo=timezone.utc)
    except ValueError:
        return None


def _read_modis_hdf(fp: Path) -> tuple[np.ndarray, np.ndarray] | None:
    """
    Read NDSI_Snow_Cover and NDSI_Snow_Cover_Basic_QA from HDF4 file.
    Returns (ndsi_array, qa_array) or None.
    """
    try:
        ds = nc4.Dataset(str(fp), "r")
    except Exception as exc:
        LOG.warning("Cannot open HDF4 %s: %s", fp.name, exc)
        return None

    try:
        # Variable may be in a group — try both root and 'MOD_Grid_Snow_500m'
        ndsi = None
        qa = None

        def _find_var(dataset, name_hints):
            for hint in name_hints:
                if hint in dataset.variables:
                    return dataset.variables[hint][:]
            for var in dataset.variables:
                if any(h.lower() in var.lower() for h in name_hints):
                    return dataset.variables[var][:]
            return None

        ndsi = _find_var(ds, ["NDSI_Snow_Cover", "NDSI", "Snow_Cover"])
        qa = _find_var(ds, ["NDSI_Snow_Cover_Basic_QA", "Basic_QA", "QA"])

        if ndsi is None:
            # Try sub-groups
            for grp_name in ds.groups:
                grp = ds.groups[grp_name]
                ndsi = _find_var(grp, ["NDSI_Snow_Cover", "NDSI", "Snow_Cover"])
                qa = _find_var(grp, ["NDSI_Snow_Cover_Basic_QA", "Basic_QA", "QA"])
                if ndsi is not None:
                    break

        if ndsi is None:
            LOG.warning("No NDSI variable found in %s", fp.name)
            return None

        ndsi = np.asarray(ndsi, dtype=np.uint8).squeeze()
        if qa is not None:
            qa = np.asarray(qa, dtype=np.uint8).squeeze()
        else:
            # If QA unavailable, assume all valid-value pixels pass QA
            qa = np.zeros_like(ndsi)

        ds.close()
        return ndsi, qa
    except Exception as exc:
        LOG.warning("Error reading %s: %s", fp.name, exc)
        try:
            ds.close()
        except Exception:
            pass
        return None


def _compute_daily_stats(ndsi: np.ndarray, qa: np.ndarray) -> dict:
    """Compute per-file statistics for the daily MODIS granule."""
    total_pixels = ndsi.size

    # Valid: data value 0-100 AND QA ≤ 2
    valid_mask = (ndsi <= 100) & (qa <= QA_MAX_VALID)
    n_valid = int(valid_mask.sum())
    valid_fraction = n_valid / total_pixels if total_pixels > 0 else np.nan

    # Snow: valid AND NDSI ≥ snow threshold
    snow_mask = valid_mask & (ndsi >= NDSI_SNOW_THRESHOLD)
    n_snow = int(snow_mask.sum())
    snow_fraction = n_snow / n_valid if n_valid > 0 else np.nan

    # Mean NDSI over snow pixels
    mean_ndsi = float(ndsi[snow_mask].mean()) if n_snow > 0 else np.nan

    # Cloud fraction (pixel value == 250)
    cloud_fraction = float((ndsi == NDSI_CLOUD).sum()) / total_pixels

    return {
        "snow_fraction": snow_fraction,
        "mean_ndsi": mean_ndsi,
        "cloud_fraction": cloud_fraction,
        "valid_fraction": valid_fraction,
        "n_valid_pixels": n_valid,
    }


def process_modis_stats() -> tuple:
    """
    Process all MOD10A1 files streaming month-by-month to avoid OOM.
    Accumulates only scalar stats (trivial memory) + computes monthly grid means
    on-the-fly (at most ~31 daily 2400×2400 arrays in RAM at once = ~700 MB peak).

    Returns (stats_parquet_path, processed_raw_paths).
    """
    out = CRYO_OUT / "modis_alps_stats.parquet"

    hdf_files = sorted(MODIS_DIR.glob("MOD10A1.A*.h18v04.*.hdf"))
    if not hdf_files:
        LOG.warning("No MODIS HDF files found in %s", MODIS_DIR)
        return out, []

    LOG.info("Processing %d MODIS MOD10A1 files …", len(hdf_files))

    records: list[dict] = []
    processed_raw: list[Path] = []
    monthly_means: list[tuple] = []   # (time_np64, mean_grid_float32)

    # Group by (year, month) using a two-pass sort so we release monthly arrays immediately
    from itertools import groupby
    def _ym_key(fp):
        d = _parse_date_from_filename(fp.name)
        return (d.year, d.month) if d else (0, 0)

    valid_files = [(fp, _parse_date_from_filename(fp.name))
                   for fp in hdf_files if _parse_date_from_filename(fp.name) is not None]
    valid_files.sort(key=lambda x: (x[1].year, x[1].month, x[1].day))

    # Group by month and process one month at a time
    for (year, month), group in groupby(valid_files, key=lambda x: (x[1].year, x[1].month)):
        month_day_arrays: list[np.ndarray] = []

        for fp, date in group:
            result = _read_modis_hdf(fp)
            if result is None:
                continue

            ndsi, qa = result
            stats = _compute_daily_stats(ndsi, qa)
            stats["date"] = pd.Timestamp(date).tz_localize(None).tz_localize("UTC")
            records.append(stats)

            # Valid NDSI for monthly mean (float32 to halve RAM vs float64)
            valid_ndsi = ndsi.astype(np.float32)
            valid_ndsi[(ndsi > 100) | (qa > QA_MAX_VALID)] = np.nan
            month_day_arrays.append(valid_ndsi)
            processed_raw.append(fp)

        # Compute and store monthly mean, then immediately release daily arrays
        if month_day_arrays:
            stacked = np.stack(month_day_arrays, axis=0)     # (ndays, 2400, 2400)
            monthly_mean = np.nanmean(stacked, axis=0).astype(np.float32)
            t = pd.Timestamp(year=year, month=month, day=1, tz="UTC").to_datetime64()
            monthly_means.append((t, monthly_mean))
            del stacked, month_day_arrays   # free RAM immediately

    if not records:
        LOG.warning("No MODIS stats computed")
        return out, processed_raw

    if not out.exists():
        df = pd.DataFrame(records).set_index("date").sort_index()
        meta = {
            "title": "MODIS MOD10A1 Snow Cover Daily Statistics — Swiss Alps (h18v04)",
            "source": "NASA NSIDC MOD10A1.061",
            "references": "https://doi.org/10.5067/MODIS/MOD10A1.061",
            "time_range": f"{df.index.min()} / {df.index.max()}",
            "units": "fractions dimensionless; mean_ndsi 0-100",
        }
        save_parquet(df, out, meta)
        register_output(MANIFEST, "modis_alps_stats", out, True, meta)

    return out, processed_raw, monthly_means


def process_modis_monthly_grid(monthly_means: list) -> None:
    """
    Save monthly mean NDSI grids as compressed NetCDF4.
    Input: list of (time_np64, mean_grid_float32) tuples — already computed per-month.
    Memory peak = one full array stack (nmonths × 2400 × 2400 float32 ≈ 6 GB for 25 yrs).
    Written with deflate=4 → output typically 200–600 MB.
    """
    out = CRYO_OUT / "modis_alps_monthly.nc"
    if out.exists():
        LOG.info("SKIP modis_alps_monthly.nc")
        return
    if not monthly_means:
        return

    time_vals = np.array([t for t, _ in monthly_means], dtype="datetime64[ns]")
    data_arr = np.stack([g for _, g in monthly_means], axis=0)   # (nmonths, 2400, 2400)
    ny, nx = data_arr.shape[1], data_arr.shape[2]

    ds = xr.Dataset(
        {"ndsi_snow_cover": (["time", "y", "x"], data_arr)},
        coords={
            "time": (["time"], time_vals),
            "y": (["y"], np.arange(ny)),
            "x": (["x"], np.arange(nx)),
        },
    )
    ds["ndsi_snow_cover"].attrs.update({
        "long_name": "Monthly mean NDSI snow cover fraction",
        "units": "0-100 (NDSI scale)",
        "valid_range": [0, 100],
        "grid_mapping": "sinusoidal",
        "comment": ("MODIS tile h18v04, 500m sinusoidal projection. "
                    "Values 0-100 are NDSI; NaN = cloud/fill/night."),
    })
    ds["time"].attrs.update({"standard_name": "time", "axis": "T"})
    ds["y"].attrs.update({"long_name": "tile row index (500m)", "units": "1"})
    ds["x"].attrs.update({"long_name": "tile column index (500m)", "units": "1"})

    meta = {
        "title": "MODIS MOD10A1 Monthly Mean Snow Cover Grid — Swiss Alps (h18v04)",
        "source": "NASA NSIDC MOD10A1.061",
        "references": "https://doi.org/10.5067/MODIS/MOD10A1.061",
        "Conventions": "CF-1.8",
    }
    save_netcdf4(ds, out, meta)
    register_output(MANIFEST, "modis_alps_monthly", out, False, meta)


def main() -> None:
    setup_logging()
    LOG.info("=== 09_process_modis.py | disk free=%.1f GB ===", disk_free_gb())
    CRYO_OUT.mkdir(parents=True, exist_ok=True)

    if not MODIS_DIR.exists():
        LOG.warning("MODIS directory not found: %s", MODIS_DIR)
        return

    # process_modis_stats returns (stats_out, processed_raw, monthly_means)
    result = process_modis_stats()
    if len(result) == 3:
        stats_out, processed_raw, monthly_means = result
    else:
        stats_out, processed_raw = result
        monthly_means = []

    process_modis_monthly_grid(monthly_means)

    if stats_out.exists() and stats_out.stat().st_size > 0:
        safe_delete(processed_raw)
    else:
        LOG.error("Stats parquet missing/empty — raw files NOT deleted")

    LOG.info("=== 09 complete | disk free=%.1f GB ===", disk_free_gb())


if __name__ == "__main__":
    main()
