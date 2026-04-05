"""
04_process_goes_xrs.py
Process GOES X-ray Sensor (XRS) NetCDF files → single Parquet.
Deletes raw NC files after verified write (~4400 files → 1 Parquet).
Output: data/processed/solar/
"""
import logging
from pathlib import Path

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).parent))
from _utils import (
    DATA_ROOT, PROCESSED_ROOT, LOG, setup_logging,
    save_parquet, register_output, safe_delete, disk_free_gb,
)

MANIFEST = PROCESSED_ROOT / "manifest.json"
SOL_OUT = PROCESSED_ROOT / "solar"

GOES_XRS_DIR = DATA_ROOT / "solar" / "goes_xrs"

# Variable name candidates for Channel A (0.05-0.4 nm) and Channel B (0.1-0.8 nm)
CHAN_A_NAMES = ["xrsa_flux", "xrsa", "a_flux", "flux_observed_0094nm", "A_FLUX",
                "xs", "XRS_A", "irradiance_a", "irr_obs_short"]
CHAN_B_NAMES = ["xrsb_flux", "xrsb", "b_flux", "flux_observed_0160nm", "B_FLUX",
                "xl", "XRS_B", "irradiance_b", "irr_obs_long"]
TIME_NAMES = ["time", "TIME", "time_coverage_start", "record_times"]
FILL_THRESHOLD = 1e30


def _find_var(ds: xr.Dataset, candidates: list[str]) -> str | None:
    """Return first candidate that exists in *ds*, or None."""
    for name in candidates:
        if name in ds:
            return name
    # Fuzzy fallback: substring match
    for var in ds.data_vars:
        vl = var.lower()
        if any(k in vl for k in ("xrsa", "flux_a", "chan_a", "short", "094")):
            if candidates is CHAN_A_NAMES:
                return var
        if any(k in vl for k in ("xrsb", "flux_b", "chan_b", "long", "160")):
            if candidates is CHAN_B_NAMES:
                return var
    return None


def _open_nc_to_dataframe(fp: Path, satellite_id: str) -> pd.DataFrame | None:
    """Open a single GOES XRS NetCDF using netCDF4 (fast) and return tidy DataFrame or None."""
    import netCDF4 as nc4
    try:
        ds = nc4.Dataset(fp, "r")
    except Exception as exc:
        LOG.warning("Could not open %s: %s", fp.name, exc)
        return None

    try:
        # --- time ---
        time_var = None
        for tname in ("time", "TIME", "record_times"):
            if tname in ds.variables:
                time_var = ds.variables[tname]
                break
        if time_var is None:
            ds.close()
            return None

        # Convert CF time to pandas DatetimeIndex
        units = getattr(time_var, "units", None)
        calendar = getattr(time_var, "calendar", "standard")
        if units:
            import cftime
            times_raw = nc4.num2date(time_var[:], units=units, calendar=calendar)
            idx = pd.DatetimeIndex(
                [pd.Timestamp(t.year, t.month, t.day, t.hour, t.minute, t.second, tzinfo=None)
                 for t in times_raw],
                tz="UTC"
            )
        else:
            idx = pd.DatetimeIndex(time_var[:].astype("datetime64[s]"), tz="UTC")

        df = pd.DataFrame(index=idx)

        # --- channel A ---
        chan_a = None
        for cname in CHAN_A_NAMES:
            if cname in ds.variables:
                chan_a = cname
                break
        if chan_a is None:
            for vname in ds.variables:
                vl = vname.lower()
                if any(k in vl for k in ("xrsa", "flux_a", "chan_a", "094")):
                    chan_a = vname
                    break

        # --- channel B ---
        chan_b = None
        for cname in CHAN_B_NAMES:
            if cname in ds.variables:
                chan_b = cname
                break
        if chan_b is None:
            for vname in ds.variables:
                vl = vname.lower()
                if any(k in vl for k in ("xrsb", "flux_b", "chan_b", "160")):
                    chan_b = vname
                    break

        if chan_a is None and chan_b is None:
            ds.close()
            return None

        if chan_a:
            vals = np.array(ds.variables[chan_a][:], dtype=float).flatten()[:len(df)]
            fill = getattr(ds.variables[chan_a], "_FillValue", FILL_THRESHOLD)
            vals[np.abs(vals) >= min(abs(fill), FILL_THRESHOLD)] = np.nan
            df["xrsa_flux"] = vals

        if chan_b:
            vals = np.array(ds.variables[chan_b][:], dtype=float).flatten()[:len(df)]
            fill = getattr(ds.variables[chan_b], "_FillValue", FILL_THRESHOLD)
            vals[np.abs(vals) >= min(abs(fill), FILL_THRESHOLD)] = np.nan
            df["xrsb_flux"] = vals

        df["satellite_id"] = satellite_id
        ds.close()
        return df

    except Exception as exc:
        LOG.warning("Error extracting %s: %s", fp.name, exc)
        try:
            ds.close()
        except Exception:
            pass
        return None


def _find_nc_files_for_satellite(sat_dir: Path) -> list[Path]:
    """Return all NetCDF files for a satellite directory (legacy or GOES-R layout)."""
    return sorted(sat_dir.rglob("*.nc"))


def process_goes_xrs() -> None:
    out = SOL_OUT / "goes_xrs.parquet"
    if out.exists():
        LOG.info("SKIP goes_xrs.parquet")
        return

    if not GOES_XRS_DIR.exists():
        LOG.warning("Missing directory: %s", GOES_XRS_DIR)
        return

    sat_parquets: list[Path] = []

    sat_dirs = sorted(d for d in GOES_XRS_DIR.iterdir() if d.is_dir())
    if not sat_dirs:
        LOG.warning("No satellite subdirectories found in %s", GOES_XRS_DIR)
        return

    for sat_dir in sat_dirs:
        satellite_id = sat_dir.name
        sat_out = SOL_OUT / f"goes_xrs_{satellite_id}.parquet"

        if sat_out.exists():
            LOG.info("SKIP %s (already processed)", sat_out.name)
            sat_parquets.append(sat_out)
            continue

        nc_files = _find_nc_files_for_satellite(sat_dir)
        if not nc_files:
            LOG.warning("No NC files for %s", satellite_id)
            continue

        LOG.info("Processing %s: %d files …", satellite_id, len(nc_files))
        sat_frames: list[pd.DataFrame] = []
        good_files: list[Path] = []

        for fp in nc_files:
            df = _open_nc_to_dataframe(fp, satellite_id)
            if df is not None and not df.empty:
                sat_frames.append(df)
                good_files.append(fp)

        if not sat_frames:
            LOG.warning("  %s: no data extracted", satellite_id)
            continue

        sat_df = pd.concat(sat_frames).sort_index()
        sat_df = sat_df[~sat_df.index.duplicated(keep="first")]

        sat_meta = {
            "title": f"GOES XRS 1-min X-ray Flux ({satellite_id})",
            "source": "NOAA/NCEI GOES XRS",
            "references": "https://doi.org/10.7289/V5BV7DSR",
            "time_range": f"{sat_df.index.min()} / {sat_df.index.max()}",
            "units": "W/m^2",
        }
        save_parquet(sat_df, sat_out, sat_meta)

        if sat_out.exists() and sat_out.stat().st_size > 0:
            safe_delete(good_files)
            LOG.info("  %s: %d rows → %s (raw deleted)", satellite_id, len(sat_df), sat_out.name)
            sat_parquets.append(sat_out)
        else:
            LOG.error("  %s: output empty — raw files NOT deleted", satellite_id)

    if not sat_parquets:
        LOG.warning("No GOES XRS data produced")
        return

    # Merge all satellite Parquets
    LOG.info("Merging %d satellite Parquets …", len(sat_parquets))
    frames = [pd.read_parquet(p) for p in sat_parquets]
    combined = pd.concat(frames).sort_index()
    numeric_cols = [c for c in combined.columns if c != "satellite_id"]
    resampled = combined[numeric_cols].resample("1min").mean()
    sat_resampled = combined["satellite_id"].resample("1min").last()
    resampled["satellite_id"] = sat_resampled

    meta = {
        "title": "GOES XRS 1-min X-ray Flux (all satellites, combined)",
        "source": "NOAA/NCEI GOES XRS",
        "references": "https://doi.org/10.7289/V5BV7DSR",
        "time_range": f"{resampled.index.min()} / {resampled.index.max()}",
        "units": "W/m^2; Channel A 0.05-0.4nm; Channel B 0.1-0.8nm",
    }
    save_parquet(resampled, out, meta)

    if out.exists() and out.stat().st_size > 0:
        # Delete per-satellite intermediate Parquets
        for p in sat_parquets:
            try:
                p.unlink()
            except Exception:
                pass
        register_output(MANIFEST, "goes_xrs", out, True, meta)
    else:
        LOG.error("Combined GOES XRS parquet missing/empty")


def main() -> None:
    setup_logging()
    LOG.info("=== 04_process_goes_xrs.py | disk free=%.1f GB ===", disk_free_gb())
    SOL_OUT.mkdir(parents=True, exist_ok=True)

    process_goes_xrs()

    LOG.info("=== 04 complete | disk free=%.1f GB ===", disk_free_gb())


if __name__ == "__main__":
    main()
