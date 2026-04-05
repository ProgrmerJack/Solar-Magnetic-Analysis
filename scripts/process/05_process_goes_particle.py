"""
05_process_goes_particle.py
Process GOES energetic particle data (legacy MAGED + GOES-R MPSH) → Parquet.

GOES-R MPSH confirmed file structure:
  {sat_dir}/mpsh_avg1m/{year}/{month}/{daily.nc}   (3 levels, use rglob)

GOES-R MPSH confirmed variables (shape per 1-min file, 1440 time steps):
  time                    (1440,)        seconds since 2000-01-01 12:00:00 UTC
  AvgDiffProtonFlux       (1440, 5, 11)  telescope × energy_channel  protons/(cm² sr keV s)
  AvgDiffElectronFlux     (1440, 5, 10)  telescope × energy_channel  electrons/(cm² sr keV s)
  AvgIntElectronFlux      (1440, 5)      telescope                    electrons/(cm² sr s)
  DiffProtonEffectiveEnergy (5, 11)      channel-centre energies in keV (use row 0)
  DiffElectronEffectiveEnergy (5, 10)    channel-centre energies in keV (use row 0)

Output columns after averaging over telescope axis (axis=1):
  proton_ch{i}_{E:.0f}keV    e.g. proton_ch0_80keV
  electron_ch{i}_{E:.0f}keV
  int_electron_flux

Deletes raw NC files after verified write per satellite.
Output: data/processed/solar/
"""
import logging
from datetime import datetime, timedelta, timezone
from pathlib import Path

import netCDF4 as nc4
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

# GOES-R MPSH time reference epoch
MPSH_EPOCH = datetime(2000, 1, 1, 12, 0, 0, tzinfo=timezone.utc)
FILL_THRESHOLD = 1e30


# ---------------------------------------------------------------------------
# Legacy MAGED (goes-13/14/15) — gzipped monthly CSV files
#
# Confirmed path structure:
#   data/solar/goes_particle/{GOESNN_EPS-MAGED_1MIN}/{YYYY}/{GOESNN_EPS-MAGED_1MIN_YYYY-MM.csv.gz}
#
# Format: 9 columns, no header row
#   col 0 : ISO datetime string  → parse with pd.to_datetime(utc=True)
#   col 1-8: E1-E8 electron differential flux channels (~40-475 keV)
#   Fill value: -1.00e+05 → NaN
# ---------------------------------------------------------------------------
MAGED_COLS  = ["time", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8"]
MAGED_FILL  = -1.0e+05
MAGED_FILL_THRESH = -9.0e+04   # anything below this is a fill


def process_goes_legacy_particle() -> None:
    out = SOL_OUT / "goes_legacy_particle.parquet"
    if out.exists():
        LOG.info("SKIP goes_legacy_particle.parquet")
        return

    src_dir = DATA_ROOT / "solar" / "goes_particle"
    if not src_dir.exists():
        LOG.warning("Missing directory: %s", src_dir)
        return

    # Confirmed structure: {SAT_EPS-MAGED_1MIN}/{YYYY}/{monthly.csv.gz}
    gz_files = sorted(src_dir.rglob("*.csv.gz"))
    if not gz_files:
        LOG.warning("No CSV.gz MAGED files found under %s", src_dir)
        return

    LOG.info("Legacy MAGED: found %d CSV.gz files under %s", len(gz_files), src_dir)

    all_frames: list[pd.DataFrame] = []
    raw_to_delete: list[Path] = []

    for fp in gz_files:
        # Satellite name from grandparent dir: "GOES13_EPS-MAGED_1MIN" → "goes13"
        sat_dir_name = fp.parent.parent.name
        sat_id = sat_dir_name.split("_")[0].lower()   # e.g. "goes13"

        try:
            df = pd.read_csv(
                fp,
                compression="gzip",
                header=None,
                names=MAGED_COLS,
            )
            df["time"] = pd.to_datetime(df["time"], utc=True, errors="coerce")
            df = df.dropna(subset=["time"]).set_index("time").sort_index()

            # Replace fill values: exact match + threshold guard for float drift
            flux_cols = [c for c in MAGED_COLS if c != "time"]
            df[flux_cols] = df[flux_cols].replace(MAGED_FILL, np.nan)
            for col in flux_cols:
                df[col] = df[col].where(df[col] > MAGED_FILL_THRESH, other=np.nan)

            df["satellite_id"] = sat_id
            all_frames.append(df)
            raw_to_delete.append(fp)

        except Exception as exc:
            LOG.warning("Error reading %s: %s", fp.name, exc)

    if not all_frames:
        LOG.warning("No MAGED data loaded")
        return

    combined = pd.concat(all_frames).sort_index()
    combined = combined[~combined.index.duplicated(keep="first")]

    meta = {
        "title": "GOES-13/14/15 EPS-MAGED 1-min Electron Differential Flux (E1-E8)",
        "source": "NOAA/NCEI GOES Particle",
        "references": "https://www.ngdc.noaa.gov/stp/satellite/goes/",
        "time_range": f"{combined.index.min()} / {combined.index.max()}",
        "units": "electrons cm-2 sr-1 s-1 keV-1; channels E1-E8 ~40-475 keV",
    }
    save_parquet(combined, out, meta)

    if out.exists() and out.stat().st_size > 0:
        safe_delete(raw_to_delete)
        register_output(MANIFEST, "goes_legacy_particle", out, True, meta)
    else:
        LOG.error("Output missing/empty — raw files NOT deleted")


# ---------------------------------------------------------------------------
# GOES-R MPSH — confirmed variable structure
# ---------------------------------------------------------------------------
def _mpsh_times(time_arr: np.ndarray) -> pd.DatetimeIndex:
    """Convert MPSH time (seconds since 2000-01-01 12:00:00 UTC) to DatetimeIndex."""
    timestamps = [MPSH_EPOCH + timedelta(seconds=float(t)) for t in time_arr]
    return pd.DatetimeIndex(timestamps, tz="UTC")


def _process_one_goes_r_sat(sat_name: str, sat_dir: Path) -> None:
    out = SOL_OUT / f"goes_r_particle_{sat_name}.parquet"
    if out.exists():
        LOG.info("SKIP %s", out.name)
        return

    # File layout: mpsh_avg1m/{year}/{month}/{daily.nc}  — use rglob
    mpsh_root = sat_dir / "mpsh_avg1m"
    if not mpsh_root.exists():
        mpsh_root = sat_dir   # fallback: search sat_dir directly

    nc_files = sorted(mpsh_root.rglob("*.nc"))
    if not nc_files:
        LOG.warning("No MPSH NC files for %s under %s", sat_name, mpsh_root)
        return

    LOG.info("Processing GOES-R %s: %d files …", sat_name, len(nc_files))

    frames: list[pd.DataFrame] = []
    raw_to_delete: list[Path] = []

    # Read energy channel centres from the first readable file
    proton_energies:   np.ndarray | None = None
    electron_energies: np.ndarray | None = None

    for fp in nc_files:
        try:
            ds = nc4.Dataset(str(fp))

            # --- Time ---
            time_raw = ds.variables["time"][:]          # (1440,) seconds
            dt_index = _mpsh_times(time_raw)

            # --- Energy channel centres (read once, same for all files) ---
            if proton_energies is None and "DiffProtonEffectiveEnergy" in ds.variables:
                # shape (5, 11) — use telescope 0 row
                proton_energies = np.asarray(
                    ds.variables["DiffProtonEffectiveEnergy"][0, :], dtype=float
                )
            if electron_energies is None and "DiffElectronEffectiveEnergy" in ds.variables:
                electron_energies = np.asarray(
                    ds.variables["DiffElectronEffectiveEnergy"][0, :], dtype=float
                )

            # --- Differential proton flux (1440, 5, 11) → mean over telescopes → (1440, 11) ---
            pf_raw = np.ma.filled(
                ds.variables["AvgDiffProtonFlux"][:].astype(float), np.nan
            )
            pf_raw[pf_raw > FILL_THRESHOLD] = np.nan
            pf_mean = np.nanmean(pf_raw, axis=1)        # (ntime, 11)

            # --- Differential electron flux (1440, 5, 10) → (1440, 10) ---
            ef_raw = np.ma.filled(
                ds.variables["AvgDiffElectronFlux"][:].astype(float), np.nan
            )
            ef_raw[ef_raw > FILL_THRESHOLD] = np.nan
            ef_mean = np.nanmean(ef_raw, axis=1)        # (ntime, 10)

            # --- Integral electron flux (1440, 5) → (1440,) ---
            ie_raw = np.ma.filled(
                ds.variables["AvgIntElectronFlux"][:].astype(float), np.nan
            )
            ie_raw[ie_raw > FILL_THRESHOLD] = np.nan
            ie_mean = np.nanmean(ie_raw, axis=1)        # (ntime,)

            ds.close()

            # --- Build column names using energy values ---
            p_labels = (
                [f"proton_ch{i}_{int(round(proton_energies[i]))}keV"
                 for i in range(pf_mean.shape[1])]
                if proton_energies is not None
                else [f"proton_ch{i}_keV" for i in range(pf_mean.shape[1])]
            )
            e_labels = (
                [f"electron_ch{i}_{int(round(electron_energies[i]))}keV"
                 for i in range(ef_mean.shape[1])]
                if electron_energies is not None
                else [f"electron_ch{i}_keV" for i in range(ef_mean.shape[1])]
            )

            data: dict[str, np.ndarray] = {"int_electron_flux": ie_mean}
            for i, lbl in enumerate(p_labels):
                data[lbl] = pf_mean[:, i]
            for i, lbl in enumerate(e_labels):
                data[lbl] = ef_mean[:, i]

            df = pd.DataFrame(data, index=dt_index)
            frames.append(df)
            raw_to_delete.append(fp)

        except Exception as exc:
            LOG.warning("Error reading %s: %s", fp.name, exc)

    if not frames:
        LOG.warning("No MPSH data loaded for %s", sat_name)
        return

    combined = pd.concat(frames).sort_index()
    combined = combined[~combined.index.duplicated(keep="first")]

    meta = {
        "title": f"GOES-R MPSH Energetic Particle Data — {sat_name}",
        "source": "NOAA/NCEI GOES-R Series",
        "references": "https://doi.org/10.7289/V5BV7DSR",
        "time_range": f"{combined.index.min()} / {combined.index.max()}",
        "units": "proton/electron differential flux protons|electrons/(cm² sr keV s)",
    }
    save_parquet(combined, out, meta)

    if out.exists() and out.stat().st_size > 0:
        safe_delete(raw_to_delete)
        register_output(MANIFEST, f"goes_r_particle_{sat_name}", out, True, meta)
    else:
        LOG.error("Output %s missing/empty — raw NOT deleted", out.name)


def process_goes_r_particle() -> None:
    out_merged = SOL_OUT / "goes_r_particle.parquet"
    if out_merged.exists():
        LOG.info("SKIP goes_r_particle.parquet")
        return

    base_dir = DATA_ROOT / "solar" / "goes_r_particle"
    if not base_dir.exists():
        LOG.warning("Missing directory: %s", base_dir)
        return

    sat_dirs = sorted(d for d in base_dir.iterdir() if d.is_dir())
    per_sat: list[Path] = []

    for sat_dir in sat_dirs:
        sat_name = sat_dir.name   # e.g. "goes16"
        _process_one_goes_r_sat(sat_name, sat_dir)
        p = SOL_OUT / f"goes_r_particle_{sat_name}.parquet"
        if p.exists() and p.stat().st_size > 0:
            per_sat.append(p)
        LOG.info("After %s | disk free=%.1f GB", sat_name, disk_free_gb())

    if not per_sat:
        LOG.warning("No per-satellite GOES-R parquets to merge")
        return

    frames = [pd.read_parquet(p) for p in per_sat]
    merged = pd.concat(frames).sort_index()

    meta = {
        "title": "GOES-R MPSH Energetic Particle Data (all satellites merged)",
        "source": "NOAA/NCEI GOES-R Series",
        "references": "https://doi.org/10.7289/V5BV7DSR",
        "time_range": f"{merged.index.min()} / {merged.index.max()}",
        "units": "proton/electron differential flux protons|electrons/(cm² sr keV s)",
    }
    save_parquet(merged, out_merged, meta)
    register_output(MANIFEST, "goes_r_particle_merged", out_merged, True, meta)


def main() -> None:
    setup_logging()
    LOG.info("=== 05_process_goes_particle.py | disk free=%.1f GB ===", disk_free_gb())
    SOL_OUT.mkdir(parents=True, exist_ok=True)

    process_goes_legacy_particle()
    LOG.info("After legacy particle | disk free=%.1f GB", disk_free_gb())

    process_goes_r_particle()
    LOG.info("=== 05 complete | disk free=%.1f GB ===", disk_free_gb())


if __name__ == "__main__":
    main()
