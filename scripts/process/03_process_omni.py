"""
03_process_omni.py
Process OMNI solar wind fixed-width text files → Parquet.
Deletes raw .dat/.asc files after verified write.
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

# OMNI2 fill values — integers stored as float after pd.read_csv(dtype=float).
# Integer sentinel columns (Kp→99, Dst→99999, AE→9999, etc.) and float columns.
# Use pd.Series.replace for exact equality; no tolerance needed for text-parsed floats.
_FILL_INTS  = [99.0, 999.0, 9999.0, 99999.0, 9999999.0]
_FILL_FLOATS = [99.9, 999.9, 9999.9, 99999.9, 9999999.0]
FILL_ALL = list(dict.fromkeys(_FILL_INTS + _FILL_FLOATS))  # unique, ordered

# OMNIweb hourly column names (55 total per omni2 format documentation)
OMNI_HOURLY_COLS = [
    "year", "doy", "hour", "bartels_rot",
    "id_IMF", "id_plasma",
    "n_IMF_pts", "n_plasma_pts", "pct_fine_scale",
    "Bx_GSE", "By_GSE", "Bz_GSE", "By_GSM", "Bz_GSM",
    "sigma_B_mag", "sigma_B_vec",
    "flow_speed", "Vx_GSE", "Vy_GSE", "Vz_GSE",
    "proton_density", "T_proton",
    "flow_pressure", "E_field", "plasma_beta", "Alfven_Mach",
    "Kp", "R_sunspot", "Dst", "AE_index",
    "proton_flux_1MeV", "proton_flux_2MeV", "proton_flux_4MeV",
    "proton_flux_10MeV", "proton_flux_30MeV", "proton_flux_60MeV",
    "flag_AE", "AL", "AU",
    "mac_number", "lyman_alpha",
    "proton_flux_gt1MeV", "proton_flux_gt2MeV", "proton_flux_gt4MeV",
    "proton_flux_gt10MeV", "proton_flux_gt30MeV", "proton_flux_gt60MeV",
    "magnetosonic_mach",
    # Padding columns to reach 55 if file has them
    "col50", "col51", "col52", "col53", "col54",
]

# OMNI 1-min column names (26 total)
OMNI_1MIN_COLS = [
    "year", "doy", "hour", "minute",
    "Bx", "By", "Bz",
    "Vx", "Vy", "Vz",
    "proton_density", "T",
    "flow_speed", "phi_v", "theta_v",
    "B_magnitude", "RMS_B", "RMS_Bsc",
    "sigma_Bx", "sigma_By", "sigma_Bz",
    "sigma_Vx", "sigma_Vy", "sigma_Vz",
    "sigma_n", "sigma_T",
]


def _fill_to_nan(df: pd.DataFrame) -> pd.DataFrame:
    """Replace OMNI fill sentinel values with NaN across all numeric columns.

    pd.Series.replace uses exact equality, which is reliable for values
    parsed directly from fixed-width text (e.g. "9999.9" → exactly 9999.9).
    """
    num_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    df[num_cols] = df[num_cols].replace(FILL_ALL, np.nan)
    return df


def _build_hourly_datetime(df: pd.DataFrame) -> pd.Series:
    """Construct UTC DatetimeIndex from year/doy/hour columns."""
    year = df["year"].fillna(0).astype(int)
    doy  = df["doy"].fillna(0).astype(int)
    hour = df["hour"].fillna(0).astype(int)
    base = pd.to_datetime(year.astype(str), format="%Y", utc=True, errors="coerce")
    return base + pd.to_timedelta(doy - 1, unit="D") + pd.to_timedelta(hour, unit="h")


def _build_1min_datetime(df: pd.DataFrame) -> pd.Series:
    """Construct UTC DatetimeIndex from year/doy/hour/minute columns."""
    year   = df["year"].fillna(0).astype(int)
    doy    = df["doy"].fillna(0).astype(int)
    hour   = df["hour"].fillna(0).astype(int)
    minute = df["minute"].fillna(0).astype(int)
    base = pd.to_datetime(year.astype(str), format="%Y", utc=True, errors="coerce")
    return (base
            + pd.to_timedelta(doy - 1, unit="D")
            + pd.to_timedelta(hour, unit="h")
            + pd.to_timedelta(minute, unit="min"))


# ---------------------------------------------------------------------------
# OMNI hourly
# ---------------------------------------------------------------------------
def process_omni_hourly() -> None:
    out_combined = SOL_OUT / "omni_hourly.parquet"
    if out_combined.exists():
        LOG.info("SKIP omni_hourly.parquet")
        return

    src_dir = DATA_ROOT / "atmospheric" / "omni_solar_wind" / "low_res"
    if not src_dir.exists():
        LOG.warning("Missing directory: %s", src_dir)
        return

    dat_files = sorted(src_dir.glob("omni2_*.dat"))
    if not dat_files:
        LOG.warning("No omni2_*.dat files found in %s", src_dir)
        return

    yearly_parquets: list[Path] = []
    raw_to_delete: list[Path] = []

    for fp in dat_files:
        year_str = fp.stem.replace("omni2_", "")
        out_year = SOL_OUT / f"omni_hourly_{year_str}.parquet"

        if not out_year.exists():
            try:
                df = pd.read_csv(
                    fp,
                    sep=r"\s+",
                    header=None,
                    names=OMNI_HOURLY_COLS[:55],
                    usecols=range(min(55, len(OMNI_HOURLY_COLS))),
                    dtype=float,
                    on_bad_lines="skip",
                )
                df = _fill_to_nan(df)
                df["time"] = _build_hourly_datetime(df)
                df = df.dropna(subset=["time"]).set_index("time").sort_index()
                drop_cols = [c for c in ("year", "doy", "hour") if c in df.columns]
                df = df.drop(columns=drop_cols, errors="ignore")
                drop_pad = [c for c in df.columns if c.startswith("col5")]
                df = df.drop(columns=drop_pad, errors="ignore")

                meta = {
                    "title": f"OMNI2 Hourly Solar Wind {year_str}",
                    "source": "NASA/GSFC OMNIweb",
                    "references": "https://omniweb.gsfc.nasa.gov/",
                    "time_range": f"{df.index.min()} / {df.index.max()}",
                    "units": "B nT; V km/s; n cm-3; T K; P nPa",
                }
                save_parquet(df, out_year, meta)
            except Exception as exc:
                LOG.warning("Failed to process %s: %s", fp.name, exc)
                continue

        if out_year.exists() and out_year.stat().st_size > 0:
            yearly_parquets.append(out_year)
            raw_to_delete.append(fp)
        else:
            LOG.warning("Output %s is zero-size, skipping delete", out_year.name)

    if not yearly_parquets:
        LOG.warning("No yearly OMNI parquets produced")
        return

    # Concatenate all years
    LOG.info("Concatenating %d yearly OMNI files …", len(yearly_parquets))
    frames = []
    for p in yearly_parquets:
        try:
            frames.append(pd.read_parquet(p))
        except Exception as exc:
            LOG.warning("Could not read %s: %s", p.name, exc)

    if not frames:
        return

    combined = pd.concat(frames).sort_index()
    combined = combined[~combined.index.duplicated(keep="first")]

    meta = {
        "title": "OMNI2 Hourly Solar Wind (all years combined)",
        "source": "NASA/GSFC OMNIweb",
        "references": "https://omniweb.gsfc.nasa.gov/",
        "time_range": f"{combined.index.min()} / {combined.index.max()}",
        "units": "B nT; V km/s; n cm-3; T K; P nPa",
    }
    save_parquet(combined, out_combined, meta)

    if out_combined.exists() and out_combined.stat().st_size > 0:
        safe_delete(raw_to_delete)
        register_output(MANIFEST, "omni_hourly", out_combined, True, meta)
    else:
        LOG.error("Combined OMNI parquet is missing/empty — raw files NOT deleted")


# ---------------------------------------------------------------------------
# OMNI 1-min
# ---------------------------------------------------------------------------
def process_omni_1min() -> None:
    out = SOL_OUT / "omni_1min.parquet"
    if out.exists():
        LOG.info("SKIP omni_1min.parquet")
        return

    src_dir = DATA_ROOT / "atmospheric" / "omni_solar_wind" / "high_res_1min"
    if not src_dir.exists():
        LOG.warning("Missing directory: %s", src_dir)
        return

    asc_files = sorted(src_dir.glob("omni_1min_*.asc"))
    if not asc_files:
        LOG.warning("No omni_1min_*.asc files found in %s", src_dir)
        return

    frames = []
    raw_to_delete: list[Path] = []
    success_files: list[Path] = []

    for fp in asc_files:
        try:
            df = pd.read_csv(
                fp,
                sep=r"\s+",
                header=None,
                names=OMNI_1MIN_COLS,
                usecols=range(len(OMNI_1MIN_COLS)),
                dtype=float,
                on_bad_lines="skip",
            )
            df = _fill_to_nan(df)
            df["time"] = _build_1min_datetime(df)
            df = df.dropna(subset=["time"]).set_index("time").sort_index()
            drop_cols = [c for c in ("year", "doy", "hour", "minute") if c in df.columns]
            df = df.drop(columns=drop_cols, errors="ignore")
            frames.append(df)
            success_files.append(fp)
        except Exception as exc:
            LOG.warning("Failed to process %s: %s", fp.name, exc)

    if not frames:
        LOG.warning("No 1-min OMNI data loaded")
        return

    combined = pd.concat(frames).sort_index()
    combined = combined[~combined.index.duplicated(keep="first")]

    meta = {
        "title": "OMNI 1-Minute Solar Wind",
        "source": "NASA/GSFC OMNIweb",
        "references": "https://omniweb.gsfc.nasa.gov/",
        "time_range": f"{combined.index.min()} / {combined.index.max()}",
        "units": "B nT; V km/s; n cm-3; T K",
    }
    save_parquet(combined, out, meta)

    if out.exists() and out.stat().st_size > 0:
        safe_delete(success_files)
        register_output(MANIFEST, "omni_1min", out, True, meta)
    else:
        LOG.error("1-min OMNI parquet is missing/empty — raw files NOT deleted")


def main() -> None:
    setup_logging()
    LOG.info("=== 03_process_omni.py | disk free=%.1f GB ===", disk_free_gb())
    SOL_OUT.mkdir(parents=True, exist_ok=True)

    process_omni_hourly()
    LOG.info("After hourly OMNI | disk free=%.1f GB", disk_free_gb())

    process_omni_1min()
    LOG.info("=== 03 complete | disk free=%.1f GB ===", disk_free_gb())


if __name__ == "__main__":
    main()
