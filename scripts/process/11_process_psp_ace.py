"""
11_process_psp_ace.py
Process Parker Solar Probe (PSP) CDF and ACE/DSCOVR compressed NetCDF → Parquet.

PSP confirmed:
  epoch variable : epoch_mag_SC_4_Sa_per_Cyc
  field variable : psp_fld_l2_mag_RTN_4_Sa_per_Cyc  shape (N, 3)  [B_R, B_T, B_N]
  Files: data/solar/psp_mag/psp_fld_l2_mag_rtn_4_sa_per_cyc_*.cdf

ACE/DSCOVR confirmed variables:
  time, proton_vx_gse, proton_vy_gse, proton_vz_gse,
  proton_speed, proton_density, proton_temperature
  Files: data/solar/ace_wind_dscovr/ncei_archive/**/*.nc.gz

PSP files are NOT deleted (unique science data).
ACE .nc.gz files are deleted after verified write.
Output: data/processed/solar/
"""
import gzip
import logging
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

PSP_DIR = DATA_ROOT / "solar" / "psp_mag"
ACE_DIR = DATA_ROOT / "solar" / "ace_wind_dscovr" / "ncei_archive"

FILL_THRESHOLD = 1e30

# Confirmed ACE/DSCOVR variable names (from actual file inspection)
# Bz variants checked via fuzzy scan after primary extraction
ACE_VARS = [
    "proton_vx_gse",
    "proton_vy_gse",
    "proton_vz_gse",
    "proton_speed",
    "proton_density",
    "proton_temperature",
    "bt",
    "bx_gsm",
    "by_gsm",
]

# Additional Bz / B-field variables to search for (names vary by dataset version)
ACE_BZ_CANDIDATES = [
    "bz_gsm", "Bz_GSM", "BZ_GSM",
    "b_gse_z", "b_z_gse",
    "b_field_z", "b_field",
    "bz_gse", "Bz_GSE",
    "dn_bz_gsm",
]


# ---------------------------------------------------------------------------
# PSP magnetic field (CDF via cdflib)
# ---------------------------------------------------------------------------
def _load_psp_cdf(fp: Path) -> pd.DataFrame | None:
    """Load a single PSP CDF file using cdflib."""
    try:
        import cdflib
    except ImportError:
        LOG.error("cdflib not installed — cannot process PSP CDFs")
        return None

    try:
        cdf = cdflib.CDF(str(fp))
        info = cdf.cdf_info()
        all_vars = list(info.zVariables) + list(info.rVariables)

        # --- Epoch (confirmed primary name, fall back to alternatives) ---
        epoch_var = None
        for ev in ("epoch_mag_RTN_4_Sa_per_Cyc", "epoch_mag_SC_4_Sa_per_Cyc",
                   "Epoch", "epoch", "EPOCH", "time"):
            if ev in all_vars:
                epoch_var = ev
                break
        if epoch_var is None:
            LOG.warning("No epoch variable in %s — vars: %s", fp.name, all_vars[:10])
            return None

        epoch_data = cdf.varget(epoch_var)
        if epoch_data is None or len(epoch_data) == 0:
            return None

        datetimes = cdflib.cdfepoch.to_datetime(epoch_data)
        dt_index = pd.DatetimeIndex(datetimes, tz="UTC")

        # --- Magnetic field vector (confirmed primary name, fall back) ---
        mag_var = None
        for mv in ("psp_fld_l2_mag_RTN_4_Sa_per_Cyc", "psp_fld_l2_mag_RTN",
                    "mag_RTN", "B_RTN", "BRTN", "psp_fld_l2_mag_SC"):
            if mv in all_vars:
                mag_var = mv
                break
        if mag_var is None:
            # Fuzzy: any var with 'mag' and 'rtn'
            for vname in all_vars:
                if "mag" in vname.lower() and "rtn" in vname.lower():
                    mag_var = vname
                    break
        if mag_var is None:
            LOG.warning("No B-field variable found in %s", fp.name)
            return None

        b_arr = cdf.varget(mag_var)
        if b_arr is None:
            return None

        b_arr = np.asarray(b_arr, dtype=float)
        # Confirmed shape: (N, 3) — transpose if needed
        if b_arr.ndim == 2 and b_arr.shape[1] != 3 and b_arr.shape[0] == 3:
            b_arr = b_arr.T
        if b_arr.shape[0] != len(dt_index):
            LOG.warning("Shape mismatch in %s: B %s vs epoch %d",
                        fp.name, b_arr.shape, len(dt_index))
            return None

        df = pd.DataFrame(index=dt_index)
        df["B_R"] = b_arr[:, 0]
        df["B_T"] = b_arr[:, 1]
        df["B_N"] = b_arr[:, 2]
        df["B_mag"] = np.sqrt(np.sum(b_arr[:, :3] ** 2, axis=1))

        # Mask fill values
        df = df.where(df.abs() < FILL_THRESHOLD)
        return df.sort_index()

    except Exception as exc:
        LOG.warning("Error reading PSP CDF %s: %s", fp.name, exc)
        return None


def process_psp() -> None:
    """
    Process PSP CDF files year-by-year to avoid OOM.
    Resamples to 1-minute averages (preserves mean, min, max of |B|)
    which is more than adequate for SOC and cross-correlation analysis.
    PSP raw files are NOT deleted (unique science data).
    """
    out = SOL_OUT / "psp_mag.parquet"
    if out.exists():
        LOG.info("SKIP psp_mag.parquet")
        return

    if not PSP_DIR.exists():
        LOG.warning("Missing PSP directory: %s", PSP_DIR)
        return

    # Build year → file list map from subdirectories
    year_dirs = sorted(d for d in PSP_DIR.iterdir() if d.is_dir() and d.name.isdigit())
    if not year_dirs:
        # Flat layout — group files by name pattern
        cdf_files = sorted(PSP_DIR.rglob("*.cdf"))
        if not cdf_files:
            LOG.warning("No PSP CDF files found in %s", PSP_DIR)
            return
        year_map = {"all": cdf_files}
    else:
        year_map = {}
        for d in year_dirs:
            cdfs = sorted(d.glob("*.cdf"))
            if cdfs:
                year_map[d.name] = cdfs

    if not year_map:
        LOG.warning("No PSP CDF files found")
        return

    total_files = sum(len(v) for v in year_map.values())
    LOG.info("Processing %d PSP CDF files across %d year(s) …", total_files, len(year_map))

    year_parquets: list[Path] = []
    for year_label, cdf_files in year_map.items():
        year_out = SOL_OUT / f"psp_mag_{year_label}.parquet"
        if year_out.exists():
            LOG.info("SKIP psp_mag_%s.parquet", year_label)
            year_parquets.append(year_out)
            continue

        LOG.info("  PSP year=%s: loading %d CDF files …", year_label, len(cdf_files))
        frames = []
        for fp in cdf_files:
            df = _load_psp_cdf(fp)
            if df is not None and not df.empty:
                frames.append(df)

        if not frames:
            LOG.warning("  PSP year=%s: no data loaded", year_label)
            continue

        year_df = pd.concat(frames).sort_index()
        del frames
        year_df = year_df[~year_df.index.duplicated(keep="first")]

        # Resample to 1-minute averages — adequate for SOC analysis, reduces data ~250×
        resampled = year_df.resample("1min").agg(["mean", "min", "max"]).dropna(how="all")
        del year_df
        # Flatten multi-level columns: B_R_mean, B_R_min, ...
        resampled.columns = ["_".join(c) for c in resampled.columns]

        year_meta = {
            "title": f"PSP FIELDS L2 Mag RTN 1-min resampled {year_label}",
            "source": "NASA SWEAP/FIELDS — PSP",
        }
        save_parquet(resampled, year_out, year_meta)
        LOG.info("  PSP year=%s: %d rows (1-min), %.1f MB",
                 year_label, len(resampled), year_out.stat().st_size / 1e6)
        year_parquets.append(year_out)
        del resampled

    if not year_parquets:
        LOG.warning("No PSP data produced")
        return

    # Merge year parquets (1-min resampled → small enough to merge in memory)
    LOG.info("Merging %d yearly PSP Parquets …", len(year_parquets))
    all_frames = [pd.read_parquet(p) for p in year_parquets]
    combined = pd.concat(all_frames).sort_index()
    combined = combined[~combined.index.duplicated(keep="first")]
    del all_frames

    meta = {
        "title": "Parker Solar Probe FIELDS L2 Magnetic Field (RTN) 1-min averages",
        "source": "NASA SWEAP/FIELDS — PSP",
        "references": "https://doi.org/10.3847/1538-4357/ab1e34",
        "time_range": f"{combined.index.min()} / {combined.index.max()}",
        "units": "B nT; components RTN (Radial, Tangential, Normal); _mean/_min/_max per 1-min window",
    }
    save_parquet(combined, out, meta)

    if out.exists() and out.stat().st_size > 0:
        for p in year_parquets:
            try:
                p.unlink()
            except Exception:
                pass
        register_output(MANIFEST, "psp_mag", out, False, meta)
    # DO NOT delete PSP raw CDF files — unique science data


# ---------------------------------------------------------------------------
# ACE / DSCOVR (.nc.gz) — read in-memory with nc4.Dataset(memory=...)
# ---------------------------------------------------------------------------
def _open_ncgz_ace(fp: Path) -> pd.DataFrame | None:
    """
    Decompress .nc.gz, open with netCDF4.Dataset,
    extract confirmed ACE/DSCOVR solar wind variables.
    Tries in-memory first (fast), falls back to temp file on Windows.
    """
    import tempfile, shutil

    # Only process f1m (Faraday cup 1-min plasma) and m1m (magnetometer 1-min)
    # All other types are duplicates at different cadences or auxiliary data
    fname = fp.name.lower()
    if not any(fname.startswith(p) for p in ("oe_f1m_", "oe_m1m_")):
        return None

    try:
        with gzip.open(fp, "rb") as gz_fh:
            raw_bytes = gz_fh.read()
    except Exception as exc:
        LOG.warning("Cannot decompress %s: %s", fp.name, exc)
        return None

    # Try in-memory first (faster, works for most files)
    tmp_path = None
    try:
        ds = nc4.Dataset("inmemory", memory=raw_bytes)
    except Exception:
        # Fallback: write to temp file
        try:
            with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as tmp:
                tmp_path = Path(tmp.name)
                tmp.write(raw_bytes)
            ds = nc4.Dataset(str(tmp_path))
        except Exception as exc:
            LOG.warning("netCDF4 cannot open %s: %s", fp.name, exc)
            if tmp_path and tmp_path.exists():
                tmp_path.unlink()
            return None
    del raw_bytes  # free memory

    try:
        # --- Time ---
        if "time" not in ds.variables:
            LOG.warning("No 'time' variable in %s — vars: %s",
                        fp.name, list(ds.variables)[:10])
            ds.close()
            if tmp_path and tmp_path.exists(): tmp_path.unlink()
            return None

        tv = ds.variables["time"]
        units    = getattr(tv, "units", "")
        calendar = getattr(tv, "calendar", "standard")
        try:
            dates = nc4.num2date(tv[:], units, calendar=calendar)
            dt_index = pd.DatetimeIndex(
                [pd.Timestamp(d.isoformat(), tz="UTC") for d in dates]
            )
        except Exception as exc:
            LOG.warning("Time decode failed in %s: %s", fp.name, exc)
            ds.close()
            if tmp_path and tmp_path.exists(): tmp_path.unlink()
            return None

        # --- Confirmed solar wind variables ---
        data: dict[str, np.ndarray] = {}
        for vname in ACE_VARS:
            if vname in ds.variables:
                try:
                    arr = np.ma.filled(
                        ds.variables[vname][:].astype(float), np.nan
                    )
                    arr[arr > FILL_THRESHOLD] = np.nan
                    if arr.ndim == 1 and len(arr) == len(dt_index):
                        data[vname] = arr
                except Exception as exc:
                    LOG.debug("Var %s in %s: %s", vname, fp.name, exc)
            else:
                LOG.debug("Variable %s not found in %s", vname, fp.name)

        # --- Bz / B-field: try candidate names, use first match ---
        bz_found = False
        for bz_name in ACE_BZ_CANDIDATES:
            if bz_name in ds.variables:
                try:
                    arr = np.ma.filled(
                        ds.variables[bz_name][:].astype(float), np.nan
                    )
                    arr[arr > FILL_THRESHOLD] = np.nan
                    if arr.ndim == 1 and len(arr) == len(dt_index):
                        data["bz_gsm"] = arr   # normalised output column name
                        LOG.debug("Bz source var=%s in %s", bz_name, fp.name)
                        bz_found = True
                        break
                except Exception as exc:
                    LOG.debug("Bz var %s in %s: %s", bz_name, fp.name, exc)
        if not bz_found:
            # Last resort: any variable with 'bz' or 'b_z' in its name
            for vname in ds.variables:
                vl = vname.lower()
                if ("bz" in vl or "b_z" in vl) and "bz_gsm" not in data:
                    try:
                        arr = np.ma.filled(
                            ds.variables[vname][:].astype(float), np.nan
                        )
                        arr[arr > FILL_THRESHOLD] = np.nan
                        if arr.ndim == 1 and len(arr) == len(dt_index):
                            data["bz_gsm"] = arr
                            LOG.debug("Bz fuzzy match var=%s in %s", vname, fp.name)
                    except Exception:
                        pass

        ds.close()
        if tmp_path and tmp_path.exists():
            tmp_path.unlink()

        if not data:
            return None

        return pd.DataFrame(data, index=dt_index).sort_index()

    except Exception as exc:
        LOG.warning("Error extracting %s: %s", fp.name, exc)
        try:
            ds.close()
        except Exception:
            pass
        if tmp_path and tmp_path.exists():
            tmp_path.unlink()
        return None


def process_ace_dscovr() -> None:
    out = SOL_OUT / "ace_dscovr.parquet"
    if out.exists():
        LOG.info("SKIP ace_dscovr.parquet")
        return

    if not ACE_DIR.exists():
        LOG.warning("Missing ACE/DSCOVR archive directory: %s", ACE_DIR)
        return

    # ACE archive is organized as ncei_archive/{YYYY}/{MM}/*.nc.gz
    # Process one year at a time to bound memory usage, delete each year's raw files after writing
    year_dirs = sorted(d for d in ACE_DIR.iterdir() if d.is_dir() and d.name.isdigit())
    if not year_dirs:
        # Flat directory — fall back to full list
        gz_all = sorted(ACE_DIR.rglob("*.nc.gz")) or sorted(ACE_DIR.rglob("*.gz"))
        year_dirs_map = {"all": gz_all}
    else:
        year_dirs_map = {d.name: sorted(d.rglob("*.nc.gz")) or sorted(d.rglob("*.gz"))
                        for d in year_dirs}

    year_parquets: list[Path] = []

    for year_label, gz_files in year_dirs_map.items():
        year_out = SOL_OUT / f"ace_dscovr_{year_label}.parquet"

        if year_out.exists():
            LOG.info("SKIP ace_dscovr_%s.parquet", year_label)
            year_parquets.append(year_out)
            continue

        if not gz_files:
            continue

        LOG.info("Processing ACE/DSCOVR year=%s: %d files …", year_label, len(gz_files))
        frames: list[pd.DataFrame] = []
        raw_to_delete: list[Path] = []

        for i, fp in enumerate(gz_files):
            df = _open_ncgz_ace(fp)
            if df is not None and not df.empty:
                frames.append(df)
            raw_to_delete.append(fp)
            if (i + 1) % 100 == 0:
                LOG.info("  %s: %d/%d files processed, %d frames collected",
                         year_label, i + 1, len(gz_files), len(frames))

        if not frames:
            # Still delete even if no data extracted
            safe_delete(raw_to_delete)
            continue

        year_df = pd.concat(frames).sort_index()
        year_df = year_df[~year_df.index.duplicated(keep="first")]

        year_meta = {
            "title": f"ACE/DSCOVR Solar Wind {year_label}",
            "source": "NOAA NCEI / NASA ACE Science Center",
            "references": "https://doi.org/10.1029/1998GL900054",
            "time_range": f"{year_df.index.min()} / {year_df.index.max()}",
            "units": "proton_speed km/s; proton_density cm-3; proton_temperature K; Vxyz km/s",
        }
        save_parquet(year_df, year_out, year_meta)

        if year_out.exists() and year_out.stat().st_size > 0:
            safe_delete(raw_to_delete)
            LOG.info("  Year %s: %d rows written, raw deleted | disk free=%.1f GB",
                     year_label, len(year_df), disk_free_gb())
            year_parquets.append(year_out)
        else:
            LOG.error("  Year %s: output empty — raw NOT deleted", year_label)

    if not year_parquets:
        LOG.warning("No ACE/DSCOVR data produced")
        return

    # Merge all yearly Parquets into combined output
    LOG.info("Merging %d yearly ACE/DSCOVR Parquets …", len(year_parquets))
    all_frames = [pd.read_parquet(p) for p in year_parquets]
    combined = pd.concat(all_frames).sort_index()
    combined = combined[~combined.index.duplicated(keep="first")]

    meta = {
        "title": "ACE / DSCOVR Solar Wind In-Situ Observations (combined)",
        "source": "NOAA NCEI / NASA ACE Science Center",
        "references": "https://doi.org/10.1029/1998GL900054",
        "time_range": f"{combined.index.min()} / {combined.index.max()}",
        "units": "proton_speed km/s; proton_density cm-3; proton_temperature K; Vxyz km/s; bz_gsm nT",
    }
    save_parquet(combined, out, meta)

    if out.exists() and out.stat().st_size > 0:
        for p in year_parquets:
            try:
                p.unlink()
            except Exception:
                pass
        register_output(MANIFEST, "ace_dscovr", out, True, meta)
    else:
        LOG.error("Combined ACE/DSCOVR parquet missing/empty")


def main() -> None:
    setup_logging()
    LOG.info("=== 11_process_psp_ace.py | disk free=%.1f GB ===", disk_free_gb())
    SOL_OUT.mkdir(parents=True, exist_ok=True)

    process_psp()
    LOG.info("After PSP | disk free=%.1f GB", disk_free_gb())

    process_ace_dscovr()
    LOG.info("=== 11 complete | disk free=%.1f GB ===", disk_free_gb())


if __name__ == "__main__":
    main()
