"""
06_process_poes.py
Process POES MEPED CDF files → Parquet using cdflib.

Processes files in yearly batches to avoid OOM:
  1. For each satellite, iterate year subdirectories.
  2. Load that year's CDFs, concat, write poes_<sat>_<year>.parquet,
     delete that year's raw CDF files.
  3. After all years, merge year parquets into poes_<sat>.parquet,
     delete the year parquets.

Output: data/processed/atmospheric/
"""
from __future__ import annotations

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
ATM_OUT = PROCESSED_ROOT / "atmospheric"

POES_DIR = DATA_ROOT / "atmospheric" / "poes_meped"

# Priority list of flux variable names to search for (lowercased for matching)
FLUX_VAR_HINTS = [
    "ted_ele_flux", "ted_pro_flux",
    "ted_ele_diff_energies", "ted_pro_diff_energies",
    "mep_ele_flux", "mep_pro_flux",
    "mep_ele_tel0_flux", "mep_pro_tel0_flux",
    "mep_ele_tel90_flux", "mep_pro_tel90_flux",
    "mep_ele_total_flux", "mep_pro_total_flux",
    "mep_0_e1", "mep_0_e2", "mep_0_e3",
    "mep_0_p1", "mep_0_p2", "mep_0_p3", "mep_0_p4", "mep_0_p5", "mep_0_p6",
]


def _is_flux_var(vname: str) -> bool:
    vl = vname.lower()
    return any(h in vl for h in ("flux", "ted_ele", "ted_pro", "mep_ele", "mep_pro",
                                  "mep_0_e", "mep_0_p", "mep_90_"))


def _load_single_cdf(fp: Path, cdflib) -> tuple[pd.DataFrame | None, bool]:
    """Load one CDF file; return (DataFrame_or_None, should_delete).

    should_delete is True for empty/no-data files and successful files.
    It is False only when an unexpected exception prevents reading the file,
    preserving the original behaviour of not deleting corrupt/unreadable CDFs.
    """
    try:
        cdf = cdflib.CDF(str(fp))
        info = cdf.cdf_info()

        # Decode Epoch (CDF_EPOCH = milliseconds since 0 AD)
        epoch_var = None
        for ev in ("Epoch", "EPOCH", "epoch", "time"):
            if ev in info.zVariables or ev in info.rVariables:
                epoch_var = ev
                break

        if epoch_var is None:
            LOG.warning("No Epoch variable in %s", fp.name)
            return None, True

        epoch_data = cdf.varget(epoch_var)
        if epoch_data is None or len(epoch_data) == 0:
            return None, True

        datetimes = cdflib.cdfepoch.to_datetime(epoch_data)
        dt_index = pd.DatetimeIndex(datetimes, tz="UTC")

        # Collect all flux variables
        all_vars = list(info.zVariables) + list(info.rVariables)
        data_dict: dict[str, np.ndarray] = {}

        for vname in all_vars:
            if not _is_flux_var(vname):
                continue
            try:
                arr = cdf.varget(vname)
                if arr is None:
                    continue
                arr = np.asarray(arr, dtype=float)
                if arr.ndim == 1 and len(arr) == len(dt_index):
                    data_dict[vname] = arr
                elif arr.ndim == 2 and arr.shape[0] == len(dt_index):
                    for i in range(arr.shape[1]):
                        data_dict[f"{vname}_ch{i}"] = arr[:, i]
            except Exception as exc:
                LOG.debug("Var %s in %s: %s", vname, fp.name, exc)

        if not data_dict:
            return None, True

        df = pd.DataFrame(data_dict, index=dt_index).sort_index()
        return df, True

    except Exception as exc:
        LOG.warning("Failed to process %s: %s", fp.name, exc)
        # Do not delete files that threw unexpected exceptions
        return None, False


def _process_year(
    sat_name: str,
    year: str,
    cdf_files: list[Path],
    cdflib,
) -> Path | None:
    """Load all CDFs for one satellite-year, write a year parquet, delete raws.

    Returns the year parquet path on success, or None if no data was produced.
    Raw CDF files are deleted only after the year parquet is verified.
    """
    year_out = ATM_OUT / f"poes_{sat_name}_{year}.parquet"

    if year_out.exists():
        LOG.info("  SKIP existing year parquet %s", year_out.name)
        return year_out

    LOG.info("  %s/%s: loading %d CDF files …", sat_name, year, len(cdf_files))

    year_frames: list[pd.DataFrame] = []
    raw_to_delete: list[Path] = []

    for fp in cdf_files:
        df, should_delete = _load_single_cdf(fp, cdflib)
        if should_delete:
            raw_to_delete.append(fp)
        if df is not None:
            year_frames.append(df)

    if not year_frames:
        LOG.warning("  No data for %s year %s — deleting %d empty/bad CDFs",
                    sat_name, year, len(raw_to_delete))
        safe_delete(raw_to_delete)
        return None

    combined = pd.concat(year_frames).sort_index()
    del year_frames  # release per-file frames before writing

    combined = combined[~combined.index.duplicated(keep="first")]

    meta = {
        "title": f"POES MEPED Particle Fluxes — {sat_name} {year}",
        "source": "NOAA/NCEI POES",
        "references": "https://doi.org/10.7289/V55H7D60",
        "time_range": f"{combined.index.min()} / {combined.index.max()}",
        "units": "electrons/protons cm-2 sr-1 s-1",
    }
    save_parquet(combined, year_out, meta)
    del combined  # release year data after write

    if year_out.exists() and year_out.stat().st_size > 0:
        safe_delete(raw_to_delete)
        return year_out

    LOG.error("  Year parquet %s missing/empty — raw CDFs NOT deleted", year_out.name)
    return None


def _build_year_file_map(sat_dir: Path) -> dict[str, list[Path]]:
    """Return a mapping of year-string → sorted list of CDF paths.

    Prefers explicit year subdirectories (e.g. sat_dir/2013/*.cdf).
    Falls back to grouping flat files by year extracted from path or filename.
    """
    year_dirs = sorted(
        d for d in sat_dir.iterdir()
        if d.is_dir() and d.name.isdigit() and len(d.name) == 4
    )

    if year_dirs:
        year_file_map: dict[str, list[Path]] = {}
        for yd in year_dirs:
            files = sorted(yd.rglob("*.cdf"))
            if files:
                year_file_map[yd.name] = files
        return year_file_map

    # Fallback: group flat files by year from path components or filename
    all_files = sorted(sat_dir.rglob("*.cdf"))
    year_file_map = {}
    for fp in all_files:
        year = None
        for part in fp.parts:
            if len(part) == 4 and part.isdigit():
                year = part
                break
        if year is None:
            stem4 = fp.stem[:4]
            year = stem4 if (len(stem4) == 4 and stem4.isdigit()) else "unknown"
        year_file_map.setdefault(year, []).append(fp)
    return year_file_map


def _process_satellite(sat_name: str, sat_dir: Path) -> None:
    out = ATM_OUT / f"poes_{sat_name}.parquet"
    if out.exists():
        LOG.info("SKIP poes_%s.parquet", sat_name)
        return

    try:
        import cdflib
    except ImportError:
        LOG.error("cdflib not installed — cannot process POES data")
        return

    year_file_map = _build_year_file_map(sat_dir)
    if not year_file_map:
        LOG.warning("No CDF files found for satellite %s in %s", sat_name, sat_dir)
        return

    total_files = sum(len(v) for v in year_file_map.values())
    LOG.info("Processing %s: %d years, %d CDF files total …",
             sat_name, len(year_file_map), total_files)

    year_parquets: list[Path] = []
    for year in sorted(year_file_map):
        yp = _process_year(sat_name, year, year_file_map[year], cdflib)
        if yp is not None:
            year_parquets.append(yp)
        LOG.info("  After %s/%s | disk free=%.1f GB", sat_name, year, disk_free_gb())

    if not year_parquets:
        LOG.warning("No year parquets produced for %s", sat_name)
        return

    # Merge all year parquets into the final per-satellite parquet
    LOG.info("Merging %d year parquets for %s …", len(year_parquets), sat_name)
    year_dfs = [pd.read_parquet(yp) for yp in year_parquets]
    combined = pd.concat(year_dfs).sort_index()
    del year_dfs  # release year DataFrames before writing

    combined = combined[~combined.index.duplicated(keep="first")]

    meta = {
        "title": f"POES MEPED Particle Fluxes — {sat_name}",
        "source": "NOAA/NCEI POES",
        "references": "https://doi.org/10.7289/V55H7D60",
        "time_range": f"{combined.index.min()} / {combined.index.max()}",
        "units": "electrons/protons cm-2 sr-1 s-1",
    }
    save_parquet(combined, out, meta)
    del combined

    if out.exists() and out.stat().st_size > 0:
        safe_delete(year_parquets)
        register_output(MANIFEST, f"poes_{sat_name}", out, True, meta)
    else:
        LOG.error("Output %s missing/empty — year parquets NOT deleted", out.name)


def main() -> None:
    setup_logging()
    LOG.info("=== 06_process_poes.py | disk free=%.1f GB ===", disk_free_gb())
    ATM_OUT.mkdir(parents=True, exist_ok=True)

    if not POES_DIR.exists():
        LOG.warning("POES directory not found: %s", POES_DIR)
        return

    sat_dirs = sorted(d for d in POES_DIR.iterdir() if d.is_dir())
    if not sat_dirs:
        LOG.warning("No satellite subdirectories found in %s", POES_DIR)
        return

    for sat_dir in sat_dirs:
        _process_satellite(sat_dir.name, sat_dir)
        LOG.info("After %s | disk free=%.1f GB", sat_dir.name, disk_free_gb())

    LOG.info("=== 06 complete | disk free=%.1f GB ===", disk_free_gb())


if __name__ == "__main__":
    main()
