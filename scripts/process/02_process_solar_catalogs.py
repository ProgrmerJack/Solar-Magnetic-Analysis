"""
02_process_solar_catalogs.py
Process solar event catalogs, SDO HMI SHARP, solar indices, and NCEP CSV → Parquet.
Outputs: data/processed/solar/, data/processed/atmospheric/

NCEP data is already pre-processed into polar-cap-mean CSV files (NOT NetCDF):
  data/atmospheric/ncep_stratosphere/{air,hgt,uwnd}/{var}_{YYYY}.csv
  CSV columns: date, level_hPa, polar_cap_mean
"""
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).parent))
from _utils import (
    DATA_ROOT, PROCESSED_ROOT, LOG, setup_logging,
    save_parquet, register_output, disk_free_gb,
)

MANIFEST = PROCESSED_ROOT / "manifest.json"
SOL_OUT = PROCESSED_ROOT / "solar"
ATM_OUT = PROCESSED_ROOT / "atmospheric"


# ---------------------------------------------------------------------------
# Flare catalog (DONKI)
# ---------------------------------------------------------------------------
def process_flare_catalog() -> None:
    out = SOL_OUT / "flares.parquet"
    if out.exists():
        LOG.info("SKIP flares.parquet")
        return

    fp = DATA_ROOT / "solar" / "flare_catalog" / "donki_flares_1995_present.json"
    if not fp.exists():
        LOG.warning("Missing: %s", fp)
        return

    try:
        raw = json.loads(fp.read_text(encoding="utf-8"))
    except Exception as exc:
        LOG.warning("Could not read flare catalog: %s", exc)
        return

    records = []
    for item in raw:
        try:
            linked = item.get("linkedEvents") or []
            records.append({
                "flareID": item.get("flareID"),
                "beginTime": pd.to_datetime(item.get("beginTime"), utc=True, errors="coerce"),
                "peakTime": pd.to_datetime(item.get("peakTime"), utc=True, errors="coerce"),
                "endTime": pd.to_datetime(item.get("endTime"), utc=True, errors="coerce"),
                "classType": item.get("classType"),
                "sourceLocation": item.get("sourceLocation"),
                "activeRegionNum": item.get("activeRegionNum"),
                "n_linked_events": len(linked),
            })
        except Exception as exc:
            LOG.warning("Flare record error: %s", exc)

    df = pd.DataFrame(records).set_index("beginTime").sort_index()
    meta = {
        "title": "DONKI Solar Flare Catalog 1995-present",
        "source": "NASA DONKI",
        "references": "https://kauai.ccmc.gsfc.nasa.gov/DONKI/",
        "time_range": f"{df.index.min()} / {df.index.max()}",
        "units": "GOES class (A/B/C/M/X)",
    }
    save_parquet(df, out, meta)
    register_output(MANIFEST, "flare_catalog", out, False, meta)


# ---------------------------------------------------------------------------
# CME catalog (DONKI)
# ---------------------------------------------------------------------------
def process_cme_catalog() -> None:
    out = SOL_OUT / "cme_catalog.parquet"
    if out.exists():
        LOG.info("SKIP cme_catalog.parquet")
        return

    fp = DATA_ROOT / "solar" / "flare_catalog" / "donki_cme_1995_present.json"
    if not fp.exists():
        LOG.warning("Missing: %s", fp)
        return

    try:
        raw = json.loads(fp.read_text(encoding="utf-8"))
    except Exception as exc:
        LOG.warning("Could not read CME catalog: %s", exc)
        return

    records = []
    for item in raw:
        try:
            analyses = item.get("cmeAnalyses") or []
            # First mostAccurate entry
            accurate = next((a for a in analyses if a.get("isMostAccurate")), analyses[0] if analyses else {})
            records.append({
                "activityID": item.get("activityID"),
                "startTime": pd.to_datetime(item.get("startTime"), utc=True, errors="coerce"),
                "speed_km_s": accurate.get("speed"),
                "halfAngle_deg": accurate.get("halfAngle"),
                "type": accurate.get("type"),
                "n_analyses": len(analyses),
            })
        except Exception as exc:
            LOG.warning("CME record error: %s", exc)

    df = pd.DataFrame(records).set_index("startTime").sort_index()
    meta = {
        "title": "DONKI CME Catalog 1995-present",
        "source": "NASA DONKI",
        "references": "https://kauai.ccmc.gsfc.nasa.gov/DONKI/",
        "time_range": f"{df.index.min()} / {df.index.max()}",
        "units": "speed km/s, halfAngle deg",
    }
    save_parquet(df, out, meta)
    register_output(MANIFEST, "cme_catalog", out, False, meta)


# ---------------------------------------------------------------------------
# SDO HMI SHARP
# ---------------------------------------------------------------------------
def process_sdo_hmi_sharp() -> None:
    out = SOL_OUT / "sdo_hmi_sharp.parquet"
    if out.exists():
        LOG.info("SKIP sdo_hmi_sharp.parquet")
        return

    src_dir = DATA_ROOT / "solar" / "sdo_hmi_aia"
    if not src_dir.exists():
        LOG.warning("Missing directory: %s", src_dir)
        return

    frames = []
    for fp in sorted(src_dir.glob("hmi_sharp_flux_*.csv")):
        try:
            df = pd.read_csv(fp, low_memory=False)
            frames.append(df)
        except Exception as exc:
            LOG.warning("Could not read %s: %s", fp.name, exc)

    if not frames:
        LOG.warning("No HMI SHARP CSV files found")
        return

    df = pd.concat(frames, ignore_index=True)
    df.columns = [c.strip() for c in df.columns]

    # Parse time column (try common names)
    for tcol in ("time", "T_REC", "DATE", "datetime"):
        if tcol in df.columns:
            df[tcol] = pd.to_datetime(df[tcol], utc=True, errors="coerce")
            df = df.dropna(subset=[tcol]).set_index(tcol).sort_index()
            break

    meta = {
        "title": "SDO/HMI SHARP Active Region Flux Parameters",
        "source": "NASA SDO / JSOC",
        "references": "https://doi.org/10.1007/s11207-014-0529-3",
        "time_range": f"{df.index.min()} / {df.index.max()}",
        "units": "USFLUX Mx, field Gauss, TOTUSJH mA/m",
    }
    save_parquet(df, out, meta)
    register_output(MANIFEST, "sdo_hmi_sharp", out, False, meta)


# ---------------------------------------------------------------------------
# Solar indices
# ---------------------------------------------------------------------------
def _try_parse_solar_index_file(fp: Path) -> pd.DataFrame:
    """Attempt multiple parse strategies for solar index text files."""
    # Strategy 1: year + 12 monthly values
    records = []
    index_name = fp.stem
    with open(fp, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            # datetime + value (2 columns)
            if len(parts) == 2:
                try:
                    t = pd.to_datetime(parts[0], utc=True)
                    v = float(parts[1])
                    records.append({"time": t, "index_name": index_name, "value": v})
                    continue
                except Exception:
                    pass
            # year + 12 monthly values
            if len(parts) >= 13:
                try:
                    year = int(parts[0])
                    for m, vs in enumerate(parts[1:13], 1):
                        v = float(vs)
                        if abs(v) >= 9990.0:
                            v = float("nan")
                        records.append({"time": pd.Timestamp(year=year, month=m, day=1, tz="UTC"),
                                        "index_name": index_name, "value": v})
                    continue
                except Exception:
                    pass
    return pd.DataFrame(records)


def process_solar_indices() -> None:
    out = SOL_OUT / "solar_indices.parquet"
    if out.exists():
        LOG.info("SKIP solar_indices.parquet")
        return

    # Source is under atmospheric/, not solar/
    src_dir = DATA_ROOT / "atmospheric" / "solar_indices"
    if not src_dir.exists():
        LOG.warning("Missing directory: %s", src_dir)
        return

    frames = []
    for fp in sorted(src_dir.iterdir()):
        if fp.is_dir():
            continue
        try:
            if fp.suffix.lower() in (".txt", ".dat", ""):
                df = _try_parse_solar_index_file(fp)
            elif fp.suffix.lower() == ".csv":
                # comment='#' handles NOAA-style header comment lines
                df = pd.read_csv(fp, comment="#", low_memory=False)
            elif fp.suffix.lower() == ".json":
                data = json.loads(fp.read_text(encoding="utf-8"))
                df = pd.DataFrame(data) if isinstance(data, list) else pd.DataFrame([data])
            else:
                continue
            if not df.empty:
                frames.append(df)
        except Exception as exc:
            LOG.warning("Could not read %s: %s", fp.name, exc)

    if not frames:
        LOG.warning("No solar index files found in %s", src_dir)
        return

    long_df = pd.concat(frames, ignore_index=True)
    long_df.columns = [str(c).strip().lower().replace(" ", "_") for c in long_df.columns]

    if {"time", "index_name", "value"}.issubset(long_df.columns):
        long_df["time"] = pd.to_datetime(long_df["time"], utc=True, errors="coerce")
        wide = long_df.dropna(subset=["time"]).pivot_table(
            index="time", columns="index_name", values="value", aggfunc="mean"
        )
        wide.columns.name = None
        df_out = wide
    else:
        df_out = long_df

    meta = {
        "title": "Solar Activity Indices (F10.7, SSN, etc.)",
        "source": "NOAA SWPC / SIDC / LASP",
        "references": "https://www.sidc.be/silso/datafiles",
        "time_range": f"{df_out.index.min()} / {df_out.index.max()}" if not df_out.empty else "unknown",
        "units": "F10.7 sfu; SSN dimensionless",
    }
    save_parquet(df_out, out, meta)
    register_output(MANIFEST, "solar_indices", out, False, meta)


# ---------------------------------------------------------------------------
# NCEP stratosphere / troposphere — polar-cap-mean CSV files
#
# Confirmed structure (NOT NetCDF):
#   {domain}/air/{var}_{YYYY}.csv   — air temperature  (K)
#   {domain}/hgt/{var}_{YYYY}.csv   — geopotential hgt (m)
#   {domain}/uwnd/{var}_{YYYY}.csv  — zonal wind        (m/s)
#   {domain}/{var}.csv              — optional pre-merged root-level summary
#
# CSV columns (confirmed): date, level_hPa, polar_cap_mean
# Output: single merged Parquet with MultiIndex (date, level_hPa),
#         columns: air_K, hgt_m, uwnd_ms
# ---------------------------------------------------------------------------
_NCEP_VARS: dict[str, tuple[str, str]] = {
    "air":  ("air_K",   "Kelvin"),
    "hgt":  ("hgt_m",   "meters"),
    "uwnd": ("uwnd_ms", "m/s"),
}


def _read_ncep_var(src_dir: Path, var_name: str) -> pd.DataFrame | None:
    """
    Concatenate all annual CSV files for one NCEP variable.
    Returns a tidy DataFrame with columns [date, level_hPa, <col_name>].
    """
    frames: list[pd.DataFrame] = []

    # Annual CSVs inside variable subdir
    subdir = src_dir / var_name
    if subdir.exists():
        for fp in sorted(subdir.glob("*.csv")):
            try:
                frames.append(pd.read_csv(fp, comment="#"))
            except Exception as exc:
                LOG.warning("Could not read %s: %s", fp.name, exc)

    # Optional root-level summary CSV (e.g., ncep_stratosphere/air.csv)
    root_csv = src_dir / f"{var_name}.csv"
    if root_csv.exists():
        try:
            frames.append(pd.read_csv(root_csv, comment="#"))
        except Exception as exc:
            LOG.warning("Could not read %s: %s", root_csv.name, exc)

    if not frames:
        LOG.warning("No CSV files for NCEP var=%s under %s", var_name, src_dir)
        return None

    df = pd.concat(frames, ignore_index=True)
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

    # Identify date, level, and value columns
    date_col  = next((c for c in df.columns if "date" in c or "time" in c), None)
    level_col = next((c for c in df.columns if "level" in c or "hpa" in c
                      or "pres" in c or "plev" in c), None)
    # Value column = anything that isn't the date or level
    skip = {date_col, level_col}
    value_col = next((c for c in df.columns if c not in skip and c), None)

    if date_col is None or value_col is None:
        LOG.warning("Cannot identify columns in NCEP %s: %s", var_name, list(df.columns))
        return None

    df[date_col] = pd.to_datetime(df[date_col], utc=True, errors="coerce")
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce")
    df = df.dropna(subset=[date_col])

    col_name = _NCEP_VARS[var_name][0]   # e.g. "air_K"
    if level_col:
        df[level_col] = pd.to_numeric(df[level_col], errors="coerce")
        out = df[[date_col, level_col, value_col]].copy()
        out.columns = ["date", "level_hPa", col_name]
    else:
        out = df[[date_col, value_col]].copy()
        out.columns = ["date", col_name]

    return out.drop_duplicates()


def _process_ncep_domain(
    src_dir: Path,
    out_name: str,
    title: str,
    key: str,
) -> None:
    out = ATM_OUT / out_name
    if out.exists():
        LOG.info("SKIP %s", out_name)
        return

    if not src_dir.exists():
        LOG.warning("Missing NCEP directory: %s", src_dir)
        return

    var_frames: dict[str, pd.DataFrame] = {}
    for var_name in _NCEP_VARS:
        df = _read_ncep_var(src_dir, var_name)
        if df is not None and not df.empty:
            var_frames[var_name] = df

    if not var_frames:
        LOG.warning("No NCEP data loaded from %s", src_dir)
        return

    # Determine merge keys — use both date+level if every frame has level
    has_level = all("level_hPa" in df.columns for df in var_frames.values())
    merge_keys = ["date", "level_hPa"] if has_level else ["date"]

    merged: pd.DataFrame | None = None
    for df in var_frames.values():
        if merged is None:
            merged = df
        else:
            common = [k for k in merge_keys if k in merged.columns and k in df.columns]
            merged = pd.merge(merged, df, on=common, how="outer")

    if merged is None or merged.empty:
        LOG.warning("Merge produced empty result for %s", key)
        return

    idx_cols = [c for c in ("date", "level_hPa") if c in merged.columns]
    merged = merged.set_index(idx_cols).sort_index()

    date_vals = merged.index.get_level_values("date") if "date" in merged.index.names else merged.index
    meta = {
        "title": title,
        "source": "NOAA/NCEP-NCAR Reanalysis I",
        "references": "https://doi.org/10.1175/1520-0477(1996)077<0437:TNYRP>2.0.CO;2",
        "units": "air_K=Kelvin; hgt_m=meters; uwnd_ms=m/s; polar cap mean ≥60°N",
        "time_range": f"{date_vals.min()} / {date_vals.max()}",
    }
    save_parquet(merged, out, meta)
    register_output(MANIFEST, key, out, False, meta)


def process_ncep() -> None:
    ATM_OUT.mkdir(parents=True, exist_ok=True)
    _process_ncep_domain(
        DATA_ROOT / "atmospheric" / "ncep_stratosphere",
        "ncep_stratosphere.parquet",
        "NCEP/NCAR Reanalysis I — Stratospheric Polar Cap Means (≥60°N)",
        "ncep_stratosphere",
    )
    _process_ncep_domain(
        DATA_ROOT / "atmospheric" / "ncep_troposphere",
        "ncep_troposphere.parquet",
        "NCEP/NCAR Reanalysis I — Tropospheric Polar Cap Means (≥60°N)",
        "ncep_troposphere",
    )


def main() -> None:
    setup_logging()
    LOG.info("=== 02_process_solar_catalogs.py | disk free=%.1f GB ===", disk_free_gb())
    SOL_OUT.mkdir(parents=True, exist_ok=True)
    ATM_OUT.mkdir(parents=True, exist_ok=True)

    process_flare_catalog()
    process_cme_catalog()
    process_sdo_hmi_sharp()
    process_solar_indices()
    process_ncep()

    LOG.info("=== 02 complete | disk free=%.1f GB ===", disk_free_gb())


if __name__ == "__main__":
    main()
