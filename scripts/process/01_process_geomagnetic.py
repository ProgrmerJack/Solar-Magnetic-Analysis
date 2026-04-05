"""
01_process_geomagnetic.py
Process geomagnetic (Dst/AE/Kp) and atmospheric catalog data → Parquet.
Outputs: data/processed/geomagnetic/, data/processed/atmospheric/
"""
import json
import logging
from pathlib import Path

import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).parent))
from _utils import (
    DATA_ROOT, PROCESSED_ROOT, LOG, setup_logging, save_parquet, register_output, disk_free_gb
)

MANIFEST = PROCESSED_ROOT / "manifest.json"

GEO_OUT = PROCESSED_ROOT / "geomagnetic"
ATM_OUT = PROCESSED_ROOT / "atmospheric"


def _read_json_files(directory: Path, pattern: str = "*.json") -> pd.DataFrame:
    """Read JSON files in *directory* matching *pattern*; each is a list of dicts."""
    frames = []
    for fp in sorted(directory.glob(pattern)):
        try:
            data = json.loads(fp.read_text(encoding="utf-8"))
            if isinstance(data, list):
                frames.append(pd.DataFrame(data))
            elif isinstance(data, dict):
                frames.append(pd.DataFrame([data]))
        except Exception as exc:
            LOG.warning("Could not read %s: %s", fp.name, exc)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def _parse_time_column(df: pd.DataFrame, col: str = "time_tag") -> pd.DataFrame:
    """Parse *col* to UTC-aware datetime and sort/deduplicate on it."""
    df[col] = pd.to_datetime(df[col], utc=True, errors="coerce")
    df = df.dropna(subset=[col]).sort_values(col).drop_duplicates(subset=[col])
    df = df.set_index(col)
    return df


# ---------------------------------------------------------------------------
# Dst / AE index
# ---------------------------------------------------------------------------
def _parse_kyoto_dst_csv(fp: Path) -> pd.DataFrame:
    """Parse a Kyoto WDC Dst CSV file (monthly, hourly rows)."""
    try:
        df = pd.read_csv(fp, comment="#", parse_dates=True)
        df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
        # Try to find a time column
        time_col = next((c for c in df.columns if "time" in c or "date" in c or "yr" in c), None)
        if time_col:
            df[time_col] = pd.to_datetime(df[time_col], utc=True, errors="coerce")
            df = df.rename(columns={time_col: "time_tag"})
        return df
    except Exception as exc:
        LOG.warning("Could not parse Kyoto CSV %s: %s", fp.name, exc)
        return pd.DataFrame()


def process_dst_ae() -> None:
    out = GEO_OUT / "dst_index.parquet"
    if out.exists():
        LOG.info("SKIP dst_index.parquet (already exists)")
        return

    src_dir = DATA_ROOT / "geomagnetic" / "dst_ae_index"
    if not src_dir.exists():
        LOG.warning("Missing directory: %s", src_dir)
        return

    # Only read Dst-specific JSON files (exclude kp files)
    dst_patterns = ["*dst*.json", "*geospace*.json", "kyoto*.json"]
    frames = []
    for pat in dst_patterns:
        df_pat = _read_json_files(src_dir, pat)
        if not df_pat.empty:
            frames.append(df_pat)

    # Also read historical Kyoto Dst CSV files
    hist_dir = src_dir / "kyoto_dst_historical"
    if hist_dir.exists():
        for fp in sorted(hist_dir.glob("*.csv")):
            df_csv = _parse_kyoto_dst_csv(fp)
            if not df_csv.empty:
                frames.append(df_csv)

    if not frames:
        LOG.warning("No Dst/AE data found in %s", src_dir)
        return

    df = pd.concat(frames, ignore_index=True)
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

    # Keep only numeric columns plus time_tag
    time_col = "time_tag"
    if time_col not in df.columns:
        LOG.warning("No time_tag column found in Dst data")
        return
    numeric_cols = [c for c in df.columns if c != time_col
                    and df[c].dtype not in ("object", "bool")]
    df = df[[time_col] + numeric_cols]
    df = _parse_time_column(df, time_col)

    meta = {
        "title": "Dst and AE Geomagnetic Indices",
        "source": "Kyoto WDC / NOAA SWPC",
        "references": "https://wdc.kugi.kyoto-u.ac.jp/",
        "units": "nT",
        "time_range": f"{df.index.min()} / {df.index.max()}",
    }
    save_parquet(df, out, meta)
    register_output(MANIFEST, "dst_ae_index", out, False, meta)


# ---------------------------------------------------------------------------
# Kp index
# ---------------------------------------------------------------------------
def _parse_gfz_kp_txt(fp: Path) -> pd.DataFrame:
    """Parse GFZ Kp file: YYYY MM DD DOY DOY.5 BART_ROT BART_INT KP1..KP8 AP1..AP8 AP SN F107 F107adj DFLAG."""
    records = []
    with open(fp, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 24:
                continue
            try:
                yr, mo, da = int(parts[0]), int(parts[1]), int(parts[2])
                kp_vals = [float(x) for x in parts[7:15]]   # 8 three-hourly Kp
                ap_vals = [float(x) for x in parts[15:23]]  # 8 three-hourly ap
                ap_daily = float(parts[23])
                for i, (kp, ap) in enumerate(zip(kp_vals, ap_vals)):
                    hour = i * 3
                    records.append({
                        "time_tag": pd.Timestamp(year=yr, month=mo, day=da, hour=hour, tz="UTC"),
                        "kp": kp,
                        "ap": ap,
                        "ap_daily": ap_daily,
                    })
            except (ValueError, IndexError):
                continue
    return pd.DataFrame(records)


def _kp_str_to_float(s: str) -> float:
    """Convert Kp string like '5Z' / '4+' / '3-' to float."""
    import re
    m = re.match(r"(\d+(?:\.\d+)?)", str(s))
    if not m:
        return float("nan")
    base = float(m.group(1))
    if "+" in str(s):
        return base + 1 / 3
    if "-" in str(s):
        return base - 1 / 3
    return base


def process_kp() -> None:
    out = GEO_OUT / "kp_index.parquet"
    if out.exists():
        LOG.info("SKIP kp_index.parquet (already exists)")
        return

    frames = []

    # 1. Historical GFZ Kp txt files (primary long record)
    kp_txt_dir = DATA_ROOT / "geomagnetic" / "kp_index"
    if kp_txt_dir.exists():
        for fp in sorted(kp_txt_dir.glob("*.txt")):
            df_txt = _parse_gfz_kp_txt(fp)
            if not df_txt.empty:
                frames.append(df_txt)

    # 2. Recent Kp JSON (from SWPC, in dst_ae_index/)
    dst_dir = DATA_ROOT / "geomagnetic" / "dst_ae_index"
    if dst_dir.exists():
        for pat in ["kp_1min.json", "noaa_planetary_kp*.json", "boulder_k_1min.json"]:
            for fp in sorted(dst_dir.glob(pat)):
                try:
                    data = json.loads(fp.read_text(encoding="utf-8"))
                    df_j = pd.DataFrame(data) if isinstance(data, list) else pd.DataFrame([data])
                    df_j.columns = [c.strip().lower().replace(" ", "_") for c in df_j.columns]
                    # Convert any string Kp column to float
                    for col in ["kp", "kp_index", "k_index"]:
                        if col in df_j.columns:
                            if df_j[col].dtype == object:
                                df_j[col] = df_j[col].apply(_kp_str_to_float)
                            df_j = df_j.rename(columns={col: "kp"})
                            break
                    # Drop duplicate string kp columns
                    for drop_col in ["estimated_kp", "observed", "noaa_scale", "a_running",
                                     "station_count"]:
                        if drop_col in df_j.columns:
                            df_j = df_j.drop(columns=[drop_col])
                    if "time_tag" in df_j.columns:
                        frames.append(df_j)
                except Exception as exc:
                    LOG.warning("Could not read Kp JSON %s: %s", fp.name, exc)

    if not frames:
        LOG.warning("No Kp data found")
        return

    df = pd.concat(frames, ignore_index=True)
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
    # Deduplicate columns (keep first occurrence of any name)
    df = df.loc[:, ~df.columns.duplicated()]
    # Ensure kp is numeric
    if "kp" in df.columns and df["kp"].dtype == object:
        df["kp"] = df["kp"].apply(_kp_str_to_float)

    df = _parse_time_column(df, "time_tag")

    meta = {
        "title": "Kp Geomagnetic Index (GFZ + NOAA SWPC)",
        "source": "GFZ Potsdam / NOAA SWPC",
        "references": "https://doi.org/10.1029/2020SW002641",
        "units": "dimensionless (0–9)",
        "time_range": f"{df.index.min()} / {df.index.max()}",
    }
    save_parquet(df, out, meta)
    register_output(MANIFEST, "kp_index", out, False, meta)


# ---------------------------------------------------------------------------
# SSW catalog
# ---------------------------------------------------------------------------
def process_ssw_catalog() -> None:
    out = ATM_OUT / "ssw_catalog.parquet"
    if out.exists():
        LOG.info("SKIP ssw_catalog.parquet (already exists)")
        return

    src_dir = DATA_ROOT / "atmospheric" / "ssw_catalog"
    if not src_dir.exists():
        LOG.warning("Missing directory: %s", src_dir)
        return

    frames = []
    for fp in sorted(src_dir.iterdir()):
        try:
            if fp.name == "metadata.json":
                continue  # skip metadata — not tabular catalog data
            if fp.suffix.lower() == ".csv":
                frames.append(pd.read_csv(fp))
            elif fp.suffix.lower() == ".json":
                data = json.loads(fp.read_text(encoding="utf-8"))
                if isinstance(data, list) and data and isinstance(data[0], dict):
                    frames.append(pd.DataFrame(data))
                # skip nested/metadata JSON structures
        except Exception as exc:
            LOG.warning("Could not read %s: %s", fp.name, exc)

    if not frames:
        LOG.warning("No SSW catalog files found in %s", src_dir)
        return

    df = pd.concat(frames, ignore_index=True)
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

    # Build datetime from year/month/day columns if present
    if {"year", "month", "day"}.issubset(df.columns):
        df["onset_date"] = pd.to_datetime(
            df[["year", "month", "day"]].astype(int).rename(
                columns={"year": "year", "month": "month", "day": "day"}),
            utc=True, errors="coerce"
        )
        df = df.drop(columns=["year", "month", "day"], errors="ignore")
        df = df.set_index("onset_date").sort_index()
    else:
        # Parse onset_date if present as a string column
        date_cols = [c for c in df.columns if "date" in c or "onset" in c or "time" in c]
        if date_cols:
            df[date_cols[0]] = pd.to_datetime(df[date_cols[0]], utc=True, errors="coerce")
            df = df.set_index(date_cols[0]).sort_index()

    meta = {
        "title": "Sudden Stratospheric Warming (SSW) Catalog",
        "source": "Various (Charlton & Polvani 2007, Butler et al. 2017)",
        "references": "https://doi.org/10.1175/2007JCLI1996.1",
        "time_range": f"{df.index.min()} / {df.index.max()}" if not df.empty else "unknown",
    }
    save_parquet(df, out, meta)
    register_output(MANIFEST, "ssw_catalog", out, False, meta)


# ---------------------------------------------------------------------------
# Climate indices
# ---------------------------------------------------------------------------
def _parse_noaa_text(fp: Path) -> pd.DataFrame:
    """Parse NOAA-style year + 12 monthly values text file."""
    records = []
    index_name = fp.stem
    with open(fp, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 13:
                continue
            try:
                year = int(parts[0])
            except ValueError:
                continue
            for month_idx, val_str in enumerate(parts[1:13], start=1):
                try:
                    val = float(val_str)
                    if abs(val) >= 999.0:
                        val = float("nan")
                    records.append({"time": pd.Timestamp(year=year, month=month_idx, day=1, tz="UTC"),
                                    "index_name": index_name, "value": val})
                except ValueError:
                    pass
    return pd.DataFrame(records)


def process_climate_indices() -> None:
    out = ATM_OUT / "climate_indices.parquet"
    if out.exists():
        LOG.info("SKIP climate_indices.parquet (already exists)")
        return

    src_dir = DATA_ROOT / "atmospheric" / "climate_indices"
    if not src_dir.exists():
        LOG.warning("Missing directory: %s", src_dir)
        return

    frames = []
    for fp in sorted(src_dir.iterdir()):
        try:
            if fp.suffix.lower() in (".txt", ".dat", ""):
                df = _parse_noaa_text(fp)
                if not df.empty:
                    frames.append(df)
            elif fp.suffix.lower() == ".json":
                data = json.loads(fp.read_text(encoding="utf-8"))
                frames.append(pd.DataFrame(data) if isinstance(data, list) else pd.DataFrame([data]))
            elif fp.suffix.lower() == ".csv":
                frames.append(pd.read_csv(fp))
        except Exception as exc:
            LOG.warning("Could not read %s: %s", fp.name, exc)

    if not frames:
        LOG.warning("No climate index files found in %s", src_dir)
        return

    long_df = pd.concat(frames, ignore_index=True)
    long_df.columns = [c.strip().lower().replace(" ", "_") for c in long_df.columns]

    # Pivot to wide if we have index_name + value columns
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
        "title": "Climate Indices (AO, NAM, QBO, ENSO, etc.)",
        "source": "NOAA CPC / ECMWF",
        "references": "https://www.cpc.ncep.noaa.gov/data/indices/",
        "time_range": f"{df_out.index.min()} / {df_out.index.max()}" if not df_out.empty else "unknown",
    }
    save_parquet(df_out, out, meta)
    register_output(MANIFEST, "climate_indices", out, False, meta)


def main() -> None:
    setup_logging()
    LOG.info("=== 01_process_geomagnetic.py | disk free=%.1f GB ===", disk_free_gb())
    GEO_OUT.mkdir(parents=True, exist_ok=True)
    ATM_OUT.mkdir(parents=True, exist_ok=True)

    process_dst_ae()
    process_kp()
    process_ssw_catalog()
    process_climate_indices()

    LOG.info("=== 01 complete | disk free=%.1f GB ===", disk_free_gb())


if __name__ == "__main__":
    main()
