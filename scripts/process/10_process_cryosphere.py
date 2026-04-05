"""
10_process_cryosphere.py
Process cryosphere tabular and spatial datasets → Parquet.
SLF avalanche, SNOTEL, NGI Norway, CAIC, IMS Snow.
Deletes raw IMS .gz files after processing (large).
Output: data/processed/cryosphere/
"""
import gzip
import logging
import re
from datetime import datetime
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
CRYO_OUT = PROCESSED_ROOT / "cryosphere"

# IMS 4km grid parameters
IMS_GRID_SIZE = 6144
IMS_CELL_M    = 4000.0
IMS_HALF_EXT  = IMS_GRID_SIZE * IMS_CELL_M / 2   # 12,288,000 m from pole

# Fallback hardcoded Alps sub-region (44-48°N, 5-11°E) in IMS LAEA NH grid
# These approximate indices are used when pyproj is unavailable.
IMS_ALPS_ROW_FALLBACK = slice(2900, 3000)
IMS_ALPS_COL_FALLBACK = slice(3100, 3200)

# IMS value legend: 1=open ocean, 2=land/no-snow, 3=sea-ice, 4=snow
IMS_SNOW_VALUE = 4
IMS_LAND_VALUES = (2, 4)   # land pixels (with or without snow)

# Source directory with YYYY sub-folders
IMS_SNOW_DIR = DATA_ROOT / "cryosphere" / "ims_snow"


# ---------------------------------------------------------------------------
# SLF Avalanche
# ---------------------------------------------------------------------------
def _normalize_df(df: pd.DataFrame) -> pd.DataFrame:
    # Replace spaces, hyphens, AND dots (SLF uses dot-separated names)
    df.columns = [
        c.strip().lower().replace(" ", "_").replace("-", "_").replace(".", "_")
        for c in df.columns
    ]
    return df


def _parse_date_col(df: pd.DataFrame) -> pd.DataFrame:
    """Try to coerce the first date-like column to UTC datetime."""
    date_hints = [c for c in df.columns if any(k in c for k in ("date", "time", "onset", "day"))]
    if date_hints:
        col = date_hints[0]
        df[col] = pd.to_datetime(df[col], utc=True, errors="coerce")
        df = df.dropna(subset=[col]).set_index(col).sort_index()
    return df


def process_slf_avalanche() -> None:
    src_dir = DATA_ROOT / "cryosphere" / "slf_avalanche"
    if not src_dir.exists():
        LOG.warning("Missing: %s", src_dir)
        return

    csv_files = sorted(src_dir.glob("*.csv"))
    if not csv_files:
        LOG.warning("No CSV files found in %s", src_dir)
        return

    for fp in csv_files:
        fname_lower = fp.stem.lower()

        # Decide output name from filename — more specific matches first
        if any(k in fname_lower for k in ("accident", "killed", "fatalities")):
            out_key = "slf_accidents"
        elif "activity_data" in fname_lower or "avalanche_activity" in fname_lower:
            out_key = "slf_activity"
        elif "snow_avalanche_data" in fname_lower or "avalanche_data" in fname_lower:
            out_key = "slf_snow_events"
        elif "wet" in fname_lower and "model" in fname_lower:
            out_key = "slf_wet_model"
        else:
            out_key = f"slf_{fp.stem}"

        out = CRYO_OUT / f"{out_key}.parquet"
        if out.exists():
            LOG.info("SKIP %s", out.name)
            continue

        try:
            # SLF avalanche-accident CSVs have 3 metadata lines before the column header:
            #   Line 0: "WSL Institute for Snow and Avalanche Research SLF"
            #   Line 1: "Avalanche accidents in Switzerland since 1970-1971"
            #   Line 2: "Update: 2025-12-15 08:58 (UTC)"
            #   Line 3: actual column header (dot-separated, e.g. avalanche.id, date.quality)
            # Activity/snow files have no metadata header but may use ; separator.
            is_accident_file = any(k in fname_lower for k in ("accident", "killed", "fatalities",
                                                               "avalanche-accidents-in-sw"))
            skip = 3 if is_accident_file else 0
            # Use sep=None + engine='python' to auto-detect comma vs semicolon
            df = pd.read_csv(fp, skiprows=skip, sep=None,
                             engine="python", encoding="utf-8",
                             encoding_errors="replace")
            df = _normalize_df(df)  # converts dot-names: avalanche.id → avalanche_id

            # For accident files: no year-column filter needed — date column handles it.
            # (hydrological_year uses "YYYY/YY" format which breaks isdigit checks)

            df = _parse_date_col(df)

            meta = {
                "title": f"SLF Avalanche Data — {fp.stem}",
                "source": "WSL Institute for Snow and Avalanche Research SLF",
                "references": "https://www.slf.ch/",
                "time_range": f"{df.index.min()} / {df.index.max()}" if not df.empty else "unknown",
            }
            save_parquet(df, out, meta)
            register_output(MANIFEST, out_key, out, False, meta)
        except Exception as exc:
            LOG.warning("Error processing %s: %s", fp.name, exc)


# ---------------------------------------------------------------------------
# SNOTEL
#
# Confirmed structure: data/cryosphere/snotel/station_data/{state}/{station_id}/{ELEMENT}_daily.json
# 6 elements per station: WTEQ, SNWD, PREC, TAVG, TMAX, TMIN
#
# JSON layout:
#   [{"stationTriplet": "0280:AK:COOP",
#     "data": [{"stationElement": {..., "storedUnitCode": "in"},
#               "values": [{"date": "2001-09-30", "value": 0.0}, ...]}]}]
#
# Unit conversions:
#   WTEQ / SNWD / PREC : inches  → mm  (× 25.4)
#   TAVG / TMAX / TMIN : °F      → °C  ((F-32) × 5/9)
# ---------------------------------------------------------------------------
_SNOTEL_ELEMENTS: dict[str, tuple[str, object]] = {
    "WTEQ": ("wteq_mm", lambda x: x * 25.4),
    "SNWD": ("snwd_mm", lambda x: x * 25.4),
    "PREC": ("prec_mm", lambda x: x * 25.4),
    "TAVG": ("tavg_c",  lambda x: (x - 32.0) * 5.0 / 9.0),
    "TMAX": ("tmax_c",  lambda x: (x - 32.0) * 5.0 / 9.0),
    "TMIN": ("tmin_c",  lambda x: (x - 32.0) * 5.0 / 9.0),
}


def process_snotel() -> None:
    import json as _json

    out = CRYO_OUT / "snotel_daily.parquet"
    if out.exists():
        LOG.info("SKIP snotel_daily.parquet")
        return

    src_dir = DATA_ROOT / "cryosphere" / "snotel" / "station_data"
    if not src_dir.exists():
        LOG.warning("Missing: %s", src_dir)
        return

    # Recurse: state_code/station_id/ELEMENT_daily.json
    json_files = sorted(src_dir.rglob("*_daily.json"))
    if not json_files:
        LOG.warning("No SNOTEL JSON files found under %s", src_dir)
        return

    LOG.info("SNOTEL: found %d element JSON files", len(json_files))

    # station_triplet → {col_name: pd.Series}
    station_series: dict[str, dict[str, pd.Series]] = {}

    for fp in json_files:
        # Derive element from filename stem: "WTEQ_daily" → "WTEQ"
        element = fp.stem.replace("_daily", "").upper()
        if element not in _SNOTEL_ELEMENTS:
            LOG.debug("Unknown SNOTEL element %s in %s — skipping", element, fp.name)
            continue

        col_name, converter = _SNOTEL_ELEMENTS[element]

        try:
            raw = _json.loads(fp.read_text(encoding="utf-8"))
            if not isinstance(raw, list):
                raw = [raw]

            for station_obj in raw:
                station_triplet = station_obj.get("stationTriplet", fp.parent.name)
                data_list = station_obj.get("data", [])
                if not data_list:
                    continue

                # Navigate confirmed JSON structure
                first_data = data_list[0] if isinstance(data_list, list) else data_list
                values_list = first_data.get("values", [])

                dates_out, vals_out = [], []
                for v in values_list:
                    d = v.get("date")
                    val = v.get("value")
                    if d is None or val is None:
                        continue
                    try:
                        dates_out.append(pd.Timestamp(d, tz="UTC"))
                        vals_out.append(converter(float(val)))
                    except (ValueError, TypeError):
                        pass

                if not dates_out:
                    continue

                series = pd.Series(
                    vals_out,
                    index=pd.DatetimeIndex(dates_out),
                    name=col_name,
                    dtype=float,
                )
                if station_triplet not in station_series:
                    station_series[station_triplet] = {}
                station_series[station_triplet][col_name] = series

        except Exception as exc:
            LOG.warning("Error reading %s: %s", fp.name, exc)

    if not station_series:
        LOG.warning("No SNOTEL data parsed")
        return

    station_frames: list[pd.DataFrame] = []
    for station_triplet, elem_dict in station_series.items():
        # Wide merge: one row per date per station
        merged = pd.concat(elem_dict.values(), axis=1)
        merged.index.name = "date"
        merged["station_id"] = station_triplet
        station_frames.append(merged.reset_index())

    result = pd.concat(station_frames, ignore_index=True)

    # Ensure all canonical columns exist
    for col in ("wteq_mm", "snwd_mm", "prec_mm", "tavg_c", "tmax_c", "tmin_c"):
        if col not in result.columns:
            result[col] = np.nan

    meta = {
        "title": "SNOTEL Daily Snow/Weather Observations",
        "source": "USDA NRCS SNOTEL Network",
        "references": "https://www.nrcs.usda.gov/wps/portal/wcc/home/",
        "units": "wteq/snwd/prec mm; tavg/tmax/tmin °C",
    }
    save_parquet(result, out, meta)
    register_output(MANIFEST, "snotel_daily", out, False, meta)


# ---------------------------------------------------------------------------
# NGI Norway — seNorge SWE dataset
#
# Confirmed file: data/cryosphere/ngi_norway/senorge_swe/norway_swe_daily_1957_2025.csv
# Columns: date, swe_national_mean_mm, swe_south_mean_mm, swe_central_mean_mm, swe_north_mean_mm
# Clean format — read directly with pd.read_csv.
# ---------------------------------------------------------------------------
def process_ngi_norway() -> None:
    out = CRYO_OUT / "norway_avalanche.parquet"
    if out.exists():
        LOG.info("SKIP norway_avalanche.parquet")
        return

    src_dir = DATA_ROOT / "cryosphere" / "ngi_norway"
    if not src_dir.exists():
        LOG.warning("Missing: %s", src_dir)
        return

    # Primary: confirmed seNorge SWE CSV
    swe_csv = src_dir / "senorge_swe" / "norway_swe_daily_1957_2025.csv"

    frames: list[pd.DataFrame] = []

    if swe_csv.exists():
        try:
            df = pd.read_csv(swe_csv, encoding="utf-8")
            df = _normalize_df(df)
            date_col = next((c for c in df.columns if "date" in c or "time" in c), None)
            if date_col:
                df[date_col] = pd.to_datetime(df[date_col], utc=True, errors="coerce")
                df = df.dropna(subset=[date_col]).set_index(date_col).sort_index()
            frames.append(df)
            LOG.info("NGI seNorge SWE: %d rows from %s", len(df), swe_csv.name)
        except Exception as exc:
            LOG.warning("Error reading seNorge SWE CSV: %s", exc)
    else:
        LOG.warning("Expected seNorge file not found: %s", swe_csv)

    # Fallback: any other CSV/JSON files in the directory tree
    for fp in sorted(src_dir.rglob("*.csv")):
        if fp == swe_csv:
            continue
        try:
            df = pd.read_csv(fp, low_memory=False, encoding="utf-8", encoding_errors="replace")
            df = _normalize_df(df)
            df = _parse_date_col(df)
            frames.append(df)
        except Exception as exc:
            LOG.warning("Error reading %s: %s", fp.name, exc)

    for fp in sorted(src_dir.rglob("*.json")):
        try:
            import json as _json
            data = _json.loads(fp.read_text(encoding="utf-8"))
            df = pd.DataFrame(data) if isinstance(data, list) else pd.DataFrame([data])
            df = _normalize_df(df)
            df = _parse_date_col(df)
            frames.append(df)
        except Exception as exc:
            LOG.warning("Error reading %s: %s", fp.name, exc)

    if not frames:
        LOG.warning("No NGI Norway data loaded")
        return

    # Only concat frames that were successfully parsed to a DatetimeIndex
    dt_frames = [f for f in frames if isinstance(f.index, pd.DatetimeIndex)]
    if not dt_frames:
        LOG.warning("No NGI Norway frames have datetime index — skipping")
        return
    combined = pd.concat(dt_frames, axis=0, sort=False) if len(dt_frames) > 1 else dt_frames[0]

    try:
        time_range = f"{combined.index.min()} / {combined.index.max()}"
    except Exception:
        time_range = "unknown"

    meta = {
        "title": "Norwegian seNorge SWE + NGI Avalanche Data",
        "source": "Norwegian Water Resources and Energy Directorate (NVE) / seNorge",
        "references": "https://senorge.no/",
        "time_range": time_range,
        "units": "SWE mm",
    }
    save_parquet(combined, out, meta)
    register_output(MANIFEST, "norway_swe", out, False, meta)


# ---------------------------------------------------------------------------
# CAIC Colorado
# ---------------------------------------------------------------------------
def process_caic() -> None:
    out = CRYO_OUT / "caic_accidents.parquet"
    if out.exists():
        LOG.info("SKIP caic_accidents.parquet")
        return

    src_dir = DATA_ROOT / "cryosphere" / "caic"
    if not src_dir.exists():
        LOG.warning("Missing: %s", src_dir)
        return

    frames = []
    for fp in sorted(src_dir.iterdir()):
        if fp.is_dir():
            continue
        try:
            if fp.suffix.lower() == ".csv":
                df = pd.read_csv(fp, low_memory=False, encoding="utf-8", encoding_errors="replace")
                df = _normalize_df(df)
                frames.append(df)
        except Exception as exc:
            LOG.warning("Error reading %s: %s", fp.name, exc)

    if not frames:
        LOG.warning("No CAIC files found")
        return

    combined = pd.concat(frames, ignore_index=True)
    combined = _parse_date_col(combined)

    meta = {
        "title": "Colorado Avalanche Information Center (CAIC) Accident Reports",
        "source": "CAIC",
        "references": "https://avalanche.state.co.us/",
    }
    save_parquet(combined, out, meta)
    register_output(MANIFEST, "caic_accidents", out, False, meta)


# ---------------------------------------------------------------------------
# IMS Snow Cover (4km gzipped ASCII)
#
# Confirmed path: data/cryosphere/ims_snow/{YYYY}/ims{YYYY}{DDD}_00UTC_4km_v1.2.asc.gz
# File format:
#   Line 1 : "Julian day of IMS data log: YYYYDDD"
#   Lines 2+: 6144 lines × 6144 space-separated single-digit integers per line
#   Values  : 1=open ocean, 2=land/no-snow, 3=sea-ice, 4=snow
# Scale    : 4 km LAEA NH centred on North Pole, ~430 KB/file × 4622 files = ~2 GB
# ---------------------------------------------------------------------------
def _ims_alps_slices() -> tuple[slice, slice]:
    """
    Compute row/col slices for Swiss Alps (44-48°N, 5-11°E) in the IMS 4km LAEA grid.
    Uses pyproj when available; falls back to hardcoded approximate indices.

    IMS LAEA grid convention:
      origin = upper-left corner  (-IMS_HALF_EXT, +IMS_HALF_EXT) in (X, Y)
      row increases → southward (−Y direction)
      col increases → eastward  (+X direction)
    """
    try:
        from pyproj import Transformer, CRS
        laea = CRS.from_proj4(
            "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 "
            "+a=6371228 +b=6371228 +units=m +no_defs"
        )
        geo = CRS.from_epsg(4326)
        tr  = Transformer.from_crs(geo, laea, always_xy=True)

        def _grid(lon: float, lat: float) -> tuple[int, int]:
            x, y = tr.transform(lon, lat)
            row = int((IMS_HALF_EXT - y) / IMS_CELL_M)
            col = int((x + IMS_HALF_EXT) / IMS_CELL_M)
            return (max(0, min(IMS_GRID_SIZE - 1, row)),
                    max(0, min(IMS_GRID_SIZE - 1, col)))

        # NW corner (lon=5°E, lat=48°N) → smallest row, smallest col
        r0, c0 = _grid(5.0, 48.0)
        # SE corner (lon=11°E, lat=44°N) → largest row, largest col
        r1, c1 = _grid(11.0, 44.0)
        LOG.debug("pyproj Alps grid: rows %d-%d  cols %d-%d", r0, r1, c0, c1)
        return slice(r0, r1 + 1), slice(c0, c1 + 1)

    except ImportError:
        LOG.debug("pyproj not available — using hardcoded Alps IMS indices")
    except Exception as exc:
        LOG.debug("pyproj projection error: %s — using hardcoded indices", exc)

    return IMS_ALPS_ROW_FALLBACK, IMS_ALPS_COL_FALLBACK


# Compute once at module load (pyproj optional)
_ALPS_ROW_SLICE, _ALPS_COL_SLICE = _ims_alps_slices()


def _parse_ims_gz(fp: Path) -> tuple:
    """
    Stream-parse a single IMS 4km .asc.gz file.
    Returns (date | None, nh_snow_fraction, alps_snow_fraction).

    Streams row-by-row to avoid loading the full 6144×6144 array into memory.
    """
    date: pd.Timestamp | None = None

    # Parse date from filename: ims{YYYY}{DDD}_00UTC_4km_v1.2.asc.gz
    m = re.search(r"ims(\d{4})(\d{3})_", fp.name)
    if m:
        try:
            date = pd.Timestamp(
                datetime.strptime(m.group(1) + m.group(2), "%Y%j")
            ).tz_localize("UTC")
        except ValueError:
            pass

    nh_snow = nh_land = 0
    alps_snow = alps_land = 0
    row_idx = 0
    header_seen = False

    try:
        with gzip.open(fp, "rt", encoding="ascii", errors="replace") as fh:
            for raw_line in fh:
                line = raw_line.strip()

                # Header line — also extract date as fallback if filename parse failed
                if not header_seen:
                    if "Julian" in raw_line or re.search(r"\d{7}", raw_line):
                        header_seen = True
                        if date is None:
                            dm = re.search(r"(\d{7})", raw_line)
                            if dm:
                                try:
                                    date = pd.Timestamp(
                                        datetime.strptime(dm.group(1), "%Y%j")
                                    ).tz_localize("UTC")
                                except ValueError:
                                    pass
                    continue

                if not line:
                    continue

                # Parse row: space-separated single-digit integers
                try:
                    row = np.fromstring(line, dtype=np.int8, sep=" ")
                    if row.size < IMS_GRID_SIZE // 2:
                        # Possibly packed (no spaces): "12342312…"
                        row = np.frombuffer(line.encode(), dtype=np.uint8) - ord("0")
                except Exception:
                    continue

                if row.size < 10:
                    continue

                # NH stats (land = {2, 4})
                nh_snow += int((row == IMS_SNOW_VALUE).sum())
                nh_land += int(np.isin(row, IMS_LAND_VALUES).sum())

                # Alps sub-region
                if _ALPS_ROW_SLICE.start <= row_idx < _ALPS_ROW_SLICE.stop:
                    alps_sub = row[_ALPS_COL_SLICE]
                    alps_snow += int((alps_sub == IMS_SNOW_VALUE).sum())
                    alps_land += int(np.isin(alps_sub, IMS_LAND_VALUES).sum())

                row_idx += 1
                if row_idx >= IMS_GRID_SIZE:
                    break

    except Exception as exc:
        LOG.warning("Stream error in %s: %s", fp.name, exc)

    if row_idx < IMS_GRID_SIZE // 2:
        return date, np.nan, np.nan

    nh_frac   = float(nh_snow   / nh_land)   if nh_land   > 0 else np.nan
    alps_frac = float(alps_snow / alps_land) if alps_land > 0 else np.nan
    return date, nh_frac, alps_frac


def process_ims_snow() -> None:
    out = CRYO_OUT / "ims_snow_daily.parquet"
    if out.exists():
        LOG.info("SKIP ims_snow_daily.parquet")
        return

    if not IMS_SNOW_DIR.exists():
        LOG.warning("Missing: %s", IMS_SNOW_DIR)
        return

    # Confirmed pattern: {YYYY}/ims{YYYY}{DDD}_00UTC_4km_v1.2.asc.gz
    gz_files = sorted(IMS_SNOW_DIR.rglob("ims*4km*.asc.gz"))
    if not gz_files:
        gz_files = sorted(IMS_SNOW_DIR.rglob("*.gz"))   # broad fallback
    if not gz_files:
        LOG.warning("No IMS gz files found under %s", IMS_SNOW_DIR)
        return

    LOG.info("Processing %d IMS 4km snow files … (Alps slice rows=%s cols=%s)",
             len(gz_files), _ALPS_ROW_SLICE, _ALPS_COL_SLICE)

    records: list[dict] = []
    processed_raw: list[Path] = []

    for fp in gz_files:
        try:
            date, nh_frac, alps_frac = _parse_ims_gz(fp)
            if date is not None:
                records.append({
                    "date": date,
                    "nh_snow_fraction":   nh_frac,
                    "alps_snow_fraction": alps_frac,
                })
            processed_raw.append(fp)
        except Exception as exc:
            LOG.warning("Error reading %s: %s", fp.name, exc)

    if not records:
        LOG.warning("No IMS records parsed")
        return

    df = pd.DataFrame(records).set_index("date").sort_index()
    meta = {
        "title": "IMS 4km Daily NH Snow Cover + Swiss Alps Fraction",
        "source": "NOAA/NESDIS IMS 4km",
        "references": "https://doi.org/10.7289/V5W9945Z",
        "time_range": f"{df.index.min()} / {df.index.max()}",
        "units": "fractions 0-1 (snow pixels / land pixels)",
    }
    save_parquet(df, out, meta)

    if out.exists() and out.stat().st_size > 0:
        safe_delete(processed_raw)
        register_output(MANIFEST, "ims_snow_daily", out, True, meta)
    else:
        LOG.error("IMS output missing/empty — raw NOT deleted")


def main() -> None:
    setup_logging()
    LOG.info("=== 10_process_cryosphere.py | disk free=%.1f GB ===", disk_free_gb())
    CRYO_OUT.mkdir(parents=True, exist_ok=True)

    process_slf_avalanche()
    process_snotel()
    process_ngi_norway()
    process_caic()

    LOG.info("Before IMS Snow | disk free=%.1f GB", disk_free_gb())
    process_ims_snow()

    LOG.info("=== 10 complete | disk free=%.1f GB ===", disk_free_gb())


if __name__ == "__main__":
    main()
