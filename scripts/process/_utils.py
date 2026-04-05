"""Shared utilities for the Solar-Magnetic-Avalanche processing pipeline."""
import json
import logging
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import xarray as xr

DATA_ROOT = Path(__file__).resolve().parents[2] / "data"
PROCESSED_ROOT = DATA_ROOT / "processed"

LOG = logging.getLogger("pipeline")


def setup_logging(level: int = logging.INFO) -> None:
    """Configure pipeline-wide logging to console and rotating file."""
    log_dir = Path(__file__).resolve().parents[2] / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "pipeline.log"

    fmt = "%(asctime)s  %(levelname)-8s  %(name)s  %(message)s"
    datefmt = "%Y-%m-%dT%H:%M:%SZ"

    logging.basicConfig(
        level=level,
        format=fmt,
        datefmt=datefmt,
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_file, encoding="utf-8"),
        ],
    )
    logging.Formatter.converter = lambda *_: datetime.now(timezone.utc).timetuple()


def save_parquet(df: pd.DataFrame, path: Path, metadata: dict[str, str]) -> None:
    """Write *df* to *path* as Parquet (snappy) with *metadata* embedded in schema."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    path.parent.mkdir(parents=True, exist_ok=True)

    # Coerce mixed-type object columns to numeric where possible (keeps strings that can't convert)
    for col in df.select_dtypes(include="object").columns:
        converted = pd.to_numeric(df[col], errors="coerce")
        if converted.notna().sum() > 0:
            df = df.copy()
            df[col] = converted

    # Drop columns with names starting with 'Unnamed' (CSV index artifacts)
    unnamed = [c for c in df.columns if str(c).startswith("Unnamed")]
    if unnamed:
        df = df.drop(columns=unnamed)

    # Encode metadata values as bytes
    encoded = {k.encode(): str(v).encode() for k, v in metadata.items()}

    table = pa.Table.from_pandas(df)
    existing_meta = table.schema.metadata or {}
    merged = {**existing_meta, **encoded}
    table = table.replace_schema_metadata(merged)

    pq.write_table(table, path, compression="snappy")
    size_mb = path.stat().st_size / 1_048_576
    LOG.info("Wrote %s  rows=%d  size=%.2f MB", path.name, len(df), size_mb)


def save_netcdf4(
    ds: xr.Dataset,
    path: Path,
    metadata: dict[str, str],
) -> None:
    """Write *ds* to CF-1.8 NetCDF4 with zlib/deflate=4 on all data variables.

    Uses a temp-file + atomic rename to prevent partial/zeroed writes, and
    verifies the HDF5 magic bytes before replacing the destination.
    """
    import tempfile

    path.parent.mkdir(parents=True, exist_ok=True)

    global_attrs = {
        "Conventions": "CF-1.8",
        "title": metadata.get("title", ""),
        "source": metadata.get("source", ""),
        "references": metadata.get("references", ""),
        "processing_date": datetime.now(timezone.utc).isoformat(),
    }
    global_attrs.update({k: v for k, v in metadata.items() if k not in global_attrs})
    ds = ds.assign_attrs(global_attrs)

    encoding = {
        var: {"zlib": True, "complevel": 4}
        for var in ds.data_vars
    }

    # Write to a temp file first; only replace destination if valid
    fd, tmp_path = tempfile.mkstemp(suffix=".nc", dir=path.parent)
    os.close(fd)
    tmp = Path(tmp_path)
    try:
        ds.to_netcdf(str(tmp), format="NETCDF4", encoding=encoding)
        # Verify HDF5 magic bytes before committing
        magic = tmp.read_bytes()[:8]
        if magic[1:4] != b"HDF":
            raise RuntimeError(
                f"save_netcdf4: temp file has bad magic {magic.hex()!r} — aborting"
            )
        # Atomic replace
        shutil.move(str(tmp), str(path))
        size_mb = path.stat().st_size / 1_048_576
        LOG.info("Wrote %s  vars=%d  size=%.2f MB", path.name, len(ds.data_vars), size_mb)
    except Exception:
        if tmp.exists():
            tmp.unlink(missing_ok=True)
        raise


def register_output(
    manifest_path: Path,
    key: str,
    out_path: Path,
    raw_deleted: bool,
    meta: dict[str, Any],
) -> None:
    """Append an entry to the pipeline manifest JSON."""
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    entries: list[dict] = []
    if manifest_path.exists():
        try:
            entries = json.loads(manifest_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            entries = []

    entry = {
        "key": key,
        "output": str(out_path),
        "raw_deleted": raw_deleted,
        "processed_at": datetime.now(timezone.utc).isoformat(),
        "size_bytes": out_path.stat().st_size if out_path.exists() else 0,
        **meta,
    }
    # Replace existing entry with same key
    entries = [e for e in entries if e.get("key") != key]
    entries.append(entry)

    manifest_path.write_text(
        json.dumps(entries, indent=2, default=str), encoding="utf-8"
    )
    LOG.info("Manifest updated: key=%s", key)


def safe_delete(path_list: list[Path]) -> int:
    """Delete files in *path_list*, log count, return number deleted."""
    deleted = 0
    for p in path_list:
        try:
            if p.is_file():
                p.unlink()
                deleted += 1
            elif p.is_dir():
                shutil.rmtree(p)
                deleted += 1
        except OSError as exc:
            LOG.warning("Could not delete %s: %s", p, exc)
    LOG.info("Deleted %d file(s)", deleted)
    return deleted


def disk_free_gb(path: Path | None = None) -> float:
    """Return free disk space in GB for the filesystem containing *path*."""
    check = str(path or DATA_ROOT.anchor)
    usage = shutil.disk_usage(check)
    return usage.free / 1_073_741_824
