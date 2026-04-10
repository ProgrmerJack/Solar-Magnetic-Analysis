"""Shared utilities for analysis scripts."""
import logging
from pathlib import Path
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

REPO_ROOT = Path(__file__).resolve().parents[2]
PROCESSED = REPO_ROOT / "data" / "processed"
RESULTS = REPO_ROOT / "data" / "results"
FIGURES = REPO_ROOT / "data" / "figures"

RESULTS.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)

LOG = logging.getLogger("analysis")


def save_analysis(df: pd.DataFrame, path: Path, meta: dict | None = None):
    """Save analysis DataFrame to Parquet with metadata."""
    path.parent.mkdir(parents=True, exist_ok=True)
    meta = meta or {}
    table = pa.Table.from_pandas(df)
    encoded = {k.encode(): str(v).encode() for k, v in meta.items()}
    existing = table.schema.metadata or {}
    table = table.replace_schema_metadata({**existing, **encoded})
    pq.write_table(table, path, compression="snappy")
    LOG.info("Saved %s: %d rows, %.1f MB", path.name, len(df),
             path.stat().st_size / 1e6)


def load_panel(winter_only: bool = True) -> pd.DataFrame:
    """Load the analysis panel, optionally filtering to winter months."""
    panel = pd.read_parquet(PROCESSED / "analysis_panel.parquet")
    if winter_only:
        panel = panel[panel["is_winter"] == 1].copy()
    return panel
