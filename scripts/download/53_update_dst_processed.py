"""Update data/processed/geomagnetic/dst_index.parquet with correct OMNI2 data."""
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path

src = pd.read_parquet("data/geomagnetic/omni_hourly_dst.parquet")
out = src.rename(columns={"kp_omni": "kp", "ae": "ae_index"})

table = pa.Table.from_pandas(out)
meta = {
    b"source": b"OMNI2 hourly (NASA GSFC SPDF)",
    b"reference": b"https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/",
    b"column_dst": b"Dst index (nT); fill masked >=99999; storms: Dst<-50 moderate, Dst<-100 intense",
    b"column_kp": b"Kp*10 planetary geomagnetic index (divide by 10 for standard Kp 0-9 scale)",
    b"column_ae_index": b"AE auroral electrojet index (nT), fill masked >=9999",
    b"temporal_coverage": b"1963-01-01 to 2025-12-31 UTC",
    b"cadence": b"1 hour",
    b"dst_min_nT": str(out["dst"].min()).encode(),
    b"dst_max_nT": str(out["dst"].max()).encode(),
    b"cf_convention": b"CF-1.8 compatible",
    b"created_by": b"scripts/download/52_download_omni_dst.py + 53_update_dst_processed.py",
}
existing_meta = table.schema.metadata or {}
existing_meta.update(meta)
table = table.replace_schema_metadata(existing_meta)

out_path = Path("data/processed/geomagnetic/dst_index.parquet")
pq.write_table(table, out_path, compression="snappy")

size_mb = out_path.stat().st_size / 1e6
n_valid = out["dst"].notna().sum()
print(f"Written: {out_path} ({size_mb:.2f} MB)")
print(f"Rows: {len(out)}, Dst valid: {n_valid} ({100*n_valid/len(out):.1f}%)")
print(f"Dst range: [{out['dst'].min():.0f}, {out['dst'].max():.0f}] nT")
print(f"Columns: {out.columns.tolist()}")
print(f"Index tz: {out.index.tz}")
