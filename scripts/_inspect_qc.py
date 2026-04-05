"""Inspect failing parquet files to understand their structure."""
import pandas as pd
from pathlib import Path

PROC = Path("data/processed")

files = [
    "cryosphere/snotel_daily.parquet",
    "solar/sdo_hmi_sharp.parquet",
    "cryosphere/slf_simulated-avalanche-probl_1999-2017_Avalanche_problem_types_from_Crocus..parquet",
    "cryosphere/slf_simulated-avalanche-probl_1999-2017_Avalanche_problem_types_from_SNOWPACK..parquet",
    "cryosphere/slf_wet_model.parquet",
    "geomagnetic/dst_index.parquet",
    "atmospheric/mls_hno3_polar.parquet",
]

for rel in files:
    df = pd.read_parquet(PROC / rel)
    date_cols = [c for c in df.columns if any(x in c.lower() for x in ("date","time","year","day","month"))]
    print(f"{rel}")
    print(f"  index type : {type(df.index).__name__}, tz={getattr(df.index,'tz','n/a')}")
    print(f"  shape      : {df.shape}")
    print(f"  columns    : {list(df.columns)[:8]}")
    print(f"  date_cols  : {date_cols}")
    if date_cols:
        print(f"  date sample: {df[date_cols[0]].head(3).tolist()}")
    # Check first numeric column
    num_cols = df.select_dtypes(include="number").columns
    if len(num_cols) > 0:
        col = num_cols[0]
        vals = df[col].dropna()
        if len(vals) > 0:
            print(f"  {col} sample: {vals.head(3).tolist()}, NaN={df[col].isna().mean():.0%}")
    print()
