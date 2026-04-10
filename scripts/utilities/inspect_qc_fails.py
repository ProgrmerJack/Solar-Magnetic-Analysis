"""Inspect files that are failing QC checks."""
import pandas as pd
from pathlib import Path

base = Path("data/processed")

# --- Check the 2 SLF timezone FAILs ---
for name in [
    "slf_simulated-avalanche-probl_1999-2017_Avalanche_problem_types_from_Crocus..parquet",
    "slf_simulated-avalanche-probl_1999-2017_Avalanche_problem_types_from_SNOWPACK..parquet",
]:
    p = base / "cryosphere" / name
    df = pd.read_parquet(p)
    tz = getattr(df.index, "tz", None)
    print(f"{name[:70]}")
    print(f"  index type: {type(df.index).__name__}, dtype: {df.index.dtype}, tz: {tz}")
    print(f"  shape: {df.shape}, cols: {list(df.columns[:6])}")
    print(f"  sample index: {list(df.index[:3])}")
    print()

# --- Check SDO HMI SHARP ---
p2 = base / "solar" / "sdo_hmi_sharp.parquet"
df2 = pd.read_parquet(p2)
print("sdo_hmi_sharp.parquet:")
print(f"  index: {type(df2.index).__name__}, dtype: {df2.index.dtype}")
print(f"  shape: {df2.shape}")
print(f"  columns: {list(df2.columns[:15])}")
date_cols = [c for c in df2.columns if any(k in c.lower() for k in ["time", "date", "t_rec", "epoch"])]
print(f"  date-like cols: {date_cols}")
if date_cols:
    print(f"  sample {date_cols[0]}: {list(df2[date_cols[0]].head(3))}")
