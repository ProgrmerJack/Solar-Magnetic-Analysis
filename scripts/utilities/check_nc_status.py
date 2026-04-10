"""Check all processed NC files for corruption."""
from pathlib import Path

for nc in sorted(Path("data/processed").rglob("*.nc")):
    d = nc.read_bytes()[:8]
    if d[1:4] == b"HDF":
        status = "HDF5-OK"
    elif all(b == 0 for b in d):
        status = "ZEROED "
    else:
        status = "OTHER  "
    size_mb = nc.stat().st_size / 1e6
    rel = str(nc.relative_to("data/processed"))
    print(f"{status}  {size_mb:8.1f} MB  {rel}")
