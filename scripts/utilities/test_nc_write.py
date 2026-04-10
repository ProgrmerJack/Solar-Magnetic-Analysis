"""Test: write same-shape MLS array to verify path works correctly."""
import xarray as xr
import numpy as np
from pathlib import Path

p = Path("data/processed/atmospheric/mls_hno3_gridded.nc")
p.unlink(missing_ok=True)

# Same shape as real MLS HNO3 gridded: 8036 x 37 x 45
arr = np.random.rand(8036, 37, 45).astype("float32")
arr[arr > 0.9] = np.nan  # add some NaNs like real data

ds = xr.Dataset(
    {"hno3": (["time", "lev", "lat"], arr)},
    coords={
        "time": np.arange(8036),
        "lev": np.linspace(0.001, 1000, 37),
        "lat": np.linspace(-88, 88, 45),
    },
)
ds.attrs["Conventions"] = "CF-1.8"
ds.to_netcdf(p, format="NETCDF4", encoding={"hno3": {"zlib": True, "complevel": 4}})

d = p.read_bytes()[:8]
valid = d[1:4] == b"HDF"
print(f"size={p.stat().st_size/1e6:.1f}MB  magic={d.hex()}  valid={valid}")
p.unlink(missing_ok=True)
print("TEST PASSED" if valid else "TEST FAILED - file zeroed!")
