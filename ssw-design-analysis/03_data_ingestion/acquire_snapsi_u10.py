#!/usr/bin/env python3
"""
acquire_snapsi_u10.py
=====================
Zonal-mean zonal wind at 10 hPa on the grid row nearest 60N, for every member of
SNAPSI's `free` arm at the four NH initialisations: the stratospheric variable
that classifies members as developing the SSW or not (natural experiment F,
07_physical_decomposition/snapsi_free_natural_experiment.py, registered 1d93ff6
before this script first ran).

WHY A SEPARATE SCRIPT  Searched first (grep "ua_6hrPt", "10 hPa", "u10"): the SNAPSI
  acquisitions read psl, tas and 100 hPa zg; none reads u. The transport is that of
  acquire_snapsi_zg.py (OPeNDAP subset; byte ranges for Meteo-France and ECCC via
  its generalised read_zg100_ranges), so only the variable, level and row change.

GATE (physical, per centre, before any member is kept): the ensemble mean of the
  first-day mean u must lie within 5 m/s of ERA5's daily mean u(10 hPa, 60N) on the
  initialisation date (era5_u10_60N_daily_cds.parquet). A wrong level, row, unit or
  sign fails it. Each member must also be finite and within -80..120 m/s.
Output cache: _snapsi_u10_free/u10_<centre>_free_<init>_<member>.parquet
  columns centre, init, member, lead_days (rebased like psl), u10_60N, lat_used
"""
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import acquire_snapsi_surface as S                  # noqa: E402
import acquire_snapsi_zg as Z                       # noqa: E402

CACHE = HERE / "_snapsi_u10_free"
NH_INITS = ["s20180125", "s20180208", "s20181213", "s20190108"]
CENTRES = ["CCCma", "CNR-ISAC", "ECCC", "ECMWF", "KMA", "Meteo-France", "NCAR", "NOAA-GFDL", "SNU", "UKMO"]
TARGET_PA = 1000.0
ERA5_U = HERE / "era5_u10_60N_daily_cds.parquet"
WORKERS = int(os.environ.get("SNAPSI_WORKERS", 8))


def nearest_60(lat):
    return [int(np.argmin(np.abs(np.asarray(lat) - 60.0)))]


def ua_path(row, tok):
    p = Z.zg_path(row, tok)
    return p.replace("/zg/", "/ua/").replace("zg_6hrPt", "ua_6hrPt")


def read_member(row, tok, sess):
    path = ua_path(row, tok)
    if row.centre in Z.RANGE_CENTRES:
        da, tv, tu = Z.read_zg100_ranges(f"{Z.FILE_ROOT}{path}", tok, var="ua", target_pa=TARGET_PA,
                                         row_sel=nearest_60)
        lat_used = float(da["lat"].values[0])
        u = da.isel(lat=0).mean("lon").values
        lead = Z.lead_days(tv, tu)
    else:
        with xr.open_dataset(f"{Z.DAP_ROOT}/{path.lstrip('/')}", engine="pydap", session=sess,
                             decode_times=False) as rem:
            lat = rem["lat"].values
            j = nearest_60(lat)[0]
            da = rem["ua"].sel(plev=TARGET_PA).isel(lat=j)
            u = Z.load_blocked(da).mean("lon").values
            lat_used = float(lat[j])
            lead = Z.lead_days(rem["time"].values, rem["time"].attrs.get("units"))
    if not (np.isfinite(u).all() and (u > -80).all() and (u < 120).all()):
        raise ValueError(f"implausible u10: {np.nanmin(u):.1f}..{np.nanmax(u):.1f}")
    return pd.DataFrame({"centre": row.centre, "init": str(row.start_date), "member": row.member,
                         "lead_days": lead + S.origin_days(row.centre, "free", str(row.start_date)),
                         "u10_60N": u, "lat_used": lat_used})


def fetch(row, tok, sess):
    cf = CACHE / f"u10_{row.centre}_free_{row.start_date}_{row.member}.parquet"
    if cf.exists():
        return pd.read_parquet(cf)
    last = None
    for attempt in range(4):
        try:
            df = read_member(row, tok, sess)
            tmp = cf.with_suffix(f".parquet.{os.getpid()}.tmp")
            df.to_parquet(tmp); tmp.replace(cf)
            return df
        except S.MissingTimeOrigin:
            raise
        except Exception as e:                      # noqa: BLE001 -- retried, then reported
            last = e
            time.sleep(3 * (attempt + 1))
    print(f"  FAILED {row.centre}/{row.start_date}/{row.member}: {type(last).__name__}: {str(last)[:150]}", flush=True)
    return None


def gate(centre, init, frames):
    era = pd.read_parquet(ERA5_U).set_index("date")["u10_60N"]; era.index = pd.to_datetime(era.index)
    obs = float(era.loc[pd.Timestamp(init[1:])])
    d0 = [f[(f.lead_days >= 0) & (f.lead_days < 1)].u10_60N.mean() for f in frames]
    m = float(np.mean(d0))
    return {"centre": centre, "init": init, "ensemble_day0": round(m, 2), "era5": round(obs, 2),
            "passes": bool(abs(m - obs) <= 5.0), "lat_used": float(frames[0].lat_used.iloc[0])}


def main():
    CACHE.mkdir(exist_ok=True)
    argv = sys.argv[1:]
    centres = argv[argv.index("--centres") + 1].split(",") if "--centres" in argv else CENTRES
    tok = S.token(); sess = Z.session(tok)
    d = pd.read_csv(S.MANIFEST)
    q = d[(d.variable == "psl") & (d.experiment == "free") & d.start_date.isin(NH_INITS) & d.centre.isin(centres)]
    gates = []
    for (c, i), g in q.groupby(["centre", "start_date"]):
        rows = list(g.itertuples())
        first = [x for x in (fetch(r, tok, sess) for r in rows[:3]) if x is not None]
        if not first:
            gates.append({"centre": c, "init": i, "passes": False, "reason": "first members failed"}); continue
        gt = gate(c, i, first); gates.append(gt)
        print(f"gate {c} {i}: day0 {gt['ensemble_day0']} vs ERA5 {gt['era5']} m/s, lat {gt['lat_used']:.2f} -> "
              f"{'PASS' if gt['passes'] else 'FAIL'}", flush=True)
        if not gt["passes"]:
            for r in rows[:3]:
                (CACHE / f"u10_{r.centre}_free_{r.start_date}_{r.member}.parquet").unlink(missing_ok=True)
            continue
        with ThreadPoolExecutor(WORKERS) as ex:
            got = list(ex.map(lambda r: fetch(r, tok, sess), rows[3:]))
        print(f"  {c} {i}: {sum(x is not None for x in got) + len(first)}/{len(rows)} members", flush=True)
    pd.DataFrame(gates).to_csv(HERE / "snapsi_u10_free_gate.csv", index=False, lineterminator="\n")
    print("Saved gate table -> snapsi_u10_free_gate.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
