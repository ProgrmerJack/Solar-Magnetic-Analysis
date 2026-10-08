#!/usr/bin/env python3
"""
acquire_snapsi_regime_precip.py
===============================
Inputs of 07_physical_decomposition/snapsi_regime_precip.py (registered 6c99ebc before
any of these values was fetched): for every SNAPSI nudged and control member of the four
Northern Hemisphere initialisations and every centre with both arms,
  pressure  daily-mean psl (6hrPt) on post-onset days 1-25, coarse-grained to 5 deg
            (cell means of the native points within +-2.5 deg) over 30-80N, 80W-40E;
  rainfall  daily-mean pr (6hr) over IBERIA 36-44N 10-1W, SCANDINAVIA 55-70N 4-20E,
            BRITISH_ISLES 49-61N 10W-1E (cos-latitude weights; kg m-2 s-1 x 86400 = mm/day).
Days are counted from the event's onset (snapsi_selection_test.ONSET) on the model's own
calendar (cftime), using the time-bounds midpoint when the file has bounds; a day is kept
only if all its 6-hourly samples are present, and a member only if all 25 days are.
pr files are found by listing <member>/6hr/pr/<grid>/<version>/; when several versions
exist the latest is taken and recorded. Whole files are fetched with the CEDA token
(acquire_snapsi_surface.token) and discarded after reduction; per-member cache in
_snapsi_regime_precip/.
AMENDED 2026-10-08 (before any rainfall outcome was computed): the 2019 initialisations
cover only post-onset days 1-24 (s20181213) and 6-25 (s20190108), so completeness is
required on days 8-24 and the analysis uses those days; days 1-25 are stored where present.
Bounds are used only when consistent with the time stamps (NRL's decode to wrong dates);
a day needs all samples finite; precipitation noise down to -0.05 mm/day is set to zero.
Assembly (--assemble, also run at the end): anomalies of the pressure maps from each
pair's control ensemble-mean field over days 8-24, as amended; writes
snapsi_regime_precip.parquet (rainfall) and snapsi_regime_psl.npz (maps, keys, lat).
"""
import io
import re
import sys
import time
import urllib.request
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "07_physical_decomposition"))
import acquire_snapsi_surface as S                    # noqa: E402  token, FILE_ROOT, resolve_archive_path

CACHE = HERE / "_snapsi_regime_precip"
OUT_RAIN = HERE / "snapsi_regime_precip.parquet"
OUT_FIELD = HERE / "snapsi_regime_psl.npz"
ONSET = {"s20180125": (2018, 2, 12), "s20180208": (2018, 2, 12), "s20181213": (2019, 1, 2), "s20190108": (2019, 1, 2)}
LATS = np.arange(30, 81, 5.0); LONS = np.arange(-80, 41, 5.0)
BOXES = {"IBERIA": (36, 44, -10, -1), "SCANDINAVIA": (55, 70, 4, 20), "BRITISH_ISLES": (49, 61, -10, 1)}
DAYS = range(1, 26)          # stored where available
NEED = range(8, 25)          # required complete: days 8-24 (amended 2026-10-08, see snapsi_regime_precip)


def get(url, tok, tries=4):
    last = None
    for a in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"Authorization": f"Bearer {tok}"}), timeout=300) as r:
                return r.read()
        except Exception as e:                       # noqa: BLE001
            last = e; time.sleep(10 * (a + 1))
    raise RuntimeError(f"{url}: {last}")


def pr_path(psl_archive_path, tok):
    """<...>/<member>/6hrPt/psl/<grid>/<ver>/<file> -> the pr file under <member>/6hr/pr/."""
    member_dir = psl_archive_path.split("/6hrPt/psl/")[0]
    base = f"{member_dir}/6hr/pr/"
    grids = re.findall(r'href="([a-z0-9]+)/"', get(S.FILE_ROOT + base, tok).decode())
    hits = []
    for gdir in grids:
        vers = sorted(set(re.findall(r'href="(v\d{8})/"', get(f"{S.FILE_ROOT}{base}{gdir}/", tok).decode())))
        for v in vers:
            files = re.findall(r'href="([^"/]+\.nc)"', get(f"{S.FILE_ROOT}{base}{gdir}/{v}/", tok).decode())
            hits += [(v, f"{base}{gdir}/{v}/{f}") for f in files if f.startswith("pr_")]
    if not hits:
        raise FileNotFoundError(f"no pr under {base}")
    latest = max(h[0] for h in hits)
    sel = [h for h in hits if h[0] == latest]
    if len(sel) != 1:
        raise FileNotFoundError(f"{base}: {len(sel)} pr files in {latest}")
    return sel[0][1], latest


def post_onset_days(ds, init):
    import cftime
    stamps = ds["time"].values
    cal = getattr(stamps[0], "calendar", "standard")
    on = cftime.datetime(*ONSET[init], calendar=cal)
    x = np.array([(v - on).total_seconds() / 86400.0 for v in stamps])
    bname = "time_bnds" if "time_bnds" in ds else "time_bounds" if "time_bounds" in ds else None
    if bname:                                  # bounds midpoint, only if consistent with the stamps (NRL's are not)
        try:
            b = ds[bname].values
            mid = np.array([((b0 + (b1 - b0) / 2) - on).total_seconds() / 86400.0 for b0, b1 in b])
            if np.all(np.abs(mid - x) <= 0.5):
                x = mid
        except Exception:                      # noqa: BLE001 -- undecodable bounds: keep the stamps
            pass
    return x


def daily(values, x):
    """values (n_time, ...) -> {day: mean} for days 1-25 with complete sampling."""
    d = np.floor(x).astype(int)
    per = np.median(np.bincount(d[(d >= 1) & (d <= 25)])) if ((d >= 1) & (d <= 25)).any() else 0
    out = {}
    for k in DAYS:
        m = d == k
        if m.sum() >= max(per, 1) and np.isfinite(values[m]).all():
            out[k] = values[m].mean(0)
    return out


def lon180(lon):
    return ((np.asarray(lon) + 180) % 360) - 180


def coarse(psl, lat, lon):
    lon = lon180(lon); out = np.full((psl.shape[0], len(LATS) * len(LONS)), np.nan); i = 0
    for la in LATS:
        mi = np.abs(lat - la) <= 2.5
        for lo in LONS:
            mj = np.abs(lon - lo) <= 2.5
            if not mi.any() or not mj.any():
                raise ValueError(f"empty 5-deg cell {la},{lo}")
            out[:, i] = psl[:, mi][:, :, mj].mean(axis=(1, 2)); i += 1
    return out


def boxes(pr, lat, lon):
    lon = lon180(lon); cols = {}
    for b, (la0, la1, lo0, lo1) in BOXES.items():
        mi = (lat >= la0) & (lat <= la1); mj = (lon >= lo0) & (lon <= lo1)
        if not mi.any() or not mj.any():
            raise ValueError(f"empty box {b}")
        w = np.cos(np.deg2rad(lat[mi]))[:, None] * np.ones(mj.sum())[None]
        v = (pr[:, mi][:, :, mj] * w).sum(axis=(1, 2)) / w.sum() * 86400.0
        if np.nanmin(v) < -0.05:
            raise ValueError(f"{b}: negative precipitation {np.nanmin(v)} mm/day")
        cols[b] = np.where(v < 0, 0.0, v)            # numerical noise (ECCC: -0.003 mm/day) set to zero
    return cols


def job(args):
    row, tok = args
    import xarray as xr
    CACHE.mkdir(exist_ok=True)
    cf = CACHE / f"{row['centre']}_{row['experiment']}_{row['start_date']}_{row['member']}.parquet"
    if cf.exists():
        return str(cf)
    try:
        psl_path = row["archive_path"] if row["verified"] else S.resolve_archive_path(row["archive_path"], tok)
        with xr.open_dataset(io.BytesIO(get(S.FILE_ROOT + psl_path + "?download=1", tok)), use_cftime=True) as ds:
            x = post_onset_days(ds, row["start_date"])
            P = daily(coarse(ds["psl"].values, ds["lat"].values, ds["lon"].values), x)
        prp, ver = pr_path(psl_path, tok)
        with xr.open_dataset(io.BytesIO(get(S.FILE_ROOT + prp + "?download=1", tok)), use_cftime=True) as ds:
            x = post_onset_days(ds, row["start_date"])
            bx = boxes(ds["pr"].values, ds["lat"].values, ds["lon"].values)
            R = {b: daily(v, x) for b, v in bx.items()}
        if any(k not in P for k in NEED) or any(k not in v for v in R.values() for k in NEED):
            return f"FAILED {cf.name}: incomplete days 8-24: psl {sorted(P)[:1]}..{sorted(P)[-1:]} pr {[len(v) for v in R.values()]}"
        rows = []
        for k in sorted(set(P) & set.intersection(*[set(v) for v in R.values()])):
            rr = {b: float(R[b][k]) for b in BOXES}
            if not all(0 <= v < 200 for v in rr.values()) or not (85000 < P[k].mean() < 110000):
                return f"FAILED {cf.name}: implausible day {k}"
            rows.append({"centre": row["centre"], "model": row["model"], "experiment": row["experiment"], "init": row["start_date"],
                         "member": row["member"], "day": k, "pr_version": ver, **rr,
                         **{f"f{i:03d}": float(v) for i, v in enumerate(P[k])}})
        pd.DataFrame(rows).to_parquet(cf)
        return str(cf)
    except Exception as e:                           # noqa: BLE001
        return f"FAILED {cf.name}: {type(e).__name__}: {str(e)[:200]}"


def assemble():
    d = pd.concat([pd.read_parquet(p) for p in sorted(CACHE.glob("*.parquet"))], ignore_index=True)
    d = d[d.day.isin(list(NEED))]                    # analysis window, days 8-24
    fcols = [c for c in d.columns if c.startswith("f") and c[1:].isdigit()]
    both = d.groupby(["centre", "init"]).experiment.nunique()
    d = d.set_index(["centre", "init"]).loc[both[both == 2].index].reset_index()
    rain = d[["centre", "model", "experiment", "init", "member", "day", *BOXES]]
    rain.to_parquet(OUT_RAIN)
    z = {"lat": np.repeat(LATS, len(LONS))}
    for c, gc in d.groupby("centre"):
        for arm in ("control", "nudged"):
            maps, keys = [], []
            for init, gi in gc.groupby("init"):
                ref = gi[gi.experiment == "control"][fcols].values.mean(0)
                ga = gi[gi.experiment == arm]
                maps.append(ga[fcols].values - ref)
                keys += [(init, m, int(dd)) for m, dd in zip(ga.member, ga.day)]
            z[f"{c}|{arm}|maps"] = np.concatenate(maps); z[f"{c}|{arm}|keys"] = np.array(keys, dtype=object)
    np.savez(OUT_FIELD, **z)
    print(rain.groupby(["centre", "init", "experiment"]).member.nunique().unstack().to_string())
    print(f"Saved -> {OUT_RAIN.name}, {OUT_FIELD.name}")


def main():
    if "--assemble" in sys.argv:
        return assemble()
    tok = S.token()
    m = pd.read_csv(HERE / "snapsi_manifest.csv")
    m = m[(m.variable == "psl") & m.start_date.isin(ONSET) & m.experiment.isin(["control", "nudged"])]
    pairs = m.groupby(["centre", "start_date"]).experiment.nunique()
    keep = pairs[pairs == 2].index
    m = m.set_index(["centre", "start_date"]).loc[keep].reset_index()
    if "--one" in sys.argv:
        m = m.head(1)
    jobs = [(r._asdict() if hasattr(r, "_asdict") else dict(r), tok) for r in m.to_dict("records")]
    print(f"{len(jobs)} members", flush=True)
    import multiprocessing as mp
    done = 0
    with ProcessPoolExecutor(8, mp_context=mp.get_context("fork")) as ex:   # netCDF from bytes: fork is safe
        for r in ex.map(job, jobs):
            done += 1
            if r.startswith("FAILED"):
                print(r, flush=True)
            if done % 100 == 0:
                print(f"  {done}/{len(jobs)}", flush=True)
    if "--one" not in sys.argv:
        assemble()
    return 0


if __name__ == "__main__":
    sys.exit(main())
