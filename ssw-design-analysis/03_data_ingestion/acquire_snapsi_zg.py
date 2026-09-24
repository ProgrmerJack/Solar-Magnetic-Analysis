#!/usr/bin/env python3
"""
acquire_snapsi_zg.py
====================
SNAPSI geopotential height at 100 hPa, reduced to a polar-cap mean, fetched by
OPeNDAP subsetting rather than by downloading whole files.

WHY A SECOND ACQUISITION SCRIPT
  `acquire_snapsi_surface.py` downloads whole `psl` files and reduces them in
  memory. That cannot work for `zg`: it is 3-D on 34 pressure levels and ~34x
  the size of psl per file, putting `nudged`+`control` for all centres at
  roughly **5 TB** against 398 GB of disk. The transport has to differ, so the
  script does. Searched before writing this: nothing in the repo speaks CEDA
  OPeNDAP -- the existing dodsC users (`acquire_heatflux.py`,
  `extend_ncep_presatellite.py`, and the superseded `scripts/acquisition/*`) all
  target NOAA PSL, a different server with no authentication.

WHAT IT IS FOR
  Loeffel et al. (2026, WCD 7, 895-913) report r = 0.85 between the week-2
  100 hPa polar-cap GPH anomaly and the surface response over weeks 3-7, across
  18 events, and conclude that SSWs differ in their capacity to couple downward.
  That correlation is computed ACROSS EVENTS on ensemble means. In a SNAPSI
  `nudged` ensemble the event is identical in every member by construction, so
  the same relation can be measured BETWEEN MEMBERS, where no event-to-event
  difference exists to explain it.
  This script only acquires the field. The test is a separate analysis.

WHY 100 hPa IS NOT NUDGED, WHICH IS WHAT MAKES THE TEST POSSIBLE
  The SNAPSI protocol (Hitchcock et al., GMD 15, 5073, 2022) nudges only the
  ZONAL-MEAN temperature and zonal wind, at full strength above 50 hPa, tapering
  to no nudging below 90 hPa. Eddies are free at every level. So 100 hPa
  polar-cap zg is not directly constrained and members genuinely differ there --
  but that has to be MEASURED, not assumed, and the analysis gates on it.

TRANSPORT, AND THE TRAP IT IS GATED AGAINST
  `https://dap.ceda.ac.uk/thredds/dodsC/<archive path>` with the same bearer JWT
  the file download uses, opened through pydap with a requests.Session carrying
  the header. Subsetting to one level and one 30-degree cap is **152x smaller**:
  204.8 MB -> 1.35 MB in 6 s.

  NCEP OPeNDAP once returned exactly 0.0 for every day of a year through this
  same class of interface -- no error, no NaN (see the project's data-access
  notes). So every centre is gated on first use: its subset is compared against
  a full-file download of the same field, and the centre is refused unless they
  agree exactly. CanESM5 passing does not license the rest.

  Two dead ends, recorded so they are not retried: hand-parsing the `.dods`
  binary returns the coordinate MAPS rather than the array, because DAP2 Grid
  payloads append them after the data; and `.nc` responses return HTTP 400 here.

SECOND TRANSPORT: BYTE RANGES, FOR CENTRES OPeNDAP REFUSES
  CEDA's OPeNDAP returns 403 for Meteo-France zg with a token that opens every
  other centre, while the same file downloads whole (1.6 GB, HTTP 200) and
  honours `Range:` (HTTP 206). Its zg is HDF5, UNCOMPRESSED, chunked
  (1 time, 17 levels, 91 lat, 180 lon), with the vertical axis named `snap34`
  -- an ordinary pressure coordinate in Pa, not a different quantity. In an
  uncompressed chunk each level is one contiguous lat x lon slab at a
  computable offset, so the 100 hPa field is read as 4 byte ranges per time
  step (~47 MB per member instead of 1.6 GB). Chunk offsets come from the HDF5
  index via h5py; the ranges are fetched with fsspec. The same full-file gate
  applies: the reader is refused unless it matches the downloaded file exactly.

Output: snapsi_polarcap_zg100.parquet
        centre, model, experiment, init, member, lead_days,
        zg100_cap_N, zg100_cap_S, hemisphere, zg100_cap
"""
from __future__ import annotations

import io
import json
import os
import pathlib
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import requests
import xarray as xr

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import acquire_snapsi_surface as S                  # noqa: E402  (time origin)

HERE = pathlib.Path(__file__).resolve().parent
MANIFEST = HERE / "snapsi_manifest.csv"
OUT = HERE / "snapsi_polarcap_zg100.parquet"
# zg gets its OWN cache directory. It must not share _snapsi_reduced: the
# analysis scripts read that directory by filename and parse the centre from
# the first underscore-delimited field, so anything not psl corrupts them.
CACHE = HERE / "_snapsi_zg100"
GATE = HERE / "snapsi_zg_opendap_gate.json"
ENV = pathlib.Path(os.environ.get("CEDA_ENV", HERE.parents[1] / ".env"))

DAP_ROOT = "https://dap.ceda.ac.uk/thredds/dodsC"
FILE_ROOT = "https://dap.ceda.ac.uk"

NH_INITS = ["s20180125", "s20180208", "s20181213", "s20190108"]
SH_INITS = ["s20190829", "s20191001"]
ALL_INITS = NH_INITS + SH_INITS
EXPERIMENTS = ["nudged", "control"]
# NRL is omitted: its submission is corrupt (nudged-minus-control effect
# identically zero, ratio 0.004 in the duplicate guard) so every downstream
# analysis excludes it. Nothing is lost by not fetching it.
ALL_CENTRES = ["CCCma", "CNR-ISAC", "ECCC", "ECMWF", "KMA", "Meteo-France",
               "NCAR", "SNU", "UKMO"]

TARGET_PA = 10000.0      # 100 hPa, the Loeffel level
CAP_LAT = 60.0
WORKERS = 8              # gentler than the psl run: each request is a server-side subset
# Centres whose OPeNDAP endpoint returns 403 (measured 2026-09-23); read by
# byte range from the file service instead. See the module docstring.
RANGE_CENTRES = {"Meteo-France", "ECCC"}
# ECCC moved from OPeNDAP to byte range on 2026-09-24 (gzip chunks inflated
# locally; 0.8 -> ~20 members/min). A centre on the range route must ALSO pass
# verify_ranges_against_cache(), recorded in the gate file as "<centre>|byte_range":
# the current reader has to reproduce members already cached by a gated route.
# OPeNDAP subsets are requested in blocks of this many time steps. CNR-ISAC zg
# is gzip-compressed with ALL 34 levels in one chunk per time step, so the
# server decompresses every level to serve one; a whole-member request took
# >340 s under load and the server cut the response ("Response ended
# prematurely") on every CNR-ISAC s20190829 member on 2026-09-23. Smaller
# requests finish inside that limit. Values are unaffected -- the same
# constraint expression, split along time. Measured, since the per-centre gates
# of 2026-09-17 exercised whole loads only: CCCma blocked == whole (array_equal),
# and CNR-ISAC control s20180125 r1i10p1f1 refetched blocked equals its
# 2026-09-17 whole-load cache exactly (max|diff| 0.0 on both caps, 181 steps).
TIME_BLOCK = 30
# The gate's reference file is fetched in byte ranges of this size, each retried
# on its own, because whole-file reads of 1.1-2.2 GB died with IncompleteRead.
GATE_PIECE = 64 * 2 ** 20
GATE_WORKERS = 8
# Concurrent range requests per member. Unbounded, one member issued ~720 at
# once; with 8 workers plus OPeNDAP jobs that was thousands of simultaneous
# requests to one server, which answered with disconnects and throughput fell
# to 0.8 members/min (2026-09-24). Total in flight <= workers x RANGE_BATCH.
RANGE_BATCH = 16
GATE_DL = HERE / "_gate_download"


# ------------------------------------------------------------------- auth ---

def token() -> str:
    txt = ENV.read_text(encoding="utf8", errors="ignore")
    m = re.search(r"^\s*CEDA\s*=\s*['\"]?([A-Za-z0-9_\-\.]+)", txt, re.M)
    if not m:
        sys.exit(f"no CEDA token in {ENV}")
    tok = m.group(1)
    import base64
    pay = json.loads(base64.urlsafe_b64decode(
        tok.split(".")[1] + "=" * (-len(tok.split(".")[1]) % 4)))
    left = pay.get("exp", 0) - time.time()
    if left <= 0:
        sys.exit("CEDA token EXPIRED -- refresh it before running")
    print(f"CEDA token valid, {left / 3600:.1f} h remaining "
          f"(user {pay.get('preferred_username', '?')})")
    return tok


def session(tok: str) -> requests.Session:
    s = requests.Session()
    s.headers.update({"Authorization": f"Bearer {tok}"})
    return s


def zg_path(row, tok=None) -> str:
    """The zg twin of a psl archive path, in the SAME data version.

    Layout is .../<freq>/<var>/<grid>/<version>/<var>_<freq>_... so only the
    variable token changes. Verified to resolve for all 10 centres by HEAD.
    For a manifest row marked verified=False the psl path is first resolved
    against the archive (see acquire_snapsi_surface.resolve_archive_path), and
    zg is then required to sit in that same version: a member whose psl and zg
    live in different versions (CNR-ISAC nudged s20190108 r1i46p1f1) is refused
    rather than paired across data vintages.
    """
    psl = row.archive_path
    if tok is not None and not bool(getattr(row, "verified", True)):
        from acquire_snapsi_surface import resolve_archive_path
        psl = resolve_archive_path(psl, tok)
        zg = psl.replace("/psl/", "/zg/").replace("psl_6hrPt", "zg_6hrPt")
        if resolve_archive_path(zg, tok) != zg:
            raise FileNotFoundError(f"zg not in psl's version for {row.member}; refusing to mix")
        return zg
    return psl.replace("/psl/", "/zg/").replace("psl_6hrPt", "zg_6hrPt")


# ------------------------------------------------------------- reduction ---

def cap_means(da, lat):
    north = da.sel(lat=lat[lat >= CAP_LAT].values)
    south = da.sel(lat=lat[lat <= -CAP_LAT].values)
    n = north.weighted(np.cos(np.deg2rad(north.lat))).mean(dim=["lat", "lon"]).values
    s = south.weighted(np.cos(np.deg2rad(south.lat))).mean(dim=["lat", "lon"]).values
    return np.asarray(n, dtype=float), np.asarray(s, dtype=float)


_UNIT_TO_DAYS = {"day": 1.0, "days": 1.0,
                 "hour": 1 / 24, "hours": 1 / 24,
                 "minute": 1 / 1440, "minutes": 1 / 1440,
                 "second": 1 / 86400, "seconds": 1 / 86400}


def lead_days(tv, units=None):
    """Forecast lead in days, without reconciling the calendar.

    Absolute dates are never needed -- lead is a difference within one file --
    so no 365-vs-366 drift from the model calendars can creep in. CCCma runs a
    365_day calendar, which is exactly why decoding is avoided.

    Opened with decode_times=False, the axis arrives as raw numbers in the
    file's own units ("days since 1850-01-01" for CanESM5, 6-hourly, so values
    step by 0.25). The first version assumed cftime objects and called
    .total_seconds() on a float, failing every member. The units attribute is
    honoured rather than assumed, and an unrecognised one RAISES: silently
    treating hours as days would scale every lead by 24 and quietly misplace
    the analysis windows.
    """
    arr = np.asarray(tv)
    t0 = arr[0]
    if np.issubdtype(arr.dtype, np.datetime64):
        return (pd.to_datetime(arr) - pd.to_datetime(t0)).total_seconds() / 86400.0
    if np.issubdtype(arr.dtype, np.number):
        if not units:
            raise ValueError("numeric time axis with no units attribute")
        unit = str(units).strip().split()[0].lower()
        if unit not in _UNIT_TO_DAYS:
            raise ValueError(f"unhandled time unit {unit!r} in {units!r}")
        return (arr.astype(float) - float(t0)) * _UNIT_TO_DAYS[unit]
    # cftime objects
    return np.array([(x - t0).total_seconds() / 86400.0 for x in arr], dtype=float)


def vertical_dim(ds) -> str:
    """The pressure dimension of `zg`, whatever the file calls it.

    Most centres use `plev`; Meteo-France uses `snap34`. Identified by
    standard_name air_pressure and units Pa, never by name alone.
    """
    for d in ds["zg"].dims:
        if d in ds.variables:
            a = ds[d].attrs
            if a.get("standard_name") == "air_pressure" and a.get("units") == "Pa":
                return d
    raise KeyError(f"no air_pressure (Pa) dimension on zg: {ds['zg'].dims}")


def _s(v):
    return v.decode() if isinstance(v, (bytes, np.bytes_)) else str(v)


def read_zg100_ranges(url: str, tok: str):
    """100 hPa zg over BOTH polar caps (time, cap-rows, lon) from a chunked HDF5
    file, fetched by byte range. Returns (DataArray, time values, time units).

    Two layouts are readable, and nothing else:
      uncompressed  only the cap rows of the 100 hPa level are fetched -- within
                    a C-order chunk they are one contiguous run per (time, chunk)
                    (Meteo-France: ~1/3 of the whole-slab bytes)
      deflate only  the whole compressed chunks covering the cap rows at 100 hPa
                    are fetched in parallel and inflated locally (ECCC, whose
                    OPeNDAP server decompresses every level and ran at
                    0.8 members/min on 2026-09-24)
    Refused: any other filter, a non-little-endian-float32 dtype, a level axis
    without an exact 100 hPa entry, an unallocated or wrong-size chunk.
    Only rows with |lat| >= CAP_LAT are returned; cap_means needs nothing else.
    """
    import zlib

    import fsspec
    import h5py
    fs = fsspec.filesystem("http", client_kwargs={
        "headers": {"Authorization": f"Bearer {tok}"}})
    # 4 KB blocks for the METADATA walk. The chunk index is scattered through
    # the file in small nodes; at fsspec's 1 MB default each node pulled a whole
    # megabyte and reading the index took 406 s per member, against 11 s here.
    with fs.open(url, "rb", block_size=2 ** 12) as fh, h5py.File(fh, "r") as f:
        v = f["zg"]
        dims = [v.dims[i][0].name.lstrip("/") if len(v.dims[i]) else None
                for i in range(v.ndim)]
        if len(dims) != 4 or dims[0] != "time" \
                or dims[2] != "lat" or dims[3] != "lon":
            raise ValueError(f"unexpected zg dims {dims}")
        lev = f[dims[1]]
        # ECCC's plev carries long_name "pressure" and units Pa but no
        # standard_name; Meteo-France's snap34 carries standard_name. Either
        # identifies a pressure axis; units must be Pa in both cases.
        is_p = (_s(lev.attrs.get("standard_name", b"")) == "air_pressure"
                or _s(lev.attrs.get("long_name", b"")).lower() == "pressure")
        if not is_p or _s(lev.attrs.get("units", b"")) != "Pa":
            raise ValueError(f"{dims[1]} is not a pressure axis in Pa")
        if v.chunks is None:
            raise ValueError("zg is contiguous; chunk offsets not available")
        plist = v.id.get_create_plist()
        filters = [plist.get_filter(i)[0] for i in range(plist.get_nfilters())]
        if filters not in ([], [h5py.h5z.FILTER_DEFLATE]):
            raise ValueError(f"unsupported filter pipeline {filters}")
        deflate = bool(filters)
        if v.dtype != np.dtype("<f4"):
            raise ValueError(f"zg dtype {v.dtype}, expected little-endian float32")
        hits = np.where(lev[:] == TARGET_PA)[0]
        if len(hits) != 1:
            raise ValueError(f"no unique exact {TARGET_PA} Pa level in {dims[1]}")
        k = int(hits[0])
        nt, _, ny, nx = v.shape
        ct, cz, cy, cx = v.chunks
        fill = v.attrs.get("_FillValue")
        lat = f["lat"][:].astype(float)
        lon = f["lon"][:].astype(float)
        tvals = f["time"][:]
        tunits = _s(f["time"].attrs.get("units", b""))
        index = {}
        for i in range(v.id.get_num_chunks()):
            ci = v.id.get_chunk_info(i)
            index[tuple(ci.chunk_offset)] = ci
    rows = np.flatnonzero(np.abs(lat) >= CAP_LAT - 1e-6)
    kc, kin = (k // cz) * cz, k % cz
    full = ct * cz * cy * cx * 4
    # (chunk y0, first row, last row) for each chunk the cap rows touch
    spans = []
    for y0 in sorted({(r // cy) * cy for r in rows}):
        rr = rows[(rows >= y0) & (rows < y0 + cy)]
        if len(rr) != rr[-1] - rr[0] + 1:
            raise ValueError("cap rows not contiguous within a chunk")
        spans.append((y0, int(rr[0]), int(rr[-1])))
    starts, ends, where = [], [], []
    for t0 in range(0, nt, ct):
        for y0, ra, rb in spans:
            for x0 in range(0, nx, cx):
                info = index.get((t0, kc, y0, x0))
                if info is None or info.byte_offset is None:
                    raise ValueError(f"chunk ({t0},{kc},{y0},{x0}) not allocated")
                if deflate:
                    starts.append(info.byte_offset)
                    ends.append(info.byte_offset + info.size)
                    where.append((t0, y0, ra, rb, x0, None))
                    continue
                if info.size != full:
                    raise ValueError(f"chunk ({t0},{kc},{y0},{x0}) size {info.size}, "
                                     f"expected {full}")
                for ti in range(min(ct, nt - t0)):
                    s0 = info.byte_offset + (((ti * cz + kin) * cy + (ra - y0)) * cx) * 4
                    starts.append(s0); ends.append(s0 + (rb - ra + 1) * cx * 4)
                    where.append((t0, y0, ra, rb, x0, ti))
    blobs = fs.cat_ranges([url] * len(starts), starts, ends, on_error="return",
                          batch_size=RANGE_BATCH)
    # Re-request only the pieces that failed or came back short; under a shared
    # link a few of several hundred parallel ranges drop (seen 2026-09-24).
    for attempt in range(5):
        bad = [i for i, b in enumerate(blobs)
               if isinstance(b, Exception) or len(b) != ends[i] - starts[i]]
        if not bad:
            break
        time.sleep(2 * (attempt + 1))
        again = fs.cat_ranges([url] * len(bad), [starts[i] for i in bad],
                              [ends[i] for i in bad], on_error="return",
                              batch_size=RANGE_BATCH)
        for i, b in zip(bad, again):
            blobs[i] = b
    pos = {r: i for i, r in enumerate(rows)}
    out = np.full((nt, len(rows), nx), np.nan, dtype=np.float32)
    for (t0, y0, ra, rb, x0, ti), b, s0, e0 in zip(where, blobs, starts, ends):
        if isinstance(b, Exception) or len(b) != e0 - s0:
            raise IOError(f"short or failed range read at t={t0} y0={y0} x0={x0}")
        x1 = min(x0 + cx, nx)
        if deflate:
            raw = zlib.decompress(b)
            if len(raw) != full:
                raise IOError(f"inflated chunk is {len(raw)} B, expected {full}")
            blk = np.frombuffer(raw, dtype="<f4").reshape(ct, cz, cy, cx)
            t1 = min(t0 + ct, nt)
            out[t0:t1, pos[ra]:pos[rb] + 1, x0:x1] = \
                blk[: t1 - t0, kin, ra - y0: rb - y0 + 1, : x1 - x0]
        else:
            blk = np.frombuffer(b, dtype="<f4").reshape(rb - ra + 1, cx)
            out[t0 + ti, pos[ra]:pos[rb] + 1, x0:x1] = blk[:, : x1 - x0]
    if np.isnan(out).any():
        raise ValueError("cap rows left unfilled by the chunk reads")
    if fill is not None:
        out = np.where(out == np.float32(np.asarray(fill).ravel()[0]), np.nan, out)
    da = xr.DataArray(out, dims=("time", "lat", "lon"),
                      coords={"lat": lat[rows], "lon": lon})
    return da, tvals, tunits


def verify_ranges_against_cache(centre, tok, n=3):
    """Gate a byte-range reader against members ALREADY cached by a gated route.

    Re-reads `n` cached members with read_zg100_ranges and requires both caps and
    the lead axis to match exactly. Used when a centre moves transport after its
    full-file gate (ECCC: OPeNDAP -> byte range) or when the range reader itself
    changes (Meteo-France: whole slab -> cap rows).
    """
    d = pd.read_csv(MANIFEST)
    q = d[(d.variable == "psl") & (d.centre == centre)
          & d.experiment.isin(EXPERIMENTS)]
    have = [r for r in q.itertuples() if (
        CACHE / f"zg100_{r.centre}_{r.experiment}_{r.start_date}_{r.member}.parquet"
    ).exists()]
    if len(have) < n:
        raise ValueError(f"{centre}: only {len(have)} cached members to verify against")
    picks = [have[i] for i in np.linspace(0, len(have) - 1, n).astype(int)]
    res = []
    for r in picks:
        old = pd.read_parquet(CACHE / f"zg100_{r.centre}_{r.experiment}_"
                                      f"{r.start_date}_{r.member}.parquet")
        da, tv, tu = read_zg100_ranges(f"{FILE_ROOT}{zg_path(r, tok)}", tok)
        cn, cs = cap_means(da, da["lat"])
        res.append({"member": f"{r.experiment}/{r.start_date}/{r.member}",
                    "max_diff_N": float(np.max(np.abs(cn - old.zg100_cap_N.values))),
                    "max_diff_S": float(np.max(np.abs(cs - old.zg100_cap_S.values))),
                    "n_steps": [int(len(cn)), int(len(old))]})
    ok = all(x["max_diff_N"] == 0.0 and x["max_diff_S"] == 0.0
             and x["n_steps"][0] == x["n_steps"][1] for x in res)
    return {"centre": centre, "transport": "byte_range", "checked": res,
            "passes": bool(ok)}


def download_ranged(url: str, tok: str) -> bytes:
    """A whole file, fetched as byte ranges in parallel, resumable across runs.

    Used only for the gate's reference copy. Each GATE_PIECE-sized range is
    written to its own file under _gate_download/<name>/ once complete, so a run
    killed partway (a session ending, a hibernation) resumes from the pieces
    already on disk instead of from zero -- on 2026-09-23 a serial download at
    ~60 KB/s died at 748 MB of 1.6 GB and lost everything. Every piece is checked
    for exact length; the assembled file is checked against Content-Length.
    """
    h = {"Authorization": f"Bearer {tok}"}
    size = int(requests.head(url, headers=h, timeout=120,
                             allow_redirects=True).headers["Content-Length"])
    pdir = GATE_DL / pathlib.Path(url).name
    pdir.mkdir(parents=True, exist_ok=True)
    spans = [(lo, min(lo + GATE_PIECE, size) - 1) for lo in range(0, size, GATE_PIECE)]

    def piece(span):
        lo, hi = span
        f = pdir / f"{lo:012d}-{hi:012d}.part"
        if f.exists() and f.stat().st_size == hi - lo + 1:
            return f
        for attempt in range(8):
            try:
                r = requests.get(url, headers={**h, "Range": f"bytes={lo}-{hi}"},
                                 timeout=900)
                if r.status_code != 206 or len(r.content) != hi - lo + 1:
                    raise IOError(f"HTTP {r.status_code}, {len(r.content)} B")
                tmp = f.with_suffix(f".{os.getpid()}.tmp")
                tmp.write_bytes(r.content)
                tmp.replace(f)
                return f
            except Exception as exc:
                if attempt == 7:
                    raise IOError(f"gate piece {lo}-{hi} failed: {exc}") from exc
                time.sleep(5 * (attempt + 1))

    with ThreadPoolExecutor(GATE_WORKERS) as ex:
        files = list(ex.map(piece, spans))
    blob = b"".join(f.read_bytes() for f in files)
    if len(blob) != size:
        raise IOError(f"gate download {len(blob)} B != Content-Length {size}")
    import shutil
    shutil.rmtree(pdir, ignore_errors=True)   # the pieces are only a resume point
    return blob


def load_blocked(da):
    """Load a lazily-indexed OPeNDAP array in time blocks of TIME_BLOCK."""
    n = da.sizes["time"]
    parts = [da.isel(time=slice(i, min(i + TIME_BLOCK, n))).load()
             for i in range(0, n, TIME_BLOCK)]
    return xr.concat(parts, dim="time")


def plausible(z) -> bool:
    """100 hPa geopotential height sits near 15-17 km. Anything else is not it."""
    if z is None or len(z) == 0 or not np.isfinite(z).any():
        return False
    return 13000.0 < float(np.nanmean(z)) < 18000.0


def cache_ok(df) -> bool:
    if df is None or len(df) == 0:
        return False
    need = {"centre", "experiment", "init", "member", "lead_days",
            "zg100_cap_N", "zg100_cap_S", "hemisphere", "zg100_cap"}
    if not need <= set(df.columns):
        return False
    return plausible(df["zg100_cap_N"].values) and plausible(df["zg100_cap_S"].values)


# ------------------------------------------------- the per-centre gate ------

def verify_centre(row, tok, sess) -> dict:
    """Compare an OPeNDAP subset against the full file for the SAME field.

    Refuses the centre unless they agree exactly. This is the guard against the
    silent-zeros failure mode, in which a subsetting server returns a perfectly
    shaped array of zeros with no error and no NaN.
    """
    path = zg_path(row, tok)
    dap = f"{DAP_ROOT}/{path.lstrip('/')}"

    # The gate downloads a WHOLE file, 0.2-2.2 GB depending on centre. On
    # 2026-09-17 Meteo-France (1.6 GB) and UKMO (2.2 GB) both died with
    # IncompleteRead partway through, which refused two good centres for a
    # transport reason rather than a data one. Retry before concluding anything
    # about the subsetting.
    t0 = time.time()
    blob = download_ranged(f"{FILE_ROOT}{path}", tok)
    dl = time.time() - t0

    with xr.open_dataset(io.BytesIO(blob)) as ds:
        lat = ds["lat"]
        truth = (ds["zg"].sel({vertical_dim(ds): TARGET_PA})
                 .sel(lat=lat[lat >= CAP_LAT].values).values)

    t0 = time.time()
    if row.centre in RANGE_CENTRES:
        da, _, _ = read_zg100_ranges(f"{FILE_ROOT}{path}", tok)
        got = da.sel(lat=da.lat[da.lat >= CAP_LAT].values).values
    else:
        with xr.open_dataset(dap, engine="pydap", session=sess, decode_times=False) as rem:
            rlat = rem["lat"]
            got = load_blocked(rem["zg"].sel(plev=TARGET_PA)
                               .sel(lat=rlat[rlat >= CAP_LAT].values)).values
    sub = time.time() - t0

    a = np.asarray(truth, float)
    b = np.asarray(got, float).reshape(a.shape)
    res = {
        "centre": row.centre,
        "transport": "byte_range" if row.centre in RANGE_CENTRES else "opendap",
        "file": pathlib.Path(path).name,
        "full_MB": round(len(blob) / 1e6, 2),
        "full_seconds": round(dl, 1),
        "subset_MB": round(b.size * 4 / 1e6, 3),
        "subset_seconds": round(sub, 1),
        "saving_x": round(len(blob) / max(b.size * 4, 1), 1),
        "all_zero": bool(np.all(b == 0)),
        "any_nan": bool(np.isnan(b).any()),
        "max_abs_diff": float(np.nanmax(np.abs(b - a))),
        "nan_pattern_equal": bool(np.array_equal(np.isnan(a), np.isnan(b))),
        "correlation": float(np.corrcoef(b.ravel(), a.ravel())[0, 1]),
        "truth_mean": float(a.mean()),
        "subset_mean": float(b.mean()),
        "plausible": plausible(b),
    }
    res["passes"] = bool(
        res["max_abs_diff"] == 0.0 and res["nan_pattern_equal"] and not res["all_zero"]
        and not res["any_nan"] and res["plausible"])
    return res


# ------------------------------------------------------------- the fetch ---

def fetch(row, sess, tok=None):
    key = f"zg100_{row.centre}_{row.experiment}_{row.start_date}_{row.member}"
    cf = CACHE / f"{key}.parquet"
    if cf.exists():
        try:
            cached = pd.read_parquet(cf)
        except Exception:
            cached = None
        if cache_ok(cached):
            return cached
        print(f"  cache REJECTED, refetching: {cf.name}", flush=True)
        cf.unlink(missing_ok=True)

    last = None
    try:
        zpath = zg_path(row, tok)
    except FileNotFoundError as exc:
        raise IOError(f"{row.centre}/{row.experiment}/{row.start_date}/{row.member}: {exc}")
    dap = f"{DAP_ROOT}/{zpath.lstrip('/')}"
    for attempt in range(4):
        try:
            if row.centre in RANGE_CENTRES:
                da, tv, tu = read_zg100_ranges(f"{FILE_ROOT}{zpath}", tok)
                cap_n, cap_s = cap_means(da, da["lat"])
                lead = lead_days(tv, tu)
            else:
                with xr.open_dataset(dap, engine="pydap", session=sess,
                                     decode_times=False) as rem:
                    lat = rem["lat"]
                    da = rem["zg"].sel(plev=TARGET_PA)
                    cap_n, cap_s = cap_means(load_blocked(da), lat)
                    lead = lead_days(rem["time"].values,
                                     rem["time"].attrs.get("units"))
            if not (plausible(cap_n) and plausible(cap_s)):
                raise ValueError(
                    f"implausible 100 hPa GPH: N={np.nanmean(cap_n):.0f} "
                    f"S={np.nanmean(cap_s):.0f} m")
            hemi = "S" if str(row.start_date) in SH_INITS else "N"
            out = pd.DataFrame({
                "centre": row.centre, "model": row.model,
                "experiment": row.experiment, "init": str(row.start_date),
                "member": row.member,
                "lead_days": lead + S.origin_days(row.centre, row.experiment,
                                                  str(row.start_date)),
                "lead_origin": S.LEAD_ORIGIN,
                "zg100_cap_N": cap_n, "zg100_cap_S": cap_s,
                "hemisphere": hemi,
                "zg100_cap": cap_s if hemi == "S" else cap_n})
            # Per-process temp name: two concurrent runs (one per transport)
            # share this cache, and each cleans up only its own leftovers.
            tmp = cf.with_suffix(f".parquet.{os.getpid()}.tmp")
            out.to_parquet(tmp)
            tmp.replace(cf)
            return out
        except (ImportError, ModuleNotFoundError) as exc:
            raise SystemExit(
                f"missing dependency, not a transfer error: {exc}\n"
                f"pydap is required; check with environment/check_environment.py") from exc
        except Exception as exc:
            last = exc
            if attempt == 0:
                print(f"  retrying {row.centre}/{row.start_date}/{row.member}: "
                      f"{type(exc).__name__}: {str(exc)[:120]}", flush=True)
            time.sleep(2 * (attempt + 1))
    raise IOError(f"{row.centre}/{row.experiment}/{row.start_date}/{row.member}: {last}")


def main() -> int:
    argv = sys.argv[1:]

    def opt(name, default=None):
        return argv[argv.index(name) + 1] if name in argv else default

    centres = opt("--centres")
    centres = [c.strip() for c in centres.split(",")] if centres else ALL_CENTRES
    inits_arg = (opt("--inits") or "all").lower()
    inits = {"all": ALL_INITS, "nh": NH_INITS, "sh": SH_INITS}.get(
        inits_arg, [i.strip() for i in inits_arg.split(",")])

    CACHE.mkdir(exist_ok=True)
    d = pd.read_csv(MANIFEST)
    q = d[(d.variable == "psl")                      # psl rows carry the paths
          & (d.experiment.isin(EXPERIMENTS))
          & (d.start_date.astype(str).isin(inits))
          & (d.centre.isin(centres))].reset_index(drop=True)

    def cached(r):
        return (CACHE / f"zg100_{r.centre}_{r.experiment}_"
                        f"{r.start_date}_{r.member}.parquet").exists()

    have = q.apply(cached, axis=1) if len(q) else pd.Series(dtype=bool)
    print(f"centres    {centres}")
    print(f"inits      {inits}")
    print(f"matched    {len(q):,} members")
    print(f"cached     {int(have.sum()):,}")
    print(f"to fetch   {len(q) - int(have.sum()):,}", flush=True)
    if "--dry-run" in argv:
        return 0

    tok = token()
    sess = session(tok)

    # ---- gate every centre before fetching any of its members ----
    gate = json.loads(GATE.read_text()) if GATE.exists() else {}
    todo_centres = [c for c in q.centre.unique()
                    if not gate.get(c, {}).get("passes")]
    if "--skip-gate" in argv and todo_centres:
        print(f"WARNING --skip-gate: {todo_centres} unverified", flush=True)
        todo_centres = []
    for c in todo_centres:
        row = q[q.centre == c].iloc[0]
        how = "byte-range read" if c in RANGE_CENTRES else "OPeNDAP subset"
        print(f"\ngating {c} -- full file vs {how} ...", flush=True)
        try:
            res = verify_centre(row, tok, sess)
        except Exception as exc:
            res = {"centre": c, "passes": False, "error": f"{type(exc).__name__}: {exc}"}
        gate[c] = res
        GATE.write_text(json.dumps(gate, indent=1))
        if res.get("passes"):
            print(f"  PASS  max|diff|={res['max_abs_diff']:g} "
                  f"corr={res['correlation']:.10f} "
                  f"mean={res['subset_mean']:.1f} m  "
                  f"{res['saving_x']}x smaller", flush=True)
        else:
            print(f"  FAIL  {res}", flush=True)

    for c in [c for c in q.centre.unique() if c in RANGE_CENTRES
              and gate.get(c, {}).get("passes")
              and not gate.get(f"{c}|byte_range", {}).get("passes")]:
        print(f"\nverifying the byte-range reader for {c} against its cache ...",
              flush=True)
        try:
            res = verify_ranges_against_cache(c, tok)
        except Exception as exc:
            res = {"centre": c, "passes": False, "error": f"{type(exc).__name__}: {exc}"}
        gate[f"{c}|byte_range"] = res
        GATE.write_text(json.dumps(gate, indent=1))
        print(f"  {'PASS' if res.get('passes') else 'FAIL'}  {res}", flush=True)
    ok_centres = [c for c in q.centre.unique() if gate.get(c, {}).get("passes")
                  and (c not in RANGE_CENTRES
                       or gate.get(f"{c}|byte_range", {}).get("passes"))]
    refused = sorted(set(q.centre.unique()) - set(ok_centres))
    if refused:
        print(f"\nREFUSED (failed the subset gate): {refused}", flush=True)
    q = q[q.centre.isin(ok_centres)].reset_index(drop=True)
    if q.empty:
        sys.exit("no centre passed the gate")

    rows, failed = [], []
    t0 = time.time()
    workers = int(opt("--workers", WORKERS))
    with ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(fetch, r, sess, tok): i for i, r in enumerate(q.itertuples())}
        for n, f in enumerate(as_completed(futs), 1):
            try:
                rows.append(f.result())
            except Exception as exc:
                failed.append(exc)
                print(f"  FAILED {exc}", flush=True)
            if n % 100 == 0:
                el = time.time() - t0
                print(f"  {n}/{len(q)}  {el/60:.1f} min  "
                      f"({n/el*60:.0f} members/min)", flush=True)

    rb = S.rebase_cache(CACHE, "zg100_*.parquet")
    if rb:
        print(f"rebased {rb} cached member(s) to lead 0 = 00 UTC on the init date")
    files = sorted(CACHE.glob("zg100_*.parquet"))
    for s in CACHE.glob(f"zg100_*.parquet.{os.getpid()}.tmp"):
        s.unlink(missing_ok=True)
    frames, dropped = [], 0
    for f in files:
        try:
            df = pd.read_parquet(f)
        except Exception:
            dropped += 1
            continue
        if cache_ok(df):
            frames.append(df)
        else:
            dropped += 1
    if dropped:
        print(f"excluded {dropped} member(s) that failed validation")
    if not frames:
        sys.exit("nothing usable")
    out = pd.concat(frames, ignore_index=True)
    out.to_parquet(OUT)
    print(f"\n{len(out):,} rows from {len(frames):,} members -> {OUT.name} "
          f"({len(failed)} failed this run)")
    print(out.groupby(["centre", "experiment"]).init.nunique()
          .unstack(fill_value=0).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
