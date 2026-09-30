#!/usr/bin/env python3
"""
inventory_datasets.py
=====================
Systematic inventory of EVERY dataset under data/: format, size, variables,
temporal coverage, row counts, and (where derivable) what it physically measures.

Reads metadata only wherever possible -- parquet footers, netCDF headers, CSV
headers plus a date-column scan -- so it runs over ~54 GB without loading it.

Written because the project's data had never been catalogued end-to-end, and
conclusions about "which observables exist" must rest on the actual contents,
not on recollection of which download script produced what.

Output: DATA_INVENTORY.md  (+ data/results/data_inventory.json)
"""
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
OUT_MD = ROOT / "docs" / "DATA_INVENTORY.md"
OUT_JSON = ROOT / "data" / "results" / "data_inventory.json"

DATE_HINTS = ("date", "datum", "time", "datetime", "timestamp", "validdate",
              "date_release", "day", "onset")
SKIP_DIRS = {"__pycache__", ".ipynb_checkpoints"}


def human(n):
    for u in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or u == "TB":
            return f"{n:.0f}{u}" if u == "B" else f"{n:.1f}{u}"
        n /= 1024


def date_range_from_series(s):
    try:
        d = pd.to_datetime(s, errors="coerce", utc=True, format="mixed")
    except Exception:
        try:
            d = pd.to_datetime(s, errors="coerce", utc=True)
        except Exception:
            return None
    d = d.dropna()
    if len(d) == 0:
        return None
    return f"{d.min().date()} .. {d.max().date()}", int(len(d))


def probe_parquet(p):
    info = {"format": "parquet"}
    try:
        import pyarrow.parquet as pq
        f = pq.ParquetFile(p)
        info["rows"] = f.metadata.num_rows
        sch = f.schema_arrow
        info["columns"] = list(sch.names)[:60]
        info["n_columns"] = len(sch.names)
        # index / time coverage: read only candidate date columns
        cand = [c for c in sch.names if any(h in c.lower() for h in DATE_HINTS)]
        got = False
        if cand:
            df = pd.read_parquet(p, columns=cand[:1])
            r = date_range_from_series(df[cand[0]])
            if r:
                info["time_coverage"], _ = r
                info["time_column"] = cand[0]
                got = True
        if not got:
            # DatetimeIndex stored as the parquet index
            df = pd.read_parquet(p, columns=[sch.names[0]])
            if isinstance(df.index, pd.DatetimeIndex) and len(df.index):
                info["time_coverage"] = f"{df.index.min().date()} .. {df.index.max().date()}"
                info["time_column"] = "(index)"
    except Exception as e:
        info["error"] = f"{type(e).__name__}: {e}"
    return info


def probe_netcdf(p):
    info = {"format": "netcdf"}
    try:
        import netCDF4
        with netCDF4.Dataset(p) as ds:
            info["dimensions"] = {k: len(v) for k, v in ds.dimensions.items()}
            info["variables"] = list(ds.variables.keys())[:40]
            info["n_variables"] = len(ds.variables)
            for tname in ("time", "TIME", "t", "date"):
                if tname in ds.variables:
                    tv = ds.variables[tname]
                    try:
                        t = netCDF4.num2date(tv[[0, -1]], tv.units,
                                             only_use_cftime_datetimes=False)
                        info["time_coverage"] = f"{t[0].date()} .. {t[-1].date()}"
                    except Exception:
                        info["time_coverage"] = f"{len(tv)} steps (units unparsed)"
                    break
    except Exception as e:
        info["error"] = f"{type(e).__name__}: {e}"
    return info


def probe_csv(p, size):
    info = {"format": "csv"}
    try:
        head = pd.read_csv(p, nrows=200, sep=None, engine="python",
                           on_bad_lines="skip")
        info["columns"] = list(head.columns)[:60]
        info["n_columns"] = len(head.columns)
        # exact row count only for files small enough to be cheap
        if size < 300 * 1024 * 1024:
            with open(p, "rb") as fh:
                info["rows"] = max(0, sum(1 for _ in fh) - 1)
        cand = [c for c in head.columns if any(h in str(c).lower() for h in DATE_HINTS)]
        if cand and size < 300 * 1024 * 1024:
            col = cand[0]
            s = pd.read_csv(p, usecols=[col], sep=None, engine="python",
                            on_bad_lines="skip")[col]
            r = date_range_from_series(s)
            if r:
                info["time_coverage"], _ = r
                info["time_column"] = col
    except Exception as e:
        info["error"] = f"{type(e).__name__}: {e}"
    return info


def probe_json(p, size):
    info = {"format": "json"}
    if size > 200 * 1024 * 1024:
        info["note"] = "too large to parse; skipped"
        return info
    try:
        d = json.loads(p.read_text(encoding="utf8", errors="replace"))
        if isinstance(d, dict):
            info["top_level_keys"] = list(d.keys())[:25]
            info["n_keys"] = len(d)
        elif isinstance(d, list):
            info["n_records"] = len(d)
            if d and isinstance(d[0], dict):
                info["record_keys"] = list(d[0].keys())[:25]
    except Exception as e:
        info["error"] = f"{type(e).__name__}: {e}"
    return info


def probe(p):
    size = p.stat().st_size
    ext = p.suffix.lower()
    if ext == ".parquet":
        info = probe_parquet(p)
    elif ext in (".nc", ".nc4", ".h5", ".hdf5"):
        info = probe_netcdf(p)
    elif ext in (".csv", ".tsv"):
        info = probe_csv(p, size)
    elif ext == ".json":
        info = probe_json(p, size)
    elif ext == ".xlsx":
        info = {"format": "xlsx"}
        try:
            x = pd.read_excel(p)
            info["rows"] = len(x)
            info["columns"] = list(x.columns)[:40]
            cand = [c for c in x.columns if any(h in str(c).lower() for h in DATE_HINTS)]
            if cand:
                r = date_range_from_series(x[cand[0]])
                if r:
                    info["time_coverage"], _ = r
                    info["time_column"] = cand[0]
        except Exception as e:
            info["error"] = str(e)
    else:
        info = {"format": ext.lstrip(".") or "none"}
    info["size_bytes"] = size
    info["size"] = human(size)
    return info


def main():
    groups = {}          # top-level dir -> list of (relpath, info)
    totals = {}
    n = 0
    for dirpath, dirnames, filenames in os.walk(DATA):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for fn in filenames:
            p = Path(dirpath) / fn
            try:
                rel = p.relative_to(DATA)
            except ValueError:
                continue
            top = rel.parts[0] if len(rel.parts) > 1 else "(root)"
            size = p.stat().st_size
            totals[top] = totals.get(top, 0) + size
            # probe everything except the very large homogeneous satellite sets,
            # which are sampled instead (one representative file per family)
            groups.setdefault(top, []).append((str(rel), p, size))
            n += 1
    print(f"walked {n} files under {DATA}")

    report = {}
    for top, items in sorted(groups.items(), key=lambda kv: -totals[kv[0]]):
        items.sort(key=lambda t: -t[2])
        # probe the 40 largest per group, plus every cryosphere/results file
        probe_all = top in ("cryosphere", "results")
        chosen = items if probe_all else items[:40]
        entries = []
        for rel, p, size in chosen:
            info = probe(p)
            info["path"] = rel
            entries.append(info)
        report[top] = {
            "total_size": human(totals[top]),
            "total_size_bytes": totals[top],
            "n_files": len(items),
            "probed": len(entries),
            "entries": entries,
        }
        print(f"  {top:14s} {human(totals[top]):>9s}  {len(items):5d} files  "
              f"probed {len(entries)}")

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(report, indent=1, default=str), encoding="utf8")

    # ---- markdown ----
    L = ["# Data inventory", "",
         f"Generated {datetime.now():%Y-%m-%d} by `scripts/utilities/inventory_datasets.py`.",
         f"Total: **{human(sum(totals.values()))}** across **{n}** files.", "",
         "| group | size | files |", "|---|---|---|"]
    for top in sorted(totals, key=lambda k: -totals[k]):
        L.append(f"| `{top}` | {human(totals[top])} | {len(groups[top])} |")
    L.append("")

    for top in sorted(totals, key=lambda k: -totals[k]):
        g = report[top]
        L += [f"## `data/{top}` — {g['total_size']}, {g['n_files']} files",
              f"(probed {g['probed']})", ""]
        for e in g["entries"]:
            bits = [f"**`{e['path']}`** ({e['size']}, {e.get('format','?')})"]
            if e.get("rows") is not None:
                bits.append(f"rows={e['rows']:,}")
            if e.get("n_records") is not None:
                bits.append(f"records={e['n_records']:,}")
            if e.get("time_coverage"):
                bits.append(f"**{e['time_coverage']}**")
            L.append("- " + " · ".join(bits))
            cols = e.get("columns") or e.get("variables") or e.get("record_keys") \
                or e.get("top_level_keys")
            if cols:
                shown = ", ".join(map(str, cols[:25]))
                more = f" … (+{e.get('n_columns', e.get('n_variables', len(cols))) - 25})" \
                    if (e.get("n_columns") or e.get("n_variables") or len(cols)) > 25 else ""
                L.append(f"  - fields: `{shown}`{more}")
            if e.get("dimensions"):
                L.append(f"  - dims: `{e['dimensions']}`")
            if e.get("error"):
                L.append(f"  - ⚠ {e['error']}")
        L.append("")
    OUT_MD.write_text("\n".join(L), encoding="utf8")
    print(f"\nwrote {OUT_MD}  and  {OUT_JSON}")


if __name__ == "__main__":
    main()
