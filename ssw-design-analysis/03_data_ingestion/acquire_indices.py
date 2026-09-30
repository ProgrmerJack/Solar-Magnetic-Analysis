#!/usr/bin/env python3
"""
acquire_indices.py
==================
Acquires the circulation indices needed for Gate 6 ("the result is not driven by
AO/NAO alone") and records a checksum for every raw file, as the plan requires.

WHAT IS ACQUIRED AND WHY

  ao   Arctic Oscillation, daily 1950-      already the primary index
  nao  North Atlantic Oscillation, daily    already a secondary index
  pna  Pacific/North American, daily        NEW. A different sector of the same
                                            hemisphere, so it tests whether the
                                            response is annular or Atlantic-specific.
  aao  Antarctic Oscillation, daily 1979-   NEW, and the important one: a
                                            NEGATIVE CONTROL. A Northern
                                            Hemisphere mid-winter SSW has no
                                            mechanism by which it should move the
                                            Southern annular mode on a 0-60 day
                                            lag. If the estimator reports an
                                            effect here, something is wrong with
                                            the estimator, not with the atmosphere.
  tele_index.nh  monthly standardised NH teleconnection indices (EA, WP, EP/NP,
                 EA/WR, SCA, TNH, POL and others) -- breadth for Gate 6 at
                 monthly resolution.

Southern-hemisphere seasonality is inverted, so the AAO control is run on the
SAME calendar days as the northern analysis (Nov-Apr). That is deliberate: the
question is whether the estimator manufactures a signal from those particular
dates, not whether the Antarctic has a winter.

Outputs (this directory):
  raw/<name>.txt            verbatim download
  index_checksums.txt       sha256 of every raw file, with source URL and date
  daily_indices.parquet     tidy daily table: date, index, value
"""
import hashlib
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
OUT_PARQUET = HERE / "daily_indices.parquet"
OUT_SUMS = HERE / "index_checksums.txt"

CWLINKS = "https://ftp.cpc.ncep.noaa.gov/cwlinks"
DAILY = {
    "ao":  f"{CWLINKS}/norm.daily.ao.index.b500101.current.ascii",
    "nao": f"{CWLINKS}/norm.daily.nao.index.b500101.current.ascii",
    "pna": f"{CWLINKS}/norm.daily.pna.index.b500101.current.ascii",
    "aao": f"{CWLINKS}/norm.daily.aao.index.b790101.current.ascii",
}
MONTHLY = {
    "tele_index_nh": "https://ftp.cpc.ncep.noaa.gov/wd52dg/data/indices/tele_index.nh",
}


def fetch(name, url):
    RAW.mkdir(parents=True, exist_ok=True)
    p = RAW / f"{name}.txt"
    if p.exists():
        print(f"  [cached] {name}")
        return p, hashlib.sha256(p.read_bytes()).hexdigest()
    r = requests.get(url, timeout=90,
                     headers={"User-Agent": "ssw-design-analysis/1.0 (research)"})
    r.raise_for_status()
    p.write_bytes(r.content)
    print(f"  [got]    {name}  {len(r.content):,} bytes")
    return p, hashlib.sha256(r.content).hexdigest()


def parse_daily(p, name):
    d = pd.read_csv(p, sep=r"\s+", header=None,
                    names=["year", "month", "day", "value"])
    d["date"] = pd.to_datetime(dict(year=d.year, month=d.month, day=d.day),
                               errors="coerce")
    d = d.dropna(subset=["date", "value"])
    d = d[d["value"].abs() < 90]          # CPC uses -99.9 for missing
    return pd.DataFrame({"date": d["date"], "index": name,
                         "value": d["value"].astype(float)})


def main():
    sums, frames = [], []
    print("daily indices:")
    for name, url in DAILY.items():
        p, h = fetch(name, url)
        sums.append((name, url, h))
        f = parse_daily(p, name)
        frames.append(f)
        print(f"           {name}: {len(f):,} days "
              f"{f['date'].min().date()}..{f['date'].max().date()}")

    print("monthly indices:")
    for name, url in MONTHLY.items():
        p, h = fetch(name, url)
        sums.append((name, url, h))

    tidy = pd.concat(frames, ignore_index=True).sort_values(["index", "date"])
    tidy.to_parquet(OUT_PARQUET, index=False)

    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    OUT_SUMS.write_text(
        "\n".join(f"{h}  {n}  {u}  retrieved {stamp}" for n, u, h in sums) + "\n",
        encoding="utf8")

    print(f"\ntidy table -> {OUT_PARQUET.name}  ({len(tidy):,} rows, "
          f"{tidy['index'].nunique()} indices)")
    print(f"checksums  -> {OUT_SUMS.name}")
    # overlap with the frozen catalogue period, reported per index
    for nm, g in tidy.groupby("index"):
        w = g[g["date"].dt.month.isin((11, 12, 1, 2, 3, 4))]
        print(f"    {nm:5s} winter days {len(w):,}  "
              f"{g['date'].min().date()}..{g['date'].max().date()}")


if __name__ == "__main__":
    main()
