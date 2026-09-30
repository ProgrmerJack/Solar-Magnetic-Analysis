#!/usr/bin/env python3
"""
acquire_snotel_stations.py
==========================
Acquires SNOTEL/COOP station metadata -- elevation, latitude, longitude -- for
the 945 stations in the daily archive.

WHY THIS IS NEEDED
  Gate 3 is closed for index outcomes and for spatially AGGREGATED outcomes, but
  not for station-level ones. The concern the gate exists to test is that
  seasonality varies with elevation, latitude and aspect; aggregating to states
  averages exactly that away, which is why the aggregated result (shared and
  unit-specific seasonality indistinguishable) cannot be extrapolated downward.

  To run the station-level calibration honestly we first need to know how
  heterogeneous the stations actually are. Without elevation and latitude there
  is no way to show that the hard case is hard.

Source: USDA NRCS Air and Water Database (AWDB) REST API
        https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1/stations
        Public, no key required.

Outputs (this directory):
  raw/snotel_stations.json           the download cache, one entry per
                                     STATE_NETWORK, each carrying the sha256 of
                                     the original API response
  snotel_stations.csv                tidy metadata joined to the daily archive
  station_checksums.txt              sha256 per cache entry

CACHE LAYOUT CHANGED 2026-08-03. This previously wrote 36 separate
raw/snotel_stations_<STATE>_<NETWORK>.json files. They were a download cache
read by nothing but this script, and they made the raw/ directory
unreadable. They are now one keyed store. The sha256 recorded per key is still
the hash of the ORIGINAL response bytes, so `station_checksums.txt` is unchanged
by the migration and provenance is preserved.
"""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RAW = HERE / "raw"
OUT_CSV = HERE / "snotel_stations.csv"
OUT_SUMS = HERE / "station_checksums.txt"
API = "https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1/stations"
TIMEOUT = 90


CACHE = RAW / "snotel_stations.json"


def load_cache():
    if CACHE.exists():
        return json.loads(CACHE.read_text(encoding="utf8"))
    return {}


def save_cache(store):
    RAW.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(store, indent=1), encoding="utf8")


def fetch_state(state, network, store):
    """Return (records, sha256_of_original_response, was_cached).

    The hash is of the bytes the API actually returned, kept verbatim in the
    cache entry, so it does not change when the cache is rewritten.
    """
    key = f"{state}_{network}"
    hit = store.get(key)
    if hit:
        return hit["stations"], hit["sha256"], True
    r = requests.get(API, params={"stationTriplets": f"*:{state}:{network}"},
                     timeout=TIMEOUT,
                     headers={"User-Agent": "ssw-design-analysis/1.0 (research)"})
    r.raise_for_status()
    h = hashlib.sha256(r.content).hexdigest()
    data = r.json()
    store[key] = {"sha256": h, "n": len(data), "stations": data}
    return data, h, False


def main():
    # states and networks actually present in the daily archive
    s = pd.read_parquet(ROOT / "data/processed/cryosphere/snotel_daily.parquet",
                        columns=["station_id"])
    ids = s["station_id"].astype(str).unique()
    parts = pd.Series(ids).str.split(":", expand=True)
    parts.columns = ["sid", "state", "network"]
    combos = sorted(set(zip(parts["state"], parts["network"])))
    print(f"{len(ids)} stations in the archive across "
          f"{parts['state'].nunique()} states, networks "
          f"{sorted(parts['network'].unique())}")

    store = load_cache()
    rows, sums, cached = [], [], 0
    for state, network in combos:
        try:
            data, h, was_cached = fetch_state(state, network, store)
        except Exception as e:
            print(f"  ! {state}:{network}: {e}")
            continue
        cached += was_cached
        sums.append((f"{state}_{network}", h, len(data)))
        for r in data:
            rows.append({
                "station_id": r.get("stationTriplet"),
                "name": r.get("name"),
                "state": r.get("stateCode"),
                "network": r.get("networkCode"),
                "elevation_ft": r.get("elevation"),
                "latitude": r.get("latitude"),
                "longitude": r.get("longitude"),
                "county": r.get("countyName"),
                "begin_date": r.get("beginDate"),
                "end_date": r.get("endDate"),
            })
    save_cache(store)
    meta = pd.DataFrame(rows).drop_duplicates(subset=["station_id"])
    meta["elevation_m"] = pd.to_numeric(meta["elevation_ft"], errors="coerce") * 0.3048

    # how many of the archive's stations did we resolve?
    have = set(meta["station_id"])
    matched = [i for i in ids if i in have]
    meta = meta[meta["station_id"].isin(ids)]
    meta.to_csv(OUT_CSV, index=False, lineterminator="\n")

    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    OUT_SUMS.write_text(
        "\n".join(f"{h}  {n}  {c} stations  retrieved {stamp}"
                  for n, h, c in sums) + "\n", encoding="utf8")

    print(f"\nresolved {len(matched)}/{len(ids)} archive stations "
          f"({cached} state files from cache)")
    print(f"metadata -> {OUT_CSV.name}   checksums -> {OUT_SUMS.name}")
    e = meta["elevation_m"].dropna()
    if len(e):
        print(f"\nelevation (m): min {e.min():.0f}  p25 {e.quantile(.25):.0f}  "
              f"median {e.median():.0f}  p75 {e.quantile(.75):.0f}  max {e.max():.0f}")
        print(f"  range spanned: {e.max() - e.min():.0f} m")
    la = meta["latitude"].dropna()
    if len(la):
        print(f"latitude: {la.min():.1f} .. {la.max():.1f} "
              f"({la.max() - la.min():.1f} deg span)")
    if len(matched) < len(ids):
        miss = [i for i in ids if i not in have][:5]
        print(f"\nunresolved examples: {miss}")


if __name__ == "__main__":
    main()
