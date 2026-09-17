#!/usr/bin/env python3
"""
download_lawis_full.py
======================
Download the COMPLETE LAWIS avalanche-incident record (all countries, all
regions, all years) from the public LAWIS API.

Why: the repo previously held only `incidents_tirol.csv` (3,060 records, Tirol
only). The full API holds ~4,860 incidents across 11 countries back to 1992.
The binding constraint on this project's inference is the number of CONTROL
WINTERS, so every additional region-winter matters.

API: https://lawis.at/lawis_api/public/swagger/  (spec at .../public/api.php)
  GET /lawis_api/public/incident?startDate=&endDate=   -> list (date, location)
  GET /lawis_api/public/incident/{id}?lang=en          -> detail (involved,
                                                          avalanche type/size,
                                                          danger rating)

Outputs (data/cryosphere/lawis_full/):
  incidents_raw.json      full list response
  details/<id>.json       per-incident detail (cached; safe to re-run)
  incidents_full.csv      tidy analysis table

Re-running is safe and incremental: cached details are not re-fetched.
"""
import csv
import json
import sys
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "cryosphere" / "lawis_full"
DETAIL_DIR = OUT / "details"
BASE = "https://lawis.at/lawis_api/public"
START, END = "1960-01-01", "2026-12-31"
PAUSE = 0.12          # polite delay between detail calls
TIMEOUT = 60

session = requests.Session()
session.headers.update({
    "User-Agent": "Solar-Magnetic-Analysis/1.0 (academic research; SSW-avalanche coupling)",
    "Accept": "application/json",
})


def fetch_list():
    OUT.mkdir(parents=True, exist_ok=True)
    raw = OUT / "incidents_raw.json"
    if raw.exists():
        print(f"[list] cached -> {raw}")
        return json.loads(raw.read_text(encoding="utf8"))
    url = f"{BASE}/incident?startDate={START}&endDate={END}"
    print(f"[list] GET {url}")
    r = session.get(url, timeout=TIMEOUT)
    r.raise_for_status()
    data = r.json()
    raw.write_text(json.dumps(data, ensure_ascii=False), encoding="utf8")
    print(f"[list] {len(data)} incidents -> {raw}")
    return data


def fetch_details(ids):
    DETAIL_DIR.mkdir(parents=True, exist_ok=True)
    todo = [i for i in ids if not (DETAIL_DIR / f"{i}.json").exists()]
    print(f"[detail] {len(ids)} total, {len(todo)} to fetch")
    ok = fail = 0
    for n, i in enumerate(todo, 1):
        try:
            r = session.get(f"{BASE}/incident/{i}?lang=en", timeout=TIMEOUT)
            r.raise_for_status()
            (DETAIL_DIR / f"{i}.json").write_text(
                json.dumps(r.json(), ensure_ascii=False), encoding="utf8")
            ok += 1
        except Exception as e:
            fail += 1
            print(f"  ! {i}: {e}", file=sys.stderr)
        if n % 250 == 0:
            print(f"  ... {n}/{len(todo)} (ok={ok} fail={fail})")
        time.sleep(PAUSE)
    print(f"[detail] done ok={ok} fail={fail}")


def g(d, *keys, default=None):
    """Nested get that tolerates nulls at any level."""
    for k in keys:
        if not isinstance(d, dict):
            return default
        d = d.get(k)
    return default if d is None else d


def build_table(lst):
    rows = []
    for it in lst:
        iid = it.get("id")
        det = {}
        f = DETAIL_DIR / f"{iid}.json"
        if f.exists():
            try:
                det = json.loads(f.read_text(encoding="utf8"))
            except Exception:
                det = {}
        src = det or it
        rows.append({
            "id": iid,
            "date": it.get("date"),
            "valid_time": it.get("valid_time"),
            "country_code": g(src, "location", "country", "code"),
            "country": g(src, "location", "country", "text"),
            "region": g(src, "location", "region", "text"),
            "region_id": g(src, "location", "region", "id"),
            "subregion": g(src, "location", "subregion", "text"),
            "longitude": g(src, "location", "longitude"),
            "latitude": g(src, "location", "latitude"),
            "aspect": g(src, "location", "aspect", "text"),
            "danger_level": g(src, "danger", "rating", "level"),
            "avalanche_type": g(det, "avalanche", "type", "text"),
            "avalanche_size": g(det, "avalanche", "size", "text"),
            "breakheight": g(det, "avalanche", "breakheight"),
            "dead": g(det, "involved", "dead"),
            "injured": g(det, "involved", "injured"),
            "uninjured": g(det, "involved", "uninjured"),
            "buried_total": g(det, "involved", "buried_total"),
            "buried_partial": g(det, "involved", "buried_partial"),
            "sweeped": g(det, "involved", "sweeped"),
        })
    out = OUT / "incidents_full.csv"
    with open(out, "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"[table] {len(rows)} rows -> {out}")
    return rows


def main():
    lst = fetch_list()
    ids = [i["id"] for i in lst if i.get("id") is not None]
    fetch_details(ids)
    rows = build_table(lst)

    # summary
    from collections import Counter
    yrs = [r["date"][:4] for r in rows if r["date"]]
    print(f"\nyears {min(yrs)}-{max(yrs)}")
    print("by country:", dict(Counter(r["country_code"] for r in rows).most_common()))
    at = Counter(r["region"] for r in rows if r["country_code"] == "AT")
    print("Austrian regions:", dict(at.most_common()))
    print("with detail (dead not null):", sum(1 for r in rows if r["dead"] is not None))


if __name__ == "__main__":
    main()
