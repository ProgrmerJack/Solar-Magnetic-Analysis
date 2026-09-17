#!/usr/bin/env python3
"""
download_albina_caaml_problems.py
=================================
Download the full ALBINA (Tyrol / South Tyrol / Trentino) avalanche-bulletin
archive in CAAML v6, which carries the EAWS-standard STRUCTURED avalanche
problems:

    problemType        new_snow | wind_slab | persistent_weak_layers |
                       wet_snow | gliding_snow
    snowpackStability  very_poor | poor | fair | good
    frequency          none | few | some | many
    avalancheSize      1..5
    aspects, elevation, validTimePeriod

WHY: R109 found that the manuscript's trigger-dissociation mechanism replicates
in the Norwegian bulletin archive (natural "spontaneous release" problems down,
human-triggerable problems up, clean placebos), while its structural /
persistent-weak-layer limb failed. ALBINA gives the same class of observable in
the ALPS -- the manuscript's own mechanism domain -- so the Norway result can be
tested where the Alpine blocking mechanism is actually claimed to operate, and
`persistent_weak_layers` gives a forecaster-assessed PWL measure independent of
the SNOWPACK model output that R108 showed to be seasonally confounded.

The repository previously held only ALBINA *danger levels*
(`data/cryosphere/albina/albina_danger.json`), not the problems.

Source: https://avalanche.report/albina_files/{date}/{date}_en_CAAMLv6.json
Licence: CC BY 4.0 (avalanche.report / EUREGIO ALBINA)

Outputs (data/cryosphere/albina_caaml/):
  raw/<date>.json          verbatim daily CAAML (cached; re-runs are incremental)
  albina_problems.csv      one row per (date, region, problem)
  albina_danger.csv        one row per (date, region) with max danger

Re-running is safe: cached days are not re-fetched.
"""
import csv
import json
import sys
import time
from datetime import date, timedelta
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "cryosphere" / "albina_caaml"
RAW = OUT / "raw"
URL = "https://avalanche.report/albina_files/{d}/{d}_en_CAAMLv6.json"

SEASON_MONTHS = (11, 12, 1, 2, 3, 4)
START = date(2018, 11, 1)
END = date(2026, 5, 31)
PAUSE = 0.15
TIMEOUT = 45

DANGER_NUM = {"low": 1, "moderate": 2, "considerable": 3, "high": 4,
              "very_high": 5, "no_rating": None, "no_snow": None}

session = requests.Session()
session.headers.update({
    "User-Agent": "Solar-Magnetic-Analysis/1.0 (academic research; SSW-avalanche coupling)",
    "Accept": "application/json",
})


def winter_days():
    d = START
    while d <= END:
        if d.month in SEASON_MONTHS:
            yield d
        d += timedelta(days=1)


def fetch_all():
    RAW.mkdir(parents=True, exist_ok=True)
    days = list(winter_days())
    todo = [d for d in days if not (RAW / f"{d}.json").exists()]
    print(f"[fetch] {len(days)} winter days in range, {len(todo)} to fetch")
    ok = miss = err = 0
    for n, d in enumerate(todo, 1):
        try:
            r = session.get(URL.format(d=d.isoformat()), timeout=TIMEOUT)
            if r.status_code == 200 and r.content.strip():
                (RAW / f"{d}.json").write_bytes(r.content)
                ok += 1
            elif r.status_code in (403, 404):
                # no bulletin that day (out of season / before launch)
                (RAW / f"{d}.json").write_text("{}", encoding="utf8")
                miss += 1
            else:
                err += 1
        except Exception as e:
            err += 1
            print(f"  ! {d}: {e}", file=sys.stderr)
        if n % 200 == 0:
            print(f"  ... {n}/{len(todo)} (ok={ok} empty={miss} err={err})")
        time.sleep(PAUSE)
    print(f"[fetch] done ok={ok} empty={miss} err={err}")


def parse():
    prob_rows, dang_rows = [], []
    for f in sorted(RAW.glob("*.json")):
        d = f.stem
        try:
            doc = json.loads(f.read_text(encoding="utf8"))
        except Exception:
            continue
        for b in (doc.get("bulletins") or []):
            regions = [r.get("regionID") for r in (b.get("regions") or [])
                       if r.get("regionID")]
            if not regions:
                continue
            drs = b.get("dangerRatings") or []
            vals = [DANGER_NUM.get(x.get("mainValue")) for x in drs]
            vals = [v for v in vals if v is not None]
            dmax = max(vals) if vals else None
            probs = b.get("avalancheProblems") or []
            for reg in regions:
                dang_rows.append({"date": d, "region": reg,
                                  "country": reg.split("-")[0],
                                  "danger_max": dmax,
                                  "n_problems": len(probs)})
                for p in probs:
                    el = p.get("elevation") or {}
                    prob_rows.append({
                        "date": d, "region": reg, "country": reg.split("-")[0],
                        "problem_type": p.get("problemType"),
                        "snowpack_stability": p.get("snowpackStability"),
                        "frequency": p.get("frequency"),
                        "avalanche_size": p.get("avalancheSize"),
                        "valid_time_period": p.get("validTimePeriod"),
                        "elev_lower": el.get("lowerBound"),
                        "elev_upper": el.get("upperBound"),
                        "aspects": "|".join(p.get("aspects") or []),
                        "danger_max": dmax,
                    })
    OUT.mkdir(parents=True, exist_ok=True)
    for rows, name in ((prob_rows, "albina_problems.csv"),
                       (dang_rows, "albina_danger.csv")):
        if not rows:
            print(f"[parse] NO ROWS for {name}")
            continue
        with open(OUT / name, "w", newline="", encoding="utf8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"[parse] {len(rows):,} rows -> {OUT / name}")
    return prob_rows, dang_rows


def main():
    fetch_all()
    prob, dang = parse()
    if not prob:
        return
    from collections import Counter
    dates = sorted({r["date"] for r in dang})
    print(f"\ncoverage {dates[0]} .. {dates[-1]}  ({len(dates)} days with bulletins)")
    print("regions:", len({r['region'] for r in dang}),
          "| countries:", dict(Counter(r["country"] for r in dang)))
    print("problem types:", dict(Counter(r["problem_type"] for r in prob)))
    print("snowpack stability:", dict(Counter(r["snowpack_stability"] for r in prob)))
    print("frequency:", dict(Counter(r["frequency"] for r in prob)))


if __name__ == "__main__":
    main()
