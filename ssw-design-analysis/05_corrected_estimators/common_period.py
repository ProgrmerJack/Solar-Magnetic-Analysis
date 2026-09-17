#!/usr/bin/env python3
"""
common_period.py
================
Separates the two things that were being conflated in "catalogue choice moves
the estimate 50x more than estimator choice".

  (a) DEFINITION  -- the catalogues genuinely disagree about which winters
                     contain an event, and about how many.
  (b) RECORD      -- they cover different periods. MERRA2 starts in 1980 and
                     contributes 27 events; NCEP-NCAR spans 1948- and
                     contributes 39. Comparing their profiles compares two
                     climates as much as two event lists.

Fixing the analysis window to the span all six reanalyses cover (1980-2019,
the ERA-Interim/MERRA2 overlap) removes (b). Whatever spread survives is (a).

If the spread collapses once the window is common, then "catalogue choice
matters" was largely "sample period matters", and the design claim has to be
restated -- reanalysis products are not really disagreeing, they are being
asked about different decades.

Output: common_period.json
"""
import json, sys, zlib
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "2_event_study"
RESULTS.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "06_simulation_validation"))
from build_catalogue import load_catalogue
import canonical_event_study as K

SETS = ["primary", "consensus_strict", "consensus_half", "union",
        "era5", "jra_55", "ncep_ncar", "merra2"]
WINDOWS = {"full_record": (None, None), "common_1980_2019": ("1980-01-01", "2019-12-31")}
WATCH = ["-30..-16", "+15..+29"]


def main():
    d0 = K.load_series("AO")
    res = {"n_boot": K.N_BOOT, "windows": {}}
    for wname, (lo, hi) in WINDOWS.items():
        d = d0 if lo is None else d0[(d0.index >= lo) & (d0.index <= hi)]
        print(f"\n=== AO, window {wname}  "
              f"({d.index.min().date()}..{d.index.max().date()}, "
              f"{d['winter'].nunique()} winters) ===")
        entry = {}
        for s in SETS:
            on = load_catalogue(s)
            on = on[(on >= d.index.min()) & (on <= d.index.max())]
            if len(on) < 10:
                print(f"  {s:18s} skipped ({len(on)} events)")
                continue
            prof, _ = K.run(d, on, seed=zlib.crc32(s.encode()) % 10000)
            entry[s] = {"n_events": int(len(on)),
                        **{w: {"effect": prof[w]["effect"],
                               "p": prof[w]["p_two_sided"]} for w in WATCH}}
            print(f"  {s:18s} n={len(on):2d}  "
                  + "  ".join(f"{w} {prof[w]['effect']:+.3f} (p={prof[w]['p_two_sided']:.3f})"
                              for w in WATCH))
        done = list(entry)          # snapshot: the loop below adds _spread keys
        for w in WATCH:
            v = [entry[s][w]["effect"] for s in done]
            entry[f"_spread{w}"] = round(float(max(v) - min(v)), 4)
            print(f"  -> spread across catalogues at {w}: {max(v)-min(v):.3f}")
        res["windows"][wname] = entry

    a = res["windows"]["full_record"]; b = res["windows"]["common_1980_2019"]
    print("\n=== attribution ===")
    for w in WATCH:
        print(f"  {w}: spread {a[f'_spread{w}']:.3f} (full record) -> "
              f"{b[f'_spread{w}']:.3f} (common window); "
              f"estimator contrast for reference: 0.002")
    (RESULTS / "common_period.json").write_text(json.dumps(res, indent=2), encoding="utf8")
    print("\nSaved -> common_period.json")


if __name__ == "__main__":
    main()
