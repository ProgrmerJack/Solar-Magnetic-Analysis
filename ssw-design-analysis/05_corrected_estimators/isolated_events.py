#!/usr/bin/env python3
"""
isolated_events.py
==================
Tests whether the pre-onset AO anomaly is a real precursor or leakage between
neighbouring events.

THE MECHANISM
  Bins are assigned by NEAREST onset. Consider a day 20 d before event B that is
  also 25 d after event A. |-20| < |25|, so it is filed in the -30..-16
  PRE-onset bin -- while the AO on that day is depressed because of A's POST-onset
  response. The estimator then reports A's response as B's precursor.

  This is the same failure mode as the three contamination traps already found in
  this project: the damage is in what the comparison contains, not the treatment.

  43 events fall in 36 winters, so 7 winters hold more than one event, and
  within-winter gaps are short enough for exactly this overlap.

THE TEST
  ISOLATED   no other catalogue event within +/-120 d, so the whole -60..+60
             window is free of any other event's -60..+60 window
  CLUSTERED  the remainder

  If the precursor is real, it survives in the isolated set -- those events have
  no neighbour to borrow a response from.
  If it is leakage, it weakens or vanishes in the isolated set and concentrates
  in the clustered set.

  Post-onset response is reported for both as a control: it should be similar in
  each, since it is not the bin at risk.

Output: isolated_events.json
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

ISOLATION_DAYS = 120
WATCH = ["-30..-16", "-15..-1", "+0..+14", "+15..+29"]


def main():
    d = K.load_series("AO")
    on = load_catalogue("primary")
    on = on[(on >= d.index.min()) & (on <= d.index.max())]
    t = np.sort(on.values.astype("datetime64[D]").astype(int))

    gap = np.full(len(t), 10**6)
    if len(t) > 1:
        gap[:-1] = np.minimum(gap[:-1], np.diff(t))
        gap[1:] = np.minimum(gap[1:], np.diff(t))
    iso = gap > ISOLATION_DAYS

    print(f"{len(t)} primary events; nearest-neighbour gap: "
          f"min {gap.min()} d, median {int(np.median(gap))} d")
    print(f"  ISOLATED  (>{ISOLATION_DAYS} d from any other event): {int(iso.sum())}")
    print(f"  CLUSTERED (<={ISOLATION_DAYS} d):                     {int((~iso).sum())}")
    print(f"  clustered gaps: {sorted(gap[~iso].tolist())}")

    res = {"isolation_days": ISOLATION_DAYS, "n_boot": K.N_BOOT,
           "n_isolated": int(iso.sum()), "n_clustered": int((~iso).sum()),
           "nearest_gap_days": gap.tolist(), "profiles": {}}

    for label, mask in (("isolated", iso), ("clustered", ~iso), ("all", np.ones(len(t), bool))):
        sub = pd.DatetimeIndex(pd.to_datetime(t[mask], unit="D"))
        if len(sub) < 8:
            print(f"  {label}: only {len(sub)} events, skipped")
            continue
        prof, _ = K.run(d, sub, seed=zlib.crc32(label.encode()) % 10000)
        res["profiles"][label] = {"n_events": int(len(sub)), "profile": prof}
        print(f"\n  {label:10s} n={len(sub):2d}  " + "  ".join(
            f"{w} {prof[w]['effect']:+.3f}(p={prof[w]['p_two_sided']:.3f})" for w in WATCH))

    p = res["profiles"]
    if "isolated" in p and "clustered" in p:
        print("\n=== verdict ===")
        for w in WATCH:
            a = p["isolated"]["profile"][w]
            b = p["clustered"]["profile"][w]
            print(f"  {w:9s} isolated {a['effect']:+.3f} (p={a['p_two_sided']:.3f})   "
                  f"clustered {b['effect']:+.3f} (p={b['p_two_sided']:.3f})")
    (RESULTS / "isolated_events.json").write_text(json.dumps(res, indent=2), encoding="utf8")
    print("\nSaved -> isolated_events.json")


if __name__ == "__main__":
    main()
