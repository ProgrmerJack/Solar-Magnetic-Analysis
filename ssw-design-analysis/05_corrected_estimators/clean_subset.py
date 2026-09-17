#!/usr/bin/env python3
"""
clean_subset.py
===============
The decisive estimate. Two INDEPENDENT design problems inflate the apparent
pre-onset AO anomaly (Fisher P=0.49 for their association, so they are not the
same confound wearing two hats):

  1. PRE-SATELLITE RECORD  pre-onset -1.74 pre-1980 vs -0.47 from 1980 (era_split);
                           the anomaly is significant in 7/8 catalogues on the
                           full record and 0/8 on 1980-2019 (common_period)
  2. EVENT CLUSTERING      nearest-onset binning files days that follow event A
                           into the pre-onset bins of event B. 14 of 43 events sit
                           43-105 d from a neighbour; those give pre-onset -2.75
                           against -0.55 for isolated events (isolated_events)

Removing both leaves events that are isolated AND in the satellite era. If the
pre-onset anomaly is a real precursor it should still be there. If it is the sum
of these two artefacts it should be gone, while the post-onset response -- which
neither artefact should touch -- stays put.

N_BOOT is 8000 here, not 1200: this subset is small, the claim is a null, and a
null asserted from a noisy interval is worthless. At p~0.05 that puts Monte Carlo
error near 0.003.

Output: clean_subset.json
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "2_event_study"
RESULTS.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "06_simulation_validation"))
from build_catalogue import load_catalogue
import canonical_event_study as K

K.N_BOOT = 8000
ISOLATION_DAYS = 120
SPLIT = 1980


def main():
    d = K.load_series("AO")
    on = load_catalogue("primary")
    on = on[(on >= d.index.min()) & (on <= d.index.max())]
    t = np.sort(on.values.astype("datetime64[D]").astype(int))
    gap = np.full(len(t), 10**6)
    gap[:-1] = np.minimum(gap[:-1], np.diff(t))
    gap[1:] = np.minimum(gap[1:], np.diff(t))
    iso = gap > ISOLATION_DAYS
    late = pd.DatetimeIndex(on).year >= SPLIT

    subsets = {
        "all": np.ones(len(t), bool),
        "isolated_only": iso,
        "satellite_only": late,
        "isolated_AND_satellite": iso & late,
    }
    res = {"n_boot": K.N_BOOT, "isolation_days": ISOLATION_DAYS,
           "split_year": SPLIT, "profiles": {}}
    for name, m in subsets.items():
        sub = pd.DatetimeIndex(pd.to_datetime(t[m], unit="D"))
        if len(sub) < 8:
            print(f"{name}: {len(sub)} events, skipped"); continue
        prof, dd = K.run(d, sub, seed=20260730)
        res["profiles"][name] = {"n_events": int(len(sub)),
                                 "n_winters": int(pd.Series(K.C.winter_of(sub)).nunique()),
                                 "profile": prof}
        star = "".join("*" if prof[l]["p_two_sided"] < 0.05 else "." for l in K.LABELS)
        vals = " ".join(f"{prof[l]['effect']:+.2f}" for l in K.LABELS)
        print(f"  {name:24s} n={len(sub):2d}  {vals}   {star}", flush=True)

    print(f"\n  bins: {' '.join(f'{l:>6s}' for l in K.LABELS)}")
    c = res["profiles"].get("isolated_AND_satellite")
    if c:
        p = c["profile"]
        print(f"\n=== cleanest subset (n={c['n_events']}, {c['n_winters']} winters) ===")
        for l in ("-30..-16", "-15..-1", "+15..+29"):
            e = p[l]
            print(f"  {l:9s} {e['effect']:+.3f}  95% CI [{e['CI95'][0]:+.3f}, "
                  f"{e['CI95'][1]:+.3f}]  p={e['p_two_sided']:.4f}")
    (RESULTS / "clean_subset.json").write_text(json.dumps(res, indent=2), encoding="utf8")
    print("\nSaved -> clean_subset.json")


if __name__ == "__main__":
    main()
