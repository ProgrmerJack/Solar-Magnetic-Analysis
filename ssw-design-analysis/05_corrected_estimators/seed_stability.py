#!/usr/bin/env python3
"""
seed_stability.py
=================
canonical_event_study.py, multi_index_event_study.py and marginal_events.py all
seeded the winter-block bootstrap with `hash(<python str>) % 10000`. Python
randomises string hashing per process (PYTHONHASHSEED), so every run drew a
DIFFERENT bootstrap sample and no reported p-value was reproducible.

This measures the consequence: the same data, the same estimator, 40 seeds.
If a conclusion moves across the 0.05 line between seeds, it was never a
finding -- it was Monte Carlo error being read as signal.

Output: seed_stability.json
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

N_SEEDS = 40
CASES = [("AO", "primary"), ("AO", "consensus_strict"), ("NAO", "primary")]
WATCH = ["-30..-16", "+15..+29"]


def main():
    res = {"n_seeds": N_SEEDS, "n_boot": K.N_BOOT, "cases": {}}
    for name, cat in CASES:
        d = K.load_series(name)
        on = load_catalogue(cat)
        on = on[(on >= d.index.min()) & (on <= d.index.max())]
        ps = {w: [] for w in WATCH}
        eff = {w: [] for w in WATCH}
        for s in range(N_SEEDS):
            prof, _ = K.run(d, on, seed=s)
            for w in WATCH:
                ps[w].append(prof[w]["p_two_sided"])
                eff[w].append(prof[w]["effect"])
        entry = {}
        print(f"\n=== {name} / {cat}  (n={len(on)}, {N_SEEDS} seeds x {K.N_BOOT} reps) ===")
        for w in WATCH:
            a = np.array(ps[w])
            flip = int((a < 0.05).sum())
            entry[w] = {"effect": round(float(eff[w][0]), 4),
                        "p_min": round(float(a.min()), 4),
                        "p_max": round(float(a.max()), 4),
                        "p_median": round(float(np.median(a)), 4),
                        "seeds_significant": flip, "n_seeds": N_SEEDS,
                        "verdict": "STABLE" if flip in (0, N_SEEDS) else "SEED-DEPENDENT"}
            print(f"  {w:9s} effect {eff[w][0]:+.3f}   p range "
                  f"{a.min():.4f}..{a.max():.4f} (median {np.median(a):.4f})   "
                  f"significant in {flip}/{N_SEEDS} seeds   {entry[w]['verdict']}")
        res["cases"][f"{name}|{cat}"] = entry
    (RESULTS / "seed_stability.json").write_text(json.dumps(res, indent=2), encoding="utf8")
    print("\nSaved -> seed_stability.json")


if __name__ == "__main__":
    main()
