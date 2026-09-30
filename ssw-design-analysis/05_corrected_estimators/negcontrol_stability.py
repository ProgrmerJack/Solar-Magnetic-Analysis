#!/usr/bin/env python3
"""
negcontrol_stability.py
=======================
The 43-event multi-index run reported the AAO negative control FAILING:
-0.66 at +45..+60, surviving BH-FDR at q=0.036 across the 32-test family.

That claim cannot be trusted as it stands, for two reasons:

  1. The seed was `hash("aao") % 10000`, which Python randomises per process.
  2. The four events the era-fix added (JAN 1963, JAN 1968, MAR 1969, JAN 1977)
     are ALL pre-1979, and the daily AAO record begins in 1979 -- so the AAO
     event set was IDENTICAL in the run that passed and the run that failed.
     The only thing that changed was the random draw.

This settles it: same data, same estimator, 60 fixed seeds, and 6000 replicates
(5x the previous 1200) so Monte Carlo error at p~0.05 falls from ~0.009 to
~0.004. If the negative control fails at most seeds it is a real estimator
defect and Gate 6 is broken. If it fails at a minority, the "failure" was
Monte Carlo error read as signal -- and so was the earlier "pass".

Output: negcontrol_stability.json
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
import multi_index_event_study as M

N_SEEDS = 60
M.N_BOOT = 6000
WATCH = "+45..+60"


def main():
    onsets_all = load_catalogue("primary")
    out = {"n_seeds": N_SEEDS, "n_boot": M.N_BOOT, "watch_bin": WATCH, "indices": {}}
    for name in ("aao", "pna"):
        d = M.load(name)
        on = onsets_all[(onsets_all >= d.index.min()) & (onsets_all <= d.index.max())]
        ps, effs = [], []
        for s in range(N_SEEDS):
            prof, _ = M.run(d, on, seed=s)
            ps.append(prof[WATCH]["p_two_sided"])
            effs.append(prof[WATCH]["effect"])
        a = np.array(ps)
        sig = int((a < 0.05).sum())
        out["indices"][name] = {
            "n_events": int(len(on)), "effect": round(float(effs[0]), 4),
            "effect_spread_across_seeds": round(float(np.ptp(effs)), 6),
            "p_min": round(float(a.min()), 4), "p_max": round(float(a.max()), 4),
            "p_median": round(float(np.median(a)), 4),
            "seeds_significant": sig, "n_seeds": N_SEEDS,
            "verdict": "STABLE" if sig in (0, N_SEEDS) else "SEED-DEPENDENT"}
        print(f"{name}: n={len(on)} effect {effs[0]:+.3f} "
              f"(point estimate varies {np.ptp(effs):.2e} across seeds -- it is not random)\n"
              f"      raw p range {a.min():.4f}..{a.max():.4f}, median {np.median(a):.4f}\n"
              f"      p<0.05 in {sig}/{N_SEEDS} seeds -> {out['indices'][name]['verdict']}",
              flush=True)
    (RESULTS / "negcontrol_stability.json").write_text(json.dumps(out, indent=2), encoding="utf8")
    print("\nSaved -> negcontrol_stability.json")


if __name__ == "__main__":
    main()
