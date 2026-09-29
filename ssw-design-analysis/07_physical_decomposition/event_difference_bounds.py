#!/usr/bin/env python3
"""
event_difference_bounds.py
==========================
HOW MUCH CAN SSW EVENTS DIFFER IN THEIR FORCED ODDS OF A DOWNWARD OUTCOME, AND
HOW MUCH OF A SINGLE EVENT'S LABEL CAN THAT EXPLAIN?

Plan approved 2026-09-29 (referee concern: "the label does not identify a kind
of event" rests on within-ensemble tests; state what is proven -- the label of a
single outcome is mostly noise plus the shift -- and what is bounded -- forced
differences between events exist but are limited).

INPUTS (read, not recomputed)
  results/current/*/forced_variance_ceiling.json -> cmip6: mean_shift, var_pseudo,
  corrected_sigma_f2 (placebo-corrected forced between-event variance) and its
  95% interval (result I); outcome in member units, days +8..+52.

MODEL
  An event's forced shift m_e ~ N(m, sigma_f^2); an outcome Y = m_e + noise,
  noise ~ N(0, s^2) with s^2 = var_pseudo (event-free variance). Forced odds of a
  negative outcome P_e = Phi(-m_e / s) (condition 1 of the criterion).
  Reported for sigma_f at the point estimate and at the upper 95% bound:
    - the mean forced odds and the central 95% range of P_e across events
    - the share of the variance of a single event's label (Bernoulli) due to
      forced differences between events: Var(P_e) / (Pbar (1 - Pbar))
  Monte Carlo with 1,000,000 draws (seeded).

Output: results/current/5_mechanism/event_difference_bounds.json
"""
import json
import zlib
from pathlib import Path

import numpy as np
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "5_mechanism"
NAME = "event_difference_bounds"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
N_MC = 1_000_000


def main():
    src = next((ROOT / "results" / "current").glob("*/forced_variance_ceiling.json"))
    c = json.loads(src.read_text())["results"]["cmip6"]
    m, s = c["mean_shift"], float(np.sqrt(c["var_pseudo"]))
    s2f, (lo2, hi2) = c["corrected_sigma_f2"], c["corrected_CI95"]
    rng = np.random.default_rng(SEED)
    z = rng.standard_normal(N_MC)
    out = {"source": str(src.relative_to(ROOT)), "seed": SEED, "n_mc": N_MC,
           "mean_shift": m, "noise_sd": round(s, 4),
           "sigma_f2_point": s2f, "sigma_f2_CI95": [lo2, hi2], "cases": {}}
    for label, v in (("zero", 0.0), ("point", max(s2f, 0.0)), ("upper95", max(hi2, 0.0))):
        sf = float(np.sqrt(v))
        pe = norm.cdf(-(m + sf * z) / s)
        pbar = float(pe.mean())
        share = float(pe.var() / (pbar * (1 - pbar)))
        out["cases"][label] = {"sigma_f": round(sf, 4), "P_mean": round(pbar, 4),
                               "P_range95": [round(float(q), 3) for q in np.quantile(pe, [0.025, 0.975])],
                               "label_variance_share_forced": round(share, 4)}
        print(f"sigma_f {label:8s} {sf:.3f}: mean P {pbar:.3f}, 95% of events "
              f"{out['cases'][label]['P_range95']}, share of label variance forced {share:.3f}")
    (RESULTS / "event_difference_bounds.json").write_text(json.dumps(out, indent=2),
                                                          encoding="utf8", newline="\n")
    print("Saved -> event_difference_bounds.json")


if __name__ == "__main__":
    main()
