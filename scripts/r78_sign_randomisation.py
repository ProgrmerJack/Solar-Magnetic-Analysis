#!/usr/bin/env python3
"""Sign-randomisation inference test for SSW-avalanche event-level rate ratios.

For each of 100,000 iterations, randomly flips the sign of each event-level
log(RR), recomputes the geometric mean RR, and tests whether the observed
magnitude+direction pattern is unlikely under the null of no directional
SSW effect.

This test is more powerful than the binary sign test because it uses the
full magnitude distribution, not just the sign.

Output: data/results/r78_sign_randomisation.json
Referenced in: Supplementary Table 63
"""

import json
import numpy as np
from pathlib import Path

# 16 event-level rate ratios (SSW / matched-control)
RRS = [0.8821, 0.4419, 0.1835, 0.1113, 0.2718, 0.1894, 0.1144, 0.2084,
       0.3162, 0.0531, 0.1078, 0.6512, 0.8531, 0.3359, 1.0397, 3.4057]

N_PERM = 100_000
RNG_SEED = 42


def main():
    rng = np.random.default_rng(RNG_SEED)
    log_rrs = np.log(RRS)
    observed_gm = float(np.exp(np.mean(log_rrs)))
    observed_sign = int(sum(1 for r in RRS if r < 1))

    gm_null = np.empty(N_PERM)
    sign_null = np.empty(N_PERM, dtype=int)

    for i in range(N_PERM):
        signs = rng.choice([-1, 1], size=len(RRS))
        perm_log = log_rrs * signs
        gm_null[i] = np.exp(np.mean(perm_log))
        sign_null[i] = int(np.sum(perm_log < 0))

    p_gm = float(np.mean(gm_null <= observed_gm))
    p_sign = float(np.mean(sign_null >= observed_sign))

    results = {
        "description": "Sign-randomisation inference on 16 event-level log(RR)",
        "method": ("For each iteration, randomly flip the sign of each event "
                    "log(RR), compute geometric mean RR. Tests whether "
                    "magnitude+direction is unlikely under no-effect null."),
        "n_events": len(RRS),
        "n_permutations": N_PERM,
        "rng_seed": RNG_SEED,
        "observed_gmRR": round(observed_gm, 4),
        "observed_sign_count": observed_sign,
        "sign_randomisation_p_gmRR": round(p_gm, 6),
        "sign_randomisation_p_sign": round(p_sign, 6),
        "null_gmRR_median": round(float(np.median(gm_null)), 4),
        "null_gmRR_95CI": [
            round(float(np.percentile(gm_null, 2.5)), 4),
            round(float(np.percentile(gm_null, 97.5)), 4),
        ],
        "bonferroni_factor": 20,
        "bonferroni_p_gmRR": round(p_gm * 20, 6),
        "bonferroni_p_sign": round(p_sign * 20, 6),
        "survives_alpha_005_twosided": p_gm * 20 < 0.05,
    }

    outpath = Path(__file__).resolve().parent.parent / "data" / "results" / "r78_sign_randomisation.json"
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "w") as f:
        json.dump(results, f, indent=2)

    print(f"gmRR = {observed_gm:.4f}")
    print(f"P(gmRR) = {p_gm:.6f}")
    print(f"P(sign) = {p_sign:.6f}")
    print(f"Bonferroni P(gmRR) = {p_gm * 20:.4f}")
    print(f"Saved to {outpath}")


if __name__ == "__main__":
    main()
