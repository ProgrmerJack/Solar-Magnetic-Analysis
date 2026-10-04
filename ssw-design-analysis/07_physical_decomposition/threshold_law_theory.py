#!/usr/bin/env python3
"""
threshold_law_theory.py
=======================
The mathematics behind the threshold-class law, computed rather than quoted
(numbers used in snapsi_free_natural_experiment F3's registered reading and in the
paper's Methods; derivation from research note theory_spread_law.md).

1. LAW. For X = mu + sigma Z with a fixed shape and a class threshold c, the contrast
   D = E[X | X < c] - E[X | X >= c] = sigma * g((c - mu) / sigma). For the Gaussian,
   g(t) = -phi(t) / (Phi(t) (1 - Phi(t))). Reported: g(0), g(+-1), the curvature
   coefficient k in |g(t)| ~ |g(0)| (1 + k t^2), and the relative change of the
   contrast when a pure mean shift moves the class rate from 0.46 to 0.85 -- the
   rate moves at first order, the contrast only at second order.
2. DISCRIMINATOR. Two response types: X = (1-w) N(0,1) + w N(Delta,1), threshold
   midway (c = Delta/2), forcing changes the weight w. Elasticity of |D| with respect
   to the total s.d. sigma(w) = sqrt(1 + w(1-w) Delta^2), fitted (log-log slope) along
   w = 0.5 -> 0.1. One population gives exactly 1.
3. RESIDUAL of the location-scale null under a mixture (forced contrast minus the
   reference contrast mapped by location and scale), for w 0.3 -> 0.7.
Output: results/current/8_experiment/threshold_law_theory.json
"""
import json
from pathlib import Path

import numpy as np
from scipy import integrate, stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "current" / "8_experiment" / "threshold_law_theory.json"
phi, Phi = stats.norm.pdf, stats.norm.cdf


def g(t):
    return -phi(t) / (Phi(t) * (1 - Phi(t)))


def contrast(pdf, c):
    m0 = integrate.quad(pdf, -40, c, limit=200)[0]
    m1 = integrate.quad(pdf, c, 40, limit=200)[0]
    return (integrate.quad(lambda x: x * pdf(x), -40, c, limit=200)[0] / m0
            - integrate.quad(lambda x: x * pdf(x), c, 40, limit=200)[0] / m1)


def mix(w, D):
    return lambda x: (1 - w) * phi(x) + w * phi(x - D)


def sd(w, D):
    return np.sqrt(1 + w * (1 - w) * D * D)


def elasticity(D):
    ws = np.linspace(0.1, 0.5, 21)
    lnD = [np.log(abs(contrast(mix(w, D), D / 2))) for w in ws]
    lns = [np.log(sd(w, D)) for w in ws]
    return float(np.polyfit(lns, lnD, 1)[0])


def ls_residual(w0, w1, D):
    c = D / 2
    mu0, s0, mu1, s1 = w0 * D, sd(w0, D), w1 * D, sd(w1, D)
    return contrast(mix(w1, D), c) - contrast(mix(w0, D), mu0 + s0 * (c - mu1) / s1) * s1 / s0


def main():
    t = np.linspace(-0.5, 0.5, 101)
    k = float(np.polyfit(t ** 2, np.abs(g(t)) / abs(g(0.0)) - 1, 1)[0])
    t0, t1 = stats.norm.ppf(0.46), stats.norm.ppf(0.85)        # rate = Phi(t)
    res = {"g0": round(float(g(0.0)), 4), "g_pm1": round(float(g(1.0)), 4),
           "curvature_k": round(k, 4),
           "contrast_change_rate_0.46_to_0.85": round(float(g(t1) / g(t0) - 1), 4),
           "rate_change_per_sigma_at_median": round(float(phi(0.0)), 4),
           "mixture_elasticity": {str(D): round(elasticity(D), 3) for D in (1.5, 2.0, 2.5, 3.0, 4.0, 6.0)},
           "one_population_elasticity": 1.0,
           "ls_residual_w0.3_to_0.7": {str(D): round(float(ls_residual(0.3, 0.7, D)), 4) for D in (1.5, 2.0, 3.0, 4.0)}}
    OUT.write_text(json.dumps(res, indent=1), encoding="utf8", newline="\n")
    print(json.dumps(res, indent=1))
    print(f"Saved -> {OUT.name}")


if __name__ == "__main__":
    main()
