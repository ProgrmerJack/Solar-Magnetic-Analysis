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
4. ASYMMETRIC MIXTURES (registered 2026-10-04 before it was computed). Fixed grid:
   separation Delta in {1.5, 2, 2.5, 3, 4}; width ratio of the second type
   r in {0.5, 0.75, 1, 1.5, 2}; threshold fixed at the starting distribution's
   quantile q in {0.3, 0.5, 0.7} or at the density minimum between the types;
   weight path w 0.5 -> 0.1, 0.5 -> 0.9, 0.2 -> 0.6. Elasticity of |D| on the total
   s.d. along each path. Reported per Delta: the maximum over the grid and the share
   of cells above 0.67 (the lower bound of the natural experiment's interval).
5. SHAPE-TEST POWER against asymmetric mixtures (seeded): reference weight 0.3,
   forced 0.6; designs (i) natural experiment: 8 ensembles of 20 forced + 30
   reference members, (ii) human forcing: 30 years of 55 + 60; Anderson-Darling
   k-sample statistic on groups standardised within unit, 200 within-unit
   permutations, 100 replicates; Delta in {1.5, 2, 3}, r in {0.5, 1, 2}.
   READING, fixed now: a (Delta, r) type pair is EXCLUDED by the natural experiment
   if every grid elasticity for it lies below 0.67, OR the shape test's power in
   design (i) is at least 0.8 (the experiment's shape test gave p = 0.62).

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


def mix2(w, D, r):
    return lambda x: (1 - w) * phi(x) + w * phi((x - D) / r) / r


def sd2(w, D, r):
    m = w * D
    return np.sqrt((1 - w) * 1 + w * (r * r + D * D) - m * m)


def contrast_any(pdf, c, lo=-40, hi=60):
    m0 = integrate.quad(pdf, lo, c, limit=300)[0]; m1 = integrate.quad(pdf, c, hi, limit=300)[0]
    return (integrate.quad(lambda x: x * pdf(x), lo, c, limit=300)[0] / m0
            - integrate.quad(lambda x: x * pdf(x), c, hi, limit=300)[0] / m1)


def quantile_of(pdf, q, lo=-15, hi=25):
    from scipy.optimize import brentq
    return brentq(lambda c: integrate.quad(pdf, lo, c, limit=300)[0] - q, lo, hi)


def density_min(pdf, D):
    from scipy.optimize import minimize_scalar
    r = minimize_scalar(pdf, bounds=(0, D), method="bounded")
    return float(r.x)


def asymmetric_grid():
    out = {}
    for D in (1.5, 2.0, 2.5, 3.0, 4.0):
        els = []
        per_r = {}
        for r in (0.5, 0.75, 1.0, 1.5, 2.0):
            n0 = len(els)
            for w0, w1 in ((0.5, 0.1), (0.5, 0.9), (0.2, 0.6)):
                start = mix2(w0, D, r)
                cs = [quantile_of(start, q) for q in (0.3, 0.5, 0.7)] + [density_min(start, D)]
                for c in cs:
                    ws = np.linspace(w0, w1, 15)
                    lnD = [np.log(abs(contrast_any(mix2(w, D, r), c))) for w in ws]
                    lns = [np.log(sd2(w, D, r)) for w in ws]
                    els.append(float(np.polyfit(lns, lnD, 1)[0]))
            per_r[str(r)] = round(float(max(els[n0:])), 3)
        els = np.asarray(els)
        out[str(D)] = {"max": round(float(els.max()), 3), "median": round(float(np.median(els)), 3),
                       "share_above_0.67": round(float((els > 0.67).mean()), 3), "n_cells": int(len(els)),
                       "max_by_width_ratio": per_r}
    return out


def shape_power(seed=20261004):
    from scipy.stats import anderson_ksamp
    import warnings
    warnings.filterwarnings("ignore")
    rng = np.random.default_rng(seed)
    def draw(n, w, D, r):
        z = rng.random(n) < w
        return np.where(z, D + r * rng.normal(size=n), rng.normal(size=n))
    def test(units):
        zs = [((a - a.mean()) / a.std(ddof=1), (b - b.mean()) / b.std(ddof=1)) for a, b in units]
        A = np.concatenate([a for a, _ in zs]); B = np.concatenate([b for _, b in zs])
        ad0 = anderson_ksamp([A, B]).statistic; ex = 0
        for _ in range(200):
            a_, b_ = [], []
            for x, y in zs:
                z = np.concatenate([x, y]); p = rng.permutation(len(z)); a_.append(z[p[:len(x)]]); b_.append(z[p[len(x):]])
            ex += anderson_ksamp([np.concatenate(a_), np.concatenate(b_)]).statistic >= ad0
        return (ex + 1) / 201
    out = {}
    for design, (nu, nf, nr) in (("natural_experiment", (8, 20, 30)), ("human_forcing", (30, 55, 60))):
        for D in (1.5, 2.0, 3.0):
            for r in (0.5, 1.0, 2.0):
                rej = [test([(draw(nf, 0.6, D, r), draw(nr, 0.3, D, r)) for _ in range(nu)]) < 0.05 for _ in range(100)]
                out[f"{design}_D{D}_r{r}"] = round(float(np.mean(rej)), 3)
    return out


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
    if True:                                     # sections 4-5 (registered 2026-10-04)
        res["asymmetric_grid"] = asymmetric_grid()
        print("asymmetric grid:", res["asymmetric_grid"], flush=True)
        res["shape_power"] = shape_power()
        print("shape power:", res["shape_power"], flush=True)
        ex = {}
        for D in ("1.5", "2.0", "3.0"):
            for r in ("0.5", "1.0", "2.0"):
                grid_ok = res["asymmetric_grid"][D]["max_by_width_ratio"][r] < 0.67
                pw = res["shape_power"][f"natural_experiment_D{D}_r{r}"]
                ex[f"D{D}_r{r}"] = {"excluded": bool(grid_ok or pw >= 0.8), "grid_max_below_0.67": bool(grid_ok),
                                    "shape_power": pw}
        res["exclusion_by_natural_experiment"] = ex
    OUT.write_text(json.dumps(res, indent=1), encoding="utf8", newline="\n")
    print(json.dumps(res, indent=1))
    print(f"Saved -> {OUT.name}")


if __name__ == "__main__":
    main()
