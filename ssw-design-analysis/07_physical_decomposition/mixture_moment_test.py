#!/usr/bin/env python3
"""
mixture_moment_test.py
======================
IS FORCED MIXING OF TWO FIXED RESPONSE TYPES INCOMPATIBLE WITH THE ENSEMBLES?
A direct test that replaces the variance-mean curvature (class_law_general T3).

Method fixed 2026-10-07, committed (with its validation gate) BEFORE it was applied to
any model data. On the CMIP6 cases of class_law_general it is a labelled REANALYSIS of
data already examined (their T3 curvature had been computed); the confirmatory use is
the fresh-ensemble replication registered separately.

THE RESTRICTION. If forcing changes only the weight p_j of two fixed distributions G, H
(any shapes, widths), F_j = (1 - p_j) G + p_j H, then every linear functional of F_j is
linear in p_j. Hence for three forcing levels j = 0, 1, 2 the points (mean, second raw
moment) are collinear, and so are the CDF vectors at any thresholds:
  PRIMARY    T = (M2_1 - M2_0)(mu_2 - mu_0) - (M2_2 - M2_0)(mu_1 - mu_0) = 0
  SECONDARY  F_1 - F_0 and F_2 - F_0 parallel at the pooled 10th..90th percentiles
No ordering of the levels and no shape of G, H is assumed; years are pooled within a
level (an average of such mixtures is again one). One population translated with
constant variance gives T = -(mu_1 - mu_0)(mu_2 - mu_0)(mu_2 - mu_1). Reported
descriptively: rho = T / that value (0 under fixed-type mixing, 1 under pure
translation with constant variance). T is location-invariant; units are s.d. of level 0.
SCOPE: two fixed components, or any number along a linear path in the weights. It does
not exclude regimes in general or heterogeneous mechanisms.

INFERENCE. Two-way (pigeonhole) bootstrap, B = 1999: whole member trajectories
resampled (one draw per member id, applied to every level containing that id, so a run
continued from historical into SSP5-8.5 stays paired), crossed with 5-year blocks of
year indices resampled jointly across levels (calendar-aligned levels share volcanic
years). PRIMARY: basic-bootstrap two-sided p for T = 0. SECONDARY: for every direction
phi in [0, pi) the vector d(phi) = cos(phi)(F_1 - F_0) - sin(phi)(F_2 - F_0) is
studentised by its bootstrap s.d.; max-|t| with its bootstrap critical value gives
p(phi); mixing is rejected only if p(phi) < 0.05 for EVERY phi (intersection over the
unknown weights).
OUTCOMES per case: "incompatible with fixed two-type mixing" (p < 0.05), otherwise
"unresolved". Family: per variable, Holm-adjusted p reported beside raw p; the
variable-level statement uses Holm.

VALIDATION GATE (fixed before any model data; --validate). 1,000 simulations per null
case, B = 999 inside each, 30 years, members (25, 40, 25) unless stated, year effects
(s.d. 0.3, shared by the members of a level-year), within-member AR(1) noise rho = 0.3:
  N1 Gaussian types, separation 3, weights 0.1/0.4/0.9   N2 widths 1 and 2, sep. 2, 0.2/0.5/0.8
  N3 gamma(2) types, hard design 0.1/0.15/0.9            N4 lognormal types
  N5 t3 types                                            N6 weak separation 0.40/0.45/0.50
  N7 non-monotone weights 0.5/0.1/0.9                    N8 weights trending +-0.05 within levels
  N9 type persistent within a member (0.8 per year)      N10 ensembles (50, 50, 50), year effects 0.6
Power cases (one population; reported, no gate): P1 translation 0/1/3; P2 translation with
s.d. 1/1.1/1.3; P3 gamma population translated 0/1/3; P4 hard design 0/0.4/3;
P5 precipitation-like 0/0.5/1.5 with constant coefficient of variation.
GATE: for each test, in every null case, the one-sided 95% Clopper-Pearson upper bound
of the false-rejection rate at nominal 0.05 must lie below 0.07. A test that fails the
gate is not applied; that is reported.

Output: results/current/8_experiment/mixture_moment_test.json
Run:    mixture_moment_test.py --validate      (writes the validation block)
        mixture_moment_test.py                 (CMIP6 reanalysis; requires a passed gate)
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "8_experiment"
OUT = RESULTS / "mixture_moment_test.json"
sys.path.insert(0, str(HERE))

NAME = "mixture_moment_test"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
B_REAL, B_VAL, N_SIM, ALPHA, GATE = 1999, 999, 1000, 0.05, 0.07
QS = np.round(np.arange(0.1, 0.91, 0.1), 2)
PHI = np.linspace(0, np.pi, 180, endpoint=False)
BLOCK = 5


# ------------------------------------------------------------------ core
def prepare(arms):
    """arms: list of (member_ids, X[members x years] with NaN). Returns standardised
    arrays (units of level 0), union of ids and per-arm index into it."""
    X0 = arms[0][1]
    m0 = np.nanmean(X0); s0 = np.sqrt(np.nanmean(np.nanvar(X0, axis=0, ddof=1)))
    ids = sorted(set().union(*[set(a[0]) for a in arms]))
    pos = {k: i for i, k in enumerate(ids)}
    out = []
    for mid, X in arms:
        Z = (X - m0) / s0
        mask = np.isfinite(Z)
        out.append({"idx": np.array([pos[k] for k in mid]), "Z": np.where(mask, Z, 0.0), "mask": mask.astype(float)})
    nY = arms[0][1].shape[1]
    allz = np.concatenate([a["Z"][a["mask"] > 0] for a in out])
    thr = np.quantile(allz, QS)
    for a in out:
        a["I"] = np.stack([(a["Z"] <= t) * a["mask"] for t in thr])          # K x M x Y
    return out, len(ids), nY


def moments(prep, Wm, Wy):
    """Weighted moments for weight rows Wm (R x n_ids), Wy (R x Y): per arm
    mean, second raw moment, CDF vector (R x K)."""
    res = []
    for a in prep:
        wm = Wm[:, a["idx"]]
        n = np.einsum("rm,my,ry->r", wm, a["mask"], Wy)
        s1 = np.einsum("rm,my,ry->r", wm, a["Z"], Wy)
        s2 = np.einsum("rm,my,ry->r", wm, a["Z"] ** 2, Wy)
        F = np.einsum("rm,kmy,ry->rk", wm, a["I"], Wy) / n[:, None]
        res.append((s1 / n, s2 / n, F))
    return res


def weights(rng, n_ids, nY, B):
    Wm = np.stack([np.bincount(rng.integers(n_ids, size=n_ids), minlength=n_ids) for _ in range(B)]).astype(float)
    nb = int(np.ceil(nY / BLOCK)); starts = np.arange(nb) * BLOCK
    Wy = np.zeros((B, nY))
    for b in range(B):
        for st in rng.choice(starts, nb):
            Wy[b, st:min(st + BLOCK, nY)] += 1
    return Wm, Wy


def tests(arms, B, rng):
    prep, n_ids, nY = prepare(arms)
    one = moments(prep, np.ones((1, n_ids)), np.ones((1, nY)))
    Wm, Wy = weights(rng, n_ids, nY, B)
    bs = moments(prep, Wm, Wy)
    def T(m):
        (mu0, a0, _), (mu1, a1, _), (mu2, a2, _) = m
        return (a1 - a0) * (mu2 - mu0) - (a2 - a0) * (mu1 - mu0)
    t0, tb = float(T(one)[0]), T(bs)
    p_T = float((1 + np.sum(np.abs(tb - t0) >= abs(t0))) / (B + 1))
    mu = [float(m[0][0]) for m in one]
    trans = -(mu[1] - mu[0]) * (mu[2] - mu[0]) * (mu[2] - mu[1])
    rho = t0 / trans if abs(trans) > 1e-9 else float("nan")
    rho_b = tb / trans if abs(trans) > 1e-9 else np.full(B, np.nan)
    D1, D2 = one[1][2][0] - one[0][2][0], one[2][2][0] - one[0][2][0]
    D1b, D2b = bs[1][2] - bs[0][2], bs[2][2] - bs[0][2]
    pmax = 0.0
    for ph in PHI:
        d = np.cos(ph) * D1 - np.sin(ph) * D2
        db = np.cos(ph) * (D1b - D1) - np.sin(ph) * (D2b - D2)
        se = db.std(axis=0, ddof=1); se = np.where(se > 1e-9, se, 1e-9)
        tobs = np.max(np.abs(d) / se); tb_ = np.max(np.abs(db) / se, axis=1)
        pmax = max(pmax, float((1 + np.sum(tb_ >= tobs)) / (B + 1)))
    return {"T": round(t0, 5), "T_ci95": [round(float(x), 5) for x in np.quantile(tb, [0.025, 0.975])], "p_T": round(p_T, 4),
            "rho": round(float(rho), 4),
            "rho_ci95": [round(float(x), 4) for x in np.nanquantile(rho_b, [0.025, 0.975])] if np.isfinite(rho) else None,
            "p_cdf": round(pmax, 4), "means": [round(x, 4) for x in mu],
            "vars": [round(float(m[1][0] - m[0][0] ** 2), 4) for m in one]}


# ------------------------------------------------------------------ simulation
def sim_case(kind, rng, nY=30):
    ns = (50, 50, 50) if kind == "N10" else (25, 40, 25)
    yeff = 0.6 if kind == "N10" else 0.3
    def ar1(n, r=0.3):
        e = rng.normal(size=(n, nY)); x = np.empty_like(e); x[:, 0] = e[:, 0]
        for y in range(1, nY):
            x[:, y] = r * x[:, y - 1] + np.sqrt(1 - r * r) * e[:, y]
        return x
    comp = {"N1": ("norm", 3.0, (1, 1)), "N2": ("norm", 2.0, (1, 2)), "N3": ("gamma", 3.0, (1, 1)), "N4": ("lognorm", 3.0, (1, 1)),
            "N5": ("t3", 3.0, (1, 1)), "N6": ("norm", 2.0, (1, 1)), "N7": ("norm", 3.0, (1, 1)), "N8": ("norm", 3.0, (1, 1)),
            "N9": ("norm", 3.0, (1, 1)), "N10": ("norm", 3.0, (1, 1))}
    wts = {"N1": (0.1, 0.4, 0.9), "N2": (0.2, 0.5, 0.8), "N3": (0.1, 0.15, 0.9), "N4": (0.2, 0.5, 0.8), "N5": (0.2, 0.5, 0.8),
           "N6": (0.40, 0.45, 0.50), "N7": (0.5, 0.1, 0.9), "N8": (0.2, 0.5, 0.8), "N9": (0.2, 0.5, 0.8), "N10": (0.1, 0.4, 0.9)}
    arms = []
    if kind.startswith("N"):
        shape, D, (wa, wb) = comp[kind]
        def draw(z, typ):
            if shape == "norm":
                base = z
            elif shape == "gamma":
                from scipy.stats import norm, gamma
                base = (gamma.ppf(norm.cdf(z), 2) - 2) / np.sqrt(2)
            elif shape == "lognorm":
                base = (np.exp(0.6 * z) - np.exp(0.18)) / 0.7
            else:
                from scipy.stats import norm, t
                base = t.ppf(norm.cdf(z), 3) / np.sqrt(3)
            return np.where(typ, D + wb * base, wa * base)
        for j, n in enumerate(ns):
            p = np.full((n, nY), wts[kind][j])
            if kind == "N8":
                p = p + np.linspace(-0.05, 0.05, nY)[None, :]
            if kind == "N9":
                typ = np.zeros((n, nY), bool); typ[:, 0] = rng.random(n) < p[:, 0]
                for y in range(1, nY):
                    keep = rng.random(n) < 0.8
                    typ[:, y] = np.where(keep, typ[:, y - 1], rng.random(n) < p[:, y])
            else:
                typ = rng.random((n, nY)) < p
            X = draw(ar1(n), typ) + yeff * rng.normal(size=nY)[None, :]
            arms.append((np.array([f"{j}_{i}" for i in range(n)]), X))
    else:
        spec = {"P1": [(0, 1), (1, 1), (3, 1)], "P2": [(0, 1), (1, 1.1), (3, 1.3)], "P3": [(0, 1), (1, 1), (3, 1)],
                "P4": [(0, 1), (0.4, 1), (3, 1)], "P5": [(0, 1), (0.5, 1), (1.5, 1)]}[kind]
        for j, (n, (mu, sd)) in enumerate(zip(ns, spec)):
            z = ar1(n)
            if kind == "P3":
                from scipy.stats import norm, gamma
                z = (gamma.ppf(norm.cdf(z), 2) - 2) / np.sqrt(2)
            if kind == "P5":
                X = (1 + 0.25 * mu) * (4 + z) / 4 * 4 - 4      # constant CV around a mean of 4 + mu
            else:
                X = mu + sd * z
            X = X + 0.3 * rng.normal(size=nY)[None, :]
            arms.append((np.array([f"m{i}" for i in range(n)]), X))   # shared ids: continued runs
    return arms


def _val_job(args):
    kind, seed, n = args
    rng = np.random.default_rng(seed)
    rT = rC = 0
    for _ in range(n):
        r = tests(sim_case(kind, rng), B_VAL, rng)
        rT += r["p_T"] < ALPHA; rC += r["p_cdf"] < ALPHA
    return kind, rT, rC


def validate():
    kinds = [f"N{i}" for i in range(1, 11)] + [f"P{i}" for i in range(1, 6)]
    chunks = 10
    jobs = [(k, SEED + 1000 * i + c, N_SIM // chunks) for i, k in enumerate(kinds) for c in range(chunks)]
    import multiprocessing as mp
    with mp.get_context("fork").Pool(8) as pool:
        outs = pool.map(_val_job, jobs)
    agg = {}
    for k, a, b in outs:
        x = agg.setdefault(k, [0, 0]); x[0] += a; x[1] += b
    ub = lambda r: float(beta.ppf(0.95, r + 1, N_SIM - r)) if r < N_SIM else 1.0
    val = {"n_sim": N_SIM, "B": B_VAL, "alpha": ALPHA, "gate_upper95": GATE, "cases": {}}
    for k, (a, b) in agg.items():
        val["cases"][k] = {"T_reject": a / N_SIM, "T_upper95": round(ub(a), 4), "cdf_reject": b / N_SIM, "cdf_upper95": round(ub(b), 4)}
    nulls = [k for k in val["cases"] if k.startswith("N")]
    val["T_passes_gate"] = all(val["cases"][k]["T_upper95"] < GATE for k in nulls)
    val["cdf_passes_gate"] = all(val["cases"][k]["cdf_upper95"] < GATE for k in nulls)
    print(json.dumps(val, indent=1), flush=True)
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    res["validation"] = val; res["method_fixed"] = "2026-10-07, before any model data"
    OUT.write_text(json.dumps(res, indent=1), encoding="utf8", newline="\n")
    print(f"Saved -> {OUT.name} (validation)")
    return 0


# ------------------------------------------------------------------ CMIP6 reanalysis
def cmip6_cases():
    import class_law_general as CL
    d = pd.read_parquet(CL.DATA).copy()
    can = d.model == "CanESM5"
    d.loc[can, "model"] = "CanESM5_" + d.loc[can, "member"].str.extract(r"i\d+(p\d+)f\d+")[0]
    for m in ("CanESM5_p1", "CanESM5_p2", "MIROC6", "ACCESS-ESM1-5"):
        dm = d[d.model == m]
        for var in ("tas", "pr"):
            for season in ("JJA", "DJF"):
                for box in CL.BOXES:
                    s = CL.seasonal(dm, var, box, season)
                    arms = []
                    for exp, y0 in (CL.LEVEL_L[m], CL.LEVEL_M, CL.LEVEL_H):
                        a = s[(s.experiment == exp) & (s.syear >= y0) & (s.syear < y0 + 30)]
                        piv = a.pivot_table(index="member", columns="syear", values="mean").reindex(columns=range(y0, y0 + 30))
                        piv.loc[:, piv.notna().sum() < CL.MIN_MEMBERS] = np.nan
                        arms.append((piv.index.values, piv.values))
                    yield f"{m}|{var}|{box}_{season}", arms


def holm(ps):
    ps = np.asarray(ps); o = np.argsort(ps); adj = np.empty_like(ps); run = 0.0
    for r, i in enumerate(o):
        run = max(run, min(1.0, (len(ps) - r) * ps[i])); adj[i] = run
    return adj


def main():
    if "--validate" in sys.argv:
        return validate()
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    v = res.get("validation")
    if not v:
        sys.exit("run --validate first")
    res["cases"] = {}
    rng = np.random.default_rng(SEED)
    for key, arms in cmip6_cases():
        r = tests(arms, B_REAL, rng)
        res["cases"][key] = r
        print(key, {k: r[k] for k in ("T", "p_T", "rho", "p_cdf")}, flush=True)
    summ = {}
    for var in ("tas", "pr"):
        ks = [k for k in res["cases"] if f"|{var}|" in k]
        for test, ok in (("p_T", v["T_passes_gate"]), ("p_cdf", v["cdf_passes_gate"])):
            if not ok:
                summ[f"{var}|{test}"] = "not applied: failed the validation gate"
                continue
            ps = [res["cases"][k][test] for k in ks]; adj = holm(ps)
            for k, a in zip(ks, adj):
                res["cases"][k][test + "_holm"] = round(float(a), 4)
                res["cases"][k][test.replace("p_", "") + "_reading"] = ("incompatible with fixed two-type mixing" if a < ALPHA else "unresolved")
            summ[f"{var}|{test}"] = {"n_cases": len(ks), "n_incompatible_raw": int(sum(p < ALPHA for p in ps)),
                                     "n_incompatible_holm": int(sum(a < ALPHA for a in adj))}
    res["summary"] = summ
    res["note"] = "reanalysis of CMIP6 cases whose T3 curvature had been examined; method and gate fixed before this run"
    res["seed"] = SEED; res["B"] = B_REAL
    print(json.dumps(summ, indent=1))
    OUT.write_text(json.dumps(res, indent=1), encoding="utf8", newline="\n")
    print(f"Saved -> {OUT.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
