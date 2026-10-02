#!/usr/bin/env python3
"""
snapsi_residual_quantiles.py
============================
IS THE IMPOSED SSW A TRANSLATION? RESIDUAL QUANTILES OF NUDGED AGAINST A SHIFTED
CONTROL, PER MODEL AND EVENT, AGAINST A DECLARED TOLERANCE.

Plan approved 2026-10-02 (user chose to act on the external audit, section 5).
Committed BEFORE it was run. The per-ensemble distributions were analysed before
only through pooled variance ratios, Kolmogorov-Smirnov tests and class rates
(snapsi_distribution_test.py, snapsi_regional_test.py); no residual quantile or
tail-probability error of a translated control has been computed.

QUESTION
  The paper calls the imposed SSW a translation of the surface distribution. A
  pooled variance ratio near 1 and a non-rejected shape test do not establish
  that. Test the approximation directly, in each nudged/control ensemble pair.

DATA
  Northern Hemisphere SNAPSI ensembles (nine centres x four initialisations, the
  set of snapsi_distribution_test.py; corrupt submissions excluded as there).
  Outcome A: polar-cap NAM proxy, days 8-25, standardised by the control mean and
  s.d. (as the distribution test). Outcome T: northern-Eurasian 2 m temperature,
  days 8-24, standardised the same way (snapsi_regional_test.region_series).

STATISTICS (per ensemble pair e)
  Delta_e = mean(nudged) - mean(control).
  R_e(q) = Q_nudged(q) - Q_control(q) - Delta_e, q = 0.05, 0.10, 0.25, 0.50, 0.75,
    0.90, 0.95 (sample quantiles, linear interpolation).
  E_e(q) = P(nudged < Q_control(q) + Delta_e) - q, q = 0.05, 0.10, 0.20: the error of
    the translated control in a tail probability.
  Uncertainty: members resampled within each arm (no matching exists), Delta
    re-estimated inside every replicate; pooled summaries (mean over pairs) with a
    two-stage bootstrap -- centres resampled, then members within arms -- 2,000
    replicates. Reported pooled, per event (February 2018, January 2019) and per
    centre.

TOLERANCE (declared now, newly chosen, before any residual was computed)
  |R(q)| <= 0.25 sigma and |E(q)| <= 0.05. Reason: moving a Gaussian by 0.25 sigma
  changes the probability below its 10th percentile by about 0.044 (0.25 x
  phi(1.28)); 0.05 is about a third of the 95% half-width of a frequency estimated
  from the 39 observed events, and a quarter of the imposed change in the downward
  rate (0.45 -> 0.84). Smaller departures would not be detectable in observations
  or matter for the class rates the paper discusses.

READING (fixed now; applied separately to A and to T)
  Approximate translation: the pooled 95% interval lies inside the tolerance at
    every q (R) and every tail q (E).
  Departure: some pooled interval excludes zero with its point estimate outside
    the tolerance.
  Unresolved: otherwise (an interval crosses a tolerance bound).
  A high p value is not evidence of equivalence; only intervals inside the
  tolerance are.

Output: results/current/8_experiment/snapsi_residual_quantiles.json
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import snapsi_distribution_test as D                 # noqa: E402
import snapsi_regional_test as SRT                   # noqa: E402

NAME = "snapsi_residual_quantiles"
RESULTS = HERE.parents[1] / "results" / "current" / "8_experiment"
SEED = 20261002
N_BOOT = 2000
QR = np.array([0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95])
QE = np.array([0.05, 0.10, 0.20])
TOL_R, TOL_E = 0.25, 0.05
EVENT = {"s20180125": "Feb-2018", "s20180208": "Feb-2018", "s20181213": "Jan-2019", "s20190108": "Jan-2019"}


def stats(N, C):
    d = N.mean() - C.mean()
    qc = np.quantile(C, QR)
    r = np.quantile(N, QR) - qc - d
    e = np.array([np.mean(N < np.quantile(C, q) + d) - q for q in QE])
    return r, e


def cases_A():
    centres = sorted({p.name.split("_")[0] for p in D.RED.glob("*.parquet")})
    out = []
    for c in centres:
        S = []
        for init in D.ONSET:
            n, ct = D.load(c, "nudged", init), D.load(c, "control", init)
            if n is not None and ct is not None:
                S.append(abs(D.window_means(n, init).mean() - D.window_means(ct, init).mean()))
        if S and max(S) < 1.0:
            continue                                   # corrupt submission, as the distribution test
        for init in sorted(D.ONSET):
            n, ct = D.load(c, "nudged", init), D.load(c, "control", init)
            if n is None or ct is None:
                continue
            nm, cm = D.window_means(n, init), D.window_means(ct, init)
            if len(nm) < 20 or len(cm) < 20:
                continue
            base, sd = float(cm.mean()), float(cm.std(ddof=1))
            out.append({"centre": c, "init": init, "N": -(nm - base) / sd, "C": -(cm - base) / sd})
    return out


def cases_T():
    centres = sorted({p.name.split("_")[0] for p in SRT.TAS.glob("*.npz")})
    out = []
    for c in centres:
        bad, _ = SRT.L.corruption_guard(c)
        if bad:
            continue
        for init in SRT.L.ONSET:
            tn, tc = SRT.region_series(c, "nudged", init), SRT.region_series(c, "control", init)
            if len(tn) < 20 or len(tc) < 20:
                continue
            x, y = tn["NEURASIA"].values, tc["NEURASIA"].values
            mu, sd = y.mean(), y.std(ddof=1)
            out.append({"centre": c, "init": init, "N": (x - mu) / sd, "C": (y - mu) / sd})
    return out


def summarise(cases, rng):
    est = [stats(k["N"], k["C"]) for k in cases]
    R, E = np.array([a for a, _ in est]), np.array([b for _, b in est])
    centres = sorted({k["centre"] for k in cases})
    byc = {c: [k for k in cases if k["centre"] == c] for c in centres}

    def boot(sub):
        bR, bE = [], []
        cs = sorted({k["centre"] for k in sub})
        by = {c: [k for k in sub if k["centre"] == c] for c in cs}
        for _ in range(N_BOOT):
            rr, ee = [], []
            for c in rng.choice(cs, len(cs)):
                for k in by[c]:
                    N = rng.choice(k["N"], len(k["N"])); C = rng.choice(k["C"], len(k["C"]))
                    a, b = stats(N, C); rr.append(a); ee.append(b)
            bR.append(np.mean(rr, axis=0)); bE.append(np.mean(ee, axis=0))
        return np.array(bR), np.array(bE)

    def block(sub):
        r = np.array([stats(k["N"], k["C"])[0] for k in sub]).mean(axis=0)
        e = np.array([stats(k["N"], k["C"])[1] for k in sub]).mean(axis=0)
        bR, bE = boot(sub)
        loR, hiR = np.percentile(bR, [2.5, 97.5], axis=0)
        loE, hiE = np.percentile(bE, [2.5, 97.5], axis=0)
        inside = bool(np.all(loR >= -TOL_R) and np.all(hiR <= TOL_R) and np.all(loE >= -TOL_E) and np.all(hiE <= TOL_E))
        dep = bool(np.any(((loR > 0) | (hiR < 0)) & (np.abs(r) > TOL_R)) or np.any(((loE > 0) | (hiE < 0)) & (np.abs(e) > TOL_E)))
        return {"n_pairs": len(sub), "n_centres": len({k["centre"] for k in sub}),
                "R": {f"{q:.2f}": [round(float(a), 4), round(float(b), 4), round(float(c), 4)] for q, a, b, c in zip(QR, r, loR, hiR)},
                "E": {f"{q:.2f}": [round(float(a), 4), round(float(b), 4), round(float(c), 4)] for q, a, b, c in zip(QE, e, loE, hiE)},
                "reading": "approximate translation" if inside else ("departure" if dep else "unresolved")}

    out = {"pooled": block(cases)}
    for ev in ("Feb-2018", "Jan-2019"):
        out[ev] = block([k for k in cases if EVENT[k["init"]] == ev])
    # registered but omitted in the first run (added 2026-10-02 after the code review)
    out["per_centre"] = {c: {k: v for k, v in block([k_ for k_ in cases if k_["centre"] == c]).items()}
                         for c in sorted({k["centre"] for k in cases})}
    out["per_pair"] = [{"centre": k["centre"], "init": k["init"], "nN": len(k["N"]), "nC": len(k["C"]),
                        "R": [round(float(x), 3) for x in stats(k["N"], k["C"])[0]],
                        "E": [round(float(x), 3) for x in stats(k["N"], k["C"])[1]]} for k in cases]
    return out


def main():
    rng = np.random.default_rng(SEED)
    res = {"plan_approved": "2026-10-02", "seed": SEED, "n_boot": N_BOOT, "q_R": QR.tolist(), "q_E": QE.tolist(),
           "tolerance": {"R_sigma": TOL_R, "E_probability": TOL_E},
           "format": "each entry [estimate, 2.5%, 97.5%]"}
    for lab, fn in (("A_polar_cap_NAM_days8_25", cases_A), ("T_northern_Eurasia_days8_24", cases_T)):
        cs = fn()
        res[lab] = summarise(cs, rng)
        print(lab, json.dumps({k: v for k, v in res[lab].items() if k != "per_pair"}), flush=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "snapsi_residual_quantiles.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print("Saved -> snapsi_residual_quantiles.json")


if __name__ == "__main__":
    main()
