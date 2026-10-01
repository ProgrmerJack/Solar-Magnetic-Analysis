#!/usr/bin/env python3
"""
response_shape_dose.py
======================
IS THE STATE-DEPENDENT SHAPE AFTER SSWs A SECOND POPULATION, OR ONE POPULATION
WHOSE FORCED DOSE VARIES?

Plan approved 2026-10-01 (user: "I do approve any new kind of analyses"); design
from research_notes/SSW paper Nature Geo strengthening/D_state_dependent_shape.md
section 4 (analyses A and C). Committed before any of the fits below was run.

WHAT HAD BEEN SEEN
  The registered comparison (class_model_comparison.py, result T-d): with
  information up to onset (tier P2), a two-component model with logistic
  membership (M2) beats a continuous Gaussian with state-dependent mean and
  log-spread (M1) after SSWs by more than on ordinary days (CRPS G = 0.00124
  against a pseudo mean of -0.0023, p 0.01), with fitted components 0.49 sigma
  apart (0.71 on ordinary days). Computed from those stored parameters (research
  note D, section 0.1, and re-checked 2026-10-01): the fitted mixture has one mode
  at every membership weight, after SSWs and on ordinary days. Nothing with
  post-onset predictors has been fitted in this framework.

QUESTION
  Published idealised and model work (White et al. 2020; Hitchcock & Simpson 2014)
  makes the surface response follow the realised lower-stratospheric anomaly after
  onset, a continuous "dose". If events differ in dose, one population with a
  linear dose response has a mean AND spread that vary with the state known at
  onset (the dose is only partly predictable then), which a two-component model can
  fit better than M1. If instead there are two kinds of response, knowing the dose
  should sharpen membership and M2's advantage should persist or grow.

DATA: exactly those of class_model_comparison.py (1,517 CMIP6 events, 20 members,
  outcome days +8..+52, 200 size-matched pseudo-onset sets, leave-one-model-out).

S1 (primary)   the registered comparison with predictors = tier P2 plus the
               realised 100 hPa zonal wind after onset (days 0-15 and 15-30, at
               60N and polar-cap mean: u100_*_p00_15, u100_*_p15_30), identically on
               the SSWs and on every pseudo set. G = score(M1) - score(M2), CRPS
               primary, ignorance secondary; p = P(G_pseudo >= G_ssw). Positive and
               negative controls as registered (planted two-to-one regimes 1 sigma
               apart, membership logistic in the first predictor), with this
               predictor set.
S2 (secondary) the same with all post-onset predictors (tier 3: 10, 50, 100 hPa).
S3 (descriptive) M2 fitted to all SSWs with P2 predictors: the component means,
               spreads and which component has the lower mean; number of modes of
               the fitted mixture over membership weights 0.01..0.99.
S4 (descriptive) quantile shift function: quantiles 0.05..0.95 of the SSW
               responses minus those of the pooled pseudo-onset responses (all 200
               sets), raw and after removing each sample's mean; model-cluster
               bootstrap (2,000) intervals.
S5 (descriptive) is the dose itself bimodal? The polar-cap 100 hPa wind days 0-30
               (mean of the two windows, standardised as the predictors) after SSWs:
               BIC of one against two Gaussian components, and the number of modes
               of the fitted two-component density.

READING (fixed now)
  S1 p > 0.05 and G_ssw - mean(G_pseudo) less than half its registered value
     (0.00354): the state-dependent shape is accounted for by the realised dose --
     one population whose forced shift varies continuously.
  S1 p <= 0.05: shape beyond the dose remains; reported as such.
  Otherwise: inconclusive. The positive-control detection rate is reported with
  it; below 0.2 the test is called underpowered.
  This is a mechanism check, not a forecast: the dose is observed after onset and
  overlaps the outcome window, on the SSWs and on the ordinary days alike.

POST-HOC POWER CHECK (--power; added 2026-10-01 AFTER the result above was seen,
labelled post hoc in the JSON). The registered controls draw outcomes as noise
around the calendar mean, with no relation to the dose; with the dose among the
predictors their G is not comparable to the real pseudo-onset G (both control
rates came out 1.0), so they do not measure power here. Instead: synthetic outcomes
= the in-sample ridge mean of the real SSW responses on the same predictors and
calendar + Gaussian noise with the residual s.d.; negative, nothing else;
positive, regimes planted exactly as in the registered control (two thirds / one
third, means 1 residual s.d. apart, membership logistic in the first predictor).
Power = share of positive datasets whose G exceeds the 95th percentile of the
negative datasets' G (100 each). Written into the existing JSON under
"posthoc_power"; registered numbers untouched.

Output: results/current/6_predictability/response_shape_dose.json
"""
import json
import multiprocessing as mp
import sys
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.stats import norm
from sklearn.linear_model import LinearRegression
from sklearn.mixture import GaussianMixture

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import class_model_comparison as CMC                 # noqa: E402

W, FV = CMC.W, CMC.FV
NAME = "response_shape_dose"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
RESULTS = CMC.RESULTS
DOSE = ["u100_60N_p00_15", "u100_60N_p15_30", "u100_cap_p00_15", "u100_cap_p15_30"]
QS = np.round(np.arange(0.05, 0.951, 0.05), 2)
N_BOOT = 2000


def n_modes(w, m1, s1, m2, s2):
    x = np.linspace(min(m1, m2) - 6 * max(s1, s2), max(m1, m2) + 6 * max(s1, s2), 40001)
    f = w * norm.pdf(x, m1, s1) + (1 - w) * norm.pdf(x, m2, s2)
    return int(np.sum((f[1:-1] > f[:-2]) & (f[1:-1] > f[2:])))


def run_set(label, ci, X, cal_ci, y, g, mdl, usable, M):
    Z = FV.std_within(X[:, ci], g)
    S, _ = CMC.evaluate(Z, X[:, cal_ci], y, g, mdl, SEED)
    real = CMC.summary(S)
    print(f"{label} SSW: {real}", flush=True)
    ok = np.isfinite(y) & np.isfinite(X[:, cal_ci]).all(1)
    f0 = LinearRegression().fit(X[ok][:, cal_ci], y[ok])
    CMC._S.clear()
    CMC._S.update(M=M, members=[f.stem for f in usable], tier_ci=ci, cal_ci=cal_ci,
                  ctrl=(X[:, ci], X[:, cal_ci], g, mdl, f0.predict(np.nan_to_num(X[:, cal_ci])),
                        float(np.std(y[ok] - f0.predict(X[ok][:, cal_ci]), ddof=1))))
    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(CMC.WORKERS, mp_context=ctx) as ex:
        pseudo = list(ex.map(CMC._pseudo_task, range(CMC.K_PSEUDO)))
        controls = list(ex.map(CMC._control_task, [(j, p) for p in (1, 0) for j in range(CMC.N_CONTROL)]))
    pos, neg = controls[:CMC.N_CONTROL], controls[CMC.N_CONTROL:]
    out = {"n_predictors": len(ci), "ssw": real}
    for sc in ("crps", "ign"):
        Gp = np.array([p[sc]["G_M1_minus_M2"] for p in pseudo])
        q95 = np.percentile(Gp, 95)
        out[sc] = {"G_ssw": real[sc]["G_M1_minus_M2"], "G_pseudo_mean": round(float(Gp.mean()), 5),
                   "G_pseudo_q025_q975": [round(float(q), 5) for q in np.percentile(Gp, [2.5, 97.5])],
                   "G_ssw_minus_pseudo": round(real[sc]["G_M1_minus_M2"] - float(Gp.mean()), 5),
                   "p_pseudo_ge_ssw": round(float(np.mean(Gp >= real[sc]["G_M1_minus_M2"])), 4),
                   "positive_detection_rate": round(float(np.mean([c[sc]["G_M1_minus_M2"] >= q95 for c in pos])), 3),
                   "negative_false_rate": round(float(np.mean([c[sc]["G_M1_minus_M2"] >= q95 for c in neg])), 3),
                   "skill_M1_vs_M0_ssw": real[sc]["skill_M1_vs_M0"], "skill_M2_vs_M0_ssw": real[sc]["skill_M2_vs_M0"],
                   "skill_M1_vs_M0_pseudo_mean": round(float(np.mean([p[sc]["skill_M1_vs_M0"] for p in pseudo])), 4),
                   "skill_M2_vs_M0_pseudo_mean": round(float(np.mean([p[sc]["skill_M2_vs_M0"] for p in pseudo])), 4)}
        print(f"{label} {sc}: {out[sc]}", flush=True)
    return out


def main():
    rng = np.random.default_rng(SEED)
    files = sorted(W.RAW.glob("*_zm.nc"))
    Xs, ys, gs, cols, usable = [], [], [], None, []
    for f in files:
        c = W.prepare_member(f)
        if c is None:
            continue
        F = W.H.features_for(c["u"], c["plev"], c["lat"], c["idx"], c["on"])
        Y = W.C6.anom(c["am"], c["on"], c["cl"], W.OUT_WIN)
        Xs.append(F.values); ys.append(Y); gs.append(np.full(len(Y), f.stem))
        cols = list(F.columns); usable.append(f)
    X = np.concatenate(Xs); y = np.concatenate(ys); g = np.concatenate(gs)
    mdl = np.array([s.split("_")[0] for s in g])
    cal_ci = [cols.index("doy_sin"), cols.index("doy_cos")]
    p2 = [cols.index(c) for c in W.P.tier_cols(cols, FV.TIERS[CMC.TIER])]
    dose = [cols.index(c) for c in DOSE]
    t3 = [cols.index(c) for c in W.P.tier_cols(cols, 3)]
    res = {"plan_approved": "2026-10-01", "seed": SEED, "n_boot": N_BOOT, "n_events": int(np.isfinite(y).sum()),
           "n_members": len(usable), "k_pseudo": CMC.K_PSEUDO, "n_control": CMC.N_CONTROL,
           "registered_reference": {"G_ssw_minus_pseudo_crps": 0.00354, "p": 0.01},
           "dose_columns": DOSE}
    W.K_DRAWS = CMC.K_PSEUDO
    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(CMC.WORKERS, mp_context=ctx) as ex:
        M = {stem: (Fk, Yk) for stem, Fk, Yk, _ in ex.map(W.matched_draws, list(enumerate(map(str, usable))))}

    # S3: component means and modes, P2, all SSWs
    Zp2 = FV.std_within(X[:, p2], g)
    okk = np.isfinite(y) & np.isfinite(Zp2).all(1)
    fit = CMC.em_mixture(Zp2[okk], X[okk][:, cal_ci], y[okk], np.random.default_rng(zlib.crc32(f"{NAME}|S3".encode())))
    _, coef, s1, s2, v, _ = fit
    cal_mean = X[okk][:, cal_ci].mean(0) @ coef[2:]
    m1, m2 = coef[0] + cal_mean, coef[1] + cal_mean
    pi = 1 / (1 + np.exp(-(np.column_stack([np.ones(okk.sum()), Zp2[okk]]) @ v)))
    lower = 1 if m1 < m2 else 2
    res["S3_components_P2"] = {
        "mean_1": round(float(m1), 4), "sd_1": round(float(s1), 4), "mean_2": round(float(m2), 4), "sd_2": round(float(s2), 4),
        "mean_weight_1": round(float(pi.mean()), 4), "lower_mean_component": lower,
        "lower_component_sd": round(float(s1 if lower == 1 else s2), 4),
        "max_modes_over_weights": max(n_modes(w, m1, s1, m2, s2) for w in np.linspace(0.01, 0.99, 99))}
    print("S3", res["S3_components_P2"], flush=True)

    # S4: quantile shift function, SSW minus pooled pseudo
    yp = {m: np.concatenate([M[m][1][k] for k in range(CMC.K_PSEUDO)]) for m in M}
    ys_m = {m: y[(g == m) & np.isfinite(y)] for m in M}
    um = sorted({s.split("_")[0] for s in M})

    def qsf(models, demean):
        a = np.concatenate([ys_m[m] for m in M if m.split("_")[0] in models])
        b = np.concatenate([yp[m][np.isfinite(yp[m])] for m in M if m.split("_")[0] in models])
        if demean:
            a, b = a - a.mean(), b - b.mean()
        return np.quantile(a, QS) - np.quantile(b, QS)
    r4 = np.random.default_rng(zlib.crc32(f"{NAME}|S4".encode()))
    res["S4_quantile_shift"] = {"q": QS.tolist()}
    for dm in (False, True):
        est = qsf(um, dm)
        bb = []
        for _ in range(N_BOOT):
            pick = list(r4.choice(um, len(um), replace=True))
            a = np.concatenate([ys_m[m] for p_ in pick for m in M if m.split("_")[0] == p_])
            b = np.concatenate([yp[m][np.isfinite(yp[m])] for p_ in pick for m in M if m.split("_")[0] == p_])
            if dm:
                a, b = a - a.mean(), b - b.mean()
            bb.append(np.quantile(a, QS) - np.quantile(b, QS))
        lo, hi = np.percentile(np.array(bb), [2.5, 97.5], axis=0)
        res["S4_quantile_shift"]["demeaned" if dm else "raw"] = {
            "estimate": [round(float(x), 4) for x in est], "ci95_lo": [round(float(x), 4) for x in lo],
            "ci95_hi": [round(float(x), 4) for x in hi]}
    print("S4", res["S4_quantile_shift"], flush=True)

    # S5: is the dose bimodal?
    cap = [cols.index("u100_cap_p00_15"), cols.index("u100_cap_p15_30")]
    d = FV.std_within(X[:, cap], g).mean(1)
    d = d[np.isfinite(d)][:, None]
    g1 = GaussianMixture(1, random_state=0).fit(d)
    g2 = GaussianMixture(2, n_init=10, random_state=0).fit(d)
    mu, sd, w = g2.means_.ravel(), np.sqrt(g2.covariances_.ravel()), g2.weights_
    res["S5_dose_unimodality"] = {"n": int(len(d)), "bic_1": round(float(g1.bic(d)), 2), "bic_2": round(float(g2.bic(d)), 2),
                                  "delta_bic_2_minus_1": round(float(g2.bic(d) - g1.bic(d)), 2),
                                  "two_comp_means": [round(float(x), 3) for x in mu],
                                  "two_comp_sds": [round(float(x), 3) for x in sd],
                                  "two_comp_weights": [round(float(x), 3) for x in w],
                                  "modes_of_fitted_two_comp": n_modes(w[0], mu[0], sd[0], mu[1], sd[1])}
    print("S5", res["S5_dose_unimodality"], flush=True)

    res["S1_P2_plus_dose"] = run_set("S1", p2 + dose, X, cal_ci, y, g, mdl, usable, M)
    res["S2_tier3"] = run_set("S2", t3, X, cal_ci, y, g, mdl, usable, M)
    s1 = res["S1_P2_plus_dose"]["crps"]
    half = s1["G_ssw_minus_pseudo"] < 0.5 * 0.00354
    res["reading"] = ("accounted for by the dose" if (s1["p_pseudo_ge_ssw"] > 0.05 and half) else
                      "shape beyond the dose remains" if s1["p_pseudo_ge_ssw"] <= 0.05 else "inconclusive")
    print("READING:", res["reading"], flush=True)
    (RESULTS / "response_shape_dose.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print(f"Saved -> {NAME}.json")


_P = {}


def _power_task(args):
    j, planted = args
    Z, cal, g, mdl, mu, sd = _P["d"]
    rng = np.random.default_rng(zlib.crc32(f"{NAME}|power|{planted}|{j}".encode()))
    y = mu + rng.normal(0, sd, len(mu))
    if planted:
        z1 = (Z[:, 0] - Z[:, 0].mean()) / Z[:, 0].std()
        lo, hi = -10.0, 10.0
        for _ in range(60):
            c = 0.5 * (lo + hi)
            lo, hi = (c, hi) if np.mean(1 / (1 + np.exp(-(c + z1)))) < 2 / 3 else (lo, c)
        k = rng.uniform(size=len(y)) < 1 / (1 + np.exp(-(c + z1)))
        y = y + np.where(k, -1.0 / 3, 2.0 / 3) * sd
    S, _ = CMC.evaluate(Z, cal, y, g, mdl, zlib.crc32(f"{NAME}|power-eval|{planted}|{j}".encode()) % (2 ** 31))
    return CMC.summary(S)["crps"]["G_M1_minus_M2"]


def power():
    from sklearn.linear_model import RidgeCV
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    out_p = RESULTS / "response_shape_dose.json"
    res = json.loads(out_p.read_text())
    files = sorted(W.RAW.glob("*_zm.nc"))
    Xs, ys, gs, cols = [], [], [], None
    for f in files:
        c = W.prepare_member(f)
        if c is None:
            continue
        F = W.H.features_for(c["u"], c["plev"], c["lat"], c["idx"], c["on"])
        Xs.append(F.values); ys.append(W.C6.anom(c["am"], c["on"], c["cl"], W.OUT_WIN))
        gs.append(np.full(len(ys[-1]), f.stem)); cols = list(F.columns)
    X = np.concatenate(Xs); y = np.concatenate(ys); g = np.concatenate(gs)
    mdl = np.array([s_.split("_")[0] for s_ in g])
    cal_ci = [cols.index("doy_sin"), cols.index("doy_cos")]
    ci = [cols.index(c) for c in W.P.tier_cols(cols, FV.TIERS[CMC.TIER])] + [cols.index(c) for c in DOSE]
    Z = FV.std_within(X[:, ci], g)
    ok = np.isfinite(y) & np.isfinite(Z).all(1) & np.isfinite(X[:, cal_ci]).all(1)
    Z, cal, y, g, mdl = Z[ok], X[ok][:, cal_ci], y[ok], g[ok], mdl[ok]
    XA = np.column_stack([Z, cal])
    m = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 4, 25))).fit(XA, y)
    mu = m.predict(XA); sd = float(np.std(y - mu, ddof=1))
    _P["d"] = (Z, cal, g, mdl, mu, sd)
    with ProcessPoolExecutor(CMC.WORKERS, mp_context=mp.get_context("fork")) as ex:
        G = list(ex.map(_power_task, [(j, p) for p in (1, 0) for j in range(CMC.N_CONTROL)]))
    pos, neg = np.array(G[:CMC.N_CONTROL]), np.array(G[CMC.N_CONTROL:])
    q95 = float(np.percentile(neg, 95))
    real = res["S1_P2_plus_dose"]["crps"]["G_ssw"]
    res["posthoc_power"] = {
        "label": "post hoc, added after the S1 result was seen",
        "residual_sd": round(sd, 4), "in_sample_r2": round(float(1 - np.var(y - mu) / np.var(y)), 4),
        "G_negative_mean": round(float(neg.mean()), 5), "G_negative_q95": round(q95, 5),
        "G_positive_mean": round(float(pos.mean()), 5),
        "power_planted_1sd_regimes": round(float(np.mean(pos > q95)), 3),
        "real_G_ssw": real, "real_G_percentile_in_negative": round(float(np.mean(neg <= real)), 3)}
    pw = res["posthoc_power"]["power_planted_1sd_regimes"]
    res["verdict"] = (f"INCONCLUSIVE: the registered criterion is met but the test is underpowered "
                      f"(post-hoc power {pw} against planted regimes, below the registered 0.2)") if pw < 0.2 else \
        f"registered reading stands (post-hoc power {pw})"
    res = {"verdict": res.pop("verdict"), **res}    # first, so the catalogue shows it
    print("POWER", res["posthoc_power"], res["verdict"], flush=True)
    out_p.write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")


if __name__ == "__main__":
    power() if "--power" in sys.argv else main()
