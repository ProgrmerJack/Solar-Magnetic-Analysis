#!/usr/bin/env python3
"""
class_model_comparison.py
=========================
ONE SHIFTED POPULATION, A CONTINUOUSLY STATE-DEPENDENT ONE, OR TWO REGIMES?

Plan approved 2026-09-30. Written before any of these models was fitted.

QUESTION
  Which statistical description of the surface response after an SSW predicts
  held-out events best: one shifted distribution, a distribution whose centre and
  spread vary continuously with the event's state, or two regimes ("downward" and
  "not") whose membership depends on the event's state? The class reading of the
  literature is the third; our reading is the first (with small continuous
  departures). Proper scores on held-out data decide, on the class model's own
  terms, rather than a threshold applied after the outcome.

A. CMIP6 (1,517 Charlton-Polvani events, 20 members, 10 models; the data, outcome
   and predictors of forecast_value_test.py: annular mode days +8..+52 in member
   sigma units, not demeaned; predictors of tier P2, information up to onset,
   standardised within member, plus day-of-year sine and cosine)
   M0 SHIFT       Gaussian; mean = calendar regression; spread = residual s.d.
                  (forecast_value_test F0).
   M1 CONTINUOUS  Gaussian; mean linear in predictors + calendar (ridge, as F1);
                  log spread linear in the same predictors (fitted by maximum
                  likelihood on the ridge residuals, L2 penalty chosen by inner
                  grouped cross-validation).
   M2 REGIMES     two Gaussian components with their own constant means (plus the
                  calendar terms, shared) and variances; membership probability
                  logistic in the predictors (EM, 20 random starts, best likelihood).
   M2b TWO POPULATIONS, state-blind: M2 with constant membership probability.
   Cross-validation: leave one MODEL out (as forecast_value_test). Scores per
   event: CRPS (closed form for Gaussians and Gaussian mixtures) and the ignorance
   (negative log predictive density). Skill of each model against M0; the class
   model's gain over the continuous one, G = score(M1) - score(M2), with a
   model-cluster bootstrap interval (2,000). The identical pipeline on 200
   size-matched sets of event-free pseudo-onsets (forecast_value_test draws).
   Primary readout: G_SSW - mean(G_pseudo); p = P(G_pseudo >= G_SSW), CRPS
   primary, ignorance secondary; and M2b against M0 at SSWs relative to pseudo.
   Positive control: synthetic outcomes with a planted two-regime structure
   (two thirds / one third, regime means 1.0 sigma apart, membership logistic in
   the first predictor with slope 1) with the real predictors; detection rate of
   G > 0 at p < 0.05 over 100 datasets. Negative control: the same with no regimes.

B. SNAPSI (36 Northern Hemisphere ensembles per arm, nine models; members' days
   +8..+25 polar-cap NAM proxy standardised by the control, as in the
   distribution test). In each ensemble the forcing is shared, so two kinds of
   response would appear as a mixture within the nudged ensemble. Five-fold
   cross-validated mean log predictive density of a two-component Gaussian mixture
   minus that of one Gaussian, per ensemble; nudged minus control on paired
   ensembles, with a 10,000-resample centre bootstrap. Positive control: control
   members with two thirds of them displaced by 1.0 sigma; detection rate over
   100 resamples.

Reading, fixed now: the class reading is supported if M2 beats M1 after SSWs by
more than on ordinary days (p < 0.05) or if the SNAPSI mixture gain is larger with
the SSW imposed (interval excluding zero). If neither, and the positive controls
detect the planted structure, the two-regime description adds nothing detectable.

IMPLEMENTATION NOTES (fixed before the run)
  EM for M2/M2b is generalised EM: each M-step does weighted least squares for the
  component intercepts and the shared calendar terms, the variance update (floor
  0.05), and one penalised Newton step (L2 1.0) for the logistic membership; at
  most 150 iterations, tolerance 1e-6 in mean log-likelihood. M1's log-spread
  model is fitted to the ridge's out-of-fold residuals inside the training
  models, its L2 penalty chosen from {1, 10, 100, 1000} by inner grouped 5-fold
  held-out likelihood. Positive/negative controls use the real predictors and
  calendar; the synthetic outcome is the M0 calendar mean plus (positive) the
  regime offset, with noise s.d. equal to the real M0 residual s.d.; "detected"
  means G at least the 95th percentile of the real pseudo-onset G distribution.
  SNAPSI: sklearn GaussianMixture (2 components, 5 initialisations) against one
  Gaussian, 5-fold cross-validation over members; the positive control's
  detection uses 1,000 centre-bootstrap resamples per planted dataset.

POST-HOC DIAGNOSTICS (added 2026-09-30 AFTER the registered result was seen: M2
beat M1 after SSWs by more than on ordinary days, p 0.01, while the negative
control's false-positive rate was 0.13-0.15 against a nominal 0.05). Labelled as
post hoc everywhere they are reported; the registered numbers are unchanged.
  PH1 M1s ONE SKEWED POPULATION: M1's mean and spread with a constant skew-normal
      shape (fitted by maximum likelihood to the training residuals), so that the
      comparison M1s vs M2 separates "two regimes" from "one skewed population".
      Scored like the others (CRPS from a tabulated standard skew-normal curve);
      same pseudo-onset reference.
  PH2 the registered G calibrated against the negative-control G distribution.
  PH4 skill of each model against M0 after SSWs with a 2,000-resample model-cluster
      bootstrap interval, and the same skill on the 200 pseudo-onset sets, so that
      "knowing the event adds skill" is judged against ordinary winter days.
  PH3 M2 fitted to all SSW events and to the pseudo-onset sets: component mean
      separation, spreads, mean weight and the spread of the membership
      probability across events; skewness of the responses.

Output: results/current/6_predictability/class_model_comparison.json
"""
import json
import multiprocessing as mp
import sys
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm
from sklearn.linear_model import LinearRegression, RidgeCV
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import GroupKFold, KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "6_predictability"
sys.path.insert(0, str(HERE))
import within_model_check as W                       # noqa: E402
import forecast_value_test as FV                     # noqa: E402
import snapsi_distribution_test as D                 # noqa: E402

NAME = "class_model_comparison"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
K_PSEUDO = 200
N_CONTROL = 100
N_STARTS = 20
N_BOOT = 2000
N_BOOT_SNAPSI = 10000
N_PLANT = 100
WORKERS = 16
TIER = "P2 at-onset"
EM_ITERS, EM_TOL, VAR_FLOOR, L2_LOGIT = 150, 1e-6, 0.05, 1.0


# ------------------------------------------------------------------ scores
def _A(m, s):
    z = m / s
    return m * (2 * norm.cdf(z) - 1) + 2 * s * norm.pdf(z)


def mix_scores(y, w, mu, sd):
    """CRPS and ignorance of Gaussian mixtures; w, mu, sd: (n, K)."""
    crps = (w * _A(y[:, None] - mu, sd)).sum(1)
    K = w.shape[1]
    for i in range(K):
        for j in range(K):
            crps -= 0.5 * w[:, i] * w[:, j] * _A(mu[:, i] - mu[:, j], np.sqrt(sd[:, i] ** 2 + sd[:, j] ** 2))
    dens = (w * norm.pdf(y[:, None], mu, sd)).sum(1)
    return crps, -np.log(np.maximum(dens, 1e-300))


# ------------------------------------------------------------------ M1 spread
def _nll_logsd(theta, Z, r, alpha):
    a, b = theta[0], theta[1:]
    ls = a + Z @ b
    e = np.exp(-2 * ls)
    f = np.sum(ls + 0.5 * r ** 2 * e) + alpha * b @ b
    g_ls = 1 - r ** 2 * e
    return f, np.concatenate([[g_ls.sum()], Z.T @ g_ls + 2 * alpha * b])


def fit_logsd(Z, r, alpha):
    th0 = np.zeros(Z.shape[1] + 1); th0[0] = np.log(np.std(r))
    return minimize(_nll_logsd, th0, args=(Z, r, alpha), jac=True, method="L-BFGS-B").x


def logsd_alpha(Z, r, g):
    best, ba = -np.inf, 10.0
    for alpha in (1.0, 10.0, 100.0, 1000.0):
        ll = 0.0
        for a, b in GroupKFold(5).split(Z, r, g):
            th = fit_logsd(Z[a], r[a], alpha)
            ls = th[0] + Z[b] @ th[1:]
            ll += np.sum(norm.logpdf(r[b], 0, np.exp(ls)))
        if ll > best:
            best, ba = ll, alpha
    return ba


# ------------------------------------------------------------------ PH1 skew-normal
from scipy.stats import skewnorm                     # noqa: E402
from scipy.optimize import minimize_scalar           # noqa: E402

_GX = np.linspace(-14, 14, 5601)


def sn_crps_table(alpha):
    """CRPS of the standard skew-normal (location 0, scale 1, shape alpha) at z,
    as an interpolating function: CRPS(z) = int_{-inf}^z F^2 + int_z^inf (1-F)^2."""
    F = skewnorm.cdf(_GX, alpha)
    dx = _GX[1] - _GX[0]
    A = np.concatenate([[0], np.cumsum(0.5 * (F[1:] ** 2 + F[:-1] ** 2) * dx)])
    Bc = np.concatenate([[0], np.cumsum(0.5 * ((1 - F[1:]) ** 2 + (1 - F[:-1]) ** 2) * dx)])
    B = Bc[-1] - Bc
    tab = A + B
    return lambda z: np.interp(z, _GX, tab)


def sn_params(mu, sd, alpha):
    """Location and scale giving mean mu and s.d. sd for shape alpha."""
    d = alpha / np.sqrt(1 + alpha ** 2)
    om = sd / np.sqrt(1 - 2 * d ** 2 / np.pi)
    return mu - om * d * np.sqrt(2 / np.pi), om


def fit_alpha(r, sd):
    def nll(a):
        xi, om = sn_params(0.0, sd, a)
        return -np.sum(skewnorm.logpdf(r, a, xi, om))
    return float(minimize_scalar(nll, bounds=(-15, 15), method="bounded").x)


# ------------------------------------------------------------------ M2 EM
def em_mixture(Z, cal, y, rng, logistic=True):
    n, p = Z.shape
    Zc = np.column_stack([np.ones(n), Z])
    best = None
    for _ in range(N_STARTS):
        r = rng.uniform(0.2, 0.8) * np.ones(n) + rng.normal(0, 0.2, n)
        r = np.clip(r + 0.5 * (y < np.median(y)) * rng.choice([-1, 1]), 0.02, 0.98)
        v = np.zeros(p + 1)
        ll_old = -np.inf
        for it in range(EM_ITERS):
            # M-step: WLS on stacked data
            Xs = np.vstack([np.column_stack([np.ones(n), np.zeros(n), cal]),
                            np.column_stack([np.zeros(n), np.ones(n), cal])])
            ws = np.concatenate([r, 1 - r]) + 1e-9
            ys = np.concatenate([y, y])
            Wsq = np.sqrt(ws)
            coef = np.linalg.lstsq(Xs * Wsq[:, None], ys * Wsq, rcond=None)[0]
            mu1 = coef[0] + cal @ coef[2:]; mu2 = coef[1] + cal @ coef[2:]
            s1 = np.sqrt(max(np.sum(r * (y - mu1) ** 2) / r.sum(), VAR_FLOOR))
            s2 = np.sqrt(max(np.sum((1 - r) * (y - mu2) ** 2) / (1 - r).sum(), VAR_FLOOR))
            if logistic:
                pi = 1 / (1 + np.exp(-(Zc @ v)))
                grad = Zc.T @ (r - pi) - L2_LOGIT * np.r_[0, v[1:]]
                H = (Zc * (pi * (1 - pi))[:, None]).T @ Zc + L2_LOGIT * np.diag(np.r_[0, np.ones(p)])
                v = v + np.linalg.solve(H + 1e-8 * np.eye(p + 1), grad)
                pi = 1 / (1 + np.exp(-(Zc @ v)))
            else:
                pi = np.full(n, np.clip(r.mean(), 1e-3, 1 - 1e-3))
            # E-step
            d1 = pi * norm.pdf(y, mu1, s1); d2 = (1 - pi) * norm.pdf(y, mu2, s2)
            tot = np.maximum(d1 + d2, 1e-300)
            r = d1 / tot
            ll = float(np.mean(np.log(tot)))
            if abs(ll - ll_old) < EM_TOL:
                break
            ll_old = ll
        if best is None or ll > best[0]:
            best = (ll, coef.copy(), s1, s2, v.copy(), float(r.mean()))
    return best


def em_predict(fit, Z, cal, logistic=True):
    _, coef, s1, s2, v, rbar = fit
    n = len(Z)
    pi = 1 / (1 + np.exp(-(np.column_stack([np.ones(n), Z]) @ v))) if logistic else np.full(n, rbar)
    mu = np.column_stack([coef[0] + cal @ coef[2:], coef[1] + cal @ coef[2:]])
    return np.column_stack([pi, 1 - pi]), mu, np.column_stack([np.full(n, s1), np.full(n, s2)])


# ------------------------------------------------------------------ one data set
def evaluate(X, cal, y, g, model, seed):
    """Leave-one-model-out scores of M0, M1, M2, M2b for every event."""
    rng = np.random.default_rng(seed)
    ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
    X, cal, y, g, model = X[ok], cal[ok], y[ok], g[ok], model[ok]
    S = {m: {"crps": np.full(len(y), np.nan), "ign": np.full(len(y), np.nan)} for m in ("M0", "M1", "M2", "M2b", "M1s")}
    for mdl in np.unique(model):
        te, tr = model == mdl, model != mdl
        f0 = LinearRegression().fit(cal[tr], y[tr])
        sd0 = np.std(y[tr] - f0.predict(cal[tr]), ddof=1)
        mu0 = f0.predict(cal[te])
        c, i_ = mix_scores(y[te], np.ones((te.sum(), 1)), mu0[:, None], np.full((te.sum(), 1), sd0))
        S["M0"]["crps"][te], S["M0"]["ign"][te] = c, i_
        # M1
        XA = np.column_stack([X, cal])
        f1 = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 4, 25)))
        oof = np.full(tr.sum(), np.nan)
        for a, b in GroupKFold(5).split(XA[tr], y[tr], g[tr]):
            oof[b] = f1.fit(XA[tr][a], y[tr][a]).predict(XA[tr][b])
        mu1 = f1.fit(XA[tr], y[tr]).predict(XA[te])
        sc = StandardScaler().fit(X[tr])
        Ztr, Zte = sc.transform(X[tr]), sc.transform(X[te])
        r = y[tr] - oof
        alpha = logsd_alpha(Ztr, r, g[tr])
        th = fit_logsd(Ztr, r, alpha)
        sd1 = np.exp(th[0] + Zte @ th[1:])
        c, i_ = mix_scores(y[te], np.ones((te.sum(), 1)), mu1[:, None], sd1[:, None])
        S["M1"]["crps"][te], S["M1"]["ign"][te] = c, i_
        # PH1 (post hoc): the same mean and spread, one constant skew-normal shape
        sd_tr = np.exp(th[0] + Ztr @ th[1:])
        al = fit_alpha(r, sd_tr)
        xi, om = sn_params(mu1, sd1, al)
        S["M1s"]["crps"][te] = om * sn_crps_table(al)((y[te] - xi) / om)
        S["M1s"]["ign"][te] = -skewnorm.logpdf(y[te], al, xi, om)
        # M2, M2b
        for nm, lg in (("M2", True), ("M2b", False)):
            fit = em_mixture(Ztr, cal[tr], y[tr], rng, logistic=lg)
            w, mu, sd = em_predict(fit, Zte, cal[te], logistic=lg)
            c, i_ = mix_scores(y[te], w, mu, sd)
            S[nm]["crps"][te], S[nm]["ign"][te] = c, i_
    return S, model


def summary(S):
    out = {}
    for sc in ("crps", "ign"):
        m = {k: float(np.nanmean(v[sc])) for k, v in S.items()}
        out[sc] = {"mean": {k: round(v, 5) for k, v in m.items()},
                   "G_M1_minus_M2": round(m["M1"] - m["M2"], 5),
                   "G_M0_minus_M2b": round(m["M0"] - m["M2b"], 5),
                   "PH1_G_M1s_minus_M2": round(m["M1s"] - m["M2"], 5),
                   "PH1_G_M1_minus_M1s": round(m["M1"] - m["M1s"], 5),
                   "skill_M1_vs_M0": round(1 - m["M1"] / m["M0"], 4) if sc == "crps" else round(m["M0"] - m["M1"], 5),
                   "skill_M2_vs_M0": round(1 - m["M2"] / m["M0"], 4) if sc == "crps" else round(m["M0"] - m["M2"], 5),
                   "skill_M1s_vs_M0": round(1 - m["M1s"] / m["M0"], 4) if sc == "crps" else round(m["M0"] - m["M1s"], 5)}
    return out


_S = {}


def _pseudo_task(k):
    ci, cal_ci = _S["tier_ci"], _S["cal_ci"]
    X = np.concatenate([_S["M"][m][0][k][:, ci] for m in _S["members"]])
    cal = np.concatenate([_S["M"][m][0][k][:, cal_ci] for m in _S["members"]])
    y = np.concatenate([_S["M"][m][1][k] for m in _S["members"]])
    g = np.concatenate([np.full(_S["M"][m][1].shape[1], m) for m in _S["members"]])
    mdl = np.array([s.split("_")[0] for s in g])
    S, _ = evaluate(FV.std_within(X, g), cal, y, g, mdl, SEED + 1 + k)
    return summary(S)


def _control_task(args):
    j, planted = args
    X, cal, g, mdl, mu0, sd0 = _S["ctrl"]
    rng = np.random.default_rng(SEED + 100000 + 1000 * planted + j)
    Z = FV.std_within(X, g)
    y = mu0 + rng.normal(0, sd0, len(mu0))
    if planted:
        z1 = (Z[:, 0] - np.nanmean(Z[:, 0])) / np.nanstd(Z[:, 0])
        # intercept so that the mean membership is two thirds
        lo, hi = -10.0, 10.0
        for _ in range(60):
            c = 0.5 * (lo + hi)
            if np.nanmean(1 / (1 + np.exp(-(c + z1)))) < 2 / 3:
                lo = c
            else:
                hi = c
        pin = 1 / (1 + np.exp(-(c + z1)))
        k = rng.uniform(size=len(y)) < pin
        y = y + np.where(k, -1.0 / 3, 2.0 / 3) * 1.0 * sd0   # means 1 sd apart, overall mean kept
    S, _ = evaluate(Z, cal, y, g, mdl, SEED + 200000 + 1000 * planted + j)
    return summary(S)


# ------------------------------------------------------------------ SNAPSI
def snapsi_cases():
    centres = sorted({p.name.split("_")[0] for p in D.RED.glob("*.parquet")})
    cases = []
    for c in centres:
        S = []
        for init in D.ONSET:
            n, ct = D.load(c, "nudged", init), D.load(c, "control", init)
            if n is None or ct is None:
                continue
            S.append(abs(D.window_means(n, init).mean() - D.window_means(ct, init).mean()))
        if S and max(S) < 1.0:
            continue
        for init in sorted(D.ONSET):
            n, ct = D.load(c, "nudged", init), D.load(c, "control", init)
            if n is None or ct is None:
                continue
            nm, cm = D.window_means(n, init), D.window_means(ct, init)
            if len(nm) < 20 or len(cm) < 20:
                continue
            base, sd = float(cm.mean()), float(cm.std(ddof=1))
            cases.append({"centre": c, "init": init, "N": -(nm - base) / sd, "C": -(cm - base) / sd})
    return cases


def cv_gain(x, seed):
    """5-fold CV mean log density: 2-component mixture minus one Gaussian."""
    x = np.asarray(x, float).reshape(-1, 1)
    g = 0.0
    for a, b in KFold(5, shuffle=True, random_state=seed).split(x):
        mu, sd = x[a].mean(), x[a].std(ddof=1)
        l1 = norm.logpdf(x[b, 0], mu, sd)
        gm = GaussianMixture(2, n_init=5, random_state=seed).fit(x[a])
        l2 = gm.score_samples(x[b])
        g += float(np.sum(l2 - l1))
    return g / len(x)


def centre_boot(diffs, centres, rng, n):
    uc = np.unique(centres)
    by = {c: diffs[centres == c] for c in uc}
    out = np.empty(n)
    for i in range(n):
        pick = rng.choice(uc, len(uc), replace=True)
        out[i] = np.mean(np.concatenate([by[c] for c in pick]))
    return out


def snapsi_arm(rng):
    cases = snapsi_cases()
    rows = []
    for k, cs in enumerate(cases):
        rows.append({"centre": cs["centre"], "init": cs["init"],
                     "gain_nudged": cv_gain(cs["N"], SEED + k), "gain_control": cv_gain(cs["C"], SEED + k)})
    d = np.array([r["gain_nudged"] - r["gain_control"] for r in rows])
    cen = np.array([r["centre"] for r in rows])
    bs = centre_boot(d, cen, rng, N_BOOT_SNAPSI)
    det = 0
    for j in range(N_PLANT):
        dp = []
        for k, cs in enumerate(cases):
            C = cs["C"].copy()
            m = rng.permutation(len(C))[: int(round(2 * len(C) / 3))]
            C[m] -= 1.0
            dp.append(cv_gain(C, SEED + 7919 * j + k) - rows[k]["gain_control"])
        b = centre_boot(np.array(dp), cen, rng, 1000)
        det += np.percentile(b, 2.5) > 0
    return {"n_ensembles": len(rows), "per_ensemble": [{**r, "gain_nudged": round(r["gain_nudged"], 4),
                                                        "gain_control": round(r["gain_control"], 4)} for r in rows],
            "mean_gain_nudged": round(float(np.mean([r["gain_nudged"] for r in rows])), 4),
            "mean_gain_control": round(float(np.mean([r["gain_control"] for r in rows])), 4),
            "nudged_minus_control": round(float(d.mean()), 4),
            "ci95_centre_bootstrap": [round(float(q), 4) for q in np.percentile(bs, [2.5, 97.5])],
            "positive_control_detection_rate": round(det / N_PLANT, 3)}


# ------------------------------------------------------------------ main
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
    tier_ci = [cols.index(c) for c in W.P.tier_cols(cols, FV.TIERS[TIER])]
    Z = FV.std_within(X[:, tier_ci], g)
    S, m_te = evaluate(Z, X[:, cal_ci], y, g, mdl, SEED)
    real = summary(S)
    print(f"CMIP6 {int(np.isfinite(y).sum())} events: {real}", flush=True)
    # model-cluster bootstrap of G
    ok = np.isfinite(S["M0"]["crps"])
    um = np.unique(m_te[ok]) if False else np.unique(mdl[np.isfinite(y) & np.isfinite(Z).all(1)])
    mm = mdl[np.isfinite(y) & np.isfinite(Z).all(1)]
    boot = {sc: [] for sc in ("crps", "ign")}
    for _ in range(N_BOOT):
        pick = rng.choice(um, len(um), replace=True)
        idx = np.concatenate([np.flatnonzero(mm == u) for u in pick])
        for sc in ("crps", "ign"):
            boot[sc].append(np.nanmean(S["M1"][sc][idx]) - np.nanmean(S["M2"][sc][idx]))
    res = {"plan_approved": "2026-09-30", "seed": SEED, "tier": TIER, "n_events": int(np.isfinite(y).sum()),
           "n_members": len(usable), "k_pseudo": K_PSEUDO, "n_starts": N_STARTS, "cv": "leave one model out",
           "ssw": real, "G_ci95_model_bootstrap": {sc: [round(float(q), 5) for q in np.percentile(v, [2.5, 97.5])]
                                                    for sc, v in boot.items()}}
    W.K_DRAWS = K_PSEUDO
    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(WORKERS, mp_context=ctx) as ex:
        M = {stem: (Fk, Yk) for stem, Fk, Yk, _ in ex.map(W.matched_draws, list(enumerate(map(str, usable))))}
    ok = np.isfinite(y) & np.isfinite(X[:, cal_ci]).all(1)
    f0 = LinearRegression().fit(X[ok][:, cal_ci], y[ok])
    _S.update(M=M, members=[f.stem for f in usable], tier_ci=tier_ci, cal_ci=cal_ci,
              ctrl=(X[:, tier_ci], X[:, cal_ci], g, mdl, f0.predict(np.nan_to_num(X[:, cal_ci])),
                    float(np.std(y[ok] - f0.predict(X[ok][:, cal_ci]), ddof=1))))
    with ProcessPoolExecutor(WORKERS, mp_context=ctx) as ex:
        pseudo = list(ex.map(_pseudo_task, range(K_PSEUDO)))
        controls = list(ex.map(_control_task, [(j, p) for p in (1, 0) for j in range(N_CONTROL)]))
    pos, neg = controls[:N_CONTROL], controls[N_CONTROL:]
    res["pseudo"], res["controls"] = {}, {}
    for sc in ("crps", "ign"):
        Gp = np.array([p[sc]["G_M1_minus_M2"] for p in pseudo])
        Bp = np.array([p[sc]["G_M0_minus_M2b"] for p in pseudo])
        q95, q95b = np.percentile(Gp, 95), np.percentile(Bp, 95)
        res["pseudo"][sc] = {
            "G_mean": round(float(Gp.mean()), 5), "G_q025_q975": [round(float(q), 5) for q in np.percentile(Gp, [2.5, 97.5])],
            "G_ssw_minus_pseudo": round(real[sc]["G_M1_minus_M2"] - float(Gp.mean()), 5),
            "p_pseudo_ge_ssw": round(float(np.mean(Gp >= real[sc]["G_M1_minus_M2"])), 4),
            "M2b_gain_pseudo_mean": round(float(Bp.mean()), 5),
            "p_M2b_pseudo_ge_ssw": round(float(np.mean(Bp >= real[sc]["G_M0_minus_M2b"])), 4)}
        res["controls"][sc] = {
            "positive_detection_rate": round(float(np.mean([c[sc]["G_M1_minus_M2"] >= q95 for c in pos])), 3),
            "negative_false_rate": round(float(np.mean([c[sc]["G_M1_minus_M2"] >= q95 for c in neg])), 3),
            "positive_M2b_detection_rate": round(float(np.mean([c[sc]["G_M0_minus_M2b"] >= q95b for c in pos])), 3),
            "positive_G_mean": round(float(np.mean([c[sc]["G_M1_minus_M2"] for c in pos])), 5),
            "negative_G_mean": round(float(np.mean([c[sc]["G_M1_minus_M2"] for c in neg])), 5)}
        print(f"{sc}: pseudo {res['pseudo'][sc]} | controls {res['controls'][sc]}", flush=True)
    # ---------------- post-hoc diagnostics (see docstring; labelled in the JSON)
    from scipy.stats import skew as _skew
    ph = {"label": "post hoc, added after the registered result was seen"}
    for sc in ("crps", "ign"):
        Gs_p = np.array([p[sc]["PH1_G_M1s_minus_M2"] for p in pseudo])
        Gneg = np.array([c[sc]["G_M1_minus_M2"] for c in neg])
        g_ssw = real[sc]["G_M1_minus_M2"]
        ph[sc] = {"PH1_G_M1s_minus_M2_ssw": real[sc]["PH1_G_M1s_minus_M2"],
                  "PH1_G_M1_minus_M1s_ssw": real[sc]["PH1_G_M1_minus_M1s"],
                  "PH1_pseudo_mean": round(float(Gs_p.mean()), 5),
                  "PH1_pseudo_q025_q975": [round(float(q), 5) for q in np.percentile(Gs_p, [2.5, 97.5])],
                  "PH1_p_pseudo_ge_ssw": round(float(np.mean(Gs_p >= real[sc]["PH1_G_M1s_minus_M2"])), 4),
                  "PH2_p_vs_negative_control": round(float(np.mean(Gneg >= g_ssw)), 4),
                  "PH2_negative_q95": round(float(np.percentile(Gneg, 95)), 5)}
    rph = np.random.default_rng(zlib.crc32(b"class_model_comparison|PH3"))

    def comp(Zk, calk, yk):
        okk = np.isfinite(yk) & np.isfinite(Zk).all(1)
        f = em_mixture(Zk[okk], calk[okk], yk[okk], rph, logistic=True)
        _, coef, s1, s2, v, rbar = f
        pi = 1 / (1 + np.exp(-(np.column_stack([np.ones(okk.sum()), Zk[okk]]) @ v)))
        return {"separation": float(abs(coef[0] - coef[1])), "sd1": float(s1), "sd2": float(s2),
                "weight1": float(pi.mean()), "sd_membership": float(pi.std()),
                "skew": float(_skew(yk[okk])), "sd_y": float(yk[okk].std())}
    ssw_c = comp(Z, X[:, cal_ci], y)
    ps_c = []
    for k in range(50):
        Xk = np.concatenate([M[m][0][k][:, tier_ci] for m in _S["members"]])
        ck = np.concatenate([M[m][0][k][:, cal_ci] for m in _S["members"]])
        yk = np.concatenate([M[m][1][k] for m in _S["members"]])
        gk = np.concatenate([np.full(M[m][1].shape[1], m) for m in _S["members"]])
        ps_c.append(comp(FV.std_within(Xk, gk), ck, yk))
    ph["PH3_ssw"] = {k: round(v, 4) for k, v in ssw_c.items()}
    ph["PH3_pseudo_mean_50sets"] = {k: round(float(np.mean([q[k] for q in ps_c])), 4) for k in ssw_c}
    ph["PH3_pseudo_q025_q975_50sets"] = {k: [round(float(x), 4) for x in np.percentile([q[k] for q in ps_c], [2.5, 97.5])]
                                         for k in ssw_c}
    rb = np.random.default_rng(zlib.crc32(b"class_model_comparison|PH4"))
    okr = np.isfinite(y) & np.isfinite(Z).all(1); mmr = mdl[okr]; umr = np.unique(mmr)
    ph["PH4_skill_vs_M0_crps"] = {}
    for nm in ("M1", "M2", "M1s"):
        bb = []
        for _ in range(N_BOOT):
            pick = rb.choice(umr, len(umr), replace=True)
            idx = np.concatenate([np.flatnonzero(mmr == u) for u in pick])
            bb.append(1 - S[nm]["crps"][idx].sum() / S["M0"]["crps"][idx].sum())
        pk = f"skill_{nm}_vs_M0"
        pv = np.array([p["crps"][pk] for p in pseudo]) if pk in pseudo[0]["crps"] else None
        ph["PH4_skill_vs_M0_crps"][nm] = {
            "ssw": round(float(1 - S[nm]["crps"].sum() / S["M0"]["crps"].sum()), 4),
            "ssw_ci95_model_bootstrap": [round(float(q), 4) for q in np.percentile(bb, [2.5, 97.5])],
            "pseudo_mean": None if pv is None else round(float(pv.mean()), 4),
            "pseudo_q025_q975": None if pv is None else [round(float(q), 4) for q in np.percentile(pv, [2.5, 97.5])],
            "p_pseudo_le_ssw": None if pv is None else round(float(np.mean(pv <= 1 - S[nm]["crps"].sum() / S["M0"]["crps"].sum())), 4)}
    res["posthoc"] = ph
    print(f"POST HOC: {ph}", flush=True)
    res["snapsi"] = snapsi_arm(rng)
    print(f"SNAPSI: { {k: v for k, v in res['snapsi'].items() if k != 'per_ensemble'} }", flush=True)
    (RESULTS / "class_model_comparison.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print("Saved -> class_model_comparison.json")


if __name__ == "__main__":
    main()
