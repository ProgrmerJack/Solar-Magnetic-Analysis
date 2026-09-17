#!/usr/bin/env python3
"""
condpois.py -- conditional Poisson regression with stratum fixed effects.

Used by the SSW/avalanche case-crossover analyses (r107+).

WHY CONDITIONAL RATHER THAN DUMMY-VARIABLE FIXED EFFECTS
  The design has one stratum per region-winter (~100-200 strata) and only a
  handful of real parameters. Fitting stratum intercepts explicitly is both
  slow (a 25,000 x 130 GLM per bootstrap replicate) and statistically awkward:
  the number of nuisance parameters grows with the data (the incidental-
  parameters problem). For Poisson the stratum intercepts have a closed form
  given the other coefficients,

      alpha_s = log( sum_{i in s} y_i / sum_{i in s} exp(x_i'b) ),

  and substituting it back gives the conditional (multinomial) log-likelihood

      l(b) = sum_i y_i x_i'b - sum_s Y_s log( sum_{i in s} exp(x_i'b) ),

  which depends only on b. This is the standard conditional-Poisson estimator
  used for case-crossover designs; it is exact, not an approximation, and it
  reduces each fit to a ~10-parameter optimisation.

Self-check: `python condpois.py` verifies the estimator against a
statsmodels dummy-variable Poisson GLM on simulated data.
"""
import numpy as np
from scipy.optimize import minimize


def fit(y, X, strata, offset=None, tol=1e-9, maxiter=500):
    """Conditional Poisson MLE.

    y       (n,)   counts
    X       (n,k)  covariates (NO intercept -- absorbed by the strata)
    strata  (n,)   integer stratum codes
    offset  (n,)   optional log-exposure (e.g. log number of reporting stations),
                   making the model a rate model: E[y] = exposure * exp(a_s + x'b)
    returns (k,) coefficients, or None if the fit fails / is unidentified.
    """
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    s = np.asarray(strata)
    off = None if offset is None else np.asarray(offset, float)

    # keep only strata containing at least one event; strata with Y_s = 0
    # contribute nothing to the conditional likelihood
    codes, s = np.unique(s, return_inverse=True)
    ns = len(codes)
    Ys = np.bincount(s, weights=y, minlength=ns)
    keep_s = Ys > 0
    if keep_s.sum() < 2:
        return None
    row = keep_s[s]
    if row.sum() < 10:
        return None
    y, X = y[row], X[row]
    if off is not None:
        off = off[row]
    _, s = np.unique(s[row], return_inverse=True)
    ns = s.max() + 1
    Ys = np.bincount(s, weights=y, minlength=ns)

    k = X.shape[1]
    suff = X.T @ y

    def negll(b):
        eta = X @ b if off is None else X @ b + off
        mx = np.full(ns, -np.inf)
        np.maximum.at(mx, s, eta)          # per-stratum max for stability
        e = np.exp(eta - mx[s])
        S = np.bincount(s, weights=e, minlength=ns)
        S = np.maximum(S, 1e-300)
        ll = y @ eta - Ys @ (np.log(S) + mx)
        w = Ys / S
        grad = suff - np.array(
            [np.bincount(s, weights=e * X[:, j], minlength=ns) @ w for j in range(k)])
        return -ll, -grad

    try:
        r = minimize(negll, np.zeros(k), jac=True, method="L-BFGS-B",
                     options={"maxiter": maxiter, "ftol": tol, "gtol": 1e-8})
    except Exception:
        return None
    if not np.all(np.isfinite(r.x)):
        return None
    return r.x


def _selfcheck():
    """Conditional estimator must reproduce dummy-variable Poisson GLM."""
    import statsmodels.api as sm
    rng = np.random.default_rng(0)
    n_s, per = 40, 60
    strata = np.repeat(np.arange(n_s), per)
    alpha = rng.normal(0, 1.0, n_s)          # stratum levels
    W = rng.integers(0, 2, n_s * per).astype(float)
    z = rng.normal(0, 1, n_s * per)
    beta = np.array([-0.35, 0.20])           # true W and z effects
    mu = np.exp(alpha[strata] + W * beta[0] + z * beta[1])
    y = rng.poisson(mu).astype(float)
    X = np.column_stack([W, z])

    b_cond = fit(y, X, strata)

    D = np.zeros((len(y), n_s))
    D[np.arange(len(y)), strata] = 1.0
    b_glm = sm.GLM(y, np.hstack([X, D]), family=sm.families.Poisson()).fit().params[:2]

    assert b_cond is not None, "conditional fit returned None"
    assert np.allclose(b_cond, b_glm, atol=1e-5), \
        f"conditional {b_cond} != dummy GLM {b_glm}"
    print(f"OK  conditional={np.round(b_cond,5)}  dummyGLM={np.round(b_glm,5)}  "
          f"true={beta}")

    # offset (rate model) must match a dummy-variable GLM with the same offset
    expo = rng.integers(1, 20, n_s * per).astype(float)
    mu2 = expo * np.exp(alpha[strata] + W * beta[0] + z * beta[1])
    y2 = rng.poisson(mu2).astype(float)
    b_cond2 = fit(y2, X, strata, offset=np.log(expo))
    b_glm2 = sm.GLM(y2, np.hstack([X, D]), family=sm.families.Poisson(),
                    offset=np.log(expo)).fit().params[:2]
    assert b_cond2 is not None and np.allclose(b_cond2, b_glm2, atol=1e-5), \
        f"offset: conditional {b_cond2} != dummy GLM {b_glm2}"
    print(f"OK  offset   conditional={np.round(b_cond2,5)}  dummyGLM={np.round(b_glm2,5)}")

    # unidentified case must return None rather than a wrong number
    assert fit(np.zeros(100), np.ones((100, 1)), np.arange(100)) is None
    print("OK  degenerate input returns None")


if __name__ == "__main__":
    _selfcheck()
