#!/usr/bin/env python3
"""
enso_regime_relative.py
=======================
IN MODELS THAT PRODUCE A STRONG EASTERN-PACIFIC EL NINO REGIME, DOES GREENHOUSE FORCING CHANGE HOW
OFTEN THAT FIXED REGIME OCCURS, OR DOES IT MOVE THE REGIMES (OR ONE POPULATION) ITSELF?

Registered 2026-10-16 (plan approved by the user the same day), committed with its validation gate and
its acquisition script (03_data_ingestion/acquire_enso_regime_fields.py) BEFORE any field was retrieved.
It answers the objection to enso_extreme_class.py that its registered rejection came from a model with
no strong regime (CanESM5) and that its null, fixed types in absolute Nino3 units, is not the regime
view: here the index is the E index of the regime literature, the frame moves with the convective
threshold, and a model is tested only if it produces a strong regime.

LITERATURE (opened 2026-10-16)
  Takahashi, Montecinos, Goubanova & Dewitte 2011, GRL 38, L10704: E and C indices from the leading
    EOFs of tropical Pacific SST; E accounts for the extreme eastern-Pacific warm events (abstract).
  Geng, Cai, Wu et al. 2022, Nat. Commun. 13, 6616 (PMC9666538, Methods "E and C indices", full text):
    EOFs of monthly SST anomalies over 15S-15N, 140E-80W, PCs normalised to unit s.d.,
    E = (PC1 - PC2)/sqrt2, C = (PC1 + PC2)/sqrt2; EOF2 cools the far east; the strong eastern-Pacific
    El Nino follows a nonlinear Bjerknes feedback that acts once warming passes a convective threshold.
  Cai et al. 2014, Nat. Clim. Change 4, 111 (full text): the "convective threshold increases with mean
    SSTs"; Takahashi & Dewitte 2016 (Clim. Dyn.): two El Nino regimes, a bimodal E at El Nino peaks.
FRAME  relative SST: ts minus the area-weighted ocean mean over 20S-20N at the same month (the
  convective threshold rises with the tropical mean), as anomalies from the L0 ensemble-mean monthly
  climatology. Forced changes of the eastern Pacific relative to the tropics therefore remain.
DATA  _enso_fields/ (2.5-degree grid; cells with OISST ocean fraction >= 0.9). Models: CESM1-LE (40),
  GFDL-CM3 (20), MIROC6, CanESM5 p1 and p2 (separately), ACCESS-ESM1-5; a model enters with >= 15
  members complete in all 1,080 months. LEVELS L0 DJF 1951-1980, L1 2021-2050, L2 2070-2099
  (Jan years; historical/RCP8.5 or SSP5-8.5).
INDEX per model: EOFs of the L0 monthly relative-SST anomalies over 15S-15N, 140E-280E (sqrt cos-lat
  weights, members pooled); EOF1 signed positive over Nino3 (5S-5N, 210-270E), EOF2 negative over
  Nino1+2 (10S-0, 270-280E); all levels projected on the L0 EOFs; PCs scaled by their L0 monthly s.d.;
  E = (PC1 - PC2)/sqrt2; DJF mean (all three months); then divided by the L0 s.d. of DJF E.
STEP 0 (before any test statistic): the model "produces a strong regime" if its L0 DJF E has
  BIC(1 component) - BIC(2) >= 10, an upper component of weight <= 0.30, and either two modes in the
  fitted density or an upper mean at least 2 lower-component s.d. above the lower mean. Models that
  fail are reported as unable to test the regime view, not as evidence against it.
P1 PRIMARY  fixed two-type mixing of DJF E over L0, L1, L2 (mixture_moment_test.tests, two-way
  bootstrap of members x 5-year blocks, B = 1999), per STEP-0 model, Holm over them; categories as in
  enso_extreme_class ("incompatible" Holm p < 0.05; "compatible" rho interval inside +-0.25; else
  "unresolved"). READING (fixed now): "the strong El Nino regime is not a fixed regime whose frequency
  forcing changes" if >= 1 model is incompatible and none compatible; "consistent with a fixed regime
  whose frequency forcing changes" if a majority is compatible; "not resolved" otherwise; "cannot test"
  if no model passes STEP 0.
SECONDARY (every STEP-0 model)
  S1 shape: L2 against L0 DJF E, each standardised by its own mean and s.d.: Anderson-Darling k-sample
     statistic, 1,999 paired permutations (each member's L0 and L2 seasons swapped with probability 1/2);
     "same shape" if p >= 0.05 and every 5th-95th percentile difference within +-0.25; "shape changes"
     if p < 0.05; else "unresolved".
  S2 regimes: two-component fits at each level, member bootstrap (1,000) of the upper weight change and
     of each component's shift L0 -> L2; and the share of the L0 -> L2 rise in the strong-event rate
     (DJF E above the pooled L0 95th percentile) that changing only the regime weights would give.
  S3 class law for that strong class, L0 -> L1 and L0 -> L2 (enso_extreme_class.law_test), +-0.25.
  S4 forced share of the strong-versus-rest contrast (enso_extreme_class.forced_share), +-0.05.
  S5 sensitivity: P1 and S1 on E from absolute SST (no tropical-mean removal).
DESCRIPTIVE  L0 skewness of DJF E and C; the nonlinearity alpha of PC2 = alpha PC1^2 + beta PC1 + gamma
  (monthly, L0); all of P1 and S1 for every model, labelled unselected.
VALIDATION GATE (--validate; before any data). 1,000 simulations per case, B = 999, three levels of 30
  seasons, year effects s.d. 0.1: V1 types N(0,1) and N(3.5,0.8), strong weight 0.05/0.10/0.20, 40
  members; V2 the same with 20 members; V3 N(0,1) and N(2.5,1.2), weights 0.08/0.15/0.25, 50 members.
  GATE: one-sided 95% Clopper-Pearson upper bound of the false-rejection rate at 0.05 below 0.07 in
  each. Power (no gate): W1 one skewed population z + 0.25 (z^2 - 1), mean 0/0.2/0.5, s.d. 1/1.1/1.25,
  40 members; W2 the V1 types both shifting by 0/0.3/0.6 with weight fixed at 0.10, 40 members.
AMENDMENT (2026-10-16, after the validation gate and BEFORE any field was retrieved; the download was
  stopped with only the land mask on disk). The gate passed, but it showed P1 almost powerless here
  (rejection 0.009 for one shifting population W1 and 0.008 for drifting regimes W2), so P1 is kept and
  reported, and a test with power is added as co-primary:
  P1b CO-PRIMARY  likelihood ratio of fixed regimes. H0: DJF E at the three levels is a two-component
     Gaussian mixture whose component means and s.d. are shared by the levels and whose weights are free
     per level (7 parameters); H1: a separate two-component mixture per level (15). LR = 2(l1 - l0),
     referred to chi-square with 8 d.f.; EM with several starts, H1 also started from the H0 solution so
     that l1 >= l0. Holm over the STEP-0 models. Category: "fixed regimes rejected" (Holm p < 0.05);
     "compatible with fixed regimes" if p >= 0.05 and the S2 bootstrap intervals of both component shifts
     L0 -> L2 lie inside +-0.25; else "unresolved". READING as for P1 with these categories. Its own gate
     (--validate-p1b): 1,000 simulations each of V1, V2, V3, upper bound of the false-rejection rate < 0.07;
     power reported for W1 and W2.
SECOND AMENDMENT (2026-10-16, after the P1b gate failed -- false rejection 0.083/0.070/0.073, upper bounds
  0.085-0.099 -- so P1b with the chi-square reference is not run; registered after the fields began
  downloading and BEFORE any statistic was computed on them). P1b-PB replaces it as co-primary: the same
  likelihood ratio (2 starts for H0, 2 for H1 plus the H0 solution, identical for data and replicates)
  referred to its parametric bootstrap under the fitted H0 (499 replicates for data). Categories and
  READING as P1b. Gate (--validate-pb): 1,000 simulations each of V1, V2, V3 with 99 replicates, upper
  bound < 0.07; power W1, W2 from 200 simulations.
THIRD AMENDMENT (2026-10-17, after the P1b-PB gate failed narrowly -- false rejection 0.062/0.063/0.062,
  upper bounds 0.076-0.077 -- and BEFORE any statistic was computed on the fields; P1b-PB at 0.05 is not
  run). P1b-PB04: the same statistic and bootstrap, with the level lowered to 0.04 (Holm family level
  0.04; "fixed regimes rejected" if the Holm-adjusted p < 0.04 -- with the gate's 99 replicates that rule is
  p <= 0.03, so p <= 0.03 is applied to the data, review 2026-10-22), the size correction taken from the failed
  gate's simulations only. Gate (--validate-pb04): 1,000 FRESH simulations each of V1, V2, V3 (new seeds),
  99 replicates, rejection at p < 0.04, upper bound < 0.07; power W1, W2 from 200. This is the third
  calibration tried; if it fails, P1b is reported as not run.
Output: results/current/8_experiment/enso_regime_relative.json
Run:    enso_regime_relative.py --validate ; --validate-p1b ; --validate-pb ; enso_regime_relative.py
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
from scipy.stats import anderson_ksamp, beta, skew

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = ROOT / "results" / "current" / "8_experiment" / "enso_regime_relative.json"
FIELDS = HERE.parents[0] / "03_data_ingestion" / "_enso_fields_v2"   # v1 had a regridding defect (2026-10-22)
sys.path.insert(0, str(HERE))
import mixture_moment_test as MM                    # noqa: E402
import enso_extreme_class as EX                      # noqa: E402

NAME = "enso_regime_relative"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
B, B_VAL, N_SIM, N_PERM, ALPHA, GATE, TOL = 1999, 999, 1000, 1999, 0.05, 0.07, 0.25
# registered 499 replicates for the data; raised to 1,999 after the first run (2026-10-17) to meet the repo's
# floor of 1,000 for any reported p value; with 499 all six p were 0.002, the floor (results kept in the JSON note)
P1B_BOOT = 1999
LAT = np.arange(-20.0, 20.01, 2.5); LON = np.arange(0.0, 360.0, 2.5)
LEVEL_JYEARS = {"L0": (1951, 1980), "L1": (2021, 2050), "L2": (2070, 2099)}
MIN_MEMBERS = 15
QS = np.round(np.arange(0.05, 0.951, 0.05), 2)


# ------------------------------------------------------------------ index
def box(lat0, lat1, lon0, lon1):
    return (LAT[:, None] >= lat0) & (LAT[:, None] <= lat1) & (LON[None, :] >= lon0) & (LON[None, :] <= lon1)


def load_model(model):
    files = sorted(FIELDS.glob(f"{model.split('-p')[0] if model.startswith('CanESM5') else model}__*.npz"))
    if model.startswith("CanESM5"):
        files = [f for f in files if f"i1{model[-2:]}" in f.stem]
    ts, ids = [], []
    for f in files:
        z = np.load(f)
        if z["ts"].shape[0] == 1080:
            ts.append(z["ts"]); ids.append(f.stem.split("__")[1]); yr, mo = z["year"], z["month"]
    return (np.array(ids), np.stack(ts), yr, mo) if ts else (np.array([]), None, None, None)


def e_index(ts, yr, mo, ocean, relative=True):
    """ts [members, 1080, lat, lon] -> DJF E and C per level (members x 30), alpha, L0 E skewness."""
    w = np.cos(np.deg2rad(LAT))[:, None] * np.ones((1, LON.size))
    band = ocean & box(-20, 20, 0, 360)
    x = ts.astype(np.float64)
    if relative:
        tm = np.nansum(x[:, :, band] * w[band], axis=2) / w[band].sum()
        x = x - tm[:, :, None, None]
    l0 = np.arange(360)
    clim = np.stack([np.nanmean(x[:, l0][:, mo[l0] == m], axis=(0, 1)) for m in range(1, 13)])
    a = x - clim[mo - 1][None]
    dom = ocean & box(-15, 15, 140, 280)
    sw = np.sqrt(w[dom])
    A = a[:, :, dom] * sw                              # members x 1080 x cells
    M = A[:, l0].reshape(-1, dom.sum())
    _, _, vt = np.linalg.svd(M - M.mean(0), full_matrices=False)
    e1, e2 = vt[0] / sw, vt[1] / sw                     # patterns in SST units
    n3 = box(-5, 5, 210, 270)[dom]; n12 = box(-10, 0, 270, 280)[dom]
    s1 = 1.0 if e1[n3].mean() > 0 else -1.0; s2 = 1.0 if e2[n12].mean() < 0 else -1.0
    pc1 = s1 * (A @ vt[0]); pc2 = s2 * (A @ vt[1])
    pc1 /= pc1[:, l0].std(); pc2 /= pc2[:, l0].std()
    Em = (pc1 - pc2) / np.sqrt(2); Cm = (pc1 + pc2) / np.sqrt(2)
    al = np.polyfit(pc1[:, l0].ravel(), pc2[:, l0].ravel(), 2)[0]
    jyear = np.where(mo == 12, yr + 1, yr)
    def djf(v):
        out = []
        for L, (y0, y1) in LEVEL_JYEARS.items():
            cols = []
            for y in range(y0, y1 + 1):
                k = (jyear == y) & np.isin(mo, (12, 1, 2))
                cols.append(v[:, k].mean(1) if k.sum() == 3 else np.full(v.shape[0], np.nan))
            out.append(np.stack(cols, 1))
        return out
    E, C = djf(Em), djf(Cm)
    s = np.nanstd(E[0]); E = [e / s for e in E]
    return E, C, float(al)


# ------------------------------------------------------------------ statistics
def gmm(x):
    from sklearn.mixture import GaussianMixture
    x = x.reshape(-1, 1)
    g1 = GaussianMixture(1, random_state=0).fit(x); g2 = GaussianMixture(2, n_init=10, random_state=0).fit(x)
    o = np.argsort(g2.means_.ravel()); mu = g2.means_.ravel()[o]; sd = np.sqrt(g2.covariances_.ravel()[o]); wt = g2.weights_[o]
    grid = np.linspace(x.min(), x.max(), 4001).reshape(-1, 1); f = np.exp(g2.score_samples(grid))
    modes = int(np.sum((f[1:-1] > f[:-2]) & (f[1:-1] > f[2:])))
    return {"dBIC_1_minus_2": round(float(g1.bic(x) - g2.bic(x)), 2), "means": [round(float(v), 3) for v in mu],
            "sd": [round(float(v), 3) for v in sd], "weights": [round(float(v), 3) for v in wt], "modes": modes}


def step0(g):
    sep = (g["means"][1] - g["means"][0]) / g["sd"][0]
    return bool(g["dBIC_1_minus_2"] >= 10 and g["weights"][1] <= 0.30 and (g["modes"] == 2 or sep >= 2.0))


def shape_test(X0, X2, rng):
    import warnings
    warnings.filterwarnings("ignore", message=".*p-value.*")
    z = lambda X: (X - X.mean()) / X.std()
    stat = lambda a, b: anderson_ksamp([z(a).ravel(), z(b).ravel()]).statistic
    t0 = stat(X0, X2); n = 0
    for _ in range(N_PERM):
        sw = rng.random(X0.shape[0]) < 0.5
        a = np.where(sw[:, None], X2, X0); b = np.where(sw[:, None], X0, X2)
        n += stat(a, b) >= t0
    p = (1 + n) / (N_PERM + 1)
    qd = np.quantile(z(X2).ravel(), QS) - np.quantile(z(X0).ravel(), QS)
    rd = ("shape changes" if p < 0.05 else "same shape" if np.all(np.abs(qd) <= TOL) else "unresolved")
    return {"AD_stat": round(float(t0), 4), "p_perm": round(float(p), 4), "max_abs_quantile_diff": round(float(np.max(np.abs(qd))), 3),
            "quantile_diff": dict(zip([str(q) for q in QS], np.round(qd, 3).tolist())), "reading": rd}


def decomposition(arms, thr):
    from sklearn.mixture import GaussianMixture
    out = {}
    for L, j in (("L0", 0), ("L2", 2)):
        x = arms[j][1].ravel()
        g = GaussianMixture(2, n_init=10, random_state=0).fit(x.reshape(-1, 1)); up = int(np.argmax(g.means_.ravel()))
        reg = g.predict(x.reshape(-1, 1)) == up
        out[L] = {"w_upper": float(reg.mean()), "p_upper": float((x[reg] > thr).mean()) if reg.any() else 0.0,
                  "p_lower": float((x[~reg] > thr).mean()), "rate": float((x > thr).mean())}
    a, b = out["L0"], out["L2"]
    w_only = b["w_upper"] * a["p_upper"] + (1 - b["w_upper"]) * a["p_lower"]
    d = b["rate"] - a["rate"]
    return {"L0": {k: round(v, 4) for k, v in a.items()}, "L2": {k: round(v, 4) for k, v in b.items()},
            "rate_weight_only": round(w_only, 4), "share_of_rate_change_from_weight": round((w_only - a["rate"]) / d, 4) if abs(d) > 1e-9 else None}


def decomposition_mixture(arms, thr):
    """POST REVIEW (2026-10-22): the weight-only counterfactual from the fitted mixtures themselves, so that
    it uses the same weights as S2's by_level fits: with the L0 components fixed, only their weights set to
    the L2 fit's; tail probabilities from the Gaussian components."""
    from scipy.stats import norm
    g0, g2 = gmm(arms[0][1].ravel()), gmm(arms[2][1].ravel())
    tail = lambda g, k: float(norm.sf(thr, g["means"][k], g["sd"][k]))
    fit_rate = lambda g, w: (1 - w) * tail(g, 0) + w * tail(g, 1)
    r0, r2 = fit_rate(g0, g0["weights"][1]), fit_rate(g2, g2["weights"][1])
    w_only = fit_rate(g0, g2["weights"][1]); d = r2 - r0
    return {"w_upper_L0_L2": [g0["weights"][1], g2["weights"][1]], "rate_fitted_L0_L2": [round(r0, 4), round(r2, 4)],
            "rate_observed_L0_L2": [round(float((arms[0][1] > thr).mean()), 4), round(float((arms[2][1] > thr).mean()), 4)],
            "rate_weight_only": round(w_only, 4), "share_of_rate_change_from_weight": round((w_only - r0) / d, 4) if abs(d) > 1e-9 else None}


def _npdf(x, m, s):
    return np.exp(-0.5 * ((x - m) / s) ** 2) / (s * np.sqrt(2 * np.pi))


def em_fixed(xs, init, n_iter=500):
    """H0: shared (m1, m2, s1, s2), weights per level. xs: list of 1-d arrays. Returns loglik, params."""
    m1, m2, s1, s2 = init[:4]; w = np.array(init[4], float)
    ll_old = -np.inf
    for _ in range(n_iter):
        rs, ll = [], 0.0
        for j, x in enumerate(xs):
            a = (1 - w[j]) * _npdf(x, m1, s1); b = w[j] * _npdf(x, m2, s2); t = a + b + 1e-300
            rs.append(b / t); ll += np.log(t).sum()
        X = np.concatenate(xs); R = np.concatenate(rs)
        w = np.array([r.mean() for r in rs]).clip(1e-4, 1 - 1e-4)
        m1 = np.sum((1 - R) * X) / np.sum(1 - R); m2 = np.sum(R * X) / np.sum(R)
        s1 = max(np.sqrt(np.sum((1 - R) * (X - m1) ** 2) / np.sum(1 - R)), 1e-3)
        s2 = max(np.sqrt(np.sum(R * (X - m2) ** 2) / np.sum(R)), 1e-3)
        if ll - ll_old < 1e-8:
            break
        ll_old = ll
    return ll, (m1, m2, s1, s2, w)


def em_level(x, init, n_iter=500):
    m1, m2, s1, s2, w = init; ll_old = -np.inf
    for _ in range(n_iter):
        a = (1 - w) * _npdf(x, m1, s1); b = w * _npdf(x, m2, s2); t = a + b + 1e-300; r = b / t; ll = np.log(t).sum()
        w = float(np.clip(r.mean(), 1e-4, 1 - 1e-4))
        m1 = np.sum((1 - r) * x) / np.sum(1 - r); m2 = np.sum(r * x) / np.sum(r)
        s1 = max(np.sqrt(np.sum((1 - r) * (x - m1) ** 2) / np.sum(1 - r)), 1e-3)
        s2 = max(np.sqrt(np.sum(r * (x - m2) ** 2) / np.sum(r)), 1e-3)
        if ll - ll_old < 1e-8:
            break
        ll_old = ll
    return ll


def lr_fixed_regimes(xs):
    """P1b: LR of H0 (shared components, free weights) against H1 (free mixture per level)."""
    from scipy.stats import chi2
    from sklearn.mixture import GaussianMixture
    X = np.concatenate(xs)
    best0 = (-np.inf, None)
    for seed in range(5):
        g = GaussianMixture(2, n_init=1, random_state=seed).fit(X.reshape(-1, 1)); o = np.argsort(g.means_.ravel())
        mu = g.means_.ravel()[o]; sd = np.sqrt(g.covariances_.ravel()[o]); wu = g.weights_[o][1]
        r = em_fixed(xs, (mu[0], mu[1], sd[0], sd[1], [wu] * len(xs)))
        if r[0] > best0[0]:
            best0 = r
    l0, (m1, m2, s1, s2, w) = best0
    l1 = 0.0
    for j, x in enumerate(xs):
        cands = [em_level(x, (m1, m2, s1, s2, w[j]))]
        g = GaussianMixture(2, n_init=5, random_state=0).fit(x.reshape(-1, 1)); cands.append(g.score(x.reshape(-1, 1)) * x.size)
        l1 += max(cands)
    lr = max(2 * (l1 - l0), 0.0); df = 3 * 5 - (4 + len(xs))
    return {"LR": round(float(lr), 3), "df": df, "p": round(float(chi2.sf(lr, df)), 5),
            "H0_components": [round(float(v), 3) for v in (m1, m2, s1, s2)], "H0_upper_weights": [round(float(v), 3) for v in w]}


def lr_stat(xs, n_seeds=2):
    """the LR used by P1b-PB, identical for data and bootstrap replicates."""
    from sklearn.mixture import GaussianMixture
    X = np.concatenate(xs); best0 = (-np.inf, None)
    for seed in range(n_seeds):
        g = GaussianMixture(2, n_init=1, random_state=seed).fit(X.reshape(-1, 1)); o = np.argsort(g.means_.ravel())
        mu = g.means_.ravel()[o]; sd = np.sqrt(g.covariances_.ravel()[o])
        r = em_fixed(xs, (mu[0], mu[1], sd[0], sd[1], [g.weights_[o][1]] * len(xs)))
        if r[0] > best0[0]:
            best0 = r
    l0, (m1, m2, s1, s2, w) = best0
    l1 = 0.0
    for j, x in enumerate(xs):
        g = GaussianMixture(2, n_init=2, random_state=0).fit(x.reshape(-1, 1))
        l1 += max(em_level(x, (m1, m2, s1, s2, w[j])), g.score(x.reshape(-1, 1)) * x.size)
    return max(2 * (l1 - l0), 0.0), (m1, m2, s1, s2, w)


def lr_pb(xs, rng, n_boot):
    """P1b-PB: the LR referred to its parametric-bootstrap distribution under the fitted H0."""
    lr, (m1, m2, s1, s2, w) = lr_stat(xs)
    n = [x.size for x in xs]; ge = 0
    for _ in range(n_boot):
        xb = []
        for j in range(len(xs)):
            t = rng.random(n[j]) < w[j]
            xb.append(np.where(t, rng.normal(m2, s2, n[j]), rng.normal(m1, s1, n[j])))
        ge += lr_stat(xb)[0] >= lr
    return {"LR": round(float(lr), 3), "n_boot": n_boot, "p": round((1 + ge) / (n_boot + 1), 5),
            "H0_components": [round(float(v), 3) for v in (m1, m2, s1, s2)], "H0_upper_weights": [round(float(v), 3) for v in w]}


def _misspec_job(args):
    """POST HOC (review 2026-10-22): false rejection of P1b under FIXED regimes that are not Gaussian."""
    kind, seed, nsim, nb, payload = args
    rng = np.random.default_rng(seed); rej = 0
    from scipy.stats import skewnorm, gamma
    for _ in range(nsim):
        xs = []
        for j in range(3):
            if kind == "empirical":
                A, Bv, w, n = payload["A"], payload["B"], payload["w"][j], payload["n"][j]
                t = rng.random(n) < w
                xs.append(np.where(t, rng.choice(Bv, n), rng.choice(A, n)))
            else:                                   # skewed fixed types, strong weight 0.05/0.10/0.20, 40 x 30
                n = 1200; w = (0.05, 0.10, 0.20)[j]; t = rng.random(n) < w
                mod = skewnorm.rvs(4, size=n, random_state=rng); mod = (mod - mod.mean()) / mod.std()
                strong = 2.0 + gamma.rvs(3, scale=0.5, size=n, random_state=rng)
                xs.append(np.where(t, strong, mod))
        rej += lr_pb(xs, rng, nb)["p"] <= 0.03
    return kind if kind != "empirical" else payload["model"], rej


def post_hoc_misspecification(arms, p1b, nsim=400, nb=99):
    """fixed regimes with each model's own L0 regime shapes (seasons assigned by the H0 fit's posterior)
    and the H0 weights per level; and a generic skewed pair. Returns false-rejection rates at p <= 0.03."""
    jobs, chunks = [], 20
    for m in arms:
        m1, m2, s1, s2 = p1b[m]["H0_components"]; w = p1b[m]["H0_upper_weights"]
        x0 = arms[m][0][1].ravel()
        post = w[0] * _npdf(x0, m2, s2) / (w[0] * _npdf(x0, m2, s2) + (1 - w[0]) * _npdf(x0, m1, s1))
        pay = {"model": m, "A": x0[post <= 0.5], "B": x0[post > 0.5], "w": w, "n": [X.size for _, X in arms[m]]}
        jobs += [("empirical", SEED + 30000 + 100 * len(jobs) + c, nsim // chunks, nb, pay) for c in range(chunks)]
    jobs += [("skewed", SEED + 40000 + c, nsim // chunks, nb, None) for c in range(chunks)]
    import multiprocessing as mp
    with mp.get_context("fork").Pool(20) as pool:
        outs = pool.map(_misspec_job, jobs)
    agg = {}
    for k, r in outs:
        agg[k] = agg.get(k, 0) + r
    ub = lambda r, n: float(beta.ppf(0.95, r + 1, n - r)) if r < n else 1.0
    return {"note": "POST HOC (review 2026-10-22): P1b false rejection (p <= 0.03, 99 replicates) when regimes are fixed but not Gaussian",
            "n_sim": nsim, "cases": {k: {"reject": v / nsim, "upper95": round(ub(v, nsim), 4)} for k, v in agg.items()}}


def _val_pb(args):
    kind, seed, nsim, nb, level = args
    rng = np.random.default_rng(seed); rej = 0
    for _ in range(nsim):
        rej += lr_pb([X.ravel() for _, X in sim(kind, rng)], rng, nb)["p"] < level
    return kind, rej


def validate_pb(level=ALPHA, seed_off=9000, key="validation_P1b_PB"):
    plan = {"V1": 1000, "V2": 1000, "V3": 1000, "W1": 200, "W2": 200}; chunks = 20
    jobs = [(k, SEED + seed_off + 1000 * i + c, n // chunks, 99, level) for i, (k, n) in enumerate(plan.items()) for c in range(chunks)]
    import multiprocessing as mp
    with mp.get_context("fork").Pool(20) as pool:
        outs = pool.map(_val_pb, jobs)
    agg = {}
    for k, r in outs:
        agg[k] = agg.get(k, 0) + r
    ub = lambda r, n: float(beta.ppf(0.95, r + 1, n - r)) if r < n else 1.0
    val = {"replicates_per_simulated_set": 99, "level": level, "cases": {k: {"n_sim": plan[k], "reject": v / plan[k], "upper95": round(ub(v, plan[k]), 4)} for k, v in agg.items()}}
    val["passes"] = all(val["cases"][k]["upper95"] < GATE for k in ("V1", "V2", "V3"))
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    res[key] = val
    res[key + "_fixed"] = ("2026-10-16, second amendment" if key == "validation_P1b_PB" else "2026-10-17, third amendment") + \
        ", after the fields were downloaded and before any statistic was computed on them"
    OUT.write_text(json.dumps(res, indent=1), encoding="utf8", newline="\n")
    print(json.dumps(val, indent=1)); return 0


def category_p1b(r, p_holm, s2boot, level=ALPHA):
    lo, up = s2boot["lower_mode_shift_K_ci95"], s2boot["upper_mode_shift_K_ci95"]
    inside = all(-TOL <= a and b <= TOL for a, b in (lo, up))
    return "fixed regimes rejected" if p_holm < level else "compatible with fixed regimes" if inside else "unresolved"


def _val_p1b(args):
    kind, seed, nsim = args
    rng = np.random.default_rng(seed); rej = 0
    for _ in range(nsim):
        rej += lr_fixed_regimes([X.ravel() for _, X in sim(kind, rng)])["p"] < ALPHA
    return kind, rej


def validate_p1b():
    kinds = ["V1", "V2", "V3", "W1", "W2"]; chunks = 10
    jobs = [(k, SEED + 5000 + 1000 * i + c, N_SIM // chunks) for i, k in enumerate(kinds) for c in range(chunks)]
    import multiprocessing as mp
    with mp.get_context("fork").Pool(10) as pool:
        outs = pool.map(_val_p1b, jobs)
    agg = {}
    for k, r in outs:
        agg[k] = agg.get(k, 0) + r
    ub = lambda r: float(beta.ppf(0.95, r + 1, N_SIM - r)) if r < N_SIM else 1.0
    val = {"n_sim": N_SIM, "cases": {k: {"reject": v / N_SIM, "upper95": round(ub(v), 4)} for k, v in agg.items()}}
    val["passes"] = all(val["cases"][k]["upper95"] < GATE for k in ("V1", "V2", "V3"))
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    res["validation_P1b"] = val; res["validation_P1b_fixed"] = "2026-10-16, amendment, before any field was retrieved"
    OUT.write_text(json.dumps(res, indent=1), encoding="utf8", newline="\n")
    print(json.dumps(val, indent=1)); return 0


# ------------------------------------------------------------------ validation
def sim(kind, rng, nY=30):
    n = {"V1": 40, "V2": 20, "V3": 50, "W1": 40, "W2": 40}[kind]
    ids = np.array([f"m{i}" for i in range(n)]); arms = []
    for j in range(3):
        ye = 0.1 * rng.normal(size=nY)[None, :]
        if kind in ("V1", "V2", "V3"):
            p = (0.05, 0.10, 0.20)[j] if kind != "V3" else (0.08, 0.15, 0.25)[j]
            hi = (3.5, 0.8) if kind != "V3" else (2.5, 1.2)
            typ = rng.random((n, nY)) < p
            X = np.where(typ, rng.normal(hi[0], hi[1], (n, nY)), rng.normal(0, 1, (n, nY)))
        elif kind == "W1":
            z = rng.normal(size=(n, nY)); X = (0, 0.2, 0.5)[j] + (1, 1.1, 1.25)[j] * (z + 0.25 * (z ** 2 - 1))
        else:
            typ = rng.random((n, nY)) < 0.10
            X = np.where(typ, rng.normal(3.5, 0.8, (n, nY)), rng.normal(0, 1, (n, nY))) + (0, 0.3, 0.6)[j]
        arms.append((ids, X + ye))
    return arms


def _val_job(args):
    kind, seed, nsim = args
    rng = np.random.default_rng(seed); rej = 0
    for _ in range(nsim):
        rej += MM.tests(sim(kind, rng), B_VAL, rng)["p_T"] < ALPHA
    return kind, rej


def validate():
    kinds = ["V1", "V2", "V3", "W1", "W2"]; chunks = 10
    jobs = [(k, SEED + 1000 * i + c, N_SIM // chunks) for i, k in enumerate(kinds) for c in range(chunks)]
    import multiprocessing as mp
    with mp.get_context("fork").Pool(8) as pool:
        outs = pool.map(_val_job, jobs)
    agg = {}
    for k, r in outs:
        agg[k] = agg.get(k, 0) + r
    ub = lambda r: float(beta.ppf(0.95, r + 1, N_SIM - r)) if r < N_SIM else 1.0
    val = {"n_sim": N_SIM, "B": B_VAL, "cases": {k: {"reject": v / N_SIM, "upper95": round(ub(v), 4)} for k, v in agg.items()}}
    val["passes"] = all(val["cases"][k]["upper95"] < GATE for k in ("V1", "V2", "V3"))
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    res["validation"] = val; res["validation_fixed"] = "2026-10-16, before any field was retrieved"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=1), encoding="utf8", newline="\n")
    print(json.dumps(val, indent=1)); return 0


# ------------------------------------------------------------------ main
def main():
    if "--validate" in sys.argv:
        return validate()
    if "--validate-p1b" in sys.argv:
        return validate_p1b()
    if "--validate-pb" in sys.argv:
        return validate_pb()
    if "--validate-pb04" in sys.argv:
        return validate_pb(level=0.04, seed_off=20000, key="validation_P1b_PB04")
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    if not res.get("validation", {}).get("passes"):
        sys.exit("validation gate not passed (run --validate)")
    import platform
    import scipy
    import sklearn
    for k in ("validation_P1b_PB", "validation_P1b_PB04"):      # inner Monte Carlo replicates of a gate, not a reported p
        if "n_boot" in res.get(k, {}):
            res[k]["replicates_per_simulated_set"] = res[k].pop("n_boot")
    ocean = np.load(FIELDS / "ocean_frac.npy") >= 0.9
    res.update({"registered": "2026-10-16, committed before any field was retrieved", "seed": SEED, "B": B, "n_perm": N_PERM,
                "inputs": [str(FIELDS.relative_to(ROOT))], "command": " ".join(sys.argv),
                "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "sklearn": sklearn.__version__}})
    models = ["CESM1-LE", "GFDL-CM3", "MIROC6", "CanESM5-p1", "CanESM5-p2", "ACCESS-ESM1-5"]
    info, arms, arms_abs = {}, {}, {}
    for m in models:
        ids, ts, yr, mo = load_model(m)
        info[m] = {"n_members": int(len(ids))}
        if len(ids) < MIN_MEMBERS:
            continue
        E, C, al = e_index(ts, yr, mo, ocean, True)
        Ea, _, _ = e_index(ts, yr, mo, ocean, False)
        ok = np.all(np.isfinite(np.concatenate(E, 1)), axis=1)
        arms[m] = [(ids[ok], X[ok]) for X in E]; arms_abs[m] = [(ids[ok], X[ok]) for X in Ea]
        g0 = gmm(E[0][ok].ravel())
        info[m].update({"n_complete": int(ok.sum()), "alpha_L0": round(al, 4), "skew_E_L0": round(float(skew(E[0][ok].ravel())), 3),
                        "skew_C_L0": round(float(skew(C[0][ok].ravel())), 3), "mean_E": [round(float(X[ok].mean()), 3) for X in E],
                        "sd_E": [round(float(X[ok].std()), 3) for X in E], "gmm_L0": g0, "step0": step0(g0)})
        print(m, json.dumps(info[m]), flush=True)
    sel = [m for m in arms if info[m]["step0"]]
    res["step0_models"] = info; res["selected"] = sel
    p1 = {m: MM.tests(arms[m], B, np.random.default_rng(SEED + 1 + i)) for i, m in enumerate(sorted(arms))}
    for fam, keys in (("selected", sel), ("all", sorted(arms))):
        if keys:
            adj = EX.holm([p1[m]["p_T"] for m in keys])
            for m, a in zip(keys, adj):
                p1[m][f"p_T_holm_{fam}"] = round(float(a), 4); p1[m][f"category_{fam}"] = EX.category(p1[m], a)
    res["P1_E_mixing"] = p1
    cats = [p1[m]["category_selected"] for m in sel]
    inc = sum(c.startswith("incompatible") for c in cats); com = sum(c == "compatible" for c in cats)
    res["P1_reading"] = ("cannot test: no model produces a strong regime" if not sel else
                         "the strong El Nino regime is not a fixed regime whose frequency forcing changes" if inc >= 1 and com == 0 else
                         "consistent with a fixed regime whose frequency forcing changes" if com > len(sel) / 2 else "not resolved")
    print("P1", res["P1_reading"], json.dumps({m: (p1[m]["rho"], p1[m]["rho_ci95"], p1[m]["p_T"]) for m in p1}), flush=True)
    res["S1_shape"] = {m: shape_test(arms[m][0][1], arms[m][2][1], np.random.default_rng(SEED + 30 + i)) for i, m in enumerate(sorted(arms))}
    S2, S3, S4 = {}, {}, {}
    for i, m in enumerate(sorted(arms)):
        thr = float(np.quantile(arms[m][0][1].ravel(), 0.95))
        S2[m] = {"by_level": {L: gmm(X.ravel()) for L, (_, X) in zip(LEVEL_JYEARS, arms[m])},
                 "bootstrap_L0_to_L2": EX.regime_boot(arms[m], np.random.default_rng(SEED + 50 + i)),
                 "strong_threshold_L0_p95": round(thr, 4), "decomposition": decomposition(arms[m], thr),
                 "decomposition_mixture_post_review": decomposition_mixture(arms[m], thr)}
        S3[m] = {"L0_to_L1": EX.law_test(arms[m], 0, 1, thr, np.random.default_rng(SEED + 70 + i)),
                 "L0_to_L2": EX.law_test(arms[m], 0, 2, thr, np.random.default_rng(SEED + 90 + i))}
        S4[m] = EX.forced_share(arms[m], thr, np.random.default_rng(SEED + 110 + i))
    res.update({"S2_regimes": S2, "S3_class_law": S3, "S4_forced_share": S4})
    pb04 = res.get("validation_P1b_PB04", {}).get("passes")
    pb_ok = res.get("validation_P1b_PB", {}).get("passes") or pb04
    # the 0.04 gate ran with 99 replicates, where p < 0.04 means p <= 0.03 (review 2026-10-22): the validated
    # rule is p <= 0.03, applied here with 1,999 replicates
    lev = 0.0301 if pb04 and not res.get("validation_P1b_PB", {}).get("passes") else ALPHA
    if res.get("validation_P1b", {}).get("passes") or pb_ok:
        p1b = ({m: lr_pb([X.ravel() for _, X in arms[m]], np.random.default_rng(SEED + 170 + i), P1B_BOOT) for i, m in enumerate(sorted(arms))}
               if pb_ok else {m: lr_fixed_regimes([X.ravel() for _, X in arms[m]]) for m in sorted(arms)})
        for fam, keys in (("selected", sel), ("all", sorted(arms))):
            if keys:
                adj = EX.holm([p1b[m]["p"] for m in keys])
                for m, a in zip(keys, adj):
                    p1b[m][f"p_holm_{fam}"] = round(float(a), 5)
                    p1b[m][f"category_{fam}"] = category_p1b(p1b[m], a, S2[m]["bootstrap_L0_to_L2"], lev)
        cb = [p1b[m]["category_selected"] for m in sel]
        rej = sum(c == "fixed regimes rejected" for c in cb); cmp_ = sum(c.startswith("compatible") for c in cb)
        res["P1b_note"] = ("registered with 499 replicates; that run (2026-10-17 01:10) gave LR 104.1/106.9/82.8/68.7/72.8/725.1 and p 0.002 (the floor) for ACCESS-ESM1-5/CESM1-LE/CanESM5-p1/CanESM5-p2/GFDL-CM3/MIROC6, Holm 0.006 for the three selected; re-run with 1,999 for the repo floor")
        res["P1b_LR_fixed_regimes"] = p1b; res["P1b_reference"] = (f"parametric bootstrap at level {lev}" if pb_ok else "chi-square")
        res["P1b_reading"] = ("cannot test: no model produces a strong regime" if not sel else
                              "the strong El Nino regime is not a fixed regime whose frequency forcing changes" if rej >= 1 and cmp_ == 0 else
                              "consistent with a fixed regime whose frequency forcing changes" if cmp_ > len(sel) / 2 else "not resolved")
        print("P1b", res["P1b_reading"], json.dumps({m: (p1b[m]["LR"], p1b[m]["p"]) for m in p1b}), flush=True)
    else:
        res["P1b_reading"] = "not run: none of the three calibrations (chi-square, bootstrap at 0.05, bootstrap at 0.04) passed its gate"
    res["S5_absolute_frame"] = {m: {"P1": MM.tests(arms_abs[m], B, np.random.default_rng(SEED + 130 + i)),
                                    "S1": shape_test(arms_abs[m][0][1], arms_abs[m][2][1], np.random.default_rng(SEED + 150 + i))}
                                for i, m in enumerate(sorted(arms_abs))}
    if "P1b_LR_fixed_regimes" in res and "--no-posthoc" not in sys.argv:
        sel_arms = {m: arms[m] for m in sel}
        res["post_hoc_P1b_misspecification"] = post_hoc_misspecification(sel_arms, res["P1b_LR_fixed_regimes"])
        print("POSTHOC", json.dumps(res["post_hoc_P1b_misspecification"]["cases"]), flush=True)
    res["S2_note"] = ("POST REVIEW: in simulations with fixed Gaussian types whose weight rises (V1, V3), per-level two-component "
                      "fits report a falling weight and a moving upper component; the per-level weights, shifts and both "
                      "decompositions are therefore descriptive of the fits and are not evidence about fixed regimes")
    for k in ("S1_shape", "S2_regimes", "S3_class_law", "S4_forced_share", "S5_absolute_frame"):
        print(k, json.dumps(res[k], default=float)[:2000], flush=True)
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf8", newline="\n")
    print(f"Saved -> {OUT.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
