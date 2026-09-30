"""
R57 — Negative Binomial & Zero-Inflated Model Comparisons for SSW–Avalanche Panel
==================================================================================

Fits Poisson, NB, ZIP, ZINB regressions on dry natural avalanche counts,
compares via AIC/BIC/Vuong, and reports GEE with winter-level clustering.

Addresses reviewer concerns:
  1. NB as natural competitor given overdispersion (var/mean ≈ 45)
  2. GEE with exchangeable correlation within winters
  3. Vuong test for ZIP vs ZINB model selection
"""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.discrete.count_model import (
    ZeroInflatedNegativeBinomialP,
    ZeroInflatedPoisson,
)
from statsmodels.discrete.discrete_model import NegativeBinomialP, Poisson
from statsmodels.genmod.generalized_estimating_equations import GEE
from statsmodels.genmod.cov_struct import Exchangeable
from statsmodels.genmod.families import Poisson as PoissonFamily, NegativeBinomial as NBFamily

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message=".*divide by zero.*")
warnings.filterwarnings("ignore", message=".*invalid value.*")

# ── paths ──────────────────────────────────────────────────────────────
DATA = Path(r"C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\analysis_panel_v2.parquet")
OUT_JSON = Path(r"C:\Users\Jack0\Solar-Magnetic-Analysis\data\results\r57_nb_models.json")

OUTCOME = "dry_natural_size_1234"
SSW_COL = "ssw_within_15d"
Z500_COL = "ncep_z500_nh"
WINTER_COL = "winter_id"


def load_data():
    """Load panel, filter to winter rows with valid outcome."""
    df = pd.read_parquet(DATA)
    # handle tz-aware index if present
    if hasattr(df.index, "tz") and df.index.tz is not None:
        df.index = df.index.tz_localize(None)

    mask = (df["is_winter"] == 1) & df[OUTCOME].notna() & df[Z500_COL].notna()
    wdf = df.loc[mask].copy()
    wdf[OUTCOME] = wdf[OUTCOME].astype(int)

    # standardise Z500 for numerical stability
    wdf["z500_std"] = (wdf[Z500_COL] - wdf[Z500_COL].mean()) / wdf[Z500_COL].std()

    print(f"Analysis sample: N={len(wdf)}, winters={wdf[WINTER_COL].nunique()}")
    y = wdf[OUTCOME]
    print(f"  Outcome mean={y.mean():.3f}, var={y.var():.2f}, "
          f"zeros={( y==0).sum()} ({(y==0).mean()*100:.1f}%)")
    print(f"  SSW prevalence: {wdf[SSW_COL].mean()*100:.1f}%")
    return wdf


def build_matrices(wdf):
    """Return endog, exog, exog_inflate, and group vector."""
    y = wdf[OUTCOME].values
    X = sm.add_constant(wdf[[SSW_COL, "z500_std"]].values)
    X_names = ["const", "ssw_within_15d", "z500_std"]
    groups = wdf[WINTER_COL].values
    return y, X, X_names, groups


def fit_poisson(y, X, X_names):
    mod = Poisson(y, X)
    res = mod.fit(disp=0, maxiter=300)
    return res


def fit_nb(y, X, X_names):
    mod = NegativeBinomialP(y, X, p=2)  # NB2 parameterisation
    res = mod.fit(disp=0, maxiter=300)
    return res


def fit_zip(y, X, X_names):
    mod = ZeroInflatedPoisson(y, X, exog_infl=X, inflation="logit")
    res = mod.fit(disp=0, maxiter=500, method="bfgs")
    return res


def fit_zinb(y, X, X_names):
    mod = ZeroInflatedNegativeBinomialP(y, X, exog_infl=X, p=2, inflation="logit")
    res = mod.fit(disp=0, maxiter=500, method="bfgs")
    return res


def vuong_test(zip_res, zinb_res, y, X):
    """
    Vuong (1989) non-nested test comparing ZIP vs ZINB.
    Positive Z → ZIP preferred; Negative Z → ZINB preferred.
    """
    ll_zip = zip_res.predict(which="prob")
    ll_zinb = zinb_res.predict(which="prob")

    # predicted probabilities for observed counts
    n = len(y)
    m_zip = np.zeros(n)
    m_zinb = np.zeros(n)
    for i in range(n):
        yi = int(y[i])
        if yi < ll_zip.shape[1]:
            m_zip[i] = ll_zip[i, yi]
        else:
            m_zip[i] = 1e-300
        if yi < ll_zinb.shape[1]:
            m_zinb[i] = ll_zinb[i, yi]
        else:
            m_zinb[i] = 1e-300

    m_zip = np.clip(m_zip, 1e-300, None)
    m_zinb = np.clip(m_zinb, 1e-300, None)

    mi = np.log(m_zip) - np.log(m_zinb)
    v_stat = np.sqrt(n) * mi.mean() / mi.std(ddof=1)
    p_val = 2 * stats.norm.sf(abs(v_stat))
    return float(v_stat), float(p_val)


def vuong_test_simple(res1, res2):
    """Simplified Vuong using per-observation log-likelihoods."""
    ll1 = res1.llf
    ll2 = res2.llf
    k1 = res1.df_model + 1
    k2 = res2.df_model + 1
    n = res1.nobs
    # AIC-corrected Vuong
    correction = (k1 - k2) / (2 * n)
    # Cannot do per-obs without predict(which='prob'), fall back to LR-based
    lr = 2 * (ll2 - ll1)
    return float(lr), k2 - k1


def fit_gee_poisson(wdf):
    """GEE Poisson with exchangeable correlation, clustered by winter."""
    y = wdf[OUTCOME].values.astype(float)
    X = sm.add_constant(wdf[[SSW_COL, "z500_std"]].values)
    groups = wdf[WINTER_COL].values

    mod = GEE(
        y, X, groups=groups,
        family=PoissonFamily(),
        cov_struct=Exchangeable(),
    )
    res = mod.fit(maxiter=200)
    return res


def fit_gee_nb(wdf):
    """GEE NB with exchangeable correlation, clustered by winter."""
    y = wdf[OUTCOME].values.astype(float)
    X = sm.add_constant(wdf[[SSW_COL, "z500_std"]].values)
    groups = wdf[WINTER_COL].values

    mod = GEE(
        y, X, groups=groups,
        family=NBFamily(alpha=1.0),
        cov_struct=Exchangeable(),
    )
    res = mod.fit(maxiter=200)
    return res


def extract_results(res, model_name, X_names, is_zi=False):
    """Pull coefficients into a dict."""
    params = res.params
    bse = res.bse
    pvals = res.pvalues

    if is_zi:
        # ZI models: first len(X_names) are inflate, next are count
        n_inflate = len(X_names)
        count_params = params[n_inflate:]
        count_se = bse[n_inflate:]
        count_pvals = pvals[n_inflate:]
        inflate_params = params[:n_inflate]
        inflate_pvals = pvals[:n_inflate]
    else:
        count_params = params
        count_se = bse
        count_pvals = pvals
        inflate_params = None
        inflate_pvals = None

    out = {
        "model": model_name,
        "ssw_coef": float(count_params[1]),
        "ssw_se": float(count_se[1]),
        "ssw_pvalue": float(count_pvals[1]),
        "ssw_irr": float(np.exp(count_params[1])),
        "z500_coef": float(count_params[2]),
        "z500_se": float(count_se[2]),
        "z500_pvalue": float(count_pvals[2]),
        "z500_irr": float(np.exp(count_params[2])),
        "intercept": float(count_params[0]),
    }

    if hasattr(res, "llf") and res.llf is not None:
        out["loglik"] = float(res.llf)
    if hasattr(res, "aic") and res.aic is not None:
        out["aic"] = float(res.aic)
    if hasattr(res, "bic") and res.bic is not None:
        out["bic"] = float(res.bic)
    if hasattr(res, "df_model"):
        out["n_params"] = int(res.df_model) + 1

    if is_zi and inflate_params is not None:
        out["inflate_ssw_coef"] = float(inflate_params[1])
        out["inflate_ssw_pvalue"] = float(inflate_pvals[1])

    return out


def extract_gee_results(res, model_name, X_names):
    params = res.params
    bse = res.bse
    pvals = res.pvalues

    return {
        "model": model_name,
        "ssw_coef": float(params[1]),
        "ssw_se": float(bse[1]),
        "ssw_pvalue": float(pvals[1]),
        "ssw_irr": float(np.exp(params[1])),
        "z500_coef": float(params[2]),
        "z500_se": float(bse[2]),
        "z500_pvalue": float(pvals[2]),
        "z500_irr": float(np.exp(params[2])),
        "intercept": float(params[0]),
        "scale": float(res.scale) if hasattr(res, "scale") else None,
    }


def print_summary(all_results, vuong_stat, vuong_p):
    hdr = f"{'Model':<22} {'SSW β':>8} {'SSW SE':>8} {'SSW P':>9} {'SSW IRR':>8} {'Z500 β':>8} {'Z500 P':>9} {'AIC':>10} {'BIC':>10}"
    sep = "─" * len(hdr)
    print(f"\n{sep}")
    print("  R57 — Count Model Comparison: SSW → Dry Natural Avalanches (Size ≥1)")
    print(sep)
    print(hdr)
    print(sep)

    for r in all_results:
        aic_s = f"{r['aic']:.1f}" if r.get("aic") else "—"
        bic_s = f"{r['bic']:.1f}" if r.get("bic") else "—"
        pstr = f"{r['ssw_pvalue']:.4f}" if r["ssw_pvalue"] >= 0.0001 else "<0.0001"
        z_pstr = f"{r['z500_pvalue']:.4f}" if r["z500_pvalue"] >= 0.0001 else "<0.0001"
        print(f"{r['model']:<22} {r['ssw_coef']:>8.4f} {r['ssw_se']:>8.4f} {pstr:>9} "
              f"{r['ssw_irr']:>8.3f} {r['z500_coef']:>8.4f} {z_pstr:>9} {aic_s:>10} {bic_s:>10}")
    print(sep)

    print(f"\n  Vuong test (ZIP vs ZINB):  Z = {vuong_stat:+.3f},  P = {vuong_p:.4f}")
    if vuong_p < 0.05:
        preferred = "ZINB" if vuong_stat < 0 else "ZIP"
        print(f"  → Significant at α=0.05: {preferred} preferred")
    else:
        print("  → Not significant: models not distinguishable")

    print()
    # Find best AIC/BIC among likelihood-based models
    lik_models = [r for r in all_results if r.get("aic")]
    if lik_models:
        best_aic = min(lik_models, key=lambda x: x["aic"])
        best_bic = min(lik_models, key=lambda x: x["bic"])
        print(f"  Best AIC: {best_aic['model']} ({best_aic['aic']:.1f})")
        print(f"  Best BIC: {best_bic['model']} ({best_bic['bic']:.1f})")

    # GEE results
    gee_models = [r for r in all_results if "GEE" in r["model"]]
    if gee_models:
        print(f"\n  GEE (cluster-robust, exchangeable within winters):")
        for g in gee_models:
            sig = "**" if g["ssw_pvalue"] < 0.05 else "ns"
            print(f"    {g['model']}: SSW P={g['ssw_pvalue']:.4f} {sig}, IRR={g['ssw_irr']:.3f}")

    print()


def main():
    wdf = load_data()
    y, X, X_names, groups = build_matrices(wdf)

    all_results = []

    # 1. Poisson
    print("\n── Fitting Poisson ──")
    pois_res = fit_poisson(y, X, X_names)
    r_pois = extract_results(pois_res, "Poisson", X_names)
    all_results.append(r_pois)
    print(f"   SSW P={r_pois['ssw_pvalue']:.4f}, AIC={r_pois['aic']:.1f}")

    # 2. Negative Binomial (NB2)
    print("\n── Fitting NB2 ──")
    nb_res = fit_nb(y, X, X_names)
    r_nb = extract_results(nb_res, "NegBin (NB2)", X_names)
    # get alpha (dispersion)
    alpha_idx = len(X_names)  # alpha is the last parameter in NB2
    if len(nb_res.params) > alpha_idx:
        r_nb["alpha"] = float(nb_res.params[alpha_idx])
    all_results.append(r_nb)
    print(f"   SSW P={r_nb['ssw_pvalue']:.4f}, AIC={r_nb['aic']:.1f}, alpha={r_nb.get('alpha', '?')}")

    # 3. ZIP
    print("\n── Fitting ZIP ──")
    zip_res = fit_zip(y, X, X_names)
    r_zip = extract_results(zip_res, "ZIP", X_names, is_zi=True)
    all_results.append(r_zip)
    print(f"   SSW P={r_zip['ssw_pvalue']:.4f}, AIC={r_zip['aic']:.1f}")

    # 4. ZINB
    print("\n── Fitting ZINB ──")
    zinb_res = fit_zinb(y, X, X_names)
    r_zinb = extract_results(zinb_res, "ZINB", X_names, is_zi=True)
    if len(zinb_res.params) > 2 * len(X_names):
        r_zinb["alpha"] = float(zinb_res.params[-1])
    all_results.append(r_zinb)
    print(f"   SSW P={r_zinb['ssw_pvalue']:.4f}, AIC={r_zinb['aic']:.1f}")

    # 5. Vuong test: ZIP vs ZINB
    print("\n── Vuong test (ZIP vs ZINB) ──")
    try:
        v_stat, v_p = vuong_test(zip_res, zinb_res, y, X)
    except Exception as e:
        print(f"   Vuong predict(which='prob') failed: {e}")
        print("   Falling back to LR-based comparison")
        lr = 2 * (zinb_res.llf - zip_res.llf)
        df_diff = 1  # ZINB has 1 extra param (alpha)
        v_p = float(stats.chi2.sf(lr, df_diff))
        v_stat = float(np.sqrt(lr) * np.sign(zinb_res.llf - zip_res.llf))
    print(f"   Vuong Z = {v_stat:+.3f}, P = {v_p:.4f}")

    # 6. GEE Poisson
    print("\n── Fitting GEE Poisson (exchangeable, winter clusters) ──")
    gee_pois = fit_gee_poisson(wdf)
    r_gee_p = extract_gee_results(gee_pois, "GEE-Poisson", X_names)
    all_results.append(r_gee_p)
    print(f"   SSW P={r_gee_p['ssw_pvalue']:.4f}, IRR={r_gee_p['ssw_irr']:.3f}")

    # 7. GEE NB
    print("\n── Fitting GEE NB (exchangeable, winter clusters) ──")
    try:
        gee_nb = fit_gee_nb(wdf)
        r_gee_nb = extract_gee_results(gee_nb, "GEE-NegBin", X_names)
        all_results.append(r_gee_nb)
        print(f"   SSW P={r_gee_nb['ssw_pvalue']:.4f}, IRR={r_gee_nb['ssw_irr']:.3f}")
    except Exception as e:
        print(f"   GEE-NB failed: {e}")

    # ── Summary ──
    print_summary(all_results, v_stat, v_p)

    # ── Overdispersion diagnostic ──
    pearson_chi2 = pois_res.pearson_chi2 if hasattr(pois_res, "pearson_chi2") else None
    n = len(y)
    k = len(X_names)
    if pearson_chi2:
        phi_hat = pearson_chi2 / (n - k)
        print(f"  Poisson overdispersion: Pearson χ²/df = {phi_hat:.1f}")
    else:
        phi_hat = y.var() / y.mean()
        print(f"  Overdispersion (var/mean): {phi_hat:.1f}")

    # ── LR test: Poisson vs NB ──
    lr_pois_nb = 2 * (nb_res.llf - pois_res.llf)
    # boundary test: 0.5 * chi2(0) + 0.5 * chi2(1)
    p_lr = 0.5 * stats.chi2.sf(lr_pois_nb, 1)
    print(f"  LR test Poisson vs NB: χ² = {lr_pois_nb:.1f}, P = {p_lr:.1e}")

    # ── Save ──
    output = {
        "description": "R57: NB/ZINB model comparison for SSW-avalanche association",
        "sample_n": int(n),
        "outcome": OUTCOME,
        "predictors": [SSW_COL, "z500_std (standardised ncep_z500_nh)"],
        "overdispersion_phi": float(phi_hat),
        "models": all_results,
        "vuong_test": {"z_statistic": v_stat, "p_value": v_p},
        "lr_test_poisson_vs_nb": {"chi2": float(lr_pois_nb), "p_value": float(p_lr)},
    }
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"  Results saved → {OUT_JSON}")


if __name__ == "__main__":
    main()
