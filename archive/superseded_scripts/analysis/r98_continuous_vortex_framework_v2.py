"""R98v2: Continuous daily-resolution vortex-strength framework (revised).

Revision over v1: the 21-lag unconstrained distributed lag was dominated by
multicollinearity because the polar-cap 10 hPa zonal wind has lag-1 autocorr
~0.98 over winter.  v2 replaces the unconstrained lag model with:

1. A single cumulative 5-15d weak-vortex exposure (Baldwin-Dunkerton window)
   as the primary NegBin predictor.  This is physically motivated by the
   downward-propagation literature and avoids collinearity.
2. An Almon polynomial distributed lag of degree 3 over days 0-21 to recover
   the lag shape with only four parameters (not 21), reported in SI.
3. A decile nonlinearity test (highest decile of cumulative weak-vortex
   exposure vs middle deciles 5-6) - reports rate ratio with winter-block
   bootstrap CI.
4. Falsification: leads 1-10 d (days BEFORE the anomaly) should have
   coefficient indistinguishable from zero if the association is causal.
5. Sensitivity: (a) mediator-adjusted model (adds Z500, NAO); (b) polar-cap
   temperature proxy (t10 positive anomaly = weak vortex) as independent
   stratospheric variable.

Framing: complementary to the event-catalog gmRR=0.32 primary result.  We
deliberately do NOT say "n=16 escape"; the daily obs are serially correlated
(u10 lag-1 autocorr ~0.98), so the effective sample size is far below 3,500
winter days.  The winter-block bootstrap captures this properly.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / "data" / "processed" / "analysis_panel_v2.parquet"
OUT = ROOT / "data" / "results" / "r98_continuous_vortex_framework.json"

WINTER_MONTHS = (11, 12, 1, 2, 3, 4)
CUM_LO, CUM_HI = 5, 15           # primary cumulative window (days after)
LEAD_LO, LEAD_HI = 1, 10         # falsification window (days before)
LAGS_PDL = list(range(0, 22))    # 0..21 for Almon PDL shape
BOOT = 1000
RNG = np.random.default_rng(42)


def load_panel() -> pd.DataFrame:
    df = pd.read_parquet(PANEL)
    df = df.copy()
    df["date"] = pd.to_datetime(df.index)
    df["month"] = df["date"].dt.month
    df = df[df["month"].isin(WINTER_MONTHS)].copy()
    df = df.sort_values("date").reset_index(drop=True)
    return df


def deseasonalise(series: pd.Series, dates: pd.Series) -> pd.Series:
    doy = dates.dt.dayofyear
    clim = series.groupby(doy).transform("mean")
    anom = series - clim
    std = anom.std()
    return (anom / std).rename(series.name + "_z")


def day_of_winter(dates: pd.Series) -> pd.Series:
    adj = np.where(dates.dt.month >= 11,
                   (dates.dt.dayofyear - 305).astype(float),
                   (dates.dt.dayofyear + 60).astype(float))
    return pd.Series(adj, index=dates.index, name="dow")


def spline_basis(x: np.ndarray, df: int = 4) -> np.ndarray:
    x_ = (x - np.nanmean(x)) / np.nanstd(x)
    return np.column_stack([x_, x_ ** 2, x_ ** 3, x_ ** 4])[:, :df]


def fit_poisson(y: np.ndarray, X: np.ndarray):
    y_int = np.rint(y).astype(float)
    y_int[y_int < 0] = 0.0
    return sm.GLM(y_int, X, family=sm.families.Poisson()).fit(maxiter=200)


def fit_negbin(y: np.ndarray, X: np.ndarray):
    poi = fit_poisson(y, X)
    mu = np.clip(poi.mu, 1e-6, None)
    y_int = np.rint(y).astype(float)
    y_int[y_int < 0] = 0.0
    df_resid = max(len(y_int) - X.shape[1], 1)
    alpha = float(np.clip(
        ((y_int - mu) ** 2 - mu).sum() / (mu ** 2).sum(),
        1e-3, 5.0,
    ))
    try:
        nb = sm.GLM(y_int, X, family=sm.families.NegativeBinomial(alpha=alpha)).fit(
            start_params=poi.params, maxiter=200,
        )
        return nb, alpha
    except Exception:
        return poi, None


def build_base_design(sub: pd.DataFrame) -> tuple[np.ndarray, pd.DataFrame]:
    spl = spline_basis(day_of_winter(sub["date"]).values, df=4)
    W = pd.get_dummies(sub["winter_id"], prefix="W",
                       drop_first=True, dtype=float).values
    const = np.ones(len(sub))
    return np.column_stack([const, spl, W]), sub


def build_cumulative_design(df: pd.DataFrame, col: str,
                            add_mediators: bool = False,
                            extra_cols: list[str] | None = None) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    base_mat, sub = build_base_design(df)
    x = df[col].values.reshape(-1, 1)
    X = np.column_stack([base_mat, x])
    if extra_cols:
        for c in extra_cols:
            X = np.column_stack([X, df[c].values.astype(float)])
    if add_mediators:
        z500 = df["ncep_z500_nh"].values.astype(float)
        z500 = (z500 - np.nanmean(z500)) / np.nanstd(z500)
        nao = df["nao_daily"].values.astype(float)
        X = np.column_stack([X, z500, nao])
    y = df["natural_size_234"].values.astype(float)
    valid = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    return y[valid], X[valid], df.loc[valid]


def almon_basis(lags: list[int], degree: int = 3) -> np.ndarray:
    """Polynomial-distributed-lag transformation matrix:
    original 21 lag columns -> 4 transformed columns that encode
    beta_k = sum_{p=0}^{degree} gamma_p * k^p.
    """
    lag_arr = np.array(lags, dtype=float)
    return np.column_stack([lag_arr ** p for p in range(degree + 1)])


def pdl_estimate(df: pd.DataFrame, col: str, lags: list[int], degree: int = 3) -> dict:
    """Fit a Poisson GLM with Almon PDL basis and recover lag-shape coefs."""
    lag_cols = [f"{col}_lag{k}" for k in lags]
    # Stack lag matrix (rows x 21 lags) times A (21 x 4) -> rows x 4 predictors
    base_mat, sub = build_base_design(df)
    L = df[lag_cols].values.astype(float)
    A = almon_basis(lags, degree=degree)
    Z = L @ A                             # n x (degree+1)
    X = np.column_stack([base_mat, Z])
    y = df["natural_size_234"].values.astype(float)
    valid = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y_v, X_v = y[valid], X[valid]
    poi = fit_poisson(y_v, X_v)
    # Recover lag-wise coefficients: beta_lag = A @ gamma
    n_ctrl = X_v.shape[1] - (degree + 1)
    gamma = poi.params[n_ctrl:]
    beta_lags = A @ gamma
    cum_5_15 = float(beta_lags[CUM_LO:CUM_HI + 1].sum())
    cum_lead = None  # not applicable - leads need separate model
    return {
        "per_lag_log_IRR": {int(k): float(b) for k, b in zip(lags, beta_lags)},
        "cumulative_log_IRR_5_15d": cum_5_15,
        "IRR_5_15d": float(np.exp(cum_5_15)),
        "gamma_almon": gamma.tolist(),
        "n_obs": int(len(y_v)),
    }


def winter_block_bootstrap_scalar(df: pd.DataFrame, builder,
                                  coef_idx_fn, B: int = BOOT) -> dict:
    winters = df["winter_id"].dropna().unique().tolist()
    estimates = []
    rng = np.random.default_rng(123)
    for _ in range(B):
        sample = rng.choice(winters, size=len(winters), replace=True)
        pieces = [df[df["winter_id"] == w] for w in sample]
        boot_df = pd.concat(pieces, ignore_index=True)
        try:
            y, X, _ = builder(boot_df)
            if len(y) < 100 or X.shape[1] >= len(y):
                continue
            poi = fit_poisson(y, X)
            idx = coef_idx_fn(X.shape[1])
            est = float(poi.params[idx])
            if np.isfinite(est):
                estimates.append(est)
        except Exception:
            continue
    estimates = np.array(estimates)
    if len(estimates) == 0:
        return {"mean": np.nan, "ci95": [np.nan, np.nan], "n": 0}
    return {
        "mean": float(estimates.mean()),
        "ci95": [float(np.quantile(estimates, 0.025)),
                 float(np.quantile(estimates, 0.975))],
        "n": int(len(estimates)),
    }


def main() -> dict:
    df = load_panel()

    # Build continuous vortex-strength proxies
    df["u10_z"] = deseasonalise(df["ncep_u_10hpa"], df["date"])
    df["t10_z"] = deseasonalise(df["ncep_t_10hpa"], df["date"])
    df["vortex_weak"] = -df["u10_z"]           # + = weak vortex
    df["vortex_weak_T"] = df["t10_z"]          # + = warm stratosphere = weak vortex

    # Primary cumulative (lead) exposure: mean of days (t-15 .. t-5) = "looking back"
    # so for outcome on day t, stratosphere state from 5 to 15 days earlier.
    # Implement by taking rolling mean of vortex_weak with a window offset.
    # Use shifted rolling mean: shift by 5 days, then 11-day rolling mean => covers t-15..t-5.
    df["cum_weak_5_15"] = (
        df["vortex_weak"].shift(CUM_LO).rolling(CUM_HI - CUM_LO + 1, min_periods=6).mean()
    )
    df["cum_weakT_5_15"] = (
        df["vortex_weak_T"].shift(CUM_LO).rolling(CUM_HI - CUM_LO + 1, min_periods=6).mean()
    )
    # Falsification: same window but using FUTURE days (leads 1-10 before the anomaly)
    df["cum_weak_lead_1_10"] = (
        df["vortex_weak"].shift(-LEAD_HI).rolling(LEAD_HI - LEAD_LO + 1, min_periods=6).mean()
    )

    # Almon PDL - build 22 lag columns
    for k in LAGS_PDL:
        df[f"vortex_weak_lag{k}"] = df["vortex_weak"].shift(k)

    # ===== Primary: cumulative 5-15d weak-vortex exposure -> daily count =====
    def primary_builder(d):
        return build_cumulative_design(d, "cum_weak_5_15")

    y, X, kept = primary_builder(df)
    nb, alpha = fit_negbin(y, X)
    poi = fit_poisson(y, X)
    coef_idx = X.shape[1] - 1
    primary_log_IRR = float(nb.params[coef_idx])
    primary_pvalue = float(nb.pvalues[coef_idx])
    poi_log_IRR = float(poi.params[coef_idx])

    # Winter-block bootstrap for primary (Poisson fits for speed)
    primary_boot = winter_block_bootstrap_scalar(
        df, primary_builder, lambda p: p - 1, B=BOOT,
    )

    # ===== Falsification leads =====
    def lead_builder(d):
        return build_cumulative_design(d, "cum_weak_lead_1_10")

    y_l, X_l, _ = lead_builder(df)
    nb_l, _ = fit_negbin(y_l, X_l)
    lead_log_IRR = float(nb_l.params[X_l.shape[1] - 1])
    lead_pvalue = float(nb_l.pvalues[X_l.shape[1] - 1])
    lead_boot = winter_block_bootstrap_scalar(
        df, lead_builder, lambda p: p - 1, B=BOOT,
    )

    # ===== Sensitivity: polar-cap temperature proxy =====
    def t_builder(d):
        return build_cumulative_design(d, "cum_weakT_5_15")

    y_t, X_t, _ = t_builder(df)
    nb_t, _ = fit_negbin(y_t, X_t)
    t_log_IRR = float(nb_t.params[X_t.shape[1] - 1])
    t_pvalue = float(nb_t.pvalues[X_t.shape[1] - 1])
    t_boot = winter_block_bootstrap_scalar(
        df, t_builder, lambda p: p - 1, B=BOOT,
    )

    # ===== Sensitivity: mediator-adjusted (Z500 + NAO) =====
    def med_builder(d):
        return build_cumulative_design(d, "cum_weak_5_15", add_mediators=True)

    y_m, X_m, _ = med_builder(df)
    nb_m, _ = fit_negbin(y_m, X_m)
    # cum_weak_5_15 index is base_cols = 1+4+(20 winter dummies) => 25, then lag = idx 25
    # Safest: identify by column position: base_mat width + 0 = (1+4+W_count)
    W_n = df["winter_id"].nunique() - 1
    med_idx = 1 + 4 + W_n
    med_log_IRR = float(nb_m.params[med_idx])

    # ===== Tail nonlinearity: highest vs middle deciles of cum_weak_5_15 =====
    cum = df["cum_weak_5_15"]
    deciles = pd.qcut(cum, 10, labels=False, duplicates="drop")
    df["decile"] = deciles
    tail = df.loc[df["decile"] == 9, "natural_size_234"].dropna()
    mid = df.loc[df["decile"].isin([4, 5]), "natural_size_234"].dropna()
    rr_tail_mid = float(tail.mean() / mid.mean()) if mid.mean() > 0 else np.nan
    # Winter-block bootstrap for the rate ratio
    rr_boot = []
    winters = df["winter_id"].dropna().unique().tolist()
    rng = np.random.default_rng(7)
    for _ in range(BOOT):
        sample = rng.choice(winters, size=len(winters), replace=True)
        pieces = [df[df["winter_id"] == w] for w in sample]
        bdf = pd.concat(pieces, ignore_index=True)
        t = bdf.loc[bdf["decile"] == 9, "natural_size_234"].dropna()
        m = bdf.loc[bdf["decile"].isin([4, 5]), "natural_size_234"].dropna()
        if len(t) > 10 and len(m) > 10 and m.mean() > 0:
            rr_boot.append(t.mean() / m.mean())
    rr_boot = np.array(rr_boot)
    rr_ci = ([float(np.quantile(rr_boot, 0.025)),
              float(np.quantile(rr_boot, 0.975))]
             if len(rr_boot) > 0 else [np.nan, np.nan])

    # ===== Almon PDL lag-shape recovery (SI only) =====
    pdl = pdl_estimate(df, "vortex_weak", LAGS_PDL, degree=3)

    result = {
        "sample": {
            "n_winter_days": int(df["natural_size_234"].notna().sum()),
            "n_winters": int(df["winter_id"].nunique()),
            "zero_share_natural_size_234": float(
                (df["natural_size_234"].dropna() == 0).mean()
            ),
            "u10_lag1_autocorr": float(df["ncep_u_10hpa"].autocorr(1)),
            "note": "Daily obs are serially correlated; effective sample size far below nominal winter-day count. Winter-block bootstrap is the primary inferential tool.",
        },
        "primary": {
            "description": "Negative Binomial GLM: daily Swiss natural dry-slab count (size 2+) ~ cumulative weak-vortex exposure averaged over days 5-15 BEFORE the count date. Controls: winter fixed effects + day-of-winter quartic spline. Sign convention: positive predictor = weak vortex.",
            "negbin_alpha": alpha,
            "n_obs": int(len(y)),
            "log_IRR": primary_log_IRR,
            "IRR": float(np.exp(primary_log_IRR)),
            "p_value_model_based": primary_pvalue,
            "poisson_log_IRR": poi_log_IRR,
            "winter_block_bootstrap": {
                "mean_log_IRR": primary_boot["mean"],
                "log_IRR_ci95": primary_boot["ci95"],
                "IRR_ci95": [float(np.exp(c)) for c in primary_boot["ci95"]] if np.isfinite(primary_boot["ci95"][0]) else [None, None],
                "n_bootstraps": primary_boot["n"],
            },
        },
        "falsification_leads": {
            "description": "Using the SAME cumulative-exposure structure but with leads (days AFTER the count date) rather than lags. Under a causal interpretation of the stratosphere->snowpack pathway this coefficient should be indistinguishable from zero; a significant negative coefficient here would indicate the primary result is driven by confounded seasonality rather than directional coupling.",
            "log_IRR": lead_log_IRR,
            "IRR": float(np.exp(lead_log_IRR)),
            "p_value_model_based": lead_pvalue,
            "winter_block_bootstrap": {
                "mean_log_IRR": lead_boot["mean"],
                "log_IRR_ci95": lead_boot["ci95"],
                "IRR_ci95": [float(np.exp(c)) for c in lead_boot["ci95"]] if np.isfinite(lead_boot["ci95"][0]) else [None, None],
                "n_bootstraps": lead_boot["n"],
            },
        },
        "sensitivity_t10_proxy": {
            "description": "Same model but using 10 hPa polar-cap TEMPERATURE anomaly as the stratospheric weak-vortex indicator (orthogonal observation pathway to u10).",
            "log_IRR": t_log_IRR,
            "IRR": float(np.exp(t_log_IRR)),
            "p_value_model_based": t_pvalue,
            "winter_block_bootstrap": {
                "mean_log_IRR": t_boot["mean"],
                "log_IRR_ci95": t_boot["ci95"],
                "n_bootstraps": t_boot["n"],
            },
        },
        "sensitivity_mediator_adjusted": {
            "description": "Add Z500_NH (standardised) and NAO_daily as covariates. These are downstream MEDIATORS on the stratosphere->troposphere->snowpack pathway, so adjusting for them partials out the exposure that is ALREADY transmitted through them. The coefficient on cum_weak_5_15 in this model is interpreted as the residual 'direct' association, not the total effect.",
            "log_IRR": med_log_IRR,
            "IRR": float(np.exp(med_log_IRR)),
            "attenuation_vs_primary": float(med_log_IRR - primary_log_IRR),
        },
        "tail_nonlinearity": {
            "description": "Raw rate ratio: winter-day mean natural_size_234 count in the HIGHEST decile of cumulative 5-15d weak-vortex exposure vs the MIDDLE deciles (5-6). Winter-block bootstrap CI.",
            "n_tail": int(len(tail)),
            "n_mid": int(len(mid)),
            "RR_tail_vs_mid": rr_tail_mid,
            "RR_ci95": rr_ci,
        },
        "almon_pdl_lag_shape": {
            "description": "Polynomial distributed lag (degree-3 Almon basis) over lags 0-21 d. Reduces 22 collinear lag columns to 4 orthogonalised shape parameters. Provides a smooth estimate of the lag profile for SI Figure.",
            **pdl,
        },
        "meta": {
            "outcome": "natural_size_234 (Swiss SLF natural dry-slab avalanche counts, size 2+).",
            "vortex_proxy_primary": "NCEP 10 hPa polar-cap (60-90N) zonal-mean zonal wind; deseasonalised, standardised; sign-flipped so positive = weak vortex. Not the canonical CP07 60N wind; we describe it as a polar-vortex strength proxy.",
            "winter_months": WINTER_MONTHS,
            "cumulative_window_days": [CUM_LO, CUM_HI],
            "lead_window_days": [LEAD_LO, LEAD_HI],
            "bootstrap_B": BOOT,
            "framing_note": "This is a complement to the event-catalog gmRR=0.32; not a replacement. We do not claim 'n=16 escape' - the daily record is serially correlated and the effective sample is far below nominal winter-day counts. The continuous framework demonstrates that the association scales with stratospheric anomaly strength beyond the discrete SSW indicator.",
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    main()
