#!/usr/bin/env python3
"""
r99_blocking_ssw_interaction.py
================================
Referee fix M2 — "Is there avalanche suppression beyond ordinary Alpine blocking?"

A single negative-binomial GLM of daily natural dry-slab avalanche counts on:
    blocking (B), SSW-window (S), and their interaction (B x S)
with winter fixed effects and a day-of-season quartic, plus a winter-block
bootstrap for cluster-robust CIs.

Parameterisation (log link):
    log E[count] = b0 + b_B*B + b_S*S + b_BS*(B*S) + quartic(dos) + C(winter_id)

Reported IRRs (exp of coefficient):
  - IRR_S   = exp(b_S)          : SSW effect on NON-blocking days
                                   ("does SSW suppress where there is no blocking?")
  - IRR_B   = exp(b_B)          : blocking effect on non-SSW days
  - IRR_BS  = exp(b_BS)         : how SSW MODIFIES blocking's potency (interaction)
  - IRR_S|B = exp(b_S + b_BS)   : total SSW effect on blocking days

Headline answer to the referee:
  IRR_S (SSW main effect controlling for blocking) and IRR_BS (interaction),
  each with a winter-block bootstrap 95% CI.

Blocking indicator follows the paper's operational definition
(scripts/analysis/56_blocking_without_ssw.py): standardized NCEP NH Z500
anomaly > 1 sigma, computed within winter-month climatology.

Outputs: data/results/r99_blocking_ssw_interaction.json
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "results" / "r99_blocking_ssw_interaction.json"
WINTER_MONTHS = [11, 12, 1, 2, 3, 4]
BLOCK_SIGMA = 1.0
N_BOOT = 1000
RNG = np.random.RandomState(20260610)


def build_frame():
    panel = pd.read_parquet(ROOT / "data/processed/analysis_panel.parquet")
    panel.index = (panel.index.tz_localize(None) if panel.index.tz is None
                   else panel.index.tz_convert("UTC").tz_localize(None))

    df = panel.copy()
    # --- blocking indicator: standardized winter-month Z500 anomaly > 1 sigma ---
    z = df["ncep_z500_nh"].copy()
    win = df.index.month.isin(WINTER_MONTHS)
    zc = z[win]
    clim = zc.groupby(zc.index.month).transform("mean")
    std = zc.groupby(zc.index.month).transform("std")
    z_anom = pd.Series(np.nan, index=df.index)
    z_anom.loc[zc.index] = (zc - clim) / std
    df["blocking"] = (z_anom > BLOCK_SIGMA).astype(float)
    df["z500anom"] = z_anom  # continuous standardized blocking strength

    # --- SSW window (+/-15 d), primary endpoint, controls ---
    df["ssw"] = df["ssw_within_15d"].astype(float)
    df["count"] = df["dry_natural_size_1234"]
    df["winter_id"] = pd.Categorical(df["winter_id"].astype(str)).codes

    # day-of-season quartic (standardised for numerical stability)
    dos = df["day_of_season"].astype(float)
    dos_z = (dos - dos.mean()) / dos.std()
    for k in range(1, 5):
        df[f"dos{k}"] = dos_z ** k

    # modelling sample: winter months, count observed, blocking defined, dos present
    # (day_of_season is NaN for early-November days; require it so all nested
    #  models share an identical sample -> coefficients are directly comparable)
    m = (win & df["count"].notna() & df["blocking"].notna()
         & df["ncep_z500_nh"].notna() & df["day_of_season"].notna())
    cols = ["count", "blocking", "z500anom", "ssw", "winter_id",
            "dos1", "dos2", "dos3", "dos4"]
    fr = df.loc[m, cols].copy()
    fr["count"] = fr["count"].round().astype(int)
    return fr


FORMULA = ("count ~ blocking * ssw + dos1 + dos2 + dos3 + dos4 + C(winter_id)")


def estimate_alpha(frame, formula=FORMULA):
    """Cameron-Trivedi auxiliary-regression estimate of NB dispersion alpha."""
    pois = smf.glm(formula, data=frame, family=sm.families.Poisson()).fit()
    mu = np.asarray(pois.fittedvalues)
    y = np.asarray(pois.model.endog)
    aux_y = ((y - mu) ** 2 - y) / mu
    alpha = float(sm.OLS(aux_y, mu).fit().params[0])
    return max(alpha, 1e-6)


def fit_nb(frame, formula=FORMULA):
    """Stable GLM negative-binomial fit with pre-estimated dispersion alpha."""
    alpha = estimate_alpha(frame, formula)
    res = smf.glm(formula, data=frame,
                  family=sm.families.NegativeBinomial(alpha=alpha)).fit()
    res._nb_alpha = alpha
    return res


def fit_poisson(frame):
    return smf.glm(FORMULA, data=frame, family=sm.families.Poisson()).fit()


def extract(res):
    p = res.params
    b_B = p.get("blocking", np.nan)
    b_S = p.get("ssw", np.nan)
    # interaction term name from patsy
    icand = [k for k in p.index if k.startswith("blocking:ssw") or k == "blocking:ssw"]
    b_BS = p[icand[0]] if icand else np.nan
    return dict(b_B=b_B, b_S=b_S, b_BS=b_BS, b_SgivenB=b_S + b_BS)


def winter_block_bootstrap(frame, n_boot=N_BOOT):
    """Resample winters with replacement; refit Poisson (fast) each draw."""
    winters = frame["winter_id"].unique()
    keys = ["b_B", "b_S", "b_BS", "b_SgivenB"]
    draws = {k: [] for k in keys}
    n_ok = 0
    for _ in range(n_boot):
        pick = RNG.choice(winters, size=len(winters), replace=True)
        parts, newid = [], 0
        for w in pick:
            sub = frame[frame["winter_id"] == w].copy()
            sub["winter_id"] = newid  # unique id so duplicated winters are distinct FE
            newid += 1
            parts.append(sub)
        bs = pd.concat(parts, ignore_index=True)
        try:
            r = fit_poisson(bs)
            e = extract(r)
            if any(not np.isfinite(e[k]) for k in keys):
                continue
            for k in keys:
                draws[k].append(e[k])
            n_ok += 1
        except Exception:
            continue
    ci = {}
    for k in keys:
        arr = np.array(draws[k])
        if arr.size > 10:
            lo, hi = np.percentile(arr, [2.5, 97.5])
            ci[k] = [float(np.exp(lo)), float(np.exp(hi))]
        else:
            ci[k] = [None, None]
    return ci, n_ok


def main():
    fr = build_frame()
    n = len(fr)
    nb_days = int(fr["blocking"].sum())
    ssw_days = int(fr["ssw"].sum())
    both = int(((fr["blocking"] == 1) & (fr["ssw"] == 1)).sum())
    print(f"Modelling sample: {n} winter-days | blocking days={nb_days} | "
          f"SSW-window days={ssw_days} | blocking&SSW={both}")

    # --- nested-model progression: shows how conditioning attenuates effects ---
    nested = {}
    nested_specs = {
        "M0_no_controls": "count ~ blocking * ssw",
        "M1_plus_season": "count ~ blocking * ssw + dos1 + dos2 + dos3 + dos4",
        "M2_plus_winterFE": FORMULA,
    }
    for name, f in nested_specs.items():
        try:
            r = fit_nb(fr, f)
            e = extract(r)
            nested[name] = {
                "IRR_S": float(np.exp(e["b_S"])),
                "IRR_B": float(np.exp(e["b_B"])),
                "IRR_BS": float(np.exp(e["b_BS"])),
                "alpha": float(r._nb_alpha),
            }
        except Exception as ex:
            nested[name] = {"error": str(ex)}
    print("\nNested-model progression (IRR_S = SSW effect | non-blocking):")
    for name in nested_specs:
        d = nested[name]
        if "error" not in d:
            print(f"  {name:18s}: IRR_S={d['IRR_S']:.3f}  IRR_B={d['IRR_B']:.3f}  IRR_BS={d['IRR_BS']:.3f}")

    # --- continuous-blocking sensitivity (uses full Z500 strength, not a 1-sigma cut) ---
    cont_f = "count ~ z500anom * ssw + dos1 + dos2 + dos3 + dos4 + C(winter_id)"
    try:
        rc = fit_nb(fr, cont_f)
        cc = rc.conf_int()
        continuous = {
            "IRR_S_ssw_at_mean_blocking": float(np.exp(rc.params["ssw"])),
            "IRR_S_95CI": [float(np.exp(cc.loc["ssw", 0])), float(np.exp(cc.loc["ssw", 1]))],
            "IRR_per_sigma_blocking": float(np.exp(rc.params["z500anom"])),
            "IRR_blocking_95CI": [float(np.exp(cc.loc["z500anom", 0])),
                                   float(np.exp(cc.loc["z500anom", 1]))],
        }
        print("\nContinuous-blocking model (blocking = standardized Z500 anomaly):")
        print(f"  IRR per +1sigma blocking = {continuous['IRR_per_sigma_blocking']:.3f} "
              f"95%CI {continuous['IRR_blocking_95CI']}")
        print(f"  IRR_S (SSW | mean blocking) = {continuous['IRR_S_ssw_at_mean_blocking']:.3f} "
              f"95%CI {continuous['IRR_S_95CI']}")
    except Exception as ex:
        continuous = {"error": str(ex)}

    # Primary = full model (M2), GLM-NB (stable, dispersion pre-estimated)
    nb = fit_nb(fr)
    nb_e = extract(nb)
    alpha = float(nb._nb_alpha)
    nb_ok = True
    cint = nb.conf_int()

    def mb(name):
        lo, hi = cint.loc[name]
        return [float(np.exp(lo)), float(np.exp(hi))]
    ibs_name = [k for k in nb.params.index if k.startswith("blocking:ssw")][0]
    model_ci = {"IRR_B": mb("blocking"), "IRR_S": mb("ssw"), "IRR_BS": mb(ibs_name)}

    # Poisson sensitivity (always)
    pois = fit_poisson(fr)
    pois_e = extract(pois)

    # Winter-block bootstrap CIs (cluster-robust)
    print(f"Running winter-block bootstrap (B={N_BOOT})...")
    ci, n_ok = winter_block_bootstrap(fr)
    print(f"  bootstrap successful draws: {n_ok}/{N_BOOT}")

    def irr(b):
        return float(np.exp(b)) if np.isfinite(b) else None

    result = {
        "description": "M2: blocking x SSW interaction NB-GLM for daily natural dry-slab counts",
        "endpoint": "dry_natural_size_1234 (natural-trigger dry slab, size 1-4)",
        "blocking_definition": f"standardized winter NH Z500 anomaly > {BLOCK_SIGMA} sigma",
        "ssw_window": "+/-15 d of SSW onset (ssw_within_15d)",
        "n_winter_days": n, "n_blocking_days": nb_days,
        "n_ssw_window_days": ssw_days, "n_blocking_and_ssw_days": both,
        "n_winters": int(fr["winter_id"].nunique()),
        "primary_model": "negative_binomial" if nb_ok else "poisson_fallback",
        "nb_dispersion_alpha": alpha,
        "IRR": {
            "IRR_B_blocking_effect_nonSSW": irr(nb_e["b_B"]),
            "IRR_S_ssw_effect_nonblocking": irr(nb_e["b_S"]),
            "IRR_BS_interaction": irr(nb_e["b_BS"]),
            "IRR_S_given_blocking": irr(nb_e["b_SgivenB"]),
        },
        "winter_block_bootstrap_95CI": {
            "IRR_B": ci["b_B"], "IRR_S": ci["b_S"],
            "IRR_BS": ci["b_BS"], "IRR_S_given_blocking": ci["b_SgivenB"],
        },
        "model_based_95CI_reference": model_ci,
        "nested_models": nested,
        "continuous_blocking_model": continuous,
        "poisson_sensitivity_IRR": {
            "IRR_B": irr(pois_e["b_B"]), "IRR_S": irr(pois_e["b_S"]),
            "IRR_BS": irr(pois_e["b_BS"]), "IRR_S_given_blocking": irr(pois_e["b_SgivenB"]),
        },
        "n_bootstrap_ok": n_ok,
    }

    # ---- interpretation ----
    irr_s = result["IRR"]["IRR_S_ssw_effect_nonblocking"]
    ci_s = ci["b_S"]
    irr_bs = result["IRR"]["IRR_BS_interaction"]
    ci_bs = ci["b_BS"]
    s_crosses = (ci_s[0] is None) or (ci_s[0] <= 1.0 <= ci_s[1])
    bs_crosses = (ci_bs[0] is None) or (ci_bs[0] <= 1.0 <= ci_bs[1])
    result["headline"] = {
        "ssw_suppression_beyond_blocking": {
            "IRR_S": irr_s, "CI": ci_s, "CI_crosses_1": bool(s_crosses),
            "reading": ("No statistically resolved SSW suppression beyond blocking "
                        "(CI includes 1)." if s_crosses else
                        "SSW retains suppression after controlling for blocking.")
        },
        "interaction": {
            "IRR_BS": irr_bs, "CI": ci_bs, "CI_crosses_1": bool(bs_crosses),
            "reading": ("SSW does not significantly modify blocking potency."
                        if bs_crosses else
                        ("SSW ATTENUATES blocking suppression (IRR_BS>1)" if irr_bs and irr_bs > 1
                         else "SSW STRENGTHENS blocking suppression (IRR_BS<1)"))
        },
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2, default=str)

    print("\n=== M2 RESULT ===")
    print(f"IRR_S  (SSW effect | non-blocking)      = {irr_s:.3f}  95%CI {ci_s}")
    print(f"IRR_B  (blocking effect | non-SSW)      = {result['IRR']['IRR_B_blocking_effect_nonSSW']:.3f}  95%CI {ci['b_B']}")
    print(f"IRR_BS (interaction, SSW x blocking)    = {irr_bs:.3f}  95%CI {ci_bs}")
    print(f"IRR_S|B(total SSW effect | blocking)    = {result['IRR']['IRR_S_given_blocking']:.3f}  95%CI {ci['b_SgivenB']}")
    print(f"NB dispersion alpha = {alpha}")
    print(f"\nHeadline: {result['headline']['ssw_suppression_beyond_blocking']['reading']}")
    print(f"          {result['headline']['interaction']['reading']}")
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    main()
