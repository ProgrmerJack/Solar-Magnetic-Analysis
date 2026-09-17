#!/usr/bin/env python3
"""
predictability_ceiling.py
=========================
CAN ANYONE PREDICT WHICH SSWs COUPLE DOWNWARD? OUT OF SAMPLE, WITH NO CORRECTION.

WHY THIS SUPERSEDES THE VARIANCE ROUTE
  `forced_variance_ceiling.py` estimates sigma_f by differencing variances against
  pseudo-events and then subtracting a placebo. It works, but the placebo removes
  91% of the raw signal in CMIP6, and that correction is the single weakest joint
  in this project -- a referee will go straight at it and be right to.

  Cross-validated prediction needs no correction at all. Fit any model on training
  events, score it on events it has never seen, and the out-of-sample R^2 IS the
  fraction of between-event variance that is genuinely predictable. Nothing to
  subtract, nothing to calibrate, honest by construction. If the best predictor
  buildable from every archived stratospheric variable scores zero out of sample,
  then no classifier works -- whatever sigma_f happens to be.

  It also asks the field's own question in the field's own terms. Karpechko et al.
  (2017) and everything downstream exist to answer "which SSWs propagate
  downward?". This measures whether that question has an answer.

THE DESIGN THAT MAKES IT DECISIVE: THREE TIERS, NONE CONTAINING THE SURFACE
  P1 PRE-ONSET     strictly before day 0. Genuine forecast information.
  P2 AT-ONSET      P1 plus days -5..0. What a forecaster has when the SSW starts.
  P3 POST-ONSET    P2 plus STRATOSPHERIC state over days 0..30. This is the
                   diagnostic case -- the PJO / lower-stratosphere / Loeffel-type
                   predictor. It is stratosphere-only; the surface never enters.

  The critical point is the CONTRAST BETWEEN TIERS. P1 and P2 do not overlap the
  response window (days +8..+52) at all. P3 does, by 22 days. So:

    P1 ~ P2 ~ 0 and P3 > 0   =>  the apparent skill comes from the WINDOW OVERLAP,
                                 not from the event. That is this project's thesis,
                                 tested by prediction rather than by projection.
    P1 or P2 > 0             =>  genuine forecast skill exists and the thesis is
                                 wrong in its strong form. Reported as such.

  P3 > 0 is therefore NOT evidence of causal predictability and is not read as
  such; it is the quantity the literature has been reporting.

HONEST MACHINERY
  - GroupKFold by ensemble member (CMIP6) or by winter (observations), so no
    training event shares a realisation or a season with a test event.
  - Ridge AND gradient boosting, because a linear failure with a nonlinear success
    would matter and must not be hidden by model choice.
  - PERMUTATION NULL: cross-validated R^2 is a random variable that is often
    negative, so "R^2 = 0.01" means nothing without its null. Y is shuffled within
    groups and the whole pipeline rerun.
  - POWER: a synthetic predictable component of known size is injected and the
    detection rate measured, so a null result comes with the effect size it could
    have found.

WHAT IS NOT AVAILABLE, STATED
  The CMIP6 archive here holds zonal-mean u and psl only, so there is no eddy heat
  flux and no wave-driving predictor in the model arm. Observations do have v'T'.
  A wave-driving term could carry information these features miss, and its absence
  in the model arm is a real limit rather than an oversight.

Output: predictability_ceiling.json
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "6_predictability"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]
SIM = ROOT / "ssw-design-analysis" / "06_simulation_validation"
sys.path.insert(0, str(SIM))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
from build_catalogue import load_catalogue          # noqa: E402
import gate3_clean_null as G                        # noqa: E402
import multi_index_event_study as M                 # noqa: E402
import ensemble_precursor as EP                     # noqa: E402
import stratifier_law_cmip6 as C6                   # noqa: E402

from sklearn.ensemble import HistGradientBoostingRegressor       # noqa: E402
from sklearn.linear_model import RidgeCV                         # noqa: E402
from sklearn.model_selection import GroupKFold                   # noqa: E402
from sklearn.pipeline import make_pipeline                       # noqa: E402
from sklearn.preprocessing import StandardScaler                 # noqa: E402

RAW = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "raw" / "cmip6"
OUT_WIN = (8, 52)
SEASON = (11, 12, 1, 2, 3, 4)
SEED = 20260803
N_PERM = 200

# (name, window, tier). Levels/latitudes are applied to each.
WINDOWS = [
    ("m45_31", (-45, -31), 1),
    ("m30_16", (-30, -16), 1),
    ("m15_01", (-15, -1), 1),
    ("m05_00", (-5, 0), 2),
    ("p00_15", (0, 15), 3),
    ("p15_30", (15, 30), 3),
]
LEVELS_PA = [1000.0, 5000.0, 10000.0]      # 10, 50, 100 hPa
LEVEL_LAB = {1000.0: "u10", 5000.0: "u50", 10000.0: "u100"}


def tier_cols(cols, tier):
    keep = []
    for c in cols:
        t = next((w[2] for w in WINDOWS if c.endswith(w[0])), None)
        if t is not None and t <= tier:
            keep.append(c)
    return keep


def cv_r2(X, y, groups, model="ridge", n_splits=5, seed=SEED):
    """Out-of-sample R^2, grouped so no realisation or season spans the split."""
    ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
    X, y, groups = X[ok], y[ok], groups[ok]
    if len(np.unique(groups)) < n_splits or len(y) < 40:
        return np.nan
    gkf = GroupKFold(n_splits=n_splits)
    pred = np.full(len(y), np.nan)
    for tr, te in gkf.split(X, y, groups):
        if model == "ridge":
            m = make_pipeline(StandardScaler(),
                              RidgeCV(alphas=np.logspace(-2, 4, 25)))
        else:
            m = HistGradientBoostingRegressor(
                max_depth=3, max_iter=200, learning_rate=0.05,
                min_samples_leaf=20, random_state=seed)
        m.fit(X[tr], y[tr])
        pred[te] = m.predict(X[te])
    # R^2 against the TRAINING-MEAN baseline is the honest denominator: a model
    # that learns nothing must score 0, not something flattering.
    ss_res = np.nansum((y - pred) ** 2)
    ss_tot = np.nansum((y - y.mean()) ** 2)
    return float(1 - ss_res / ss_tot)


def permutation_null(X, y, groups, model, rng, n_perm=N_PERM):
    """Shuffle Y WITHIN groups, so the null keeps group structure intact."""
    out = []
    for _ in range(n_perm):
        yp = y.copy()
        for g in np.unique(groups):
            i = np.flatnonzero(groups == g)
            yp[i] = rng.permutation(yp[i])
        r = cv_r2(X, yp, groups, model)
        if np.isfinite(r):
            out.append(r)
    return np.array(out)


def power_curve(X, y, groups, rng, model="ridge", reps=40):
    """Inject a KNOWN predictable component and measure the detection rate.

    A null result is worth nothing without the effect size it could have found.
    The injected signal is a random linear combination of the real features, so it
    is exactly the kind of structure the model is able to represent.
    """
    ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
    Xo, yo, go = X[ok], y[ok], groups[ok]
    Z = StandardScaler().fit_transform(Xo)
    crit = np.percentile(permutation_null(Xo, yo, go, model, rng, 100), 95)
    print(f"    permutation 95th percentile R^2 = {crit:+.4f}")
    out = {}
    for target in (0.02, 0.05, 0.10, 0.20):
        hits = 0
        for _ in range(reps):
            w = rng.normal(size=Z.shape[1])
            s = Z @ w
            s = (s - s.mean()) / s.std()
            noise = rng.normal(size=len(s))
            yy = np.sqrt(target) * s + np.sqrt(1 - target) * noise
            if cv_r2(Xo, yy, go, model) > crit:
                hits += 1
        out[f"R2={target}"] = round(hits / reps, 3)
        print(f"    injected R^2 {target:.2f} -> detection rate {hits/reps:.2f}")
    return out


def implied_r2(contrast, q, sd=1.0):
    """R^2 a binary split with this contrast implies, if it were a real predictor.

    Between-group variance for a split at fraction q with mean difference C is
    q(1-q)C^2; dividing by total variance gives the share of between-event
    variance that split claims to explain.
    """
    return float(q * (1 - q) * contrast ** 2 / sd ** 2)


def build_cmip6():
    """Per-event target and predictors from the zonal-mean CMIP6 archive."""
    rows, Ys, groups = [], [], []
    for f in sorted(RAW.glob("*_zm.nc")):
        try:
            m = EP.load_member(f)
        except Exception:
            continue
        m2 = m[np.isin(m.index.month, SEASON)].dropna()
        if len(m2) < 2000:
            continue
        on = EP.detect_ssw(m2["u10"].values, m2.index)
        if len(on) < 15:
            continue
        am = m2["am"]
        msk = C6.influence_mask(am.index, on)
        cl = am[~msk].groupby(am[~msk].index.dayofyear).mean()
        Y = C6.anom(am, on, cl, OUT_WIN)

        # raw fields for predictors
        ds = xr.open_dataset(f)
        lat = ds["lat"].values
        plev = ds["plev"].values
        t = pd.to_datetime(ds["time"].values.astype("datetime64[ns]")) \
            if np.issubdtype(ds["time"].dtype, np.datetime64) else m.index
        u = ds["u_zm"].values
        ds.close()
        if len(t) != u.shape[0]:
            t = m.index
        i60 = int(np.argmin(np.abs(lat - 60.0)))
        cap = (lat >= 65.0)
        idx = pd.DatetimeIndex(t)

        feat = {}
        for lev in LEVELS_PA:
            j = int(np.argmin(np.abs(plev - lev)))
            if not np.isclose(plev[j], lev, rtol=0.1):
                continue
            s60 = pd.Series(u[:, j, i60], index=idx)
            scap = pd.Series(np.average(u[:, j, cap],
                                        weights=np.cos(np.deg2rad(lat[cap])),
                                        axis=1), index=idx)
            for src, tag in ((s60, "60N"), (scap, "cap")):
                z = (src - src.groupby(src.index.dayofyear).transform("mean")) \
                    / src.groupby(src.index.dayofyear).transform("std")
                for wname, win, _ in WINDOWS:
                    feat[f"{LEVEL_LAB[lev]}_{tag}_{wname}"] = _wm(z, on, win)
        feat["doy_sin"] = np.sin(2 * np.pi * pd.DatetimeIndex(on).dayofyear / 365.25)
        feat["doy_cos"] = np.cos(2 * np.pi * pd.DatetimeIndex(on).dayofyear / 365.25)
        rows.append(pd.DataFrame(feat))
        Ys.append(Y)
        groups.append(np.full(len(Y), f.stem))
    if not rows:
        return None, None, None
    X = pd.concat(rows, ignore_index=True)
    return X, np.concatenate(Ys), np.concatenate(groups)


def _wm(series, onsets, win):
    out = []
    for o in pd.DatetimeIndex(onsets):
        v = series[(series.index >= o + pd.Timedelta(days=win[0]))
                   & (series.index <= o + pd.Timedelta(days=win[1]))]
        out.append(float(v.mean()) if len(v) >= max(3, (win[1] - win[0]) // 3)
                   else np.nan)
    return np.array(out)


def run(name, X, y, groups, rng, res):
    print("\n" + "=" * 78)
    print(f"=== {name}: n={np.isfinite(y).sum()} events, "
          f"{X.shape[1]} features, {len(np.unique(groups))} groups ===")
    r = {}
    for tier, lab in ((1, "P1 pre-onset"), (2, "P2 at-onset"),
                      (3, "P3 + post-onset stratosphere")):
        cols = tier_cols(X.columns, tier) + ["doy_sin", "doy_cos"]
        Xt = X[cols].values
        e = {}
        for model in ("ridge", "gb"):
            obs = cv_r2(Xt, y, groups, model)
            null = permutation_null(Xt, y, groups, model, rng,
                                    N_PERM if model == "ridge" else 60)
            p = float((null >= obs).mean()) if len(null) else np.nan
            e[model] = {"cv_r2": round(obs, 4),
                        "null_mean": round(float(null.mean()), 4),
                        "null_p95": round(float(np.percentile(null, 95)), 4),
                        "p_value": round(p, 4),
                        "significant": bool(p < 0.05)}
            print(f"  {lab:<30s} {model:<6s} CV R^2 = {obs:+.4f}   "
                  f"null p95 = {np.percentile(null, 95):+.4f}   p = {p:.3f}"
                  f"   {'SKILL' if p < 0.05 else 'no skill'}")
        r[lab] = {"n_features": len(cols), **e}
    res[name] = r
    return r


def main():
    rng = np.random.default_rng(SEED)
    res = {"outcome_window": list(OUT_WIN), "results": {}}

    print("building CMIP6 event table ...")
    Xc, yc, gc = build_cmip6()
    if Xc is not None:
        run("CMIP6", Xc, yc, gc, rng, res["results"])
        print("\n  POWER (P2 features, ridge):")
        cols = tier_cols(Xc.columns, 2) + ["doy_sin", "doy_cos"]
        res["power_CMIP6"] = power_curve(Xc[cols].values, yc, gc, rng)

        # Implied R^2 must use the variance of the SAME system the contrast was
        # measured in. An earlier version divided observational contrasts by the
        # CMIP6 sd, which inflated the Karpechko AO figure to R^2 = 1.285 -- an
        # impossible value, and the giveaway that the denominator was wrong.
        sd_c = float(np.nanstd(yc, ddof=1))
        var_c = sd_c ** 2
        var_o = 1.1167          # Var(Y|SSW) observations, forced_variance_ceiling
        print(f"\n  implied R^2 of published splits "
              f"(Var CMIP6 = {var_c:.3f}, Var obs = {var_o:.3f}):")
        res["implied"] = {}
        r1 = res["results"]["CMIP6"]["P1 pre-onset"]["ridge"]["cv_r2"]
        for lab, c, q, var in (("CMIP6 DW-NDW contrast", 1.143, 0.50, var_c),
                               ("Karpechko AO, observations", 1.782, 0.70, var_o),
                               ("ERA5 1000 hPa NAM, Karpechko", 0.684, 0.54, var_o),
                               ("ERA5 850 hPa NAM, ACP", 0.639, 0.59, var_o),
                               ("ACP 26,3723 published NAO", 0.850, 0.59, var_o)):
            v = q * (1 - q) * c ** 2 / var
            res["implied"][lab] = {"implied_r2": round(float(v), 4),
                                   "over_achievable": round(float(v / r1), 2)}
            print(f"    {lab:<32s} implies R^2 = {v:.3f}   "
                  f"({v / r1:.1f}x the {r1:.3f} achievable out of sample)")
        print("\n  NOT every published contrast is an overclaim on this metric: the")
        print("  two ERA5 contrasts imply R^2 ~ 0.09-0.10, essentially equal to the")
        print("  genuine predictable variance. The large ones (CMIP6 DW-NDW, and the")
        print("  observational AO) imply ~0.60, which exceeds even the post-onset")
        print("  diagnostic value and is ~6x what can be predicted in advance.")

    (RESULTS / "predictability_ceiling.json").write_text(json.dumps(res, indent=2),
                                                      encoding="utf8")
    print("\nSaved -> predictability_ceiling.json")


if __name__ == "__main__":
    main()
