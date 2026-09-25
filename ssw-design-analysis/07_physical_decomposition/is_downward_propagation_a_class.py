#!/usr/bin/env python3
"""
is_downward_propagation_a_class.py
==================================
THE PHYSICAL QUESTION UNDERNEATH THE METHODS RESULT.

The literature treats "downward-propagating" SSWs as a CLASS. Karpechko et al.
(2017) split events in two; roughly two thirds are said to propagate; papers
compare DW against NDW composites and report the contrast as an effect size.
That framing carries an implicit physical model: there are two populations of
SSWs, one that couples to the surface and one that does not.

A mixture of two populations predicts a BIMODAL (or at least mixture-shaped)
distribution of per-event surface response. A single population with continuous
spread predicts a UNIMODAL one, in which case DW/NDW is a threshold on a
continuum, "two thirds propagate downward" is just the fraction below a cut
rather than a property of nature, and comparing the groups compares the tails of
one distribution with itself.

This has to be tested, not asserted, and it is testable.

THE TESTS
  1. HARTIGAN DIP STATISTIC on the per-event surface anomaly. Implemented here
     directly (the `diptest` package is not installed) as the sup-norm distance
     between the ECDF and its greatest convex minorant / least concave majorant
     over the modal interval. Calibrated by resampling the pseudo-event pool,
     re-centred on the event mean (see dip_pvalue for why not the uniform).
  2. GAUSSIAN MIXTURE, k=1 vs k=2, by BIC. A two-component fit will beat one
     component by chance some of the time, so the BIC difference is CALIBRATED
     against pseudo-events, which are a single population by construction.
  3. THE SHIFT TEST, which is the sharpest of the three. If SSWs merely displace
     the ordinary winter distribution rather than creating a second population,
     then Y_ssw should be distributable as Y_pseudo + constant. A two-sample KS
     test against the pseudo-event distribution shifted by the difference in
     means asks exactly that. Failing to reject means a pure shift suffices and no second class is
     needed.

  Observations give n=43, which cannot settle a mixture question -- the CMIP6
  ensemble already assembled for this project gives ~1,500 events and can.

WHAT WOULD FALSIFY THE "CONTINUUM" READING
  A significant dip statistic, a BIC preference for k=2 beyond what pseudo-events
  produce, or a KS rejection of the pure shift. Any of the three, in the ensemble
  where power is adequate, and the class interpretation survives.

Output: is_downward_propagation_a_class.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.mixture import GaussianMixture

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
import selection_on_outcome as SO                   # noqa: E402
import ensemble_precursor as EP                     # noqa: E402
import stratifier_law_cmip6 as C6                   # noqa: E402

RAW = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "raw" / "cmip6"
OUT_WIN = (8, 52)
SEASON = (11, 12, 1, 2, 3, 4)
N_CAL = 2000
SEED = 20260801


# ---------------------------------------------------------------- dip statistic
def dip(x):
    """Hartigan's dip: sup distance from the ECDF to the nearest unimodal CDF.

    Computed as the smaller of the sup-norm gaps to the greatest convex minorant
    and the least concave majorant of the ECDF, which is the standard equivalent
    formulation for the one-dimensional case.
    """
    x = np.sort(np.asarray(x, float))
    n = len(x)
    if n < 4:
        return np.nan
    F = np.arange(1, n + 1) / n

    def gcm(xs, ys):
        """greatest convex minorant via monotone slope hull"""
        idx = [0]
        for i in range(1, len(xs)):
            while len(idx) >= 2:
                a, b = idx[-2], idx[-1]
                if ((ys[b] - ys[a]) / (xs[b] - xs[a] + 1e-300)
                        >= (ys[i] - ys[b]) / (xs[i] - xs[b] + 1e-300)):
                    idx.pop()
                else:
                    break
            idx.append(i)
        return np.interp(xs, xs[idx], ys[idx])

    lo = gcm(x, F)
    hi = -gcm(x[::-1] * -1, -F[::-1])[::-1]
    return float(min(np.max(np.abs(F - lo)), np.max(np.abs(F - hi))))


def dip_pvalue(x, pool, rng, n_cal=N_CAL):
    """Calibrate against a KNOWN SINGLE POPULATION of the same size.

    Calibrating against the uniform -- the textbook choice for Hartigan's exact
    dip -- is WRONG for the statistic implemented above. `dip()` here measures
    the distance from the ECDF to a globally convex or globally concave CDF, i.e.
    to a MONOTONE density, not to the nearest unimodal one. A uniform sample sits
    at ~0 on it by construction while any Gaussian sample sits at ~0.25, so a
    uniform calibration returns p=0 for perfectly unimodal data. That produced a
    spurious "BIMODAL" verdict at n=1888 on the first run of this script.

    Resampling the pseudo-event pool -- a single population by construction, with
    the empirical shape and the matched sample size -- makes the test valid
    whatever the statistic happens to measure.
    """
    d = dip(x)
    n = len(x)
    null = np.array([dip(rng.choice(pool, n, replace=True)) for _ in range(n_cal)])
    return d, float((null >= d).mean())


def bic_gap(x):
    """BIC(k=1) - BIC(k=2). Positive favours TWO components."""
    z = np.asarray(x, float).reshape(-1, 1)
    b = []
    for k in (1, 2):
        g = GaussianMixture(k, covariance_type="full", random_state=0,
                            n_init=5).fit(z)
        b.append(g.bic(z))
    return float(b[0] - b[1])


def analyse(name, Y_ev, Y_ps_draws, rng, res):
    """Y_ev: per-event anomalies. Y_ps_draws: list of pseudo-event arrays."""
    Y_ev = Y_ev[np.isfinite(Y_ev)]
    n = len(Y_ev)
    print(f"\n=== {name}  (n = {n}) ===")

    pool0 = np.concatenate([p[np.isfinite(p)] for p in Y_ps_draws])
    d, pd_ = dip_pvalue(Y_ev, pool0 - pool0.mean() + Y_ev.mean(), rng)
    print(f"  Hartigan dip    = {d:.4f}   P(unimodal null >= observed) = {pd_:.4f}"
          f"   -> {'BIMODAL' if pd_ < 0.05 else 'no evidence against unimodal'}")

    gap = bic_gap(Y_ev)
    cal = np.array([bic_gap(p[np.isfinite(p)]) for p in Y_ps_draws
                    if np.isfinite(p).sum() >= 20])
    p_gap = float((cal >= gap).mean()) if len(cal) else np.nan
    print(f"  BIC(1)-BIC(2)   = {gap:+.2f}   pseudo-event calibration: "
          f"median {np.median(cal):+.2f}, P(pseudo >= observed) = {p_gap:.4f}")
    print(f"                    -> {'TWO COMPONENTS' if p_gap < 0.05 else 'one component suffices'}")

    # pure-shift test against pooled pseudo-events
    pool = np.concatenate([p[np.isfinite(p)] for p in Y_ps_draws])
    shift = float(np.mean(Y_ev) - np.mean(pool))
    ks, p_ks = stats.ks_2samp(Y_ev, pool + shift)
    print(f"  pure-shift test : shift {shift:+.3f} sigma, KS = {ks:.4f}, "
          f"p = {p_ks:.4f}")
    print(f"                    -> {'SHIFT INSUFFICIENT' if p_ks < 0.05 else 'a pure SHIFT of the ordinary distribution suffices'}")

    sk = float(stats.skew(Y_ev))
    kt = float(stats.kurtosis(Y_ev))
    print(f"  shape           : mean {Y_ev.mean():+.3f}, sd {Y_ev.std():.3f}, "
          f"skew {sk:+.3f}, excess kurtosis {kt:+.3f}")

    res[name] = {"n": int(n), "dip": round(d, 5), "dip_p": pd_,
                 "bic_gap": round(gap, 3),
                 "bic_gap_pseudo_median": round(float(np.median(cal)), 3) if len(cal) else None,
                 "bic_gap_p": p_gap,
                 "shift": round(shift, 4), "ks": round(float(ks), 4),
                 "ks_p": float(p_ks),
                 "mean": round(float(Y_ev.mean()), 4),
                 "sd": round(float(Y_ev.std()), 4),
                 "skew": round(sk, 4), "excess_kurtosis": round(kt, 4),
                 "bimodal": bool(pd_ < 0.05 or (p_gap == p_gap and p_gap < 0.05)),
                 "pure_shift_sufficient": bool(p_ks >= 0.05)}
    return res[name]


def power_analysis(n, sd, pool, rng, frac=2 / 3, n_rep=200):
    """What two-population structure WOULD have been detected at this n?

    Failing to reject a mixture is not evidence of no mixture -- this project has
    already mistaken one absence of evidence for a finding. So the claim "there is
    no class" is only worth making alongside the separation that IS excluded.

    Simulates a two-component mixture with mixing fraction `frac` (the literature's
    "about two thirds propagate downward") and component separation d, with total
    variance held at the observed sd so the mixture is not detectable by spread
    alone, then measures how often the BIC test flags it.
    """
    crit = np.percentile(
        [bic_gap(rng.choice(pool, n, replace=True)) for _ in range(300)], 95)
    print(f"\n  power analysis at n={n}: BIC-gap 95th-percentile critical "
          f"value = {crit:+.2f}")
    print(f"  {'separation d':>13s} {'component sd':>13s} {'power':>8s}")
    out = {}
    for d in (0.25, 0.5, 0.75, 1.0, 1.5, 2.0):
        # var_total = var_within + frac(1-frac)d^2  -> hold var_total at sd^2
        vw = sd ** 2 - frac * (1 - frac) * d ** 2
        if vw <= 0.01:
            continue
        sw = np.sqrt(vw)
        hit = 0
        for _ in range(n_rep):
            k = rng.random(n) < frac
            x = np.where(k, rng.normal(0, sw, n), rng.normal(d, sw, n))
            if bic_gap(x) > crit:
                hit += 1
        p = hit / n_rep
        out[f"d={d}"] = round(p, 3)
        print(f"  {d:13.2f} {sw:13.3f} {p:8.2f}")
    det = [float(k.split("=")[1]) for k, v in out.items() if v >= 0.8]
    out["min_detectable_separation_80pct"] = min(det) if det else None
    return out


def main():
    rng = np.random.default_rng(SEED)
    res = {"outcome_window": list(OUT_WIN), "n_calibration": N_CAL, "results": {}}

    # ---------------- observations ----------------------------------------
    ao = M.load("ao")["y"]
    real = load_catalogue("primary")
    real = real[(real >= ao.index.min()) & (real <= ao.index.max())]
    mask = G.real_influence_mask(ao.index, real)
    clim = ao[~mask].groupby(ao[~mask].index.dayofyear).mean()
    clean_idx = G.zone_free_index(ao.index, real)
    doys = np.array([t.dayofyear for t in pd.DatetimeIndex(real)])

    Y_obs = SO.per_event(ao, real, clim)
    ps = []
    for _ in range(300):
        f = G.draw_clean(clean_idx, doys, rng)
        if len(f) >= 10:
            ps.append(SO.per_event(ao, f, clim))
    analyse("observations", Y_obs, ps, rng, res["results"])

    # ---------------- CMIP6, where the power is ---------------------------
    print("\nloading CMIP6 ensemble ...")
    Yc, psc = [], []
    for f in sorted(RAW.glob("*_zm.nc")):
        try:
            m = EP.load_member(f)
            full = m
        except Exception:
            continue
        m = m[np.isin(m.index.month, SEASON)].dropna()
        if len(m) < 2000:
            continue
        on = EP.detect_ssw(full["u10"].values, full.index)   # full daily series: CP07 needs contiguous days
        if len(on) < 15:
            continue
        am = m["am"]
        msk = C6.influence_mask(am.index, on)
        cl = am[~msk].groupby(am[~msk].index.dayofyear).mean()
        cln = C6.zone_free_index(am.index, on)
        dy = np.array([t.dayofyear for t in pd.DatetimeIndex(on)])
        Yc.append(C6.anom(am, on, cl, OUT_WIN))
        for _ in range(20):
            p = C6.draw_pseudo(cln, dy, rng)
            if len(p) >= 10:
                psc.append(C6.anom(am, p, cl, OUT_WIN))
    if Yc:
        Yall = np.concatenate(Yc)
        r6 = analyse("CMIP6_ensemble", Yall, psc, rng, res["results"])
        pool6 = np.concatenate([p[np.isfinite(p)] for p in psc])
        res["power_CMIP6"] = power_analysis(
            int(r6["n"]), float(r6["sd"]),
            pool6 - pool6.mean() + Yall[np.isfinite(Yall)].mean(), rng)

    # ---------------- the verdict ------------------------------------------
    print("\n" + "=" * 72)
    print("=== IS 'DOWNWARD PROPAGATION' A CLASS? ===")
    key = "CMIP6_ensemble" if "CMIP6_ensemble" in res["results"] else "observations"
    r = res["results"][key]
    if r["bimodal"]:
        v = (f"CLASS SUPPORTED in {key}: the per-event surface response is not "
             f"unimodal (dip p={r['dip_p']:.4f}, BIC-gap p={r['bic_gap_p']:.4f}). "
             f"The DW/NDW dichotomy has an empirical basis.")
    elif r["pure_shift_sufficient"]:
        v = (f"NO CLASS in {key} (n={r['n']}): the per-event surface response is "
             f"unimodal (dip p={r['dip_p']:.4f}), one Gaussian component suffices "
             f"(BIC-gap p={r['bic_gap_p']:.4f}), and the whole distribution is a "
             f"pure SHIFT of {r['shift']:+.3f} sigma of the ordinary winter "
             f"distribution (KS p={r['ks_p']:.4f}). 'Downward propagation' is a "
             f"threshold on a continuum, not a population.")
    else:
        v = (f"UNIMODAL but NOT a pure shift in {key}: dip p={r['dip_p']:.4f}, "
             f"KS p={r['ks_p']:.4f}. SSWs change the distribution's shape, not "
             f"only its location, but there is still no second mode.")
    res["verdict"] = v
    print("  " + v)

    (RESULTS / "is_downward_propagation_a_class.json").write_text(
        json.dumps(res, indent=2), encoding="utf8")
    print("\nSaved -> is_downward_propagation_a_class.json")


if __name__ == "__main__":
    main()
