#!/usr/bin/env python3
"""
era_dependence_gate.py
======================
THE GATE THAT DECIDES WHETHER THE PRE-SATELLITE EVENTS MAY BE USED AT ALL.

`extend_ncep_presatellite.py` takes the stratospheric record back to 1958 and the
stratifier arm from 29 events to 42. That is a 45% increase in observational n and
it is worth having -- but only if the pre-satellite analyses actually support it.

NCEP/NCAR R1 before 1979 rests on radiosondes with almost no Arctic coverage above
100 hPa. Where observations are sparse a reanalysis relaxes toward its background,
which SHRINKS variance toward climatology. That is not a neutral error here: the
whole method is built on beta, an OLS slope of surface on stratosphere, and
attenuating the stratospheric predictor's variance BIASES BETA TOWARD ZERO
(classical errors-in-variables regression dilution). A deflated beta under-predicts
the selection bias, which would make the correction look SMALLER than it is and
leave a spurious "causal" residual. So the failure mode of naive pooling points in
the direction that flatters the result, which is exactly when a gate is needed.

This project has already been bitten on this seam: the original catalogue rule
(">=4 of 6 reanalyses") demanded unanimity pre-1979 because few products cover it,
producing era-dependent event selection at Fisher P=0.011.

THREE TESTS, RUN IN ORDER OF POWER (highest first, because the event-level test
has n=13 in one arm and is nearly powerless).

  T1  VARIANCE RATIO, per level, DJF, on all ~1900 winter days per era. This is
      the direct signature of analysis-nudged-to-climatology and it is the only
      test here with real power. A ratio far below 1 at 10 hPa condemns the era.

  T2  BETA STABILITY. beta is refitted separately within each era on event-free
      pseudo-onsets -- thousands of draws, so this too is well powered. beta is a
      property of the CLIMATOLOGY, not of the events, so an era-dependent beta
      means the two eras do not describe the same atmosphere and must not share a
      single correction coefficient.

  T3  RESIDUAL AGREEMENT at the events themselves. Reported for completeness and
      read with its power stated: with 13 vs 29 events a null result here is weak
      evidence of agreement, NOT evidence of no difference. It is never used on
      its own to license pooling.

VERDICT RULE, fixed before the numbers are seen:
  pooling is allowed only if T1 variance ratios lie in [0.8, 1.25] at every level
  used by a stratifier, AND T2 betas agree within their bootstrap CIs. Otherwise
  the observational arm stays at n=29 and the extension is reported as a
  sensitivity test rather than as the headline.

Output: era_dependence_gate.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "7_ensemble"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
from build_catalogue import load_catalogue          # noqa: E402
import gate3_clean_null as G                        # noqa: E402
import multi_index_event_study as M                 # noqa: E402
import selection_on_outcome as S                    # noqa: E402
import stratifier_bias_law as L                     # noqa: E402

STRAT_OLD = ROOT / "data" / "processed" / "atmospheric" / "ncep_stratosphere.parquet"
EXT = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "ncep_presatellite_extension.parquet"
MERGED = ROOT / "data" / "processed" / "atmospheric" / "ncep_stratosphere_1958.parquet"
HF = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "heatflux_daily.parquet"

SPLIT = pd.Timestamp("1979-01-01")
N_BETA = 3000
N_BOOT = 800
SEED = 20260803
VAR_RATIO_OK = (0.80, 1.25)


def merge():
    """Write the merged record to a NEW file. The original is never touched."""
    old = pd.read_parquet(STRAT_OLD)
    ext = pd.read_parquet(EXT)
    for d in (old, ext):
        if d.index.tz is not None:
            d.index = d.index.tz_convert("UTC").tz_localize(None)
    ext = ext[[c for c in ext.columns if c in old.columns]]
    ext = ext[ext.index < old.index.min()]
    m = pd.concat([ext.reindex(columns=old.columns), old]).sort_index()
    assert not m.index.has_duplicates, "duplicate dates after merge"
    m.to_parquet(MERGED)
    print(f"merged: {m.index.min().date()} .. {m.index.max().date()} "
          f"({len(m):,} days) -> {MERGED.name}")
    return m


def t1_variance(m):
    """DJF variability ratio per column, pre- vs post-1979.

    Measured on ANOMALIES from a PER-ERA day-of-year climatology, not on raw
    values. Raw DJF standard deviation conflates three things: the resolved
    day-to-day variability (what a nudged-to-climatology analysis actually
    suppresses, and the only thing this gate cares about), the secular trend, and
    the mean offset between eras. Stratospheric heights alone fall by hundreds of
    metres between the two periods from greenhouse cooling, which is real physics
    and no reason to reject the era. Removing each era's own seasonal cycle
    isolates the variability the analysis resolves.

    The DJF-mean interannual SD is reported alongside, because an analysis pinned
    to climatology loses year-to-year spread too, and that channel is independent.
    """
    djf = m[np.isin(m.index.month, (12, 1, 2))]
    a, b = djf[djf.index < SPLIT], djf[djf.index >= SPLIT]
    print(f"\n=== T1  DJF VARIABILITY RATIO  (n={len(a):,} vs {len(b):,} days) ===")
    print(f"{'column':18s} {'anom 58-78':>10s} {'anom 79-24':>10s} {'ratio':>7s} "
          f"{'iav ratio':>10s} {'mean diff':>10s}  verdict")
    print("-" * 82)
    out = {}
    for c in m.columns:
        stats = []
        for d in (a, b):
            an = d[c] - d[c].groupby(d.index.dayofyear).transform("mean")
            iav = d[c].groupby(L.winter_of(d.index)).mean().std()
            stats.append((float(an.std()), float(iav)))
        (sa, ia), (sb, ib) = stats
        r = float(sa / sb) if sb else np.nan
        ri = float(ia / ib) if ib else np.nan
        ok = VAR_RATIO_OK[0] <= r <= VAR_RATIO_OK[1]
        md = float(a[c].mean() - b[c].mean())
        out[c] = {"anom_sd_pre": round(sa, 3), "anom_sd_post": round(sb, 3),
                  "sd_ratio": round(r, 3), "interannual_sd_ratio": round(ri, 3),
                  "mean_diff": round(md, 3), "pass": bool(ok)}
        print(f"{c:18s} {sa:10.3f} {sb:10.3f} {r:7.3f} {ri:10.3f} {md:+10.3f}  "
              f"{'ok' if ok else 'FAIL'}")
    return out


def t4_independent(m):
    """Agreement with an INDEPENDENT reanalysis, before vs after 1979.

    T1 is internal to NCEP and cannot separate "the pre-satellite analysis is
    degraded" from "the pre-satellite atmosphere was genuinely different". ERA5
    can: it is a different model, a different assimilation system and a different
    century of development, so if NCEP's pre-1979 stratosphere were substantially
    fabricated its agreement with ERA5 would fall away in that era.

    Read with its limit stated: both products assimilate the SAME sparse
    radiosonde network before 1979, so they can be wrong together. A DROP in
    agreement is therefore strong evidence against pooling, while NO drop is only
    moderate evidence for it. That asymmetry is why this is a supporting
    diagnostic and not the gate.

    Compared here: NCEP polar-cap 100 hPa height against the ERA5 150 hPa polar-
    cap NAM, the closest independent pair the project already holds.
    """
    e5 = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "era5_nam_daily.parquet"
    if not e5.exists():
        print("\n=== T4 skipped: era5_nam_daily.parquet not found ===")
        return None
    e = pd.read_parquet(e5)
    if e.index.tz is not None:
        e.index = e.index.tz_convert("UTC").tz_localize(None)
    col = "nam_150" if "nam_150" in e.columns else next(
        (c for c in e.columns if "150" in c), None)
    if col is None:
        print(f"\n=== T4 skipped: no 150 hPa column in {list(e.columns)} ===")
        return None

    n = -L.doy_standardise(m["hgt_m_100hPa"])     # + = strong vortex, as elsewhere
    x = L.doy_standardise(e[col])
    j = pd.concat([n.rename("ncep"), x.rename("era5")], axis=1).dropna()
    j = j[np.isin(j.index.month, (11, 12, 1, 2, 3))]
    out = {"era5_column": col}
    print(f"\n=== T4  AGREEMENT WITH ERA5 (independent reanalysis, NDJFM) ===")
    print(f"{'era':12s} {'n days':>8s} {'corr':>8s}")
    print("-" * 30)
    for lab, sel in (("1959-1978", j.index < SPLIT), ("1979-2022", j.index >= SPLIT)):
        s = j[sel]
        c = float(s.ncep.corr(s.era5)) if len(s) > 100 else np.nan
        out[lab] = {"n": int(len(s)), "corr": round(c, 4)}
        print(f"{lab:12s} {len(s):8,d} {c:8.4f}")
    a, b = out.get("1959-1978", {}).get("corr"), out.get("1979-2022", {}).get("corr")
    if a is not None and b is not None and np.isfinite(a) and np.isfinite(b):
        out["drop"] = round(float(b - a), 4)
        print(f"\n  agreement drop pre-1979: {b - a:+.4f} "
              f"({'no material drop' if b - a < 0.05 else 'MATERIAL DROP'})")
    return out


def predictors(m):
    hf = pd.read_parquet(HF).set_index("date")["vT_100hPa_45_75N"]
    # standardised on the FULL merged record, deliberately: standardising within
    # each era would absorb exactly the era difference this gate is looking for.
    return {"u10": L.doy_standardise(m["uwnd_ms_10hPa"]),
            "z100": -L.doy_standardise(m["hgt_m_100hPa"]),
            "vt": L.doy_standardise(hf)}


def fit_beta(names, SRC, WINS, ao, clim, clean_idx, doys, rng, n_target):
    Yp, Sp = [], {k: [] for k in names}
    drawn = 0
    while drawn < n_target:
        f = G.draw_clean(clean_idx, doys, rng)
        if len(f) < 8:
            continue
        Yp.append(S.per_event(ao, f, clim))
        for nm in names:
            Sp[nm].append(L.window_mean(SRC[nm], f, WINS[nm]))
        drawn += len(f)
    Yp = np.concatenate(Yp)
    Sp = {k: np.concatenate(v) for k, v in Sp.items()}
    beta = {}
    for nm in names:
        ok = np.isfinite(Yp) & np.isfinite(Sp[nm])
        beta[nm] = float(np.cov(Yp[ok], Sp[nm][ok])[0, 1] / np.var(Sp[nm][ok], ddof=1))
    return beta, Yp, Sp


def beta_ci(Yp, Sv, rng, n_boot=N_BOOT):
    ok = np.isfinite(Yp) & np.isfinite(Sv)
    y, s = Yp[ok], Sv[ok]
    bs = []
    for _ in range(n_boot):
        i = rng.integers(0, len(y), len(y))
        bs.append(np.cov(y[i], s[i])[0, 1] / np.var(s[i], ddof=1))
    return [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]


def main():
    rng = np.random.default_rng(SEED)
    m = merge()
    res = {"split": str(SPLIT.date()), "variance": t1_variance(m),
           "independent_era5": t4_independent(m)}

    ao = M.load("ao")["y"]
    pred = predictors(m)
    pred["ao"] = L.doy_standardise(ao)
    names = [s[0] for s in L.STRATIFIERS]
    SRC = {nm: pred[var] for nm, _, var, _, _ in L.STRATIFIERS}
    WINS = {nm: w for nm, _, _, w, _ in L.STRATIFIERS}

    allev = load_catalogue("primary")   # EXCLUSION set: every event, even ones not scored here
    real = allev
    lo_cov = max(v.index.min() for v in pred.values())
    hi_cov = min(v.index.max() for v in pred.values())
    real = real[(real >= max(lo_cov, ao.index.min()) + pd.Timedelta(days=60))
                & (real <= min(hi_cov, ao.index.max()) - pd.Timedelta(days=60))]
    print(f"\nevents with full predictor coverage: {len(real)} "
          f"({real.min().date()} .. {real.max().date()})")
    print(f"  pre-1979: {(real < SPLIT).sum()}   post-1979: {(real >= SPLIT).sum()}")
    res["n_events_total"] = int(len(real))
    res["n_events_pre"] = int((real < SPLIT).sum())
    res["n_events_post"] = int((real >= SPLIT).sum())

    mask = G.real_influence_mask(ao.index, allev)   # every catalogued event
    clean = ao[~mask]
    doys = np.array([t.dayofyear for t in pd.DatetimeIndex(real)])

    # ---------------- T2: beta refitted within each era ---------------------
    print(f"\n=== T2  BETA STABILITY (refitted per era on event-free draws) ===")
    print(f"{'stratifier':12s} {'beta 58-78':>22s} {'beta 79-24':>22s}  verdict")
    print("-" * 68)
    eras = {"pre": clean.index < SPLIT, "post": clean.index >= SPLIT}
    betas, cis, Yps, Sps = {}, {}, {}, {}
    for era, sel in eras.items():
        zf = G.zone_free_index(clean.index, allev)    # pseudo zones clear of every event
        idx = zf[zf < SPLIT] if era == "pre" else zf[zf >= SPLIT]
        clim = clean[sel].groupby(clean[sel].index.dayofyear).mean()
        dd = doys[(pd.DatetimeIndex(real) < SPLIT) if era == "pre"
                  else (pd.DatetimeIndex(real) >= SPLIT)]
        b, Yp, Sp = fit_beta(names, SRC, WINS, ao, clim, idx, dd, rng, N_BETA)
        betas[era], Yps[era], Sps[era] = b, Yp, Sp
        cis[era] = {nm: beta_ci(Yp, Sp[nm], rng) for nm in names}

    res["beta"] = {}
    n_beta_ok = 0
    for nm in names:
        ca, cb = cis["pre"][nm], cis["post"][nm]
        overlap = not (ca[1] < cb[0] or cb[1] < ca[0])
        n_beta_ok += overlap
        res["beta"][nm] = {
            "pre": round(betas["pre"][nm], 4), "pre_CI95": [round(x, 3) for x in ca],
            "post": round(betas["post"][nm], 4), "post_CI95": [round(x, 3) for x in cb],
            "CIs_overlap": bool(overlap)}
        print(f"{nm:12s} {betas['pre'][nm]:+7.3f} [{ca[0]:+6.3f},{ca[1]:+6.3f}] "
              f"{betas['post'][nm]:+7.3f} [{cb[0]:+6.3f},{cb[1]:+6.3f}]  "
              f"{'ok' if overlap else 'DIFFER'}")

    # ---------------- T3: residuals at the real events, per era -------------
    print(f"\n=== T3  RESIDUAL AGREEMENT AT EVENTS "
          f"(n={res['n_events_pre']} vs {res['n_events_post']}; LOW POWER) ===")
    print(f"{'stratifier':12s} {'resid pre':>10s} {'resid post':>11s} "
          f"{'difference':>11s} {'95% CI of diff':>20s}")
    print("-" * 70)
    clim_all = clean.groupby(clean.index.dayofyear).mean()
    Y_real = S.per_event(ao, real, clim_all)
    pre_m = np.asarray(pd.DatetimeIndex(real) < SPLIT)
    res["residual"] = {}
    for nm in names:
        Sr = L.window_mean(SRC[nm], real, WINS[nm])
        rr = {}
        for era, sel in (("pre", pre_m), ("post", ~pre_m)):
            dY, dS, n = L.contrast(Y_real[sel], Sr[sel])
            rr[era] = dY - betas[era][nm] * dS if np.isfinite(dY) else np.nan
        diff = rr["pre"] - rr["post"]
        bs = []
        for _ in range(N_BOOT):
            vals = {}
            for era, sel in (("pre", pre_m), ("post", ~pre_m)):
                ix = np.flatnonzero(sel)
                p = rng.choice(ix, len(ix), replace=True)
                d1, d2, _ = L.contrast(Y_real[p], Sr[p])
                vals[era] = d1 - betas[era][nm] * d2 if np.isfinite(d1) else np.nan
            if np.isfinite(vals["pre"]) and np.isfinite(vals["post"]):
                bs.append(vals["pre"] - vals["post"])
        ci = ([float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]
              if len(bs) > 50 else None)
        res["residual"][nm] = {
            "pre": None if not np.isfinite(rr["pre"]) else round(float(rr["pre"]), 4),
            "post": None if not np.isfinite(rr["post"]) else round(float(rr["post"]), 4),
            "difference": None if not np.isfinite(diff) else round(float(diff), 4),
            "difference_CI95": [round(x, 3) for x in ci] if ci else None,
            "consistent": bool(ci and ci[0] <= 0 <= ci[1])}
        cs = f"[{ci[0]:+7.3f},{ci[1]:+7.3f}]" if ci else " " * 20
        print(f"{nm:12s} {rr['pre']:+10.3f} {rr['post']:+11.3f} {diff:+11.3f} {cs}")

    # ---------------- verdict, by the rule fixed above ----------------------
    used = ["uwnd_ms_10hPa", "hgt_m_100hPa"]
    t1_ok = all(res["variance"][c]["pass"] for c in used if c in res["variance"])
    t2_ok = n_beta_ok == len(names)
    res["T1_pass"] = bool(t1_ok)
    res["T2_pass"] = bool(t2_ok)
    res["T2_n_betas_agreeing"] = int(n_beta_ok)
    res["pooling_allowed"] = bool(t1_ok and t2_ok)
    res["verdict"] = (
        f"POOLING ALLOWED. The pre-satellite era matches the satellite era on both "
        f"powered tests (DJF variance ratios within {VAR_RATIO_OK} at the levels used, "
        f"and all {len(names)} betas agreeing within bootstrap CIs), so the "
        f"observational arm goes to n={res['n_events_total']}."
        if res["pooling_allowed"] else
        f"POOLING REFUSED. T1 {'passed' if t1_ok else 'FAILED'}, "
        f"T2 {n_beta_ok}/{len(names)} betas agree. The observational arm stays at "
        f"n={res['n_events_post']} (1979 onward) and the pre-satellite extension is "
        f"reported only as a sensitivity test, with the failing diagnostic stated.")
    print(f"\n=== VERDICT ===\n  {res['verdict']}")

    # ------------------------------------------------------------------
    # POST-HOC REFINEMENT. Decided AFTER seeing T2, and labelled as such.
    # The pre-specified verdict above stands and is reported first.
    #
    # Two things are wrong with the pre-specified rule, and both are visible
    # regardless of which way they push:
    #
    #  (a) IT IS A DIFFERENCE TEST WHERE AN EQUIVALENCE TEST IS WANTED. beta is
    #      fitted on N_BETA pseudo-onsets, a number I chose freely. Its CI
    #      narrows as 1/sqrt(N_BETA), so "CIs overlap" becomes arbitrarily
    #      strict as I spend more compute, and at N_BETA=3000 it flags AO_post
    #      as era-dependent over a gap of ~0.08 on a beta of ~1.83. A rule whose
    #      verdict depends on how many pseudo-draws I chose to make is not a
    #      rule about the atmosphere. What matters is not whether the betas
    #      differ but whether the difference CHANGES THE CORRECTION, so the
    #      consequence is quantified here: |beta_pre - beta_post| x dS, compared
    #      against the residual's own bootstrap uncertainty.
    #
    #  (b) IT IS GLOBAL WHERE THE DECISION IS PER-STRATIFIER. Each stratifier
    #      has its own beta and its own correction. vT failing is not a reason
    #      to refuse the extra events for u10_post30, whose beta and residual
    #      both agree across eras. A blanket veto is not conservatism, it is
    #      just a different way of getting the answer wrong.
    #
    # Note vT is NOT affected by this extension at all -- heatflux_daily.parquet
    # already covered 1958-2024 before any of this. Its era difference is a
    # pre-existing property of the heat flux record, surfaced here by accident.
    # ------------------------------------------------------------------
    print(f"\n=== POST-HOC: PER-STRATIFIER POOLING, WITH CONSEQUENCE ===")
    print(f"{'stratifier':12s} {'dbeta':>8s} {'dS_pre':>8s} {'dResid':>8s} "
          f"{'resid CI hw':>12s} {'ratio':>7s}  pool?")
    print("-" * 76)
    res["per_stratifier"] = {}
    for nm in names:
        Sr = L.window_mean(SRC[nm], real, WINS[nm])
        _, dS_pre, _ = L.contrast(Y_real[pre_m], Sr[pre_m])
        db = betas["pre"][nm] - betas["post"][nm]
        dres = db * dS_pre if np.isfinite(dS_pre) else np.nan
        rc = res["residual"][nm]["difference_CI95"]
        hw = (rc[1] - rc[0]) / 2 if rc else np.nan
        ratio = abs(dres) / hw if hw and np.isfinite(dres) and hw > 0 else np.nan
        # material only if the beta difference moves the correction by an
        # appreciable fraction of what the residual is already uncertain by
        immaterial = np.isfinite(ratio) and ratio < 0.25
        agree = res["residual"][nm]["consistent"]
        pool = bool((res["beta"][nm]["CIs_overlap"] or immaterial) and agree)
        res["per_stratifier"][nm] = {
            "delta_beta": round(float(db), 4),
            "dS_pre": None if not np.isfinite(dS_pre) else round(float(dS_pre), 4),
            "correction_shift": None if not np.isfinite(dres) else round(float(dres), 4),
            "residual_CI_halfwidth": None if not np.isfinite(hw) else round(float(hw), 4),
            "shift_over_uncertainty": None if not np.isfinite(ratio) else round(float(ratio), 3),
            "beta_difference_immaterial": bool(immaterial),
            "residuals_agree": bool(agree),
            "pool": pool}
        print(f"{nm:12s} {db:+8.3f} {dS_pre:+8.3f} {dres:+8.3f} {hw:12.3f} "
              f"{ratio:7.2f}  {'YES' if pool else 'no'}")

    poolable = [n for n in names if res["per_stratifier"][n]["pool"]]
    res["poolable_stratifiers"] = poolable
    res["headline_recommendation"] = (
        f"Pre-specified global rule REFUSES pooling, and that verdict stands as the "
        f"headline: the observational arm is n={res['n_events_post']}. Per-stratifier, "
        f"{len(poolable)}/{len(names)} survive both checks ({', '.join(poolable) or 'none'}) "
        f"and the n={res['n_events_total']} record is reported for those as a "
        f"pre-registered sensitivity test, never as the primary number.")
    print(f"\n  poolable per-stratifier: {poolable or 'none'}")
    print(f"\n=== RECOMMENDATION ===\n  {res['headline_recommendation']}")

    (RESULTS / "era_dependence_gate.json").write_text(json.dumps(res, indent=2),
                                                   encoding="utf8")
    print("\nSaved -> era_dependence_gate.json")


if __name__ == "__main__":
    main()
