#!/usr/bin/env python3
"""
era5_recompute_and_two_thirds.py
================================
Two things the audit still needed, in the literature's own variables.

(A) THE BODY COUNT IN ERA5, NOT IN CPC PROXIES
    `recompute_published_criterion.py` applied Karpechko et al. (2017) using the
    CPC AO as a stand-in for the 1000 hPa NAM. That establishes the decomposition
    but cannot speak to a published number. ERA5 polar-cap NAM at 1000, 850 and
    150 hPa is now available (`03_data_ingestion/era5_nam_daily.parquet`,
    validated at corr = +0.839 against the CPC AO over 11,609 shared days), so the
    criterion runs on the actual fields:
      Karpechko original : conditions on NAM_1000 and NAM_150
      ACP 26, 3723 (2026): the same with 850 hPa substituted for 1000 hPa,
                           reporting NAO -0.620 (all DW; -0.762 is the BOTH
                           subtype only) vs +0.088 (NDW).

(B) DOES THE SHIFT PREDICT "TWO THIRDS"?
    Baldwin et al. (2021), Reviews of Geophysics 59, e2020RG000708, section 7.2 --
    the field's consensus statement -- says verbatim:

      "Most studies agree that about two thirds (Charlton-Perez et al., 2018;
       Domeisen, 2019; White et al., 2019) of SSW events are characterized as
       having a visible downward impact"

    quoted with no uncertainty range, no null-model comparison, and no statement
    of what fraction of random winter dates would pass the same test.

    `is_downward_propagation_a_class.py` found the per-event surface response
    distribution to be a PURE SHIFT of the ordinary winter distribution (-0.892
    sigma in observations, -0.511 in CMIP6; unimodal, one Gaussian component,
    KS p = 0.19-0.33). If that is right it makes a hard, falsifiable prediction:

      applying the criterion to pseudo-events whose surface NAM is displaced by
      the measured shift must reproduce the observed classification rate.

    Get two thirds out, and the whole phenomenology -- "about two thirds of SSWs
    propagate downward" -- follows from one shifted population with no classes in
    it. Get something else out, and the pure-shift model is wrong.

    Tested on the SURFACE conditions (1 and 2), which are what the shift model
    speaks to. Condition 3 is stratospheric and an SSW plainly does displace the
    stratosphere; it is reported separately rather than folded in.

Output: era5_recompute_and_two_thirds.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "9_literature"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]
SIM = ROOT / "ssw-design-analysis" / "06_simulation_validation"
sys.path.insert(0, str(SIM))
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
from build_catalogue import load_catalogue          # noqa: E402
import gate3_clean_null as G                        # noqa: E402
import selection_on_outcome as SO                   # noqa: E402

ERA5 = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "era5_nam_daily.parquet"
SEL_WIN = (8, 52)
OUT_WIN = (1, 60)
N_NULL = 2000
SEED = 20260801
# Lu & Rao (2026, ACP 26, 3723) give ERA5 60-day NAO means by DW subtype: BOTH -0.762 (n=13), EA -0.567 (14), NA -0.435 (6); NDW +0.088 (19). All-DW mean -0.620, contrast -0.708, DW fraction 33/52 = 0.635. The -0.850 used until 2026-09-25 subtracted the BOTH subtype alone.
PUBLISHED = {"ACP_2026_NAO_DW": -0.620, "ACP_2026_NAO_NDW": 0.088,
             "ACP_2026_contrast": -0.708,
             "Baldwin_2021_fraction": "about two thirds"}


def win_days(series, o, win):
    lo, hi = o + pd.Timedelta(days=win[0]), o + pd.Timedelta(days=win[1])
    return series[(series.index >= lo) & (series.index <= hi)]


def classify(onsets, surf, strat, shift=0.0):
    """Karpechko conditions 1-3. `shift` displaces the SURFACE field only."""
    lab, c1s, c2s, c3s = [], [], [], []
    for o in onsets:
        a = win_days(surf, o, SEL_WIN) + shift
        if len(a) < 20:
            lab.append(False); c1s.append(False); c2s.append(False); c3s.append(False)
            continue
        c1 = bool(a.mean() < 0)
        c2 = bool((a < 0).mean() > 0.5)
        n = win_days(strat, o, SEL_WIN)
        c3 = bool((n < 0).mean() > 0.7) if len(n) >= 20 else True
        lab.append(c1 and c2 and c3); c1s.append(c1); c2s.append(c2); c3s.append(c3)
    return (np.array(lab), np.array(c1s), np.array(c2s), np.array(c3s))


def main():
    rng = np.random.default_rng(SEED)
    e = pd.read_parquet(ERA5)
    if e.index.tz is not None:
        e.index = e.index.tz_convert("UTC").tz_localize(None)
    print(f"ERA5 NAM {e.index.min().date()} .. {e.index.max().date()} ({len(e):,} d)")

    real = load_catalogue("primary")
    real = real[(real >= e.index.min() + pd.Timedelta(days=70))
                & (real <= e.index.max() - pd.Timedelta(days=70))]
    print(f"events inside ERA5 coverage: {len(real)}")

    mask = G.real_influence_mask(e.index, real)
    clean_idx = e[~mask].index
    doys = np.array([t.dayofyear for t in pd.DatetimeIndex(real)])
    clim = {c: e[c][~mask].groupby(e[c][~mask].index.dayofyear).mean()
            for c in ("nam_1000", "nam_850")}

    res = {"published": PUBLISHED, "n_events": int(len(real)),
           "selection_window": list(SEL_WIN), "outcome_window": list(OUT_WIN),
           "variants": {}}

    for name, selcol in (("Karpechko_1000hPa", "nam_1000"),
                         ("ACP2026_850hPa", "nam_850")):
        surf, strat = e[selcol], e["nam_150"]
        lab, c1, c2, c3 = classify(real, surf, strat)
        Y = SO.per_event(e["nam_1000"], real, clim["nam_1000"])
        S = np.array([float(win_days(surf, o, SEL_WIN).mean()) for o in real])
        ok = np.isfinite(Y) & np.isfinite(S)

        dw, ndw = float(Y[lab & ok].mean()), float(Y[(~lab) & ok].mean())
        dY = dw - ndw
        dS = float(S[lab & ok].mean() - S[(~lab) & ok].mean())

        # beta and the null, from pseudo-onsets on event-cleaned days
        nl, nY, nS, rates = [], [], [], []
        for _ in range(N_NULL):
            f = G.draw_clean(clean_idx, doys, rng)
            if len(f) < 10:
                continue
            fl, *_ = classify(f, surf, strat)
            Yf = SO.per_event(e["nam_1000"], f, clim["nam_1000"])
            Sf = np.array([float(win_days(surf, o, SEL_WIN).mean()) for o in f])
            k = np.isfinite(Yf) & np.isfinite(Sf)
            rates.append(float(fl.mean()))
            nY.append(Yf[k]); nS.append(Sf[k])
            if (fl & k).sum() >= 3 and ((~fl) & k).sum() >= 3:
                nl.append((float(Yf[fl & k].mean() - Yf[(~fl) & k].mean()),
                           float(Sf[fl & k].mean() - Sf[(~fl) & k].mean())))
        nY_all, nS_all = np.concatenate(nY), np.concatenate(nS)
        beta = float(np.cov(nY_all, nS_all)[0, 1] / np.var(nS_all, ddof=1))
        nl = np.array(nl)
        nres = nl[:, 0] - beta * nl[:, 1]
        resid = dY - beta * dS
        p = float((np.abs(nres) >= abs(resid)).mean())

        res["variants"][name] = {
            "selection_field": selcol,
            "n_DW": int(lab.sum()), "n_NDW": int((~lab).sum()),
            "rate_DW": round(float(lab.mean()), 3),
            "rate_DW_null": round(float(np.mean(rates)), 3),
            "DW_composite": round(dw, 4), "NDW_composite": round(ndw, 4),
            "contrast": round(dY, 4), "dS": round(dS, 4),
            "beta": round(beta, 4),
            "selection_term": round(float(beta * dS), 4),
            "residual": round(float(resid), 4), "p_residual": p,
            "pct_from_selection": round(float(100 * beta * dS / dY), 1)
            if abs(dY) > 1e-9 else None,
            "null_contrast": round(float(nl[:, 0].mean()), 4),
            "pct_reproduced_under_null": round(
                float(100 * nl[:, 0].mean() / dY), 1) if abs(dY) > 1e-9 else None,
            "cond_pass_rates": {"c1_mean_negative": round(float(c1.mean()), 3),
                                "c2_frac_negative": round(float(c2.mean()), 3),
                                "c3_stratospheric": round(float(c3.mean()), 3)}}
        r = res["variants"][name]
        print(f"\n=== {name}  (selection on {selcol}) ===")
        print(f"  DW {r['n_DW']}/{len(real)} = {100*r['rate_DW']:.0f}%   "
              f"(null rate {100*r['rate_DW_null']:.0f}%)")
        print(f"  DW {dw:+.3f}   NDW {ndw:+.3f}   contrast {dY:+.3f}")
        print(f"  beta {beta:+.3f} x dS {dS:+.3f} = {beta*dS:+.3f}  "
              f"-> {r['pct_from_selection']:.0f}% of the contrast")
        print(f"  residual {resid:+.3f} (p={p:.3f})")
        print(f"  same criterion on pseudo-events reproduces "
              f"{r['pct_reproduced_under_null']:.0f}% of it")

    # ---------------- (B) does the shift predict two thirds? ---------------
    print("\n" + "=" * 72)
    print("=== DOES A PURE SHIFT PREDICT 'ABOUT TWO THIRDS'? ===")
    surf, strat = e["nam_1000"], e["nam_150"]
    lab, c1, c2, c3 = classify(real, surf, strat)
    Yr = SO.per_event(e["nam_1000"], real, clim["nam_1000"])
    pseudo_sets = []
    for _ in range(400):
        f = G.draw_clean(clean_idx, doys, rng)
        if len(f) >= 10:
            pseudo_sets.append(f)
    Yp = np.concatenate([SO.per_event(e["nam_1000"], f, clim["nam_1000"])
                         for f in pseudo_sets])
    delta = float(np.nanmean(Yr) - np.nanmean(Yp))
    print(f"  measured surface shift, real minus pseudo: {delta:+.3f} sigma")

    obs_c12 = float((c1 & c2).mean())
    null_c12, pred_c12 = [], []
    for f in pseudo_sets[:200]:
        _, a1, a2, _ = classify(f, surf, strat)
        null_c12.append(float((a1 & a2).mean()))
        _, b1, b2, _ = classify(f, surf, strat, shift=delta)
        pred_c12.append(float((b1 & b2).mean()))
    nm, pm = float(np.mean(null_c12)), float(np.mean(pred_c12))
    lo, hi = np.percentile(pred_c12, [2.5, 97.5])
    print(f"  surface conditions (1 AND 2) pass rate:")
    print(f"    observed at real SSWs        {100*obs_c12:5.1f}%")
    print(f"    pseudo-events, no shift      {100*nm:5.1f}%")
    print(f"    pseudo-events SHIFTED by delta {100*pm:5.1f}%  "
          f"[{100*lo:.1f}, {100*hi:.1f}]")
    agree = bool(lo <= obs_c12 <= hi)
    print(f"    -> shift model {'REPRODUCES' if agree else 'FAILS TO REPRODUCE'} "
          f"the observed rate")
    res["two_thirds_test"] = {
        "measured_shift_sigma": round(delta, 4),
        "observed_pass_rate_c1c2": round(obs_c12, 4),
        "null_pass_rate_c1c2": round(nm, 4),
        "shifted_pseudo_pass_rate_c1c2": round(pm, 4),
        "shifted_pseudo_CI95": [round(float(lo), 4), round(float(hi), 4)],
        "shift_model_reproduces_observed_rate": agree,
        "full_criterion_rate_observed": round(float(lab.mean()), 4)}

    res["verdict"] = (
        f"In ERA5 the published criterion gives a DW-minus-NDW contrast of "
        f"{res['variants']['Karpechko_1000hPa']['contrast']:+.3f} of which "
        f"{res['variants']['Karpechko_1000hPa']['pct_from_selection']:.0f}% is the "
        f"selection term. A pure surface shift of {delta:+.3f} sigma applied to "
        f"event-free dates gives a {100*pm:.0f}% pass rate against {100*obs_c12:.0f}% "
        f"observed -- so 'about two thirds propagate downward' is what ONE shifted "
        f"population produces, with no classes in it."
        if agree else
        f"The pure-shift model does NOT reproduce the observed classification rate "
        f"({100*pm:.0f}% predicted vs {100*obs_c12:.0f}% observed); the one-population "
        f"reading is incomplete.")
    print(f"\n=== VERDICT ===\n  {res['verdict']}")

    (RESULTS / "era5_recompute_and_two_thirds.json").write_text(
        json.dumps(res, indent=2), encoding="utf8")
    print("\nSaved -> era5_recompute_and_two_thirds.json")


if __name__ == "__main__":
    main()
