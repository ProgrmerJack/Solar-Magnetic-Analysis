#!/usr/bin/env python3
"""
selection_on_outcome.py
=======================
THE DESIGN CHOICE THAT ACTUALLY MATTERS.

Two design axes have now been tested on the frozen catalogue and returned null:
estimator choice moves the AO estimate 0.002 sigma, catalogue choice 0.115 sigma
on a common analysis window. This tests a third axis, which is standard practice:
SELECTING EVENTS ON THE OUTCOME.

THE PRACTICE, AS PUBLISHED
  Karpechko et al. (2017) classify an SSW as "downward-propagating" (dSSW) if:
    1. mean NAM at 1000/850 hPa over days +8..+52 is NEGATIVE
    2. fraction of those days with negative surface NAM > 0.5
    3. fraction with negative NAM at 150/100 hPa > 0.7
  Still in active use (ACP 24, 1389, 2024; ACP 26, 3723, 2026).

  Criteria 1 and 2 ARE the surface response. This is NOT a criticism of the
  classification, which is sound for what it was built for -- precursors,
  predictability, and what distinguishes coupling from non-coupling events. A
  case-control design may legitimately select on outcome and look back at
  predictors. The problem is narrower: when the SELECTED GROUP'S MEAN SURFACE
  ANOMALY is reported as an impact, that number is partly manufactured by the
  selection, and the size has not been quantified.

TWO CORRECTIONS TO THE FIRST VERSION OF THIS SCRIPT, both of which changed the
conclusion and are documented here so the error is on the record:

  (a) CONTAMINATED REFERENCE. The day-of-year climatology was computed on the
      FULL record, which includes real-event days. Those days are depressed, so
      clean pseudo-onset days scored +0.23 sigma against it and the null did not
      centre on zero (observed +0.279, predicted +0.2296). The climatology is now
      built from CLEAN days only. This is the same contamination trap already
      recorded three times in this project: the damage is in the reference, not
      the treatment.

  (b) UNMATCHED SELECTION RATE. The criterion selects 70% of real events but only
      38% of pseudo-events, and selecting a SMALLER fraction is a STRONGER
      selection. Comparing the two directly overstated the bias. Both the
      identical-criterion comparison and a rate-matched comparison are now
      reported.

WHAT THIS MEASURES
  Pseudo-onsets are drawn on clean data (every day within -60..+75 of a real
  onset deleted first), so the true effect is ZERO by construction. Applying the
  published classification to them and reporting the selected group's mean gives
  the number a study would report if SSWs had no surface effect at all.

  Reported for three outcomes, in decreasing dependence on the selection variable:
    AO   -- the selection variable itself; circularity total, the benchmark
    NAO  -- a different index, strongly correlated with AO
    PNA  -- a different index, weakly correlated with AO

Output: selection_on_outcome.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "4_selection_bias"
RESULTS.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
from build_catalogue import load_catalogue          # noqa: E402
import gate3_clean_null as G                        # noqa: E402
import multi_index_event_study as M                 # noqa: E402

ROOT = HERE.parents[1]
STRAT = ROOT / "data" / "processed" / "atmospheric" / "ncep_stratosphere.parquet"
WIN = (8, 52)
N_DRAW = 1500
OUTCOMES = ("ao", "nao", "pna")


def strat_nam():
    s = pd.read_parquet(STRAT)
    if s.index.tz is not None:
        s.index = s.index.tz_convert("UTC").tz_localize(None)
    z = s["hgt_m_100hPa"]
    doy = z.index.dayofyear
    return -((z - z.groupby(doy).transform("mean")) / z.groupby(doy).transform("std"))


def classify(onsets, ao, snam):
    lab = []
    for o in onsets:
        lo, hi = o + pd.Timedelta(days=WIN[0]), o + pd.Timedelta(days=WIN[1])
        a = ao[(ao.index >= lo) & (ao.index <= hi)]
        if len(a) < 20:
            lab.append(False)
            continue
        c1 = a.mean() < 0
        c2 = (a < 0).mean() > 0.5
        n = snam[(snam.index >= lo) & (snam.index <= hi)]
        c3 = ((n < 0).mean() > 0.7) if len(n) >= 20 else None
        lab.append(bool(c1 and c2 and (c3 is not False)))
    return np.array(lab)


def per_event(series, onsets, clim):
    """Day-of-year-adjusted window anomaly for each onset (NaN if unusable)."""
    out = []
    for o in onsets:
        lo, hi = o + pd.Timedelta(days=WIN[0]), o + pd.Timedelta(days=WIN[1])
        s = series[(series.index >= lo) & (series.index <= hi)]
        out.append(float(s.mean() - clim.reindex(s.index.dayofyear).mean())
                   if len(s) >= 20 else np.nan)
    return np.array(out)


def main():
    idx = {k: M.load(k)["y"] for k in OUTCOMES}
    ao = idx["ao"]
    snam = strat_nam()
    real = load_catalogue("primary")
    real = real[(real >= ao.index.min()) & (real <= ao.index.max())]

    # CLEAN climatology -- correction (a). Real-event days must not define the
    # reference the null is scored against.
    mask = G.real_influence_mask(ao.index, real)
    clims = {}
    for k, v in idx.items():
        cv = v[~G.real_influence_mask(v.index, real)]
        clims[k] = cv.groupby(cv.index.dayofyear).mean()

    lab = classify(real, ao, snam)
    n_real_d = int(lab.sum())
    rate_real = float(lab.mean())
    print(f"REAL events: {len(real)}  ->  dSSW {n_real_d}, nSSW {len(real) - n_real_d} "
          f"({100 * rate_real:.0f}% classified downward)")

    real_all, real_d = {}, {}
    print(f"\n{'outcome':6s} {'ALL events':>12s} {'dSSW only':>12s}")
    print("-" * 32)
    for k in OUTCOMES:
        v = per_event(idx[k], real, clims[k])
        real_all[k] = float(np.nanmean(v))
        real_d[k] = float(np.nanmean(v[lab]))
        print(f"{k:6s} {real_all[k]:+12.3f} {real_d[k]:+12.3f}")

    # ---------------- null: pseudo-onsets, true effect exactly zero ----------
    clean_idx = G.zone_free_index(ao.index, real)
    doys = np.array([t.dayofyear for t in pd.DatetimeIndex(real)])
    rng = np.random.default_rng(20260731)
    null = {k: {"all": [], "crit": [], "matched": []} for k in OUTCOMES}
    rates = []
    for i in range(N_DRAW):
        fake = G.draw_clean(clean_idx, doys, rng)
        if len(fake) < 10:
            continue
        fl = classify(fake, ao, snam)
        if fl.sum() < 5:
            continue
        rates.append(float(fl.mean()))
        sel_ao = per_event(ao, fake, clims["ao"])
        # rate-matched selection -- correction (b): take the same FRACTION of
        # pseudo-events as the criterion took of real ones, most negative first
        k_match = max(1, int(round(rate_real * len(fake))))
        order = np.argsort(np.where(np.isnan(sel_ao), np.inf, sel_ao))
        matched = np.zeros(len(fake), bool)
        matched[order[:k_match]] = True
        for k in OUTCOMES:
            v = per_event(idx[k], fake, clims[k])
            if np.all(np.isnan(v)):
                continue
            null[k]["all"].append(np.nanmean(v))
            if fl.sum():
                null[k]["crit"].append(np.nanmean(v[fl]))
            null[k]["matched"].append(np.nanmean(v[matched]))
        if (i + 1) % 300 == 0:
            print(f"  {i + 1}/{N_DRAW} null draws", flush=True)

    nd = len(rates)
    rate_null = float(np.mean(rates))
    # is the real classification RATE distinguishable from the null rate?
    p_rate = float((np.array(rates) >= rate_real).mean())
    print(f"\n=== NULL ({nd} draws, true effect = 0) ===")
    print(f"  null centres on zero? all-pseudo-event mean AO = "
          f"{np.mean(null['ao']['all']):+.3f}  (was +0.279 before correction (a))")
    print(f"  classification RATE: real {100 * rate_real:.0f}% vs null "
          f"{100 * rate_null:.0f}%   P(null >= real) = {p_rate:.4f}")

    res = {"window": list(WIN), "n_draw": nd, "n_real": int(len(real)),
           "n_real_dSSW": n_real_d, "rate_real": round(rate_real, 3),
           "rate_null": round(rate_null, 3), "rate_p_value": p_rate,
           "null_centres_on_zero": round(float(np.mean(null["ao"]["all"])), 4),
           "outcomes": {}}

    print(f"\n{'outcome':6s} {'real dSSW':>11s} {'null crit':>11s} {'null matched':>13s} "
          f"{'reproduced':>11s}")
    print("-" * 60)
    for k in OUTCOMES:
        c = np.array(null[k]["crit"])
        m = np.array(null[k]["matched"])
        rep = 100 * m.mean() / real_d[k] if real_d[k] else np.nan
        res["outcomes"][k] = {
            "real_all_events": round(real_all[k], 4),
            "real_dSSW_only": round(real_d[k], 4),
            "null_all_events": round(float(np.mean(null[k]["all"])), 4),
            "null_same_criterion": round(float(c.mean()), 4),
            "null_rate_matched": round(float(m.mean()), 4),
            "null_rate_matched_CI95": [round(float(np.percentile(m, 2.5)), 4),
                                       round(float(np.percentile(m, 97.5)), 4)],
            "pct_of_real_reproduced_by_selection": round(float(rep), 1),
            "p_real_vs_matched_null": float((m <= real_d[k]).mean())}
        print(f"{k:6s} {real_d[k]:+11.3f} {c.mean():+11.3f} {m.mean():+13.3f} "
              f"{rep:10.0f}%")

    (RESULTS / "selection_on_outcome.json").write_text(json.dumps(res, indent=2),
                                                    encoding="utf8")
    print("\n  'null matched' selects the same FRACTION of pseudo-events as the")
    print("  criterion selects of real ones, so it is the fair comparison.")
    print("  'reproduced' = how much of the published-style dSSW number that")
    print("  selection alone accounts for, with a true effect of exactly zero.")
    print("\nSaved -> selection_on_outcome.json")


if __name__ == "__main__":
    main()
