#!/usr/bin/env python3
"""
selection_projection_law.py
===========================
Turns the selection-bias observation into a quantitative, falsifiable law.

WHAT PROMPTED THIS
  selection_on_outcome.py found the Karpechko dSSW classification manufactures
  30-44% of the reported anomaly for AO, NAO and PNA under a true null.
  selection_snotel.py then found the bias for western-US SNOTEL temperature and
  snow depth has the OPPOSITE SIGN -- it attenuates rather than inflates.

  So "selecting on the outcome inflates impacts" is wrong as a general claim.
  The pattern instead looks like a PROJECTION: the selection is made on the
  surface NAM, and each outcome inherits the bias in proportion to how strongly
  it co-varies with that selection variable.

THE LAW BEING TESTED, AND WHY IT SHOULD HOLD
      bias(outcome)  =  beta(Y_window on A_window)  x  bias(selection variable)

  Derivation: the selection keeps events whose SELECTION-VARIABLE window mean
  A_w falls in the lowest k. For jointly normal (Y_w, A_w),
      E[Y_w | sel] - E[Y_w] = rho * (sigma_Y / sigma_A) * (E[A_w | sel] - E[A_w])
  The bracket is exactly the bias measured for the selection variable itself,
  and rho * sigma_Y/sigma_A is the ORDINARY REGRESSION SLOPE of Y_w on A_w. So
  the law is not an empirical curiosity; it is what a linear projection must give.

  TWO SPECIFICATION ERRORS WERE MADE BEFORE ARRIVING HERE, both recorded:

   1. DAILY instead of WINDOW correlation. The selection acts on the +8..+52
      window mean, not on daily values. Daily correlation gave R2=0.868 with NAO
      off by 66%.
   2. CORRELATION instead of REGRESSION SLOPE -- the sigma_Y/sigma_A factor was
      dropped. Window means of persistent variables (AO) have much larger spread
      than those of less persistent ones (NAO), so rho alone over-predicts the
      bias for every non-selection outcome. With window correlation but no sigma
      ratio: R2=0.823, NAO predicted -0.412 against -0.133 measured, a factor of
      three.

  Both were my errors, not evidence against the law. This version uses the
  regression slope, which is the quantity the derivation actually specifies.

  Slopes are computed on CLEAN pseudo-event windows only, so the real events
  never define the relationship used to predict their own bias.

  If this holds, it is genuinely useful: the bias for any outcome can be
  predicted from a single correlation without rerunning the null, its SIGN is
  known in advance, and "selection is conservative" can never be assumed.

  If it fails, the honest claim collapses back to "measure the bias separately
  for every outcome", which is still actionable but far weaker.

DESIGN
  All outcomes on ONE common event set and ONE set of null draws, so the biases
  are directly comparable. Previous runs used 43 events for indices and 29 for
  SNOTEL at different selection rates, which are not comparable.

Output: selection_projection_law.json
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
import selection_on_outcome as S                    # noqa: E402
import selection_snotel as SN                       # noqa: E402

N_DRAW = 1200


def main():
    ao = M.load("ao")["y"]
    snam = S.strat_nam()
    outcomes = {k: M.load(k)["y"] for k in ("ao", "nao", "pna")}
    outcomes.update(SN.load_snotel())
    names = list(outcomes)

    real = load_catalogue("primary")
    # common window: every outcome must cover it
    lo = max(max(v.index.min() for v in outcomes.values()), ao.index.min())
    hi = min(min(v.index.max() for v in outcomes.values()), ao.index.max())
    real = real[(real >= lo) & (real <= hi)]
    print(f"common period {lo.date()}..{hi.date()}, {len(real)} events")

    # clean climatologies and the clean-data correlation with the selection variable
    clims, clean_series = {}, {}
    for k, v in outcomes.items():
        cv = v[~G.real_influence_mask(v.index, real)]
        clims[k] = cv.groupby(cv.index.dayofyear).mean()
        clean_series[k] = cv
    ao_clean = clean_series["ao"]

    # window-level correlation: accumulated below from the same clean pseudo-event
    # windows the null uses, because that is the level the selection acts on
    corr = {}

    lab = S.classify(real, ao, snam)
    rate = float(lab.mean())
    print(f"dSSW {int(lab.sum())}/{len(real)} ({100 * rate:.0f}%)\n")

    clean_idx = ao[~G.real_influence_mask(ao.index, real)].index
    clean_idx = clean_idx[(clean_idx >= lo) & (clean_idx <= hi)]
    doys = np.array([t.dayofyear for t in pd.DatetimeIndex(real)])
    rng = np.random.default_rng(20260731)
    null = {k: [] for k in names}
    win_pairs = {k: [] for k in names}      # per-window values, for the correlation
    win_sel = []
    for i in range(N_DRAW):
        fake = G.draw_clean(clean_idx, doys, rng)
        if len(fake) < 8:
            continue
        sel = S.per_event(ao, fake, clims["ao"])
        k_match = max(1, int(round(rate * len(fake))))
        order = np.argsort(np.where(np.isnan(sel), np.inf, sel))
        m = np.zeros(len(fake), bool)
        m[order[:k_match]] = True
        win_sel.append(sel)
        for k in names:
            pe = S.per_event(outcomes[k], fake, clims[k])
            win_pairs[k].append(pe)
            if not np.all(np.isnan(pe)):
                null[k].append(np.nanmean(pe) - np.nanmean(pe[m]))  # all minus selected
        if (i + 1) % 300 == 0:
            print(f"  {i + 1}/{N_DRAW} draws", flush=True)

    A = np.concatenate(win_sel)
    rho = {}
    for k in names:
        Y = np.concatenate(win_pairs[k])
        ok = np.isfinite(Y) & np.isfinite(A)
        if ok.sum() > 100:
            # regression slope beta = cov(Y,A)/var(A) = rho * sigma_Y/sigma_A
            corr[k] = float(np.cov(Y[ok], A[ok])[0, 1] / np.var(A[ok], ddof=1))
            rho[k] = float(np.corrcoef(Y[ok], A[ok])[0, 1])
        else:
            corr[k] = rho[k] = np.nan
    print(f"\nwindow-level regression slopes from {len(A):,} clean pseudo-event windows")

    # bias defined as (selected - all), sign consistent with earlier scripts
    bias = {k: -float(np.mean(null[k])) for k in names}
    b_ao = bias["ao"]

    print(f"\n{'outcome':10s} {'rho':>7s} {'beta':>8s} {'predicted':>11s} "
          f"{'measured':>10s} {'error':>9s}")
    print("-" * 60)
    rows = []
    for k in names:
        pred = b_ao * corr[k]
        rows.append((k, corr[k], pred, bias[k], bias[k] - pred))
        print(f"{k:10s} {rho[k]:+7.3f} {corr[k]:+8.3f} {pred:+11.3f} {bias[k]:+10.3f} "
              f"{bias[k] - pred:+9.3f}")

    c = np.array([r[1] for r in rows])
    b = np.array([r[3] for r in rows])
    ok = np.isfinite(c) & np.isfinite(b)
    sl, ic, r, p, se = stats.linregress(c[ok], b[ok])
    print(f"\n  regression of measured bias on correlation:")
    print(f"    slope {sl:+.3f} (bias at corr=1 should equal bias_AO = {b_ao:+.3f})")
    print(f"    intercept {ic:+.3f} (should be ~0)   r = {r:+.3f}   R2 = {r**2:.3f}   P = {p:.4f}")
    verdict = ("LAW HOLDS" if (r ** 2 > 0.9 and abs(ic) < 0.05) else
               "LAW APPROXIMATE" if r ** 2 > 0.7 else "LAW FAILS")
    print(f"    -> {verdict}")

    res = {"n_events": int(len(real)), "n_draw": N_DRAW, "selection_rate": round(rate, 3),
           "bias_of_selection_variable": round(b_ao, 4),
           "outcomes": {k: {"corr_with_AO": round(cc, 4),
                            "predicted_bias": round(pp, 4),
                            "measured_bias": round(bb, 4),
                            "error": round(ee, 4)}
                        for k, cc, pp, bb, ee in rows},
           "regression": {"slope": round(float(sl), 4), "intercept": round(float(ic), 4),
                          "r": round(float(r), 4), "r2": round(float(r ** 2), 4),
                          "p": float(p)},
           "verdict": verdict}
    (RESULTS / "selection_projection_law.json").write_text(json.dumps(res, indent=2),
                                                        encoding="utf8")
    print("\nSaved -> selection_projection_law.json")


if __name__ == "__main__":
    main()
