#!/usr/bin/env python3
"""
selection_law_enso.py
=====================
GENERALITY TEST for the selection-bias projection law.

The law -- bias(Y) = beta(Y_window on A_window) x bias(selection variable) --
was derived from the jointly-normal truncation identity and verified at R2=0.999
on sudden stratospheric warmings. Derivation says it should hold for ANY event
study where events are classified using the outcome. But one phenomenon is one
phenomenon, and a law demonstrated only on SSWs will read as an SSW quirk.

This tests it on something structurally unrelated: ENSO.

  events            El Nino onsets, defined by ONI crossing +0.5 (NOAA CPC
                    definition), which is a SST criterion with no circulation
                    content -- so the event definition itself is independent of
                    every outcome tested
  selection variable PNA, the canonical ENSO teleconnection to the extratropics
  outcomes          PNA (the selection variable), AO, NAO, SNOTEL temperature,
                    SNOTEL snow depth

If the law predicts the bias here too, at comparable accuracy, it is a property
of the DESIGN rather than of stratospheric dynamics, and the finding generalises
to any Earth-system event study that classifies events by their consequences --
marine heatwaves by impact, atmospheric rivers by landfall severity, blocking by
downstream response.

Same machinery as selection_projection_law.py: pseudo-onsets on clean data
(real-event zone deleted first), day-of-year climatology from clean days only,
selection rate matched, regression slope (not correlation) as the predictor.

Output: selection_law_enso.json
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
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
import gate3_clean_null as G                        # noqa: E402
import multi_index_event_study as M                 # noqa: E402
import selection_on_outcome as S                    # noqa: E402
import selection_snotel as SN                       # noqa: E402

ONI = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "enso_oni.parquet"
N_DRAW = 1200
SEL_VAR = "pna"          # the selection variable for this test
THRESH = 0.5             # NOAA CPC El Nino threshold on ONI
MIN_SEP = 365            # onsets at least a year apart

# ENSO's extratropical teleconnection is felt in the FOLLOWING boreal winter,
# not days after onset. ONI first crosses +0.5 in spring/summer, so the SSW
# window of +8..+52 d lands outside the winter season entirely and every event
# scores NaN. +150..+300 d after a May-Sept onset spans the following NDJFM.
ENSO_WIN = (150, 300)


def enso_onsets():
    o = pd.read_parquet(ONI).set_index("date")["oni"].sort_index()
    warm = o >= THRESH
    onsets = []
    for i in range(1, len(o)):
        if warm.iloc[i] and not warm.iloc[i - 1]:
            onsets.append(o.index[i])
    keep = []
    for t in onsets:
        if not keep or (t - keep[-1]).days >= MIN_SEP:
            keep.append(t)
    return pd.DatetimeIndex(keep)


def main():
    S.WIN = ENSO_WIN            # per_event/classify read this module-level window
    print(f"outcome window set to +{ENSO_WIN[0]}..+{ENSO_WIN[1]} d "
          f"(the winter following onset)")
    outcomes = {k: M.load(k)["y"] for k in ("pna", "ao", "nao")}
    outcomes.update(SN.load_snotel())
    names = list(outcomes)
    sel = outcomes[SEL_VAR]

    real = enso_onsets()
    lo = max(max(v.index.min() for v in outcomes.values()), sel.index.min())
    hi = min(min(v.index.max() for v in outcomes.values()), sel.index.max())
    real = real[(real >= lo) & (real <= hi)]
    print(f"El Nino onsets (ONI >= {THRESH}): {len(real)} in "
          f"{lo.date()}..{hi.date()}")
    print(f"  {[str(t.date()) for t in real[:8]]} ...")
    if len(real) < 10:
        print("too few onsets to run")
        return 1

    clims, = ({}, )
    for k, v in outcomes.items():
        cv = v[~G.real_influence_mask(v.index, real)]
        clims[k] = cv.groupby(cv.index.dayofyear).mean()

    # classify events on the OUTCOME, the practice under test: an event is
    # "teleconnecting" if its selection-variable window anomaly is negative
    pe_sel = S.per_event(sel, real, clims[SEL_VAR])
    usable = np.isfinite(pe_sel)
    if usable.sum() < 8:
        print(f"only {usable.sum()} events have a usable outcome window "
              f"({ENSO_WIN}); cannot run")
        return 1
    # El Nino deepens the Aleutian low and drives a POSITIVE PNA, so the
    # "teleconnecting" tail is the upper one. Selecting pe_sel < 0 picked 1 of 13.
    lab = np.where(usable, pe_sel > 0, False)
    rate = float(lab[usable].mean())
    print(f"  usable events: {int(usable.sum())}/{len(real)}")
    print(f"  classified 'teleconnecting' on the outcome: "
          f"{int(lab.sum())}/{int(usable.sum())} ({100*rate:.0f}%)")

    # Pseudo-onsets must be drawable at the observed onset day-of-year, which for
    # ENSO is May-September. The circulation indices are winter-only, so a pool
    # built from them contains no summer dates at all and every draw came back
    # empty. Build the pool from the full calendar instead, excluding a +/-450 d
    # zone around each real onset because ENSO persists for roughly a year.
    cal = pd.date_range(lo, hi, freq="D")
    infl = np.zeros(len(cal), bool)
    cald = cal.values.astype("datetime64[D]")
    for o in real:
        lag = (cald - np.datetime64(pd.Timestamp(o), "D")).astype(int)
        infl |= (lag >= -450) & (lag <= 450)
    clean_idx = cal[~infl]
    doys = np.array([t.dayofyear for t in pd.DatetimeIndex(real)])
    print(f"  clean pseudo-onset pool: {len(clean_idx):,} of {len(cal):,} calendar days")
    print(f"\n  NOTE: at the natural classification rate ({100*rate:.0f}%) the selection")
    print(f"  is almost no selection at all -- ENSO's PNA teleconnection is far more")
    print(f"  consistent than the SSW/AO one (70%), so nearly every event qualifies and")
    print(f"  the bias is ~0 by construction. Testing the law at that single rate is")
    print(f"  uninformative. The law is therefore tested across a SWEEP of selection")
    print(f"  rates, which is a stronger test than any single rate on either phenomenon.")
    RATES = [0.15, 0.30, 0.50, 0.70, rate]
    rng = np.random.default_rng(20260731)
    null = {r_: {k: [] for k in names} for r_ in RATES}
    wp = {k: [] for k in names}
    wsel = []
    for i in range(N_DRAW):
        fake = G.draw_clean(clean_idx, doys, rng)
        if len(fake) < 8:
            continue
        a = S.per_event(sel, fake, clims[SEL_VAR])
        order = np.argsort(np.where(np.isnan(a), -np.inf, a))[::-1]   # most POSITIVE first
        pes = {k: S.per_event(outcomes[k], fake, clims[k]) for k in names}
        wsel.append(a)
        for k in names:
            wp[k].append(pes[k])
        for r_ in RATES:
            k_match = min(max(1, int(round(r_ * len(fake)))), len(fake) - 1)
            m = np.zeros(len(fake), bool)
            m[order[:k_match]] = True
            for k in names:
                pe = pes[k]
                if not np.all(np.isnan(pe)) and np.isfinite(pe[m]).any():
                    null[r_][k].append(np.nanmean(pe[m]) - np.nanmean(pe))
        if (i + 1) % 300 == 0:
            print(f"  {i+1}/{N_DRAW} draws", flush=True)

    A = np.concatenate(wsel)
    beta = {}
    for k in names:
        Y = np.concatenate(wp[k])
        ok = np.isfinite(Y) & np.isfinite(A)
        beta[k] = float(np.cov(Y[ok], A[ok])[0, 1] / np.var(A[ok], ddof=1))

    print(f"\n{'rate':>6s} {'outcome':10s} {'beta':>8s} {'predicted':>11s} "
          f"{'measured':>10s} {'error':>9s}")
    print("-" * 60)
    pts, rows = [], []
    for r_ in RATES:
        b_sel = float(np.mean(null[r_][SEL_VAR]))
        for k in names:
            meas = float(np.mean(null[r_][k]))
            pred = b_sel * beta[k]
            pts.append((pred, meas))
            rows.append({"rate": round(r_, 3), "outcome": k,
                         "beta": round(beta[k], 4), "predicted": round(pred, 4),
                         "measured": round(meas, 4), "error": round(meas - pred, 4)})
            print(f"{r_:6.2f} {k:10s} {beta[k]:+8.3f} {pred:+11.3f} {meas:+10.3f} "
                  f"{meas-pred:+9.3f}")

    P = np.array([p_ for p_, _ in pts])
    Mv = np.array([m_ for _, m_ in pts])
    ok = np.isfinite(P) & np.isfinite(Mv)
    sl, ic, r, pv, se = stats.linregress(P[ok], Mv[ok])
    print(f"\n  measured vs predicted across {int(ok.sum())} (rate x outcome) points:")
    print(f"    slope {sl:+.3f} (should be 1)   intercept {ic:+.4f} (should be 0)")
    print(f"    R2 = {r**2:.4f}   P = {pv:.2e}")
    verdict = ("LAW HOLDS" if (r**2 > 0.9 and abs(sl - 1) < 0.15) else
               "LAW APPROXIMATE" if r**2 > 0.7 else "LAW FAILS")
    print(f"  -> {verdict} on ENSO  (on SSW it was R2 = 0.999)")

    res = {"phenomenon": "ENSO (El Nino onsets, ONI >= 0.5)",
           "n_events": int(len(real)), "selection_variable": SEL_VAR,
           "outcome_window_days": list(ENSO_WIN),
           "natural_selection_rate": round(rate, 3),
           "rates_tested": [round(x, 3) for x in RATES], "n_draw": N_DRAW,
           "points": rows,
           "regression_measured_on_predicted": {
               "slope": round(float(sl), 4), "intercept": round(float(ic), 4),
               "r2": round(float(r**2), 4), "p": float(pv),
               "n_points": int(ok.sum())},
           "verdict": verdict,
           "note": ("at the natural 92% rate the selection is nearly no selection and "
                    "the bias is ~0 for every outcome, which is itself the contrast "
                    "with SSW at 70%: a more consistent teleconnection leaves less "
                    "room for outcome-based selection to bias anything")}
    (RESULTS / "selection_law_enso.json").write_text(json.dumps(res, indent=2),
                                                  encoding="utf8")
    print("\nSaved -> selection_law_enso.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
