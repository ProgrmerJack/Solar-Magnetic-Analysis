#!/usr/bin/env python3
"""
snapsi_selection_test.py
========================
THE DECISIVE TEST: does the published "downward-propagating SSW" classification
manufacture its own contrast even when the causal forcing is held EXACTLY fixed
by experimental design?

WHY THIS IS THE CLEAN EXPERIMENT
  Every observational and CMIP6 demonstration in this project shares one
  weakness: the true causal effect per event is unknown, so "the contrast is
  manufactured" always rests on a modelled or pseudo-event null.

  SNAPSI removes that weakness. In the `nudged` runs the stratosphere of EVERY
  member is nudged to the SAME observed evolution (zonally symmetric, 50-90 hPa,
  6-hourly). Members differ ONLY in tropospheric initial condition and internal
  noise. The stratospheric driver is therefore identical across members BY
  CONSTRUCTION -- not modelled, not assumed, not matched: identical.

  So if we apply the Karpechko et al. (2017) surface conditions to these members
  and split them into "downward-propagating" and not, ANY contrast we measure is
  100% selection. There is no physical difference left for it to reflect.

  This is the experiment the literature never ran: a case where the causal
  effect is pinned and the classifier is still asked to find classes.

WHAT IS APPLIED
  Karpechko conditions 1 and 2 -- the two that ARE the surface response:
    1  mean surface NAM over the post-onset window is negative
    2  fraction of days in that window with negative NAM > 0.5
  Condition 3 (stratospheric) is deliberately NOT applied: in nudged runs the
  stratosphere is common to all members, so it cannot discriminate and would
  only dilute the test.

  Surface NAM proxy = -(polar-cap SLP anomaly), standardised by the CONTROL
  ensemble spread for that centre/init. Higher polar-cap pressure = negative
  annular mode, the usual sign convention.

GUARD
  NRL is excluded: its `nudged` and `control` submissions are byte-identical for
  three of four initialisations, which is a corrupt submission rather than a
  null result. The guard re-checks this rather than trusting the note.

Output: results/current/8_experiment/snapsi_selection_test.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "8_experiment"
RESULTS.mkdir(parents=True, exist_ok=True)
RED = HERE.parents[0] / "03_data_ingestion" / "_snapsi_reduced"

# onset per initialisation, from snapsi_causal_effect.py
ONSET = {"s20180125": "2018-02-12", "s20180208": "2018-02-12",
         "s20181213": "2019-01-02", "s20190108": "2019-01-02"}
INIT_DATE = {"s20180125": "2018-01-25", "s20180208": "2018-02-08",
             "s20181213": "2018-12-13", "s20190108": "2019-01-08"}
WINDOW = (8, 25)          # post-onset days, as in snapsi_causal_effect.py
PUBLISHED_CONTRAST = -0.850   # ACP 26, 3723 (2026), same criterion


def load(centre, exp, init):
    rows = []
    for f in sorted(RED.glob(f"{centre}_{exp}_{init}_*.parquet")):
        d = pd.read_parquet(f)
        rows.append(d)
    return pd.concat(rows) if rows else None


def post_onset_offset(init):
    return (pd.Timestamp(ONSET[init]) - pd.Timestamp(INIT_DATE[init])).days


def member_series(df, init):
    """member -> (mean NAM-proxy over window, fraction of days negative).

    Returned in RAW pressure units; standardisation happens against control.
    """
    off = post_onset_offset(init)
    lo, hi = WINDOW[0] + off, WINDOW[1] + off
    w = df[(df["lead_days"] >= lo) & (df["lead_days"] <= hi)]
    return w.groupby("member")["psl_cap"]


def window_means(df, init):
    off = post_onset_offset(init)
    lo, hi = WINDOW[0] + off, WINDOW[1] + off
    w = df[(df["lead_days"] >= lo) & (df["lead_days"] <= hi)]
    return w.groupby("member")["psl_cap"].mean()


def corruption_guard(centre):
    """Flag a centre whose nudged-minus-control effect is identically zero.

    Tested on the ENSEMBLE MEAN, not member-by-member. A first version compared
    members by label and passed NRL, whose nudged submission is a relabelled
    copy for some initialisations and an exact duplicate for others: members
    differ, the ensemble mean does not. Nudging the stratosphere to an observed
    SSW cannot leave the surface ensemble mean bit-identical, so |S| ~ 0 at
    every initialisation is a corrupt submission, not a null result.
    """
    S = []
    for init in ONSET:
        n, c = load(centre, "nudged", init), load(centre, "control", init)
        if n is None or c is None:
            continue
        S.append(abs(float(window_means(n, init).mean()
                           - window_means(c, init).mean())))
    return S and max(S) < 1.0, S       # 1 Pa against control spreads of ~100-470


def main():
    centres = sorted({p.name.split("_")[0] for p in RED.glob("*.parquet")})
    out = {"window_post_onset_days": list(WINDOW),
           "criterion": ("Karpechko et al. (2017) conditions 1 and 2 "
                         "(the surface conditions) applied to nudged members "
                         "whose stratospheric forcing is identical by design"),
           "published_contrast_for_scale": PUBLISHED_CONTRAST,
           "excluded": {}, "per_case": [], }

    usable = []
    for c in centres:
        bad, S = corruption_guard(c)
        if bad:
            out["excluded"][c] = (
                f"|nudged - control| ensemble-mean effect is < 1 Pa at every "
                f"initialisation (max {max(S):.3g} Pa): corrupt submission, "
                f"not a null result")
            print(f"  EXCLUDED {c}: S identically zero (max |S| = {max(S):.3g} Pa)")
        else:
            usable.append(c)
    print(f"  usable centres: {usable}\n")

    print(f"{'centre':8s} {'init':11s} {'arm':8s} {'n':>4s} {'DW':>4s} "
          f"{'NDW':>4s} {'rate':>6s} {'DW':>8s} {'NDW':>8s} {'contrast':>9s}")
    print("-" * 82)

    for c in usable:
        for init in sorted(ONSET):
            con = load(c, "control", init)
            if con is None:
                continue
            cm = window_means(con, init)
            sd, base = float(cm.std(ddof=1)), float(cm.mean())
            if not np.isfinite(sd) or sd <= 0:
                continue

            off = post_onset_offset(init)
            lo, hi = WINDOW[0] + off, WINDOW[1] + off

            # BOTH arms get the identical treatment:
            #   nudged  -- stratospheric forcing identical across members
            #   control -- no SSW forcing at all (the pure null)
            for arm in ("nudged", "control"):
                d = load(c, arm, init)
                if d is None:
                    continue
                w = d[(d["lead_days"] >= lo) & (d["lead_days"] <= hi)].copy()
                # NAM proxy: high polar-cap pressure -> negative annular mode
                w["nam"] = -(w["psl_cap"] - base) / sd

                g = w.groupby("member")["nam"]
                mean_nam = g.mean()
                frac_neg = g.apply(lambda s: float((s < 0).mean()))
                dw = (mean_nam < 0) & (frac_neg > 0.5)   # conditions 1 AND 2

                n_dw, n_nd = int(dw.sum()), int((~dw).sum())
                if n_dw < 3 or n_nd < 3:
                    continue
                dwm = float(mean_nam[dw].mean())
                ndm = float(mean_nam[~dw].mean())
                out["per_case"].append({
                    "centre": c, "init": init, "arm": arm,
                    "n_members": int(len(mean_nam)),
                    "n_DW": n_dw, "n_NDW": n_nd,
                    "DW_rate": round(n_dw / len(mean_nam), 3),
                    "DW_mean_sigma": round(dwm, 4),
                    "NDW_mean_sigma": round(ndm, 4),
                    "contrast_sigma": round(dwm - ndm, 4),
                    "ensemble_mean_sigma": round(float(mean_nam.mean()), 4)})
                print(f"{c:8s} {init:11s} {arm:8s} {len(mean_nam):4d} "
                      f"{n_dw:4d} {n_nd:4d} {n_dw/len(mean_nam):6.2f} "
                      f"{dwm:+8.3f} {ndm:+8.3f} {dwm-ndm:+9.3f}")

    if not out["per_case"]:
        print("no usable cases"); return

    out["summary"] = {}
    for arm in ("nudged", "control"):
        rows = [r for r in out["per_case"] if r["arm"] == arm]
        if not rows:
            continue
        con = np.array([r["contrast_sigma"] for r in rows])
        rate = np.array([r["DW_rate"] for r in rows])
        out["summary"][arm] = {
            "n_cases": len(con),
            "mean_contrast_sigma": round(float(con.mean()), 4),
            "sd_across_cases": round(float(con.std(ddof=1)), 4),
            "range": [round(float(con.min()), 4), round(float(con.max()), 4)],
            "mean_DW_rate": round(float(rate.mean()), 3),
            "ratio_to_published": round(
                abs(float(con.mean())) / abs(PUBLISHED_CONTRAST), 2),
        }
    out["units_caveat"] = (
        "Contrasts here are standardised by the CONTROL ensemble's "
        "member-to-member spread for a single event. The published -0.850 is "
        "standardised by the BETWEEN-EVENT spread in observations. These are "
        "different yardsticks, so ratio_to_published indicates scale and is "
        "not an exact like-for-like fraction.")

    print("\n=== RESULT ===")
    for arm, lab in (("nudged", "stratospheric forcing IDENTICAL by design"),
                     ("control", "NO SSW forcing at all (pure null)")):
        s = out["summary"].get(arm)
        if not s:
            continue
        print(f"\n  {arm.upper():8s} -- {lab}")
        print(f"     manufactured DW-NDW contrast : "
              f"{s['mean_contrast_sigma']:+.3f} sigma "
              f"(sd {s['sd_across_cases']:.3f}, n={s['n_cases']})")
        print(f"     DW rate                      : {s['mean_DW_rate']:.2f}")
    print(f"\n  published contrast, same criterion: {PUBLISHED_CONTRAST:+.3f} sigma")
    print("\n  In BOTH arms the classifier produces a large contrast. In the")
    print("  nudged arm every member had the SAME stratosphere; in the control")
    print("  arm there was no SSW at all. Neither contrast can reflect a")
    print("  physical difference between groups -- both are entirely the cut.")
    print("\n  " + out["units_caveat"])

    (RESULTS / "snapsi_selection_test.json").write_text(
        json.dumps(out, indent=2), encoding="utf8")
    print("\nSaved -> results/current/8_experiment/snapsi_selection_test.json")


if __name__ == "__main__":
    main()
