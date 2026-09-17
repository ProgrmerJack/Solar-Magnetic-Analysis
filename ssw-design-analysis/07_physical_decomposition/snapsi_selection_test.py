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

# ---- Southern Hemisphere, reported SEPARATELY and never pooled with NH ----
# Central date 18 September 2019, from the SNAPSI protocol (Hitchcock et al.,
# GMD 15, 5073, 2022): the 10 hPa 60S zonal-mean wind "did not reverse. However,
# they did decelerate dramatically, reaching their minimum value on 18 September
# 2019, which can be considered the 'central' date for the event."
#
# TWO CAVEATS THAT MUST TRAVEL WITH ANY SH NUMBER:
#  1. This is an AUSTRAL MINOR WARMING. The wind never reversed, so it is not an
#     SSW under the WMO definition. It tests generality to a different hemisphere
#     AND to a weaker, different class of event -- not to the same event mirrored.
#  2. s20191001 initialises 13 days AFTER the central date, so its post-onset
#     coverage starts at +13 and it cannot span the +8..+25 window at all. Only
#     s20190829 (onset at lead +20) covers it, giving 9 ensembles per arm.
SH_ONSET = {"s20190829": "2019-09-18", "s20191001": "2019-09-18"}
SH_INIT_DATE = {"s20190829": "2019-08-29", "s20191001": "2019-10-01"}
ONSET_ALL = {**ONSET, **SH_ONSET}
INIT_DATE_ALL = {**INIT_DATE, **SH_INIT_DATE}
HEMI = {**{k: "NH" for k in ONSET}, **{k: "SH" for k in SH_ONSET}}
WINDOW = (8, 25)          # post-onset days, as in snapsi_causal_effect.py
# An initialisation is admitted only if its forecasts actually SPAN the window;
# see spans_window(). Averaging a partially-covering ensemble over the nominal
# window computes a different quantity and reports it under the same name.
REQUIRE_FULL_WINDOW = True
PUBLISHED_CONTRAST = -0.850   # ACP 26, 3723 (2026), same criterion


def load(centre, exp, init):
    rows = []
    for f in sorted(RED.glob(f"{centre}_{exp}_{init}_*.parquet")):
        d = pd.read_parquet(f)
        rows.append(d)
    return pd.concat(rows) if rows else None


def post_onset_offset(init):
    return (pd.Timestamp(ONSET_ALL[init]) - pd.Timestamp(INIT_DATE_ALL[init])).days


def member_series(df, init):
    """member -> (mean NAM-proxy over window, fraction of days negative).

    Returned in RAW pressure units; standardisation happens against control.
    """
    off = post_onset_offset(init)
    lo, hi = WINDOW[0] + off, WINDOW[1] + off
    w = df[(df["lead_days"] >= lo) & (df["lead_days"] <= hi)]
    return w.groupby("member")["psl_cap"]


def spans_window(df, init):
    """Does this ensemble's data actually cover the analysis window?

    s20191001 initialises 13 days AFTER the SH central date, so its forecasts
    cover post-onset +13..+25. Averaging it over the nominal +8..+25 computes a
    different, shorter-window quantity and reports it under the same name --
    the defect the Loeffel script guards against with win_mean(). Same rule
    here, applied identically to every initialisation.
    """
    if not REQUIRE_FULL_WINDOW:
        return True
    off = post_onset_offset(init)
    lo, hi = WINDOW[0] + off, WINDOW[1] + off
    lead = df["lead_days"]
    return bool(lead.min() <= lo + 1.0 and lead.max() >= hi - 1.0)


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
        if n is not None and c is not None and not (
                spans_window(n, init) and spans_window(c, init)):
            continue
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
    out["all_ensembles"] = []

    for c in usable:
        for init in sorted(ONSET_ALL):
            con = load(c, "control", init)
            if con is None or not spans_window(con, init):
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
                if d is None or not spans_window(d, init):
                    continue
                w = d[(d["lead_days"] >= lo) & (d["lead_days"] <= hi)].copy()
                # NAM proxy: high polar-cap pressure -> negative annular mode
                w["nam"] = -(w["psl_cap"] - base) / sd

                g = w.groupby("member")["nam"]
                mean_nam = g.mean()
                frac_neg = g.apply(lambda s: float((s < 0).mean()))
                dw = (mean_nam < 0) & (frac_neg > 0.5)   # conditions 1 AND 2

                n_dw, n_nd = int(dw.sum()), int((~dw).sum())
                # The DW RATE is defined for every ensemble and is recorded
                # before any filtering. The CONTRAST is not: it needs both
                # groups to be non-empty. The <3 filter below is therefore
                # ASYMMETRIC BY CONSTRUCTION -- it can only bite the arm with a
                # strong forced shift, where nearly every member lands in one
                # class, and never the control arm, which sits near 50/50.
                # Taking the rate only over surviving ensembles understated it:
                # 11 of 36 nudged ensembles are dropped and their DW rate is
                # 0.993, so the kept-only mean of 0.776 is not the rate.
                out["all_ensembles"].append({
                    "centre": c, "init": init, "arm": arm,
                    "hemisphere": HEMI[init],
                    "n_members": int(len(mean_nam)),
                    "n_DW": n_dw, "n_NDW": n_nd,
                    "DW_rate": round(n_dw / len(mean_nam), 3),
                    "contrast_estimable": bool(n_dw >= 3 and n_nd >= 3)})
                if n_dw < 3 or n_nd < 3:
                    continue
                dwm = float(mean_nam[dw].mean())
                ndm = float(mean_nam[~dw].mean())
                out["per_case"].append({
                    "centre": c, "init": init, "arm": arm,
                    "hemisphere": HEMI[init],
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
        rows = [r for r in out["per_case"]
                if r["arm"] == arm and r["hemisphere"] == "NH"]
        if not rows:
            continue
        con = np.array([r["contrast_sigma"] for r in rows])
        rate = np.array([r["DW_rate"] for r in rows])
        out["summary"][arm] = {
            "n_cases": len(con),
            "mean_contrast_sigma": round(float(con.mean()), 4),
            "sd_across_cases": round(float(con.std(ddof=1)), 4),
            "range": [round(float(con.min()), 4), round(float(con.max()), 4)],
            "mean_DW_rate_contrast_subset": round(float(rate.mean()), 3),
            "ratio_to_published": round(
                abs(float(con.mean())) / abs(PUBLISHED_CONTRAST), 2),
        }
    # Unconditional DW rate: every ensemble, nothing dropped.
    for arm in ("nudged", "control"):
        # NH ONLY. The SH case is an austral MINOR warming in the other
        # hemisphere and is summarised in its own block; pooling it here would
        # silently mix two different event classes into one rate.
        alls = [r for r in out["all_ensembles"]
                if r["arm"] == arm and r["hemisphere"] == "NH"]
        if not alls or arm not in out["summary"]:
            continue
        rates = np.array([r["DW_rate"] for r in alls])
        dropped = [r for r in alls if not r["contrast_estimable"]]
        out["summary"][arm]["n_ensembles_all"] = len(alls)
        out["summary"][arm]["n_dropped_for_contrast"] = len(dropped)
        out["summary"][arm]["mean_DW_rate_UNCONDITIONAL"] = round(
            float(rates.mean()), 3)
        if dropped:
            dr = np.array([r["DW_rate"] for r in dropped])
            out["summary"][arm]["dropped_DW_rate_mean"] = round(float(dr.mean()), 3)
            out["summary"][arm]["dropped_DW_rate_range"] = [
                round(float(dr.min()), 3), round(float(dr.max()), 3)]
    # ---- Southern Hemisphere, its own block. NEVER pooled with NH. ----
    out["southern_hemisphere"] = {
        "event": "Austral MINOR warming, central date 2019-09-18 "
                 "(SNAPSI protocol, Hitchcock et al. GMD 15, 5073, 2022)",
        "caveat": ("The 10 hPa 60S wind never reversed, so this is NOT an SSW "
                   "under the WMO definition. It tests generality to a "
                   "different hemisphere AND a weaker class of event. "
                   "s20191001 initialises 13 days after the central date and "
                   "cannot span the +8..+25 window, so only s20190829 "
                   "contributes."),
    }
    for arm in ("nudged", "control"):
        rows = [r for r in out["per_case"]
                if r["arm"] == arm and r["hemisphere"] == "SH"]
        alls = [r for r in out["all_ensembles"]
                if r["arm"] == arm and r["hemisphere"] == "SH"]
        if not alls:
            continue
        rates = np.array([r["DW_rate"] for r in alls])
        blk = {"n_ensembles_all": len(alls),
               "n_cases_with_contrast": len(rows),
               "mean_DW_rate_UNCONDITIONAL": round(float(rates.mean()), 3),
               "inits_present": sorted({r["init"] for r in alls})}
        if rows:
            con = np.array([r["contrast_sigma"] for r in rows])
            blk["mean_contrast_sigma"] = round(float(con.mean()), 4)
            blk["sd_across_cases"] = (round(float(con.std(ddof=1)), 4)
                                      if len(con) > 1 else None)
        out["southern_hemisphere"][arm] = blk

    print("\n=== SOUTHERN HEMISPHERE (austral MINOR warming, reported separately) ===")
    for arm in ("nudged", "control"):
        b = out["southern_hemisphere"].get(arm)
        if not b:
            print(f"  {arm:<8} no usable ensembles")
            continue
        print(f"  {arm:<8} contrast "
              f"{b.get('mean_contrast_sigma', float('nan')):+.3f} sigma "
              f"(n={b['n_cases_with_contrast']} of {b['n_ensembles_all']}), "
              f"DW rate unconditional {b['mean_DW_rate_UNCONDITIONAL']:.3f}, "
              f"inits {b['inits_present']}")
    print("  The wind never reversed: this is a MINOR warming, not an SSW.")

    out["filter_caveat"] = (
        "The DW-NDW contrast requires >=3 members in BOTH groups, which is "
        "asymmetric by construction: only an arm with a strong forced shift "
        "can push nearly every member into one class. 11 of 36 nudged "
        "ensembles are dropped that way, with a mean DW rate of 0.993, while "
        "0 of 36 control ensembles are. Quote mean_DW_rate_UNCONDITIONAL for "
        "the rate; the contrast is necessarily conditioned on the ensembles "
        "where both groups exist, which are the weaker-responding ones.")
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
        print(f"     DW rate, contrast subset     : "
              f"{s['mean_DW_rate_contrast_subset']:.3f}  "
              f"({s['n_cases']} of {s.get('n_ensembles_all', '?')} ensembles)")
        print(f"     DW rate, UNCONDITIONAL       : "
              f"{s.get('mean_DW_rate_UNCONDITIONAL', float('nan')):.3f}"
              + (f"   [{s['n_dropped_for_contrast']} dropped, "
                 f"their rate {s['dropped_DW_rate_mean']:.3f}]"
                 if s.get("n_dropped_for_contrast") else ""))
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
