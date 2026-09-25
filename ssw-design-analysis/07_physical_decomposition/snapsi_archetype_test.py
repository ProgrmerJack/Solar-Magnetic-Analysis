#!/usr/bin/env python3
"""
snapsi_archetype_test.py
========================
WERE FEB-2018 AND JAN-2019 DIFFERENT EVENTS, OR TWO DRAWS FROM ONE DISTRIBUTION?

THE ARCHETYPES
  The two NH SNAPSI events are the field's standard pair of opposites. Rao,
  Garfinkel & White (2020, JGR-Atmos 125, e2019JD031919) describe Feb-2018 as a
  strong, downward-propagating event (north Eurasia ~4 C colder) and Jan-2019 as
  one whose surface response was weak or opposite (north Eurasia anomalously
  warm), and attribute the difference chiefly to SSW strength. Nebel et al.
  (2024, GRL 51, e2024GL110529) frame such non-propagating outcomes as potential
  forecast busts.

WHAT SNAPSI ALREADY SAYS, AND WHAT THIS ADDS
  With the stratosphere nudged to each observed event, the models give the two
  events comparable forced surface shifts (result K) and the same probability
  of a downward-propagating outcome (result L). This script
  asks where the OBSERVED outcome of each event falls inside its own nudged
  ensemble, and whether the observed difference between the two events is
  unusual for two members drawn from the two nudged ensembles.

DEFINITIONS -- IDENTICAL TO RESULT L, IMPORTED FROM IT
  Window post-onset days +8..+25. NAM proxy A = -(psl_cap - base)/sd, where base
  and sd are the mean and sd of the CONTROL ensemble's window means for that
  centre and initialisation. DW = mean A < 0 AND fraction of 6-hourly values
  with A < 0 > 0.5 (Karpechko et al. 2017 conditions 1 and 2). The observed
  value is ERA5 polar-cap (60-90N, cos-lat) psl, 6-hourly, at exactly the
  forecast times the window covers (acquire_era5_psl_cap.py).

TWO REFERENCES FOR THE OBSERVATION, BECAUSE THE MODELS HAVE BIASES
  primary   obs placed against the model's own control base, as members are
  adjusted  obs shifted by the model-minus-ERA5 offset at forecast leads 0-1 d
            (control ensemble mean minus ERA5 at the SAME forecast times; UKMO
            and Meteo-France start at 06 UTC, so their first day has four
            steps, not five), which removes a representation
            offset (psl extrapolation over Greenland differs between models) but
            not lead-dependent drift
  The EVENT-PAIR test is the robust one: a model bias common to both winters
  enters both events' anomalies and largely cancels in their difference.

PAIRS, DECLARED BEFORE THE DATA WERE EXAMINED
  primary     s20180125 (onset at lead 18 d) with s20181213 (lead 20 d) --
              the two initialisations with comparable lead to onset
  short-lead  s20180208 (lead 4 d) with s20190108 (lead -6 d), reported
              separately, as every SNAPSI result here reports s20190108

WHAT WOULD SHOW THE CLAIM WRONG
  A < 0 is the DOWNWARD sign, so the archetype reading -- Jan-2019 weak or
  opposite -- predicts observed A in the TOP tail of its nudged ensemble. The
  claim fails if observed Jan-2019 lies above the 97.5th percentile in most
  models (>= 5 of 9), or if the observed 2018-minus-2019 difference lies outside
  the central 95% of the model pair distribution in most models. (The first
  version of this text named the 2.5th percentile -- the wrong tail, which
  could never catch the archetype being right; found in review 2026-09-24.
  Both tails are reported.)

THE OBSERVED LABELS THEMSELVES
  The archetype reading presumes Jan-2019 IS non-propagating. The observed class
  is recomputed here with the project's own ERA5 implementation of Karpechko et
  al. (2017) (`08_literature_audit/era5_recompute_and_two_thirds.classify`,
  conditions 1-3), at 1000 hPa (Karpechko) and 850 hPa (the ACP 26, 3723 2026
  variant), over the published days +8..+52 and over this script's +8..+25.
  A label that flips between those choices is a threshold decision, not a
  property of the event.

LIMITS, STATED UP FRONT
  Two events. This tests the two events the field treats as canonical; it says
  nothing general about events it does not contain. One observation per event,
  so per-model percentiles share the same observed value and are not
  independent across models.

Output: results/current/8_experiment/snapsi_archetype_test.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "8_experiment"
ERA5 = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "era5_psl_cap_6h.parquet"
sys.path.insert(0, str(HERE))
import snapsi_selection_test as L                     # noqa: E402
sys.path.insert(0, str(HERE.parents[0] / "08_literature_audit"))
import era5_recompute_and_two_thirds as K17            # noqa: E402

OBS_ONSET = {"Feb-2018": "2018-02-12", "Jan-2019": "2019-01-02"}

EVENT = {"s20180125": "Feb-2018", "s20180208": "Feb-2018",
         "s20181213": "Jan-2019", "s20190108": "Jan-2019"}
PAIRS = {"primary": ("s20180125", "s20181213"),
         "short_lead": ("s20180208", "s20190108")}
OFFSET_LEADS = (0.0, 1.0)


def nam(values, base, sd):
    return -(np.asarray(values, float) - base) / sd


def classify(a_series):
    """(mean A, fraction negative, DW) for one 6-hourly A series."""
    m = float(np.mean(a_series))
    f = float(np.mean(np.asarray(a_series) < 0))
    return m, f, bool(m < 0 and f > 0.5)


def percentile(x, dist):
    dist = np.asarray(dist, float)
    return float((np.sum(dist < x) + 0.5 * np.sum(dist == x)) / len(dist))


def obs_window(obs, init, leads):
    """ERA5 values at EXACTLY the forecast leads the ensemble has in the window.

    Not the nominal window: an ensemble can be admitted by spans_window() while
    ending one 6-hourly step short, and the observation must be averaged over
    the same times as the members it is compared with.
    """
    lead = ((obs["time"] - pd.Timestamp(L.INIT_DATE[init])).dt.total_seconds()
            / 86400).round(4)
    want = np.round(np.asarray(leads, float), 4)
    w = obs[lead.isin(want)]
    if len(w) != len(want):
        raise ValueError(f"{init}: ERA5 has {len(w)} of {len(want)} forecast times")
    return w["psl_cap_N"].values


def main() -> int:
    obs = pd.read_parquet(ERA5)
    obs["time"] = pd.to_datetime(obs["time"])
    centres = sorted({p.name.split("_")[0] for p in L.RED.glob("*.parquet")})
    usable = [c for c in centres if not L.corruption_guard(c)[0]]
    print(f"usable centres: {usable}")

    rows, members, member_A = [], {}, {}
    for c in usable:
        for init in L.ONSET:                       # NH only
            con, nud = L.load(c, "control", init), L.load(c, "nudged", init)
            if con is None or nud is None or not (
                    L.spans_window(con, init) and L.spans_window(nud, init)):
                continue
            cm = L.window_means(con, init)
            base, sd = float(cm.mean()), float(cm.std(ddof=1))
            off = L.post_onset_offset(init)
            lo, hi = L.WINDOW[0] + off, L.WINDOW[1] + off

            def arm_stats(d):
                w = d[(d["lead_days"] >= lo) & (d["lead_days"] <= hi)]
                return {m: classify(nam(g["psl_cap"].values, base, sd))
                        for m, g in w.groupby("member")}
            ns, cs = arm_stats(nud), arm_stats(con)
            a_n = np.array([v[0] for v in ns.values()])
            a_c = np.array([v[0] for v in cs.values()])
            members[(c, init)] = a_n
            member_A[f"{c}|{init}"] = {"nudged": [round(float(x), 4) for x in a_n],
                                       "control": [round(float(x), 4) for x in a_c]}

            wl = nud[(nud["lead_days"] >= lo) & (nud["lead_days"] <= hi)]
            per_member = wl.groupby("member")["lead_days"].apply(lambda s: tuple(sorted(s)))
            if per_member.nunique() != 1:
                raise ValueError(f"{c} {init}: members sample different leads")
            ow = obs_window(obs, init, per_member.iloc[0])
            # representation offset at the start of the forecast
            # at the control ensemble's own forecast times: averaging ERA5 over
            # leads 0-1 d while a 06 UTC-start model has only 0.25-1 d moved the
            # UKMO offset by up to 25 Pa (methods audit 2026-09-25)
            c0 = con[(con["lead_days"] >= OFFSET_LEADS[0])
                     & (con["lead_days"] <= OFFSET_LEADS[1])]
            e0 = c0["psl_cap"].mean()
            o0 = float(obs_window(obs, init, np.unique(c0["lead_days"].values)).mean())
            offset = float(e0 - o0)

            a_obs, f_obs, dw_obs = classify(nam(ow, base, sd))
            a_adj, f_adj, dw_adj = classify(nam(ow + offset, base, sd))
            rows.append({
                "centre": c, "init": init, "event": EVENT[init],
                "n_nudged": int(len(a_n)), "n_control": int(len(a_c)),
                "control_sd_Pa": round(sd, 1),
                "init_offset_model_minus_era5_Pa": round(offset, 1),
                "nudged_DW_rate": round(float(np.mean([v[2] for v in ns.values()])), 3),
                "control_DW_rate": round(float(np.mean([v[2] for v in cs.values()])), 3),
                "nudged_mean_A": round(float(a_n.mean()), 3),
                "obs_A": round(a_obs, 3), "obs_frac_neg": round(f_obs, 3), "obs_DW": dw_obs,
                "obs_pct_in_nudged": round(percentile(a_obs, a_n), 3),
                "obs_pct_in_control": round(percentile(a_obs, a_c), 3),
                "obs_A_adjusted": round(a_adj, 3), "obs_DW_adjusted": dw_adj,
                "obs_pct_in_nudged_adjusted": round(percentile(a_adj, a_n), 3),
                "_a_obs": a_obs, "_a_adj": a_adj})

    df = pd.DataFrame(rows)
    out = {"window_post_onset_days": list(L.WINDOW), "pairs": PAIRS,
           "centres": usable, "era5_input": str(ERA5.relative_to(ROOT)),
           "per_ensemble": [{k: v for k, v in r.items() if not k.startswith("_")}
                            for r in rows]}

    print("\nper ensemble (A < 0 is the downward-propagating sign):")
    print(df[["centre", "init", "event", "nudged_DW_rate", "nudged_mean_A", "obs_A",
              "obs_DW", "obs_pct_in_nudged", "obs_pct_in_nudged_adjusted",
              "init_offset_model_minus_era5_Pa"]].to_string(index=False))

    summ = {}
    for init, g in df.groupby("init"):
        summ[init] = {
            "event": EVENT[init], "n_centres": int(len(g)),
            "obs_DW_in_centres": int(g.obs_DW.sum()),
            "nudged_DW_rate_mean": round(float(g.nudged_DW_rate.mean()), 3),
            "obs_pct_median": round(float(g.obs_pct_in_nudged.median()), 3),
            "obs_pct_range": [float(g.obs_pct_in_nudged.min()),
                              float(g.obs_pct_in_nudged.max())],
            "centres_obs_below_2p5pct": int((g.obs_pct_in_nudged < 0.025).sum()),
            "centres_obs_above_97p5pct": int((g.obs_pct_in_nudged > 0.975).sum()),
            "obs_pct_median_adjusted": round(float(g.obs_pct_in_nudged_adjusted.median()), 3),
            "centres_below_2p5pct_adjusted": int((g.obs_pct_in_nudged_adjusted < 0.025).sum()),
            "centres_above_97p5pct_adjusted": int((g.obs_pct_in_nudged_adjusted > 0.975).sum())}
    out["per_init"] = summ
    print("\nper initialisation:")
    for k, v in summ.items():
        print(f"  {k} {v['event']}: obs DW in {v['obs_DW_in_centres']}/{v['n_centres']} "
              f"centres' frames; nudged DW rate {v['nudged_DW_rate_mean']}; obs percentile "
              f"median {v['obs_pct_median']} (range {v['obs_pct_range']}); "
              f"<2.5%: {v['centres_obs_below_2p5pct']}, >97.5%: {v['centres_obs_above_97p5pct']} | "
              f"adjusted median {v['obs_pct_median_adjusted']}")

    # ------------------------------------------------ the event-pair test
    out["pair_test"] = {}
    for tag, (i18, i19) in PAIRS.items():
        per = []
        for c in usable:
            r18 = df[(df.centre == c) & (df.init == i18)]
            r19 = df[(df.centre == c) & (df.init == i19)]
            if r18.empty or r19.empty:
                continue
            a18, a19 = members[(c, i18)], members[(c, i19)]
            diffs = (a18[:, None] - a19[None, :]).ravel()   # every member pair
            d_obs = float(r18._a_obs.iloc[0] - r19._a_obs.iloc[0])
            d_adj = float(r18._a_adj.iloc[0] - r19._a_adj.iloc[0])
            pct = percentile(d_obs, diffs)
            pct_adj = percentile(d_adj, diffs)
            per.append({"centre": c, "n_pairs": int(len(diffs)),
                        "model_diff_mean": round(float(diffs.mean()), 3),
                        "model_diff_sd": round(float(diffs.std(ddof=1)), 3),
                        "obs_diff": round(d_obs, 3), "obs_diff_pct": round(pct, 3),
                        "two_sided_p": round(2 * min(pct, 1 - pct), 3),
                        "obs_diff_adjusted": round(d_adj, 3),
                        "obs_diff_pct_adjusted": round(pct_adj, 3)})
        p = pd.DataFrame(per)
        out["pair_test"][tag] = {
            "inits": [i18, i19], "per_centre": per,
            "median_obs_diff_pct": round(float(p.obs_diff_pct.median()), 3),
            "centres_outside_central_95": int(((p.obs_diff_pct < 0.025)
                                               | (p.obs_diff_pct > 0.975)).sum()),
            "centres_outside_central_95_adjusted": int(((p.obs_diff_pct_adjusted < 0.025)
                                                        | (p.obs_diff_pct_adjusted > 0.975)).sum()),
            "n_centres": int(len(p))}
        print(f"\npair test [{tag}] {i18} minus {i19} "
              f"(A_2018 - A_2019; negative = 2018 more downward):")
        print(p[["centre", "model_diff_mean", "model_diff_sd", "obs_diff",
                 "obs_diff_pct", "two_sided_p", "obs_diff_pct_adjusted"]].to_string(index=False))
        t = out["pair_test"][tag]
        print(f"  outside the central 95%: {t['centres_outside_central_95']} of "
              f"{t['n_centres']} centres (adjusted: "
              f"{t['centres_outside_central_95_adjusted']})")

    # ------------------------------------------ the observed labels, recomputed
    e = pd.read_parquet(K17.ERA5)
    if e.index.tz is not None:
        e.index = e.index.tz_convert("UTC").tz_localize(None)
    on = pd.DatetimeIndex(list(OBS_ONSET.values()))
    lab_rows = []
    print("\nobserved class, Karpechko conditions 1-3 on ERA5 (project implementation):")
    # The published window goes through the project's classify() unchanged.
    # The +8..+25 window CANNOT: classify() marks any window with < 20 days as
    # unclassifiable (False), and +8..+25 has 18 -- patching its window gave
    # Feb-2018 "NDW" at -0.98 sigma with 100% of days negative. So for that
    # window the same three conditions are evaluated explicitly.
    for level, col in (("1000hPa", "nam_1000"), ("850hPa", "nam_850")):
        lab, c1, c2, c3 = K17.classify(on, e[col], e["nam_150"])
        per = {(K17.SEL_WIN, ev): (bool(l), bool(a), bool(b), bool(c))
               for ev, l, a, b, c in zip(OBS_ONSET, lab, c1, c2, c3)}
        for ev, o in OBS_ONSET.items():
            sw = K17.win_days(e[col], pd.Timestamp(o), tuple(L.WINDOW))
            st = K17.win_days(e["nam_150"], pd.Timestamp(o), tuple(L.WINDOW))
            a, b = bool(sw.mean() < 0), bool((sw < 0).mean() > 0.5)
            c = bool((st < 0).mean() > 0.7)
            per[(tuple(L.WINDOW), ev)] = (a and b and c, a, b, c)
        for (win, ev), (l, a, b, c) in per.items():
            sw = K17.win_days(e[col], pd.Timestamp(OBS_ONSET[ev]), win)
            lab_rows.append({"event": ev, "level": level, "window": list(win),
                             "nam_mean": round(float(sw.mean()), 3),
                             "frac_negative": round(float((sw < 0).mean()), 3),
                             "c1": a, "c2": b, "c3": c,
                             "class": "DW" if l else "NDW",
                             "method": ("classify() as published" if win == K17.SEL_WIN
                                        else "conditions 1-3 evaluated explicitly")})
            print(f"  {ev} {level} days {win}: NAM mean {sw.mean():+.3f}, "
                  f"{(sw < 0).mean():.0%} negative, c1={a} c2={b} c3={c} "
                  f"-> {'DW' if l else 'NDW'}")
    out["observed_labels"] = lab_rows
    out["member_A"] = member_A          # per-member NAM proxy, for the display item

    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "snapsi_archetype_test.json").write_text(
        json.dumps(out, indent=1), encoding="utf8", newline="\n")
    print("\nSaved -> snapsi_archetype_test.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
