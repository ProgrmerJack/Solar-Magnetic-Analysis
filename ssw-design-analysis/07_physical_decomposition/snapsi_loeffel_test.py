#!/usr/bin/env python3
"""
snapsi_loeffel_test.py
======================
IS DOWNWARD COUPLING A PROPERTY OF THE EVENT, OR OF WHAT FOLLOWS IT?

THE PAPER, IN ITS PUBLISHED FORM
  Loeffel, Rupp, Kiefer, Pinto, Birner, Garny (2026), "Quantifying the
  tropospheric response to individual sudden stratospheric warmings revealed by
  an ensemble simulation strategy", Weather and Climate Dynamics 7, 895-913.
  (The earlier `loeffel2025_missing_control.py` in 08_literature_audit was
  written against the EGUsphere preprint and tests their THRESHOLD diagnostic on
  observations. This tests the published CORRELATION, in a designed experiment.
  Different claim, different data, different question.)

  Verified from the published text:
    predictor  100 hPa polar-cap (60-90N) standardised GPH anomaly, week 2,
               days 8-14 after onset
    response   1000 hPa polar-cap GPH anomaly, weeks 3-7
    result     r = 0.85, p < 0.001, n = 18 events
    computed   ACROSS EVENTS, on ENSEMBLE MEANS -- one point per event

  Their conclusion: "pronounced and robust event-to-event differences in the
  tropospheric response to SSWs".

WHY THAT CONCLUSION DOES NOT FOLLOW FROM THAT CORRELATION, AND WHAT IS TESTED
  A correlation of ensemble means across events is SCALE-FREE. r = 0.85 says the
  forced lower-stratospheric and forced surface responses move TOGETHER across
  events. It says nothing about how far apart those events are, which is the
  quantity `snapsi_distribution_test.py` bounds (variance ratio 0.939
  [0.781, 1.130], sigma_f^2 <= +0.13).

  So both can be true at once, and the question this script asks is which:
  is the 100 hPa -> surface relation a property that DISTINGUISHES events, or a
  general lower-stratosphere-to-surface coupling that operates just as strongly
  BETWEEN MEMBERS of a single event, where there is no event-to-event difference
  to explain it?

  In a SNAPSI `nudged` ensemble every member's stratosphere is nudged to the
  same observed evolution, so the SSW is identical BY CONSTRUCTION. Members
  differ only in tropospheric noise and in whatever the nudging does not fix.

WHY 100 hPa IS FREE, WHICH IS WHAT MAKES THIS POSSIBLE
  The protocol (Hitchcock et al., GMD 15, 5073, 2022) nudges only the ZONAL-MEAN
  temperature and zonal wind, at full strength above 50 hPa, tapering to NO
  nudging below 90 hPa. Eddies are free at every level. So 100 hPa polar-cap zg
  is not directly constrained -- but that is an argument, not a measurement, and
  the DESIGN GATE below measures it before anything else is computed.

THE DESIGN GATE, PRE-SPECIFIED
  If nudging above 50 hPa effectively pins the 100 hPa level, members do not
  differ there and the test is empty. Proceed only if

      sd(week-2 100 hPa zg across NUDGED members)
      ---------------------------------------------  >  0.5
      sd(week-2 100 hPa zg across CONTROL members)

  This threshold was fixed before the data were examined. If it fails, the
  reported result is that SNAPSI cannot address the question.

WINDOWS, AND THE COVERAGE CONSTRAINT THAT FORCES THREE OF THEM
  Post-onset coverage differs by initialisation AND by centre, because run
  lengths differ, so the exact published response window -- weeks 3-7, days
  15-49 -- is covered in full by only some ensembles. Three windows are
  reported, all declared here rather than chosen after seeing an answer:
    primary        days 15-25, the widest every NH initialisation covers
    fidelity       days 15-49, Loeffel's own window, where coverage allows
    max_available  days 15 to the last day the shortest member supports,
                   resolved from the data at run time
  The convention follows snapsi_causal_effect.py: "each initialisation's own
  maximum window is reported alongside as a sensitivity, so the choice is
  visible rather than load-bearing". An ensemble is only admitted to a window
  if its members actually span it; a partial ensemble is dropped, not averaged
  over a shorter stretch, because that would silently change the estimand.

SURFACE VARIABLE
  `psl` is used where Loeffel use 1000 hPa GPH. For a CORRELATION this is
  immaterial: over a polar cap the two are related by the hydrostatic relation
  to very good approximation, and a Pearson correlation is invariant under any
  affine transform of either variable. It also keeps this consistent with
  results K, L and M, which are all computed on psl.

WHAT WOULD SHOW THE APPROACH IS WRONG
  - the design gate fails -> SNAPSI cannot address this, reported as such;
  - the CONTROL arm (no SSW at all) reproduces the nudged correlation at full
    strength -> the relation is pure window overlap and adds nothing beyond the
    baseline `loeffel2025_missing_control.py` already established;
  - the relation fails to appear anywhere -> the windows or the standardisation
    are wrong here, not the published result.

SAMPLE SIZE, STATED UP FRONT
  SNAPSI has TWO NH events at two initialisations each. This script CANNOT test
  Loeffel's between-event correlation -- n=2 is not a sample. It tests the
  INFERENCE drawn from r = 0.85, not the coefficient.

Output: results/current/8_experiment/snapsi_loeffel_test.json
"""
from __future__ import annotations

import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "8_experiment"
ING = ROOT / "ssw-design-analysis" / "03_data_ingestion"
PSL = ING / "snapsi_polarcap_psl.parquet"
ZG = ING / "snapsi_polarcap_zg100.parquet"

# WMO 10 hPa 60N wind-reversal central dates for the two NH SNAPSI events.
ONSET = {"s20180125": "2018-02-12", "s20180208": "2018-02-12",
         "s20181213": "2019-01-02", "s20190108": "2019-01-02"}
EVENT = {"s20180125": "Feb-2018", "s20180208": "Feb-2018",
         "s20181213": "Jan-2019", "s20190108": "Jan-2019"}

PRED_WIN = (8, 14)          # Loeffel week 2
RESP_PRIMARY = (15, 25)     # widest window every NH initialisation covers
RESP_FIDELITY = (15, 49)    # Loeffel weeks 3-7, where coverage allows
GATE_MIN_RATIO = 0.5        # pre-specified, fixed before looking at the data
N_PERM = 5000
SEED = zlib.crc32(b"snapsi_loeffel_test") % (2 ** 32)


def post_onset(d):
    off = {k: (pd.Timestamp(ONSET[k]) - pd.Timestamp(k[1:])).days for k in ONSET}
    d = d[d.init.isin(ONSET)].copy()
    d["post_onset"] = d["lead_days"] - d["init"].map(off)
    return d


def win_mean(d, col, lo, hi):
    """Per-member mean of `col` over post-onset days [lo, hi]; NaN if not covered."""
    m = d[(d.post_onset >= lo) & (d.post_onset <= hi)]
    if m.empty:
        return pd.Series(dtype=float)
    covered = m.groupby(["centre", "experiment", "init", "member"]).post_onset.agg(
        ["min", "max"])
    ok = covered[(covered["min"] <= lo + 1.0) & (covered["max"] >= hi - 1.0)].index
    g = m.groupby(["centre", "experiment", "init", "member"])[col].mean()
    return g.reindex(ok)


def duplicate_guard(psl):
    """Refuse a centre whose nudged and control ensembles are the same data.

    NRL/NAVGEM's submission is corrupt: paired member differences are exactly
    zero at three of four initialisations. Tests the EFFECT, not the labels --
    an earlier version compared members by name and passed it.
    """
    bad = []
    for c, s in psl.groupby("centre"):
        eff = []
        for init, si in s.groupby("init"):
            n = si[si.experiment == "nudged"].groupby("lead_days").psl_cap.mean()
            k = si[si.experiment == "control"].groupby("lead_days").psl_cap.mean()
            j = n.index.intersection(k.index)
            if len(j):
                eff.append(abs(float((n[j] - k[j]).mean())))
        if eff and max(eff) < 1.0:          # Pa; real effects are O(100 Pa)
            bad.append(c)
    return bad


def fisher_pool(rs, ns):
    """Fisher-z pooled correlation with a 95% interval."""
    rs, ns = np.asarray(rs, float), np.asarray(ns, float)
    keep = np.isfinite(rs) & (ns > 3)
    if not keep.any():
        return None
    z = np.arctanh(np.clip(rs[keep], -0.999999, 0.999999))
    w = ns[keep] - 3
    zbar = float(np.sum(w * z) / np.sum(w))
    se = float(np.sqrt(1.0 / np.sum(w)))
    return {"r": float(np.tanh(zbar)),
            "CI95": [float(np.tanh(zbar - 1.96 * se)),
                     float(np.tanh(zbar + 1.96 * se))],
            "n_ensembles": int(keep.sum()),
            "n_members_total": int(ns[keep].sum())}


def perm_p(x, y, r_obs, rng, n=N_PERM):
    """P(|r| >= |r_obs|) when members are paired at random within the ensemble."""
    cnt = 0
    for _ in range(n):
        if abs(np.corrcoef(x, rng.permutation(y))[0, 1]) >= abs(r_obs):
            cnt += 1
    return (cnt + 1) / (n + 1)


def main() -> int:
    if not ZG.exists():
        sys.exit(f"missing {ZG.name} -- run 03_data_ingestion/acquire_snapsi_zg.py")
    zg = post_onset(pd.read_parquet(ZG))
    psl = post_onset(pd.read_parquet(PSL))

    dup = duplicate_guard(psl)
    if dup:
        print(f"EXCLUDED by duplicate guard (corrupt submission): {dup}")
        zg, psl = zg[~zg.centre.isin(dup)], psl[~psl.centre.isin(dup)]

    # Per CENTRE and initialisation. Aggregating over centres reports the
    # longest run as if every centre reached it, which is how the fidelity
    # window can look available when no single ensemble covers it.
    print("\npost-onset coverage per centre x initialisation (days), from the data:")
    cov = (psl.groupby(["centre", "init"]).post_onset.agg(["min", "max"])
           .round(1).reset_index())
    cov["covers_15_49"] = (cov["min"] <= 16) & (cov["max"] >= 48)
    print(cov.to_string(index=False))
    n_fid = int(cov["covers_15_49"].sum())
    print(f"  ensembles covering Loeffel's exact weeks 3-7 window: {n_fid}"
          f" of {len(cov)}")
    print(f"\ncentres with zg: {sorted(zg.centre.unique())}")
    print(f"centres with psl: {sorted(psl.centre.unique())}")

    rng = np.random.default_rng(SEED)
    out = {
        "paper": {"citation": "Loeffel et al. 2026, WCD 7, 895-913",
                  "reported_r": 0.85, "n_events": 18,
                  "computed": "across events, on ensemble means"},
        "windows": {"predictor_days": list(PRED_WIN),
                    "response_primary_days": list(RESP_PRIMARY),
                    "response_fidelity_days": list(RESP_FIDELITY)},
        "seed": SEED, "n_perm": N_PERM,
        "excluded_centres": dup,
    }

    # ---------------------------------------------------------- design gate
    pred = win_mean(zg, "zg100_cap", *PRED_WIN)
    gate_rows = []
    for (centre, init), _ in pred.groupby(level=[0, 2]):
        try:
            sd_n = float(pred.loc[(centre, "nudged", init)].std(ddof=1))
            sd_k = float(pred.loc[(centre, "control", init)].std(ddof=1))
        except KeyError:
            continue
        gate_rows.append({"centre": centre, "init": init,
                          "sd_nudged_m": round(sd_n, 2),
                          "sd_control_m": round(sd_k, 2),
                          "ratio": round(sd_n / sd_k, 3) if sd_k else None})
    out["design_gate"] = {"threshold_ratio": GATE_MIN_RATIO, "per_ensemble": gate_rows}
    ratios = [g["ratio"] for g in gate_rows if g["ratio"] is not None]
    out["design_gate"]["median_ratio"] = round(float(np.median(ratios)), 3) if ratios else None
    out["design_gate"]["passes"] = bool(ratios and np.median(ratios) > GATE_MIN_RATIO)

    print("\n=== DESIGN GATE: is 100 hPa free enough for the test to mean anything? ===")
    for g in gate_rows:
        print(f"  {g['centre']:<12} {g['init']:<11} sd nudged {g['sd_nudged_m']:7.2f} m   "
              f"control {g['sd_control_m']:7.2f} m   ratio {g['ratio']}")
    print(f"  median ratio {out['design_gate']['median_ratio']} "
          f"vs threshold {GATE_MIN_RATIO} -> "
          f"{'PASS' if out['design_gate']['passes'] else 'FAIL'}")
    if not out["design_gate"]["passes"]:
        print("\n  Nudging pins the 100 hPa level: members do not differ there, so this"
              "\n  experiment cannot separate an event property from a general coupling."
              "\n  Reported as a null result about the METHOD, not about the physics.")
        RESULTS.mkdir(parents=True, exist_ok=True)
        (RESULTS / "snapsi_loeffel_test.json").write_text(
            json.dumps(out, indent=1), encoding="utf8")
        return 0

    # --------------------------------------------------------- the test
    # RESP_MAXAVAIL is the widest response window the shortest-running member of
    # an ensemble supports. It is declared here, not chosen after seeing a
    # result, and follows the convention snapsi_causal_effect.py already uses:
    # "each initialisation's own maximum window is reported alongside as a
    # sensitivity, so the choice is visible rather than load-bearing".
    for tag, win in (("primary", RESP_PRIMARY),
                     ("fidelity", RESP_FIDELITY),
                     ("max_available", None)):
        if win is None:
            hi = float(psl.groupby(
                ["centre", "experiment", "init", "member"]).post_onset.max().min())
            win = (RESP_PRIMARY[0], int(np.floor(hi)))
            out["windows"]["response_max_available_days"] = list(win)
            print(f"\n(max_available response window resolved to days "
                  f"{win[0]}-{win[1]} from the shortest member)")
        resp = win_mean(psl, "psl_cap", *win)
        per = []
        for arm in ("nudged", "control"):
            rs, ns = [], []
            for (centre, init) in sorted({(c, i) for c, e, i, _ in pred.index if e == arm}):
                try:
                    x = pred.loc[(centre, arm, init)]
                    y = resp.loc[(centre, arm, init)]
                except KeyError:
                    continue
                j = x.index.intersection(y.index)
                if len(j) < 10:
                    continue
                xv, yv = x[j].values, y[j].values
                r = float(np.corrcoef(xv, yv)[0, 1])
                p = perm_p(xv, yv, r, rng)
                per.append({"centre": centre, "init": init, "event": EVENT[init],
                            "arm": arm, "n_members": int(len(j)),
                            "r": round(r, 4), "perm_p": round(p, 5)})
                rs.append(r); ns.append(len(j))
            out.setdefault(tag, {})[arm] = {"per_ensemble": [q for q in per if q["arm"] == arm],
                                            "pooled": fisher_pool(rs, ns)}
        print(f"\n=== {tag.upper()} response window: post-onset days {win[0]}-{win[1]} ===")
        for arm in ("nudged", "control"):
            pooled = out[tag][arm]["pooled"]
            print(f"  {arm:<8} " + (
                f"pooled r = {pooled['r']:+.3f} "
                f"[{pooled['CI95'][0]:+.3f}, {pooled['CI95'][1]:+.3f}]  "
                f"({pooled['n_ensembles']} ensembles, {pooled['n_members_total']} members)"
                if pooled else "no ensemble with enough members"))
            for q in out[tag][arm]["per_ensemble"]:
                print(f"      {q['centre']:<12} {q['init']:<11} {q['event']:<9} "
                      f"n={q['n_members']:<3} r={q['r']:+.3f}  perm p={q['perm_p']:.4f}")

    # ---------------------------------------------------- power, as promised
    # What within-ensemble correlation could this design have detected? Fisher-z,
    # two-sided alpha = 0.05, at the median ensemble size actually used.
    from scipy import stats as _st
    ns_used = [q["n_members"] for q in out["primary"]["nudged"]["per_ensemble"]]
    n_med = int(np.median(ns_used)) if ns_used else 0
    power = {}
    if n_med > 3:
        zc = _st.norm.ppf(0.975)
        se = 1.0 / np.sqrt(n_med - 3)
        for r_true in (0.1, 0.2, 0.3, 0.4, 0.5, 0.85):
            zt = np.arctanh(r_true)
            power[str(r_true)] = round(float(
                _st.norm.sf(zc - zt / se) + _st.norm.cdf(-zc - zt / se)), 4)
    out["power_within_ensemble"] = {
        "n_members_median": n_med,
        "detection_prob_by_true_r": power,
        "note": ("Power to detect a WITHIN-ensemble correlation at n members. "
                 "The 0.85 row is included only for scale: it is Loeffel's "
                 "BETWEEN-EVENT value on ensemble means and is not the same "
                 "estimand as anything measured here.")}
    if power:
        print(f"\n=== POWER (within-ensemble, n={n_med}, two-sided alpha=0.05) ===")
        for k, v in power.items():
            print(f"  true r = {k:<5} detected with probability {v:.3f}")

    # ------------------------------ the closest SNAPSI analogue of THEIR estimator
    # Loeffel correlate ENSEMBLE MEANS across 18 events. SNAPSI has two NH events,
    # so the between-EVENT correlation cannot be estimated here -- n=2. What can be
    # computed is the correlation across ENSEMBLES (centre x initialisation) of the
    # same two quantities. It is reported as descriptive only: its spread is
    # dominated by differences BETWEEN MODELS, not between events, so it is NOT a
    # test of their claim and must never be quoted as one.
    resp_p = win_mean(psl, "psl_cap", *RESP_PRIMARY)
    ens = []
    for (centre, arm, init) in sorted({k[:3] for k in pred.index}):
        try:
            x = float(pred.loc[(centre, arm, init)].mean())
            y = float(resp_p.loc[(centre, arm, init)].mean())
        except KeyError:
            continue
        ens.append({"centre": centre, "arm": arm, "init": init,
                    "event": EVENT[init], "zg100_week2_m": round(x, 2),
                    "psl_resp_Pa": round(y, 1)})
    out["between_ensemble_descriptive"] = {"points": ens}
    for arm in ("nudged", "control"):
        pts = [e for e in ens if e["arm"] == arm]
        n_ev = len({e["event"] for e in pts})
        if len(pts) >= 4:
            xv = np.array([e["zg100_week2_m"] for e in pts])
            yv = np.array([e["psl_resp_Pa"] for e in pts])
            r = float(np.corrcoef(xv, yv)[0, 1])
            out["between_ensemble_descriptive"][arm] = {
                "n_ensembles": len(pts), "n_distinct_events": n_ev, "r": round(r, 4)}
            print(f"\nbetween-ensemble (DESCRIPTIVE, {arm}): r = {r:+.3f} over "
                  f"{len(pts)} ensembles spanning only {n_ev} distinct event(s)")
            print("  Not a test of Loeffel: the spread is between MODELS, not events.")
        else:
            out["between_ensemble_descriptive"][arm] = {
                "n_ensembles": len(pts), "n_distinct_events": n_ev,
                "r": None, "note": "too few ensembles"}

    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "snapsi_loeffel_test.json").write_text(
        json.dumps(out, indent=1), encoding="utf8")
    print(f"\nSaved -> snapsi_loeffel_test.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
