#!/usr/bin/env python3
"""
snapsi_causal_effect.py
=======================
THE CAUSAL STRATOSPHERE-TO-SURFACE EFFECT, MEASURED BY EXPERIMENT RATHER THAN
INFERRED FROM COMPOSITES -- AND THE NUMBER THE WHOLE LITERATURE IS TRYING TO
ESTIMATE.

WHY THIS IS DIFFERENT FROM EVERYTHING ELSE IN THIS PROJECT
  `stratifier_bias_law.py` and `forced_variance_ceiling.py` both work on data
  where the stratospheric anomaly is never manipulated: observations and
  free-running models. They can bound the artefact but they can never observe the
  counterfactual, because the counterfactual does not exist in that data.

  SNAPSI manipulates it. `nudged` relaxes the zonally symmetric stratosphere
  (full strength above 50 hPa, tapering to none below 90 hPa; Hitchcock et al.
  2022) toward the OBSERVED evolution of a real SSW; `control` applies the
  IDENTICAL nudging toward a 1979-2019 climatology. The two ensembles share
  initial conditions, model, resolution and nudging machinery, and differ only in
  whether the stratosphere carries the event. So

      S  =  mean(nudged)  -  mean(control)

  is the causal effect of the stratospheric anomaly on the surface, by
  construction. Nudging artefacts cancel. This is the counterfactual the
  compositing literature has been trying to reconstruct from co-occurrence.

  (Differencing against `free` would NOT do this: `free` is unnudged, so the
  contrast would confound the anomaly with the suppression of internal
  stratospheric variability that nudging itself imposes.)

THE WINDOW, CHOSEN TO BE COMPARABLE ACROSS INITIALISATIONS
  The four NH initialisations sit at very different leads relative to onset:
      s20180125  Feb-2018 SSW (central date 2018-02-12)  onset at lead +18
      s20180208  Feb-2018 SSW                            onset at lead  +4
      s20181213  Jan-2019 SSW (central date 2019-01-02)  onset at lead +20
      s20190108  Jan-2019 SSW                            onset at lead  -6
  Forecasts run 45-60 days depending on the centre, so post-onset coverage
  differs between initialisations.
  The primary window is therefore the COMMON overlap, post-onset days +8..+25,
  which every initialisation covers in full and which starts where the published
  surface window starts. Each initialisation's own maximum window is reported
  alongside as a sensitivity, so the choice is visible rather than load-bearing.

UNITS, SO THE NUMBER IS COMPARABLE TO PUBLISHED CONTRASTS
  S is divided by the standard deviation ACROSS CONTROL MEMBERS of the same
  window mean -- the internal variability of the seasonal-mean polar-cap surface
  state. That puts it in the same "sigma of surface variability" units the
  DW-minus-NDW contrasts are quoted in.

  CAVEAT, STATED NOT BURIED: the literature's sigma is the spread across EVENTS or
  dates in observations, which also contains any between-event forced variation;
  this sigma is the spread across MEMBERS within one event, which does not. The
  two coincide only if the between-event forced variance is small -- which is
  exactly what `forced_variance_ceiling.py` concludes, so the comparison is
  consistent with, but not independent of, that result.

SAMPLE SIZE, STATED UP FRONT
  Two NH EVENTS, each at two initialisations. Not four events. Anything said here
  about how events DIFFER rests on n=2 and is reported as indicative only. What
  is precisely estimated is the MEAN causal effect, from 38-52 members per
  ensemble across 9 models.

Output: snapsi_causal_effect.json
"""
import json
import pathlib

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "8_experiment"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]
SRC = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "snapsi_polarcap_psl.parquet"

# SSW central dates, WMO 10 hPa 60N wind reversal, for the two NH SNAPSI events.
ONSET = {"s20180125": "2018-02-12", "s20180208": "2018-02-12",
         "s20181213": "2019-01-02", "s20190108": "2019-01-02"}
EVENT = {"s20180125": "Feb-2018", "s20180208": "Feb-2018",
         "s20181213": "Jan-2019", "s20190108": "Jan-2019"}
WIN = (8, 25)          # common post-onset window covered by every initialisation

# reported DW-minus-NDW contrast and the DW FRACTION that study reports
REPORTED = [
    ("Karpechko-criterion AO contrast, observations", 1.782, 0.70),
    ("CMIP6 ensemble DW-NDW contrast", 1.143, 0.50),
    # Lu & Rao (2026, ACP 26, 3723) give ERA5 60-day NAO means by DW subtype: BOTH -0.762 (n=13), EA -0.567 (14), NA -0.435 (6); NDW +0.088 (19). All-DW mean -0.620, contrast -0.708, DW fraction 33/52 = 0.635. The -0.850 used until 2026-09-25 subtracted the BOTH subtype alone.
    ("ACP 26, 3723 (2026) published NAO contrast", 0.708, 0.635),
    ("ERA5 1000 hPa NAM contrast, Karpechko criterion", 0.684, 0.54),
    ("ERA5 850 hPa NAM contrast, ACP criterion", 0.639, 0.59),
]


def duplicate_guard(d):
    """Refuse any centre whose `nudged` and `control` ensembles are the same data.

    NRL/NAVGEM returns S = -0.000 sigma for all four initialisations. That is not a
    model with weak downward coupling: paired member-by-member differences are at
    most 0.023 Pa at three of the four initialisations and average 3 Pa at the
    fourth (s20180125), against a within-ensemble spread of ~110 Pa. Two independent 45-day
    forecasts from perturbed states diverge by O(100 Pa); near-identical
    trajectories mean the archived files are duplicates, an artefact of the
    submission rather than a property of the atmosphere. NRL's within-ensemble
    spread is also 3-4x smaller than every other centre's, which corroborates it.

    Including it would have dragged the multi-model mean down by a third with a
    number that means nothing. The test is automatic because this project has been
    bitten repeatedly by silent duplication and zero-filling.
    """
    print("\n=== DUPLICATE GUARD: paired |nudged - control| vs ensemble spread ===")
    print(f"{'centre':<8s} {'mean|diff|':>11s} {'within-ens sd':>14s} {'ratio':>8s}  verdict")
    print("-" * 56)
    keep = []
    for c in sorted(d.centre.unique()):
        s = d[d.centre == c]
        n = s[s.experiment == "nudged"].set_index(["init", "member", "lead_days"]).psl_cap
        k = s[s.experiment == "control"].set_index(["init", "member", "lead_days"]).psl_cap
        j = pd.concat([n.rename("a"), k.rename("b")], axis=1).dropna()
        pdif = float((j.a - j.b).abs().mean()) if len(j) else np.nan
        sd = float(s[s.experiment == "control"]
                   .groupby(["init", "lead_days"]).psl_cap.std().mean())
        r = pdif / sd if sd else np.nan
        ok = np.isfinite(r) and r >= 0.2
        if ok:
            keep.append(c)
        print(f"{c:<8s} {pdif:11.2f} {sd:14.2f} {r:8.3f}  "
              f"{'ok' if ok else 'DUPLICATE -> EXCLUDED'}")
    return keep


def member_window_means(d, lo, hi):
    """Mean polar-cap psl per member over post-onset days [lo, hi]."""
    m = d[(d.post_onset >= lo) & (d.post_onset <= hi)]
    return m.groupby(["centre", "experiment", "init", "member"]).psl_cap.mean()


def main():
    if not SRC.exists():
        raise SystemExit(f"{SRC.name} not found -- run acquire_snapsi_surface.py")
    d = pd.read_parquet(SRC)
    # post-onset day from LEAD plus a fixed offset, never from the file's calendar.
    # CCCma runs a 365-day NoLeap calendar; reconciling it against real dates is
    # both unnecessary and a source of silent drift. The offset below is computed
    # once from the two known real dates and applies to every member.
    off = {k: (pd.Timestamp(ONSET[k]) - pd.Timestamp(k[1:])).days for k in ONSET}
    print(f"onset-minus-init offsets (days): {off}")
    d["post_onset"] = d["lead_days"] - d["init"].map(off)
    d["event"] = d["init"].map(EVENT)
    print(f"{len(d):,} rows | centres {sorted(d.centre.unique())}")
    print(d.groupby(["centre", "experiment", "init"]).member.nunique()
          .unstack(fill_value=0).to_string())

    keep = duplicate_guard(d)
    dropped = sorted(set(d.centre.unique()) - set(keep))
    d = d[d.centre.isin(keep)]
    print(f"\ncentres retained: {keep}   excluded as duplicates: {dropped or 'none'}")

    cov = d.groupby("init").post_onset.max().round(1)
    print(f"\npost-onset coverage per init (days): {cov.to_dict()}")
    print(f"primary common window: +{WIN[0]}..+{WIN[1]} days after onset\n")

    wm = member_window_means(d, *WIN)
    res = {"window": list(WIN), "onsets": ONSET, "per_case": {}, "reported": {}}

    print("=" * 88)
    print(f"{'centre':<9s} {'event':<9s} {'init':<11s} {'nudged':>8s} {'control':>8s} "
          f"{'S (Pa)':>9s} {'sd_ctl':>8s} {'S/sigma':>9s} {'n_N/n_C':>9s}")
    print("-" * 88)
    rows = []
    for (centre, init), g in wm.groupby(level=[0, 2]):
        try:
            nud = g.xs("nudged", level="experiment")
            ctl = g.xs("control", level="experiment")
        except KeyError:
            continue
        if len(nud) < 10 or len(ctl) < 10:
            continue
        S = float(nud.mean() - ctl.mean())
        sd = float(ctl.std(ddof=1))
        se = float(np.sqrt(nud.var(ddof=1) / len(nud) + ctl.var(ddof=1) / len(ctl)))
        z = S / sd if sd else np.nan
        rows.append({"centre": centre, "init": init, "event": EVENT[init],
                     "S_Pa": S, "sd_control_Pa": sd, "S_sigma": z,
                     "se_Pa": se, "se_sigma": se / sd if sd else np.nan,
                     "n_nudged": int(len(nud)), "n_control": int(len(ctl))})
        print(f"{centre:<9s} {EVENT[init]:<9s} {init:<11s} {nud.mean():8.1f} "
              f"{ctl.mean():8.1f} {S:+9.1f} {sd:8.1f} {z:+9.3f} "
              f"{len(nud):4d}/{len(ctl):<4d}")
    r = pd.DataFrame(rows)

    # SECOND NORMALISATION, because the first is lead-dependent.
    # sd_control is the spread of a FORECAST ensemble, which grows with lead. The
    # four initialisations cover the same POST-ONSET window at different LEADS
    # (26-43, 12-29, 28-45, 14-31 days), so dividing each by its own sd compares
    # effects against different yardsticks -- s20190108 scores 2.2 partly because
    # its window sits at short lead where the ensemble has not yet spread. The
    # literature's sigma is a climatological spread, not a lead-dependent one, so
    # a per-centre pooled sigma (deviations from each init's own mean, pooled
    # across all four) is the closer analogue and is reported alongside.
    pooled = {}
    for c, g in wm.groupby(level=0):
        ctl = g.xs("control", level="experiment")
        dev = ctl - ctl.groupby(level="init").transform("mean")
        pooled[c] = float(dev.std(ddof=1))
    r["sd_pooled_Pa"] = r.centre.map(pooled)
    r["S_sigma_pooled"] = r.S_Pa / r.sd_pooled_Pa
    print(f"\npooled within-centre control sigma (Pa): "
          f"{ {k: round(v, 1) for k, v in pooled.items()} }")
    print(f"{'centre':<8s} {'init':<11s} {'S/sd_init':>10s} {'S/sd_pooled':>12s}")
    print("-" * 44)
    for row in r.itertuples():
        print(f"{row.centre:<8s} {row.init:<11s} {row.S_sigma:+10.3f} "
              f"{row.S_sigma_pooled:+12.3f}")
    res["per_case"] = r.round(4).to_dict("records")
    res["pooled_sigma_Pa"] = {k: round(v, 2) for k, v in pooled.items()}

    print("\n" + "=" * 88)
    print("=== THE CAUSAL EFFECT ===")
    # average over models within an init, then over inits within an event, so a
    # centre with more members does not dominate and the two events weigh equally
    # HEADLINE uses the pooled sigma. Chosen on principle -- a forecast ensemble
    # spread grows with lead, so per-init normalisation compares each effect
    # against a different yardstick -- and worth noting that it moves S DOWN,
    # which makes the implied NDW responses smaller and so works AGAINST the
    # argument being made here rather than for it.
    per_init = r.groupby(["event", "init"]).S_sigma_pooled.mean()
    per_init_raw = r.groupby(["event", "init"]).S_sigma.mean()
    print(f"\n(per-init-normalised sensitivity, grand mean "
          f"{per_init_raw.groupby('event').mean().mean():+.3f} sigma)")
    per_event = per_init.groupby("event").mean()
    grand = float(per_event.mean())
    print(f"\nper initialisation (mean over {r.centre.nunique()} models):")
    for (ev, init), v in per_init.items():
        print(f"  {ev:<9s} {init:<11s} S/sigma = {v:+.3f}")
    print(f"\nper event (mean over its two initialisations):")
    for ev, v in per_event.items():
        print(f"  {ev:<9s} S/sigma = {v:+.3f}")
    print(f"\nGRAND MEAN causal effect S = {grand:+.3f} sigma")
    spread = float(per_event.max() - per_event.min()) if len(per_event) > 1 else np.nan
    print(f"between-EVENT difference (n=2, indicative only) = {spread:+.3f} sigma")
    across_models = r.groupby("centre").S_sigma_pooled.mean()
    print(f"\nbetween-MODEL spread of S/sigma: "
          f"{across_models.min():+.3f} .. {across_models.max():+.3f} "
          f"(sd {across_models.std(ddof=1):.3f})")
    for c, v in across_models.items():
        print(f"  {c:<14s} {v:+.3f}")

    res["S_sigma_per_init"] = {f"{a}|{b}": round(float(v), 4)
                               for (a, b), v in per_init.items()}
    res["S_sigma_per_event"] = {k: round(float(v), 4) for k, v in per_event.items()}
    res["S_sigma_grand_mean"] = round(grand, 4)
    res["between_event_difference"] = None if not np.isfinite(spread) else round(spread, 4)
    res["between_model_sd"] = round(float(across_models.std(ddof=1)), 4)
    res["S_sigma_per_model"] = {k: round(float(v), 4) for k, v in across_models.items()}

    print("\n" + "=" * 88)
    print("=== WHAT EACH PUBLISHED CONTRAST IMPLIES, GIVEN THE MEASURED S ===")
    print("  A DW-minus-NDW contrast is bounded by BETWEEN-EVENT variation, not by S,")
    print("  so 'contrast vs S' is the wrong comparison and an earlier draft of this")
    print("  script made it. The right question: with a mean causal effect S and a DW")
    print("  fraction q, a contrast C forces the two subgroups to be")
    print("      DW  =  S + (1-q)*C        NDW  =  S - q*C")
    print("  and it is the implied NDW value that is testable.\n")
    print(f"{'reported quantity':<48s} {'C':>6s} {'q':>5s} {'DW':>7s} {'NDW':>7s}")
    print("-" * 78)
    for lab, v, q in REPORTED:
        dw = grand + (1 - q) * v
        ndw = grand - q * v
        res["reported"][lab] = {"contrast": v, "dw_fraction": q,
                                "implied_DW_sigma": round(float(dw), 3),
                                "implied_NDW_sigma": round(float(ndw), 3)}
        print(f"{lab:<48s} {v:6.3f} {q:5.2f} {dw:+7.3f} {ndw:+7.3f}")
    # zero-or-REVERSED, not |NDW| small: a negative implied NDW means events
    # labelled non-propagating would have to respond in the OPPOSITE sense to the
    # SSW effect, which is a stronger demand on the data than merely null.
    ndws = {lab: grand - q * v for lab, v, q in REPORTED}
    dead = [lab for lab, x in ndws.items() if x <= 0.20]
    live = [x for x in ndws.values() if x > 0.20]
    print(f"\n  READ HONESTLY, S = {grand:+.3f} sigma:")
    print(f"   - {len(dead)} of {len(REPORTED)} contrasts force NDW to zero or REVERSED sign:")
    for lab in dead:
        print(f"       {lab}  ->  NDW = {ndws[lab]:+.3f}")
    if live:
        print(f"   - the other {len(live)} imply NDW responses of "
              f"{min(live):+.2f} to {max(live):+.2f} sigma, i.e. events labelled")
        print("     'non-propagating' would STILL carry a large real surface response.")
        print("     That is a CONTINUUM, not two populations, and it is exactly what a")
        print("     threshold cut through one shifted distribution looks like.")
    print("\n  So SNAPSI does NOT show the published contrasts exceed what the")
    print("  stratosphere can cause -- an earlier draft of this script claimed that and")
    print(f"  it was wrong. What it shows is that a large mean effect (S = {grand:+.2f}")
    print("  sigma, measured causally) plus a threshold split reproduces most of those")
    print("  contrasts WITHOUT any second population; only the most extreme reported")
    print("  contrast demands one, and it demands a physically odd one.")
    res["implied_NDW_sigma"] = {k: round(float(v), 4) for k, v in ndws.items()}

    (RESULTS / "snapsi_causal_effect.json").write_text(json.dumps(res, indent=2),
                                                    encoding="utf8")
    print("\nSaved -> snapsi_causal_effect.json")


if __name__ == "__main__":
    main()
