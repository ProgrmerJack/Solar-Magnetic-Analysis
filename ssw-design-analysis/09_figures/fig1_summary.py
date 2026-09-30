#!/usr/bin/env python3
"""
fig1_summary.py -- Figure 1: one shifted population, cut by a threshold.

Reads (recomputes nothing but a histogram):
  results/current/8_experiment/snapsi_distribution_test.json   member outcomes (M)
  results/current/8_experiment/snapsi_selection_test.json       DW rates (L)
  results/current/*/era5_recompute_and_two_thirds.json          "two thirds" (H)
  results/current/8_experiment/snapsi_archetype_test.json       archetype pair (O)
  results/current/8_experiment/snapsi_regional_test.json        cold risk (R)
  results/current/2_event_study/era5_regional_test.json         observed cold risk (R)

  a  SNAPSI, all 36 NH ensembles: control members (grey) and SSW-imposed members,
     each re-centred on its own ensemble mean and placed at the mean shift
     (vermilion); the threshold of the downward criterion at zero; the share of
     members below it in each arm
  b  observations: share of the observed SSWs (n from the result) meeting the surface conditions, against
     event-free dates as they are and displaced by the measured shift
  c  the archetype pair: forced probability of a downward outcome in each model
     (nudged members) for February 2018 and January 2019, primary initialisations
  d  probability of a cold northern-Eurasian fortnight: SNAPSI without and with
     the imposed SSW and as predicted from the circulation shift; ERA5 on
     event-free dates and after 39 observed SSWs (95% interval)
"""
import glob
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nature_style as S  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402


def main():
    S.apply()
    R = S.RESULTS
    dist = json.loads((R / "8_experiment" / "snapsi_distribution_test.json").read_text())
    sel = json.loads((R / "8_experiment" / "snapsi_selection_test.json").read_text())
    two = json.loads(Path(glob.glob(str(R / "*" / "era5_recompute_and_two_thirds.json"))[0]).read_text())
    arc = json.loads((R / "8_experiment" / "snapsi_archetype_test.json").read_text())
    reg = json.loads((R / "8_experiment" / "snapsi_regional_test.json").read_text())
    ere = json.loads((R / "2_event_study" / "era5_regional_test.json").read_text())

    fig = plt.figure(figsize=(S.DOUBLE, 118 * S.MM))
    gs = fig.add_gridspec(2, 2, wspace=0.32, hspace=0.55)
    grey, verm = "0.6", S.C["vermillion"]

    # a -- one population, shifted, cut at zero
    ax = fig.add_subplot(gs[0, 0])
    MA = dist["member_A"]
    ctl = np.concatenate([np.array(v["control"]) - np.mean(v["control"]) for v in MA.values()])
    nud = np.concatenate([np.array(v["nudged"]) - np.mean(v["nudged"]) for v in MA.values()])
    shift = float(np.mean([np.mean(v["nudged"]) - np.mean(v["control"]) for v in MA.values()]))
    bins = np.linspace(-5, 4, 55)
    ax.hist(ctl, bins=bins, density=True, color=grey, alpha=0.55, label="no SSW (control)")
    ax.hist(nud + shift, bins=bins, density=True, histtype="step", color=verm, lw=1.2,
            label="SSW imposed (nudged)")
    ax.axvline(0, color="k", lw=0.8)
    nh = [e for e in sel["all_ensembles"] if e["hemisphere"] == "NH"]
    rate = {arm: np.mean([e["DW_rate"] for e in nh if e["arm"] == arm]) for arm in ("nudged", "control")}
    ax.text(-4.9, 0.49, f"downward:\n{rate['nudged']:.0%} with the SSW\n{rate['control']:.0%} without",
            fontsize=5.6, va="top", color="k")
    vr = dist["pooled_all"]["variance_ratio"]; ci = dist["pooled_all"]["variance_ratio_CI95"]
    ax.text(2.0, 0.30, f"same shape:\nvariance ratio\n{vr:.2f} [{ci[0]:.2f}, {ci[1]:.2f}]", fontsize=5.6, va="top")
    ax.set_xlabel("surface NAM proxy, days 8-25 (σ of control)"); ax.set_ylabel("density")
    ax.legend(fontsize=5.2, loc="upper right")
    ax.set_title("an imposed SSW shifts one population; a threshold cuts it", fontsize=6)
    S.panel_label(ax, "a", x=-0.14, y=1.04)

    # b -- the field's "two thirds"
    bx = fig.add_subplot(gs[0, 1])
    t = two["two_thirds_test"]
    vals = [t["null_pass_rate_c1c2"], t["shifted_pseudo_pass_rate_c1c2"], t["observed_pass_rate_c1c2"]]
    labs = ["event-free\ndates", "event-free dates\n+ measured shift", f"observed\nSSWs ({two['n_events']})"]
    cols = ["0.85", grey, verm]
    bx.bar(range(3), vals, color=cols, width=0.6)
    lo, hi = t["shifted_pseudo_CI95"]
    bx.errorbar(1, vals[1], yerr=[[vals[1] - lo], [hi - vals[1]]], color="k", lw=0.7, capsize=2)
    for i, v in enumerate(vals):
        bx.text(i - (0.2 if i == 1 else 0), v + 0.02, f"{v:.0%}", ha="center", fontsize=5.8)
    bx.set_xticks(range(3)); bx.set_xticklabels(labs, fontsize=5.6)
    bx.set_ylim(0, 1); bx.set_ylabel("share meeting the surface conditions")
    bx.set_title("'about two thirds' is one shifted population (ERA5)", fontsize=6)
    S.panel_label(bx, "b", x=-0.16, y=1.04)

    # c -- the archetype pair
    cx = fig.add_subplot(gs[1, 0])
    pe = [e for e in arc["per_ensemble"] if e["init"] in ("s20180125", "s20181213")]
    cen = sorted({e["centre"] for e in pe})
    for c in cen:
        a = [e for e in pe if e["centre"] == c and e["init"] == "s20180125"]
        b = [e for e in pe if e["centre"] == c and e["init"] == "s20181213"]
        if a and b:
            cx.plot([0, 1], [a[0]["nudged_DW_rate"], b[0]["nudged_DW_rate"]], color="0.7", lw=0.6, zorder=1)
            cx.scatter([0, 1], [a[0]["nudged_DW_rate"], b[0]["nudged_DW_rate"]], s=8, color=grey, zorder=2)
    for x, init, lab in ((0, "s20180125", "Feb 2018\n(observed: downward)"),
                         (1, "s20181213", "Jan 2019\n(observed: labelled\nnon-downward)")):
        m = np.mean([e["nudged_DW_rate"] for e in pe if e["init"] == init])
        cx.plot([x - 0.18, x + 0.18], [m, m], color=verm, lw=1.6)
        cx.text(x + 0.21, m, f"{m:.2f}", va="center", fontsize=5.6, color=verm)
    cx.set_xticks([0, 1]); cx.set_xticklabels(["Feb 2018\n(observed: downward)",
                                               "Jan 2019\n(labelled non-downward)"], fontsize=5.6)
    cx.set_xlim(-0.5, 1.6); cx.set_ylim(0, 1.05)
    cx.set_ylabel("forced P(downward), per model")
    cx.set_title("the field's archetype pair: similar forced odds, two draws", fontsize=6)
    S.panel_label(cx, "c", x=-0.14, y=1.04)

    # d -- cold risk
    dx = fig.add_subplot(gs[1, 1])
    r, rx = reg["results"]["NEURASIA"], reg["sensitivity_excl_ECCC"]["NEURASIA"]
    e = ere["regions"]["NEURASIA"]
    bars = [("SNAPSI\nno SSW", r["cold_control"]["mean"], "0.85", None),
            ("SNAPSI\nSSW imposed", r["cold_obs"]["mean"], verm, rx["cold_obs"]["mean"]),
            ("predicted\nfrom shift", r["cold_pred"]["mean"], grey, rx["cold_pred"]["mean"]),
            ("ERA5\nevent-free", e["cold_event_free"], "0.85", None),
            (f"ERA5 after\n{ere['n_events']} SSWs", e["cold_obs"], verm, None)]
    for i, (lab, v, col, alt) in enumerate(bars):
        dx.bar(i, v, color=col, width=0.62)
        if alt is not None:
            dx.scatter([i], [alt], marker="_", s=60, color="k", zorder=3)
    lo, hi = e["cold_obs_ci95"]
    dx.errorbar(4, e["cold_obs"], yerr=[[e["cold_obs"] - lo], [hi - e["cold_obs"]]], color="k", lw=0.7, capsize=2)
    dx.set_xticks(range(len(bars))); dx.set_xticklabels([b[0] for b in bars], fontsize=5.4)
    dx.set_ylabel("P(cold fortnight), northern Eurasia")
    dx.set_title("regional cold risk follows the shift", fontsize=6)
    S.panel_label(dx, "d", x=-0.16, y=1.04)
    S.save(fig, "fig1_summary")


if __name__ == "__main__":
    main()
