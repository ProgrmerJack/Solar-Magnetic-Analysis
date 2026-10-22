#!/usr/bin/env python3
"""
fig3_scope.py -- Figure 3: across events and systems, where forcing differs -- the forced
share of class contrasts, and how much events differ.

Reads (recomputes nothing):
  results/current/8_experiment/snapsi_residual_quantiles.json   residual quantiles (Z1)
  results/current/5_mechanism/icon_event_heterogeneity.json     ICON event ensembles (Z3)
  results/current/8_experiment/snapsi_archetype_test.json       archetype pair (O)

  results/current/8_experiment/class_forcing_content.json      forced shares (Z20)
  results/current/8_experiment/enso_extreme_class.json         extreme El Nino (Z22)
  results/current/8_experiment/mjo_class_law.json              MJO crossing under the QBO (Z26)

  a  forced share R of the class contrast (observed, filled; 95% interval where stored) and
     its prediction (open: the law's R_LS where forcing levels are known, else h2), per system
     (residual-quantile panels of earlier versions are Extended Data Fig. 8)
  b  18 ICON SSWs re-run as 40-member ensembles (Loeffel et al. 2026): ensemble-mean
     surface response, days 8-25, against the week-2 100 hPa anomaly (both NAM sign);
     open, the six ensembles that start on their onset day (already easterly at t = 0)
  c  the field's archetype pair: forced probability of a downward outcome in each
     model (nudged members), February 2018 and January 2019, primary initialisations
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nature_style as S  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

M = lambda s: s.replace("-", "−")                     # noqa: E731


def main():
    S.apply()
    R = S.RESULTS
    cf = json.loads((R / "8_experiment" / "class_forcing_content.json").read_text())
    en = json.loads((R / "8_experiment" / "enso_extreme_class.json").read_text())
    mj = json.loads((R / "8_experiment" / "mjo_class_law.json").read_text())["primary"]["M3_forced_share"]
    ic = json.loads((R / "5_mechanism" / "icon_event_heterogeneity.json").read_text())["primary_registered"]
    arc = json.loads((R / "8_experiment" / "snapsi_archetype_test.json").read_text())
    grey, verm = "0.6", S.C["vermillion"]
    fig = plt.figure(figsize=(S.DOUBLE, 140 * S.MM))
    gs = fig.add_gridspec(2, 2, wspace=0.34, hspace=0.45, height_ratios=[1.35, 1])
    rows = []          # (label, observed R, ci or None, predicted, kind of prediction)
    p2 = cf["P2_obs_ncep_episodes"]; s2 = cf["S2_obs_era5_episodes"]; s1 = cf["S1_free_natural_experiment"]
    rows.append(("SSW, 114 observed weakenings (NCEP)", p2["R_obs"]["mean"], p2["R_obs"]["ci95"], p2["R_h2"]["mean"], "h2"))
    rows.append(("SSW, 143 observed weakenings (ERA5)", s2["R_obs"]["mean"], s2["R_obs"]["ci95"], s2["R_h2"]["mean"], "h2"))
    rows.append(("SSW, spontaneous in SNAPSI", s1["R_obs"]["mean"], s1["R_obs"]["ci95"], s1["R_LS"]["mean"], "law"))
    s5 = cf["S5_cmip6_events_dose"]
    rows.append(("SSW, CMIP6 realised 100 hPa dose", s5["R_obs"], s5["R_obs_summary"]["ci95"], s5["R_h2"], "h2"))
    for k, lab in (("MAJOR", "cyclones, major (genesis SST)"), ("RI", "cyclones, rapid intensification")):
        t = cf["S3_tropical_cyclones"]["global_terciles"][k]; rows.append((lab, t["R_obs"], None, t["R_LS"], "law"))
    rows.append(("MJO crossing the Maritime Continent (QBO)", mj["R_obs"], None, mj["R_LS"], "law"))
    for k, lab in (("CEUR_JJA", "hot summers, central Europe"), ("WRUS_JJA", "hot summers, western Russia"),
                   ("CNAM_JJA", "hot summers, central N. America"), ("CEAS_DJF", "cold winters, central Eurasia")):
        t = cf["S4_human_forcing_hadgem3"][k]; rows.append((lab + " (HadGEM3)", t["R_obs"], None, t["R_LS"], "law"))
    for m, t in en["T4_forced_share"].items():
        if "R_obs" in t:
            rows.append((f"extreme El Niño, {m}", t["R_obs"], None, t["R_LS"], "law"))
    ax = fig.add_subplot(gs[0, :])
    yy = np.arange(len(rows))[::-1]
    for y_, (lab, r, ci, pr, kind) in zip(yy, rows):
        if ci:
            ax.plot(ci, [y_, y_], color=verm, lw=0.9)
        ax.scatter([r], [y_], s=14, color=verm, zorder=3)
        ax.scatter([pr], [y_], s=14, facecolor="white", edgecolor="k", lw=0.6, zorder=4)
        ax.text(0.465, y_, f"×{1 / r:.0f}", va="center", fontsize=5.2, color="0.3")
    ax.set_yticks(yy); ax.set_yticklabels([r[0] for r in rows], fontsize=5.2)
    ax.set_xlim(0, 0.5); ax.set_xlabel("share of the class contrast that is a difference in forcing")
    ax.text(0.465, len(rows) - 0.2, "contrast /\nforcing", ha="center", va="bottom", fontsize=5.0, color="0.3")
    ax.scatter([], [], s=14, color=verm, label="observed (95% interval where stored)")
    ax.scatter([], [], s=14, facecolor="white", edgecolor="k", lw=0.6, label="predicted (law, or variance share h²)")
    ax.legend(loc="lower right", fontsize=4.8, frameon=False, bbox_to_anchor=(0.9, 0.0))
    ax.set_title("classes cut from the outcome overstate their difference in forcing", fontsize=6)
    S.panel_label(ax, "a", x=-0.28, y=1.02)

    cx = fig.add_subplot(gs[1, 0])
    ev = ic["events"]
    x = np.array([e["z100_wk2"] for e in ev]); y = np.array([e["y"] for e in ev])
    start = np.array([e["onset_day"] == 0.0 for e in ev])
    cx.scatter(x[~start], y[~start], s=12, color=verm, zorder=3, label="onset within the run")
    cx.scatter(x[start], y[start], s=12, facecolor="white", edgecolor=verm, zorder=3, label="run starts on onset day")
    b = np.polyfit(x, y, 1); xx = np.linspace(x.min(), x.max(), 10)
    cx.plot(xx, np.polyval(b, xx), color="k", lw=0.6)
    m4, m2, m3 = ic["M4_corr_week2_100hPa_vs_outcome"], ic["M2_between_event_variance"], ic["M3_share_of_within_ensemble_variance"]
    cx.text(0.03, 0.97, (f"r = {m4['r']:.2f} [{m4['ci95'][0]:.2f}, {m4['ci95'][1]:.2f}], {ic['n_events']} events\n"
                          f"between-event variance {m2['lower_bound_noise_sd_equals_daily_sd']:.2f}–{m2['upper_bound_noise_zero']:.2f}"
                          f"\n({m3['share_lower']:.0%}–{m3['share_upper']:.0%} of within-ensemble)"),
            transform=cx.transAxes, va="top", fontsize=5.4)
    cx.set_xlabel("week-2 100 hPa anomaly (NAM sign, s.d.)")
    cx.set_ylabel("ensemble-mean surface response,\ndays 8-25 (NAM sign, s.d.)")
    cx.legend(loc="lower right", fontsize=4.8)
    cx.set_title("events differ, with the lower stratosphere (ICON)", fontsize=6)
    S.panel_label(cx, "b", x=-0.2, y=1.04)

    dx = fig.add_subplot(gs[1, 1])
    pe = [e for e in arc["per_ensemble"] if e["init"] in ("s20180125", "s20181213")]
    for c in sorted({e["centre"] for e in pe}):
        a = [e for e in pe if e["centre"] == c and e["init"] == "s20180125"]
        b_ = [e for e in pe if e["centre"] == c and e["init"] == "s20181213"]
        if a and b_:
            dx.plot([0, 1], [a[0]["nudged_DW_rate"], b_[0]["nudged_DW_rate"]], color="0.7", lw=0.6, zorder=1)
            dx.scatter([0, 1], [a[0]["nudged_DW_rate"], b_[0]["nudged_DW_rate"]], s=8, color=grey, zorder=2)
    for xpos, init in ((0, "s20180125"), (1, "s20181213")):
        m = np.mean([e["nudged_DW_rate"] for e in pe if e["init"] == init])
        dx.plot([xpos - 0.18, xpos + 0.18], [m, m], color=verm, lw=1.6)
        dx.text(xpos + 0.21, m, f"{m:.2f}", va="center", fontsize=5.6, color=verm)
    dx.set_xticks([0, 1]); dx.set_xticklabels(["Feb 2018\n(observed: downward)", "Jan 2019\n(labelled non-downward)"], fontsize=5.6)
    dx.set_xlim(-0.5, 1.6); dx.set_ylim(0, 1.05)
    dx.set_ylabel("forced P(downward), per model")
    dx.set_title("the archetype pair: similar forced odds", fontsize=6)
    S.panel_label(dx, "c", x=-0.16, y=1.04)
    S.save(fig, "fig3_scope")


if __name__ == "__main__":
    main()
