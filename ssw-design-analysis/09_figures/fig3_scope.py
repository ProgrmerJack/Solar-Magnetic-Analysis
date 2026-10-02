#!/usr/bin/env python3
"""
fig3_scope.py -- Figure 3: what the result permits -- distributional departures and
differences between events.

Reads (recomputes nothing):
  results/current/8_experiment/snapsi_residual_quantiles.json   residual quantiles (Z1)
  results/current/5_mechanism/icon_event_heterogeneity.json     ICON event ensembles (Z3)
  results/current/8_experiment/snapsi_archetype_test.json       archetype pair (O)

  a  polar-cap NAM, days 8-25: residual quantile of the SSW-imposed members against
     the control translated by the mean shift, pooled over 36 ensembles and per event
     (95% two-stage bootstrap); grey, the tolerance declared before the test
  b  the same for northern-Eurasian temperature, days 8-24
  c  18 ICON SSWs re-run as 40-member ensembles (Loeffel et al. 2026): ensemble-mean
     surface response, days 8-25, against the week-2 100 hPa anomaly (both NAM sign);
     open, the six ensembles that start on their onset day (already easterly at t = 0)
  d  the field's archetype pair: forced probability of a downward outcome in each
     model (nudged members), February 2018 and January 2019, primary initialisations
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nature_style as S  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from figED8_translation import panel  # noqa: E402

M = lambda s: s.replace("-", "−")                     # noqa: E731


def main():
    S.apply()
    R = S.RESULTS
    rq = json.loads((R / "8_experiment" / "snapsi_residual_quantiles.json").read_text())
    ic = json.loads((R / "5_mechanism" / "icon_event_heterogeneity.json").read_text())["primary_registered"]
    arc = json.loads((R / "8_experiment" / "snapsi_archetype_test.json").read_text())
    grey, verm = "0.6", S.C["vermillion"]
    fig = plt.figure(figsize=(S.DOUBLE, 112 * S.MM))
    gs = fig.add_gridspec(2, 2, wspace=0.34, hspace=0.62)
    tR = rq["tolerance"]["R_sigma"]

    ax = fig.add_subplot(gs[0, 0])
    A = rq["A_polar_cap_NAM_days8_25"]
    panel(ax, A, tR, "R", "residual quantile (σ of control)")
    ax.set_ylim(-0.6, 0.6); ax.legend(loc="lower left", fontsize=4.8)
    ax.set_title(f"polar-cap NAM: pooled within tolerance; single events {A['Feb-2018']['reading']}", fontsize=5.8)
    S.panel_label(ax, "a", x=-0.18, y=1.04)

    bx = fig.add_subplot(gs[0, 1])
    T = rq["T_northern_Eurasia_days8_24"]
    panel(bx, T, tR, "R", "residual quantile (σ of control)")
    bx.set_ylim(-0.6, 0.6)
    bx.set_title(f"northern-Eurasian temperature: {T['pooled']['reading']}", fontsize=5.8)
    S.panel_label(bx, "b", x=-0.18, y=1.04)

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
    S.panel_label(cx, "c", x=-0.2, y=1.04)

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
    S.panel_label(dx, "d", x=-0.16, y=1.04)
    S.save(fig, "fig3_scope")


if __name__ == "__main__":
    main()
