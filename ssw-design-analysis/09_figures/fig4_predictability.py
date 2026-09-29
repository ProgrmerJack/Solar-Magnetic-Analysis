#!/usr/bin/env python3
"""
fig4_predictability.py -- Figure 4: knowing the event adds no forecast
information after an SSW.

Reads results/current/6_predictability/forecast_value_test.json (Fig. 3a) and
within_model_check.json (result J; Fig. 3b,c). Recomputes nothing but histograms.

  a  CRPS skill of the event-aware forecast over the shifted distribution, on
     held-out models: after SSWs (point, 95% model-bootstrap interval) and on 200
     size-matched sets of ordinary winter days (mean, 2.5-97.5% range), for
     predictors before onset (P1) and up to onset (P2)
  b,c within-model out-of-sample R^2 at real SSWs (line) against 1,000
     event-free sets with the same onset counts (histogram), before onset and
     including the post-onset stratosphere

The implied-R^2 comparison of earlier versions is gone: it set diagnostic
variance explained beside predictive R^2 (review 2026-09-25).
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nature_style as S  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

TIERS_J = [("P1 pre-onset", "before onset"),
           ("P3 + post-onset stratosphere", "adds post-onset stratosphere")]
TIERS_C = [("P1 pre-onset", "P1: before onset"), ("P2 at-onset", "P2: up to onset")]


def main():
    S.apply()
    fv = json.loads((S.RESULTS / "6_predictability" / "forecast_value_test.json").read_text())
    j = json.loads((S.RESULTS / "6_predictability" / "within_model_check.json").read_text())
    fig = plt.figure(figsize=(S.DOUBLE, 60 * S.MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 1], wspace=0.45)

    ax = fig.add_subplot(gs[0, 0])
    for y, (key, lab) in enumerate(TIERS_C):
        t = fv["tiers"][key]
        lo, hi = t["CRPSS_pseudo_q025_q975"]
        ax.plot([100 * lo, 100 * hi], [y - 0.12, y - 0.12], color="0.6", lw=3, solid_capstyle="butt")
        ax.scatter([100 * t["CRPSS_pseudo_mean"]], [y - 0.12], s=14, color="0.45", zorder=3,
                   label="ordinary winter days" if y == 0 else None)
        lo, hi = t["CRPSS_ssw_CI95_model_bootstrap"]
        ax.plot([100 * lo, 100 * hi], [y + 0.12, y + 0.12], color=S.C["vermillion"], lw=1.2)
        ax.scatter([100 * t["CRPSS_ssw"]], [y + 0.12], s=14, color=S.C["vermillion"], zorder=3,
                   label="after SSWs" if y == 0 else None)
    ax.axvline(0, color="k", lw=0.5)
    ax.set_yticks(range(len(TIERS_C)))
    ax.set_yticklabels([lab for _, lab in TIERS_C])
    ax.set_ylim(-0.6, len(TIERS_C) - 0.4)
    ax.set_xlabel("skill of event-aware over shift-only\nforecast (% of CRPS)")
    ax.legend(loc="lower right", fontsize=5.3)
    ax.set_title(f"{fv['n_events']:,} CMIP6 SSWs, held-out models", fontsize=6)
    S.panel_label(ax, "a", x=-0.42, y=1.08)

    for i, (key, sub) in enumerate(TIERS_J):
        bx = fig.add_subplot(gs[0, i + 1])
        sm = j[key]["size_matched"]
        null = np.array(sm["null_draws"])
        bx.hist(null, bins=30, color="0.78", density=True)
        bx.axvline(sm["within_real_cv_r2"], color=S.C["vermillion"], lw=1.2)
        bx.set_title(sub, fontsize=6)
        bx.set_xlabel("out-of-sample R²")
        bx.set_yticks([])
        if i == 0:
            bx.set_ylabel("density (event-free sets)")
        bx.text(0.03, 0.97, f"SSW adds\n{sm['event_specific']:+.3f}\np = {sm['p_null_ge_real']:.2f}",
                transform=bx.transAxes, va="top", fontsize=5.5)
        S.panel_label(bx, "bc"[i], x=-0.12, y=1.08)
    # colours mean the same in all panels (vermilion: SSWs; grey: no SSW), so
    # panel a's legend serves b and c and nothing covers the histograms
    S.save(fig, "fig4_predictability")


if __name__ == "__main__":
    main()
