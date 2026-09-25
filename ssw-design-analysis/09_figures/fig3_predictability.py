#!/usr/bin/env python3
"""
fig3_predictability.py -- Figure 3: an SSW adds no out-of-sample predictability.

Reads results/current/6_predictability/within_model_check.json (result J) and
the R^2 each DW/NDW split implies (eta^2: q(1-q)C^2 over the variance of the
same outcome) from the scripts that apply the criterion:
predictability_ceiling.json (CMIP6), and 9_literature/
recompute_published_criterion.json (AO) and era5_recompute_and_two_thirds.json
(ERA5 NAM). Recomputes nothing but histograms.

  a-c  within-member cross-validated R^2 at real SSWs (line) against 1,000
       event-free pseudo-onset draws of the SAME size (histogram), per
       predictor tier
  d    event-specific R^2 (real minus null mean) with its null-based and
       member-bootstrap 95% intervals, against the R^2 implied by the DW/NDW
       split when the published criterion is applied to this project's data
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nature_style as S  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

TIERS = [("P1 pre-onset", "predictors before day 0"),
         ("P2 at-onset", "adds days −5 to 0"),
         ("P3 + post-onset stratosphere", "adds stratosphere, days 0 to +30")]


def main():
    S.apply()
    j = json.loads((S.RESULTS / "6_predictability" / "within_model_check.json").read_text())
    pc = json.loads((S.RESULTS / "6_predictability" / "predictability_ceiling.json").read_text())
    fig = plt.figure(figsize=(S.DOUBLE, 60 * S.MM))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 1.35], wspace=0.42)
    for i, (key, sub) in enumerate(TIERS):
        ax = fig.add_subplot(gs[0, i])
        sm = j[key]["size_matched"]
        null = np.array(sm["null_draws"])
        ax.hist(null, bins=30, color="0.78", density=True)
        ax.axvline(sm["within_real_cv_r2"], color=S.C["vermillion"], lw=1.2)
        ax.set_title(f"{key.split(' ', 1)[0]}: {sub}", fontsize=6)
        ax.set_xlabel("out-of-sample R²")
        if i == 0:
            ax.set_ylabel("density (event-free draws)")
        ax.set_yticks([])
        ax.text(0.03, 0.97, f"SSW adds\n{sm['event_specific']:+.3f}\np = {sm['p_null_ge_real']:.2f}",
                transform=ax.transAxes, va="top", fontsize=5.5)
        S.panel_label(ax, "abc"[i], x=-0.12, y=1.08)
    ax.plot([], [], color=S.C["vermillion"], lw=1.2, label="real SSWs")
    ax.plot([], [], color="0.78", lw=4, label="no SSW, same onset counts")
    ax.legend(loc="upper right", fontsize=5.3, bbox_to_anchor=(1.02, 0.8))

    dx = fig.add_subplot(gs[0, 3])
    y = 0
    labels = []
    for key, short in (("P1 pre-onset", "SSW-specific R², pre-onset"),
                       ("P3 + post-onset stratosphere", "SSW-specific R², post-onset")):
        sm = j[key]["size_matched"]
        lo_b, hi_b = sm["bootstrap_CI95"]
        lo_n, hi_n = sm["event_specific_CI95_from_null"]
        dx.plot([lo_b, hi_b], [y, y], color=S.C["blue"], lw=0.8)
        dx.plot([lo_n, hi_n], [y, y], color=S.C["blue"], lw=2.4, solid_capstyle="butt")
        dx.scatter([sm["event_specific"]], [y], s=12, color="k", zorder=3)
        labels.append(short)
        y += 1
    y += 0.4
    lit = S.RESULTS / "9_literature"
    rc = json.loads((lit / "recompute_published_criterion.json").read_text())
    er = json.loads((lit / "era5_recompute_and_two_thirds.json").read_text())
    imp = {"CMIP6, Karpechko cond. 1-2":
               pc["implied"]["CMIP6 Karpechko conditions 1-2"]["implied_r2"],
           "obs. AO, Karpechko": rc["outcomes"]["ao"]["implied_r2_eta2"],
           "ERA5 1000 hPa, Karpechko": er["variants"]["Karpechko_1000hPa"]["implied_r2_eta2"],
           "ERA5 850 hPa, ACP 2026": er["variants"]["ACP2026_850hPa"]["implied_r2_eta2"]}
    order = sorted(imp, key=imp.get)
    for name in order:
        dx.scatter([imp[name]], [y], s=14, marker="D", color=S.C["vermillion"], zorder=3)
        labels.append("DW/NDW split implies: " + name)
        y += 1
    dx.axvline(0, color="k", lw=0.4)
    dx.set_yticks(range(2))
    ticks = [0, 1] + list(np.arange(2.4, 2.4 + len(order)))
    dx.set_yticks(ticks)
    dx.set_yticklabels(labels, fontsize=5.2)
    dx.yaxis.tick_right()
    dx.set_xlabel("R² of the between-event surface response")
    dx.set_xlim(-0.08, 0.7)
    dx.plot([], [], color=S.C["blue"], lw=2.4, label="95%, event-free null")
    dx.plot([], [], color=S.C["blue"], lw=0.8, label="95%, member bootstrap")
    dx.legend(loc="lower right", fontsize=5.2, bbox_to_anchor=(1.0, 0.22))
    S.panel_label(dx, "d", x=-0.05, y=1.08)
    S.save(fig, "fig3_predictability")


if __name__ == "__main__":
    main()
