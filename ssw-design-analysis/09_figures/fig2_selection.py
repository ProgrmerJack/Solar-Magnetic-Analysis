#!/usr/bin/env python3
"""
fig2_selection.py -- Figure 2: the DW/NDW contrast needs no SSW.

Reads results/current/8_experiment/snapsi_selection_test.json (result L).
Recomputes nothing except one closed-form reference line:

  splitting a unit-variance Gaussian at its mean gives a difference of group
  means of 2*sqrt(2/pi) = 1.596 sigma. Control members are standardised by
  their own ensemble (mean 0, sd 1), so this is what the Karpechko surface
  conditions produce when applied to pure noise.

  a  DW-minus-NDW surface contrast per ensemble, nudged vs control, NH, with
     the SH minor warming and the published ACP 26, 3723 (2026) value for scale
  b  the fraction of members classified DW, per ensemble (unconditional)
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nature_style as S  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

REF = -2 * math.sqrt(2 / math.pi)


def main():
    S.apply()
    d = json.loads((S.RESULTS / "8_experiment" / "snapsi_selection_test.json").read_text())
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(S.DOUBLE, 62 * S.MM),
                                 gridspec_kw={"width_ratios": [1.35, 1], "wspace": 0.3})
    rng = np.random.default_rng(0)          # jitter only
    arms = ("nudged", "control")
    lab = {"nudged": "SSW imposed\n(nudged, identical in every member)",
           "control": "no SSW\n(control)"}
    for x, arm in enumerate(arms):
        nh = [r["contrast_sigma"] for r in d["per_case"]
              if r["arm"] == arm and r["hemisphere"] == "NH"]
        sh = [r["contrast_sigma"] for r in d["per_case"]
              if r["arm"] == arm and r["hemisphere"] == "SH"]
        ax.scatter(x + rng.uniform(-0.12, 0.12, len(nh)), nh, s=6,
                   color=S.ARM[arm], alpha=0.8, linewidths=0, label=None)
        ax.scatter(x + 0.3 + rng.uniform(-0.04, 0.04, len(sh)), sh, s=9, marker="^",
                   facecolor="white", edgecolor=S.ARM[arm], linewidths=0.6)
        m = d["summary"][arm]["mean_contrast_sigma"]
        ax.plot([x - 0.2, x + 0.2], [m, m], color="k", lw=1.2)
        ax.text(x - 0.24, m, f"{m:+.2f}", ha="right", va="center", fontsize=6)
    ax.axhline(REF, color=S.C["green"], lw=0.8, ls=(0, (3, 2)))
    ax.text(1.5, REF - 0.05, "pure noise split at its mean:\n$-2\\sqrt{2/\\pi}$ = −1.60σ",
            color=S.C["green"], fontsize=5.5, va="top", ha="left")
    pub = d["published_contrast_for_scale"]
    ax.axhline(pub, color="0.4", lw=0.6, ls=":")
    ax.text(1.5, pub - 0.04, "published contrast\n(ACP 26, 3723, 2026)", color="0.4",
            fontsize=5.5, va="top", ha="left")
    # Caption, not panel: 11 of 36 nudged ensembles have < 3 NDW members (DW
    # rate 0.99) and cannot form a contrast, so panel a shows 25 of 36.
    ax.set_xticks([0, 1])
    ax.set_xticklabels([lab[a] for a in arms])
    ax.set_xlim(-0.5, 2.45)
    ax.set_ylim(-2.6, -0.75)
    ax.set_ylabel("DW − NDW surface contrast (σ)")
    ax.scatter([], [], s=6, color="0.4", label="NH ensemble (centre × initialisation)")
    ax.scatter([], [], s=9, marker="^", facecolor="white", edgecolor="0.4",
               linewidths=0.6, label="SH minor warming, Sep 2019")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2)
    S.panel_label(ax, "a", x=-0.14, y=1.08)

    for x, arm in enumerate(arms):
        r = [e["DW_rate"] for e in d["all_ensembles"]
             if e["arm"] == arm and e["hemisphere"] == "NH"]
        bx.scatter(x + rng.uniform(-0.12, 0.12, len(r)), r, s=6, color=S.ARM[arm],
                   alpha=0.8, linewidths=0)
        m = d["summary"][arm]["mean_DW_rate_UNCONDITIONAL"]
        bx.plot([x - 0.2, x + 0.2], [m, m], color="k", lw=1.2)
        bx.text(x + 0.24, m, f"{m:.2f}", ha="left", va="center", fontsize=6)
    bx.set_xticks([0, 1])
    bx.set_xticklabels(["SSW imposed", "no SSW"])
    bx.set_xlim(-0.5, 1.6)
    bx.set_ylim(0, 1.05)
    bx.set_ylabel("fraction of members classified DW")
    bx.set_title("the RATE separates the arms; the contrast does not", fontsize=6.5)
    S.panel_label(bx, "b", x=-0.2, y=1.08)
    S.save(fig, "fig2_selection")


if __name__ == "__main__":
    main()
