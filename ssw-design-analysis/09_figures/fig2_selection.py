#!/usr/bin/env python3
"""
fig2_selection.py -- Figure 2: the downward label is a threshold on one population.

Reads results/current/8_experiment/snapsi_selection_test.json (result L).
Recomputes nothing except one closed-form reference line:

  splitting a unit-variance Gaussian at its mean gives a difference of group
  means of 2*sqrt(2/pi) = 1.596 sigma. Control members are standardised by
  their own ensemble (mean 0, sd 1), so this is what the Karpechko surface
  conditions produce when applied to pure noise.

  a  PAIRED: DW-minus-NDW contrast of the same centre x initialisation with the
     SSW imposed (nudged) against no SSW (control), NH, with the paired
     difference and its centre-bootstrap interval from paired_NH
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
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(S.DOUBLE * 0.82, 66 * S.MM),
                                 gridspec_kw={"width_ratios": [1, 1], "wspace": 0.45})
    rng = np.random.default_rng(0)          # jitter only (panel b)
    arms = ("nudged", "control")
    # a: PAIRED -- the same centre x initialisation with and without the SSW
    rows = {}
    for r in d["per_case"]:
        if r["hemisphere"] == "NH":
            rows.setdefault((r["centre"], r["init"]), {})[r["arm"]] = r["contrast_sigma"]
    pairs = [v for v in rows.values() if {"nudged", "control"} <= set(v)]
    xs = np.array([v["control"] for v in pairs]); ys = np.array([v["nudged"] for v in pairs])
    lim = (-2.45, -0.95)
    ax.plot(lim, lim, color="k", lw=0.6, ls=(0, (3, 2)))
    ax.axvline(REF, color=S.C["green"], lw=0.6, ls=":")
    ax.axhline(REF, color=S.C["green"], lw=0.6, ls=":")
    ax.scatter(xs, ys, s=10, color=S.ARM["nudged"], edgecolor="white", linewidths=0.3, zorder=3)
    pr = d["paired_NH"]
    lo, hi = pr["paired_difference_CI95_centre_bootstrap"]
    ax.text(0.03, 0.97, f"{pr['n_pairs']} pairs, {pr['n_centres']} models\n"
                        f"SSW imposed − no SSW:\n{pr['paired_difference']:+.2f}σ [{lo:+.2f}, {hi:+.2f}]",
            transform=ax.transAxes, va="top", fontsize=5.8)
    ax.text(lim[0] + 0.03, REF - 0.03, "dotted: −2√(2/π)\n(one Gaussian cut)",
            color=S.C["green"], fontsize=5.2, va="top", ha="left")
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
    ax.set_xlabel("DW − NDW contrast, no SSW (control, σ)")
    ax.set_ylabel("DW − NDW contrast, SSW imposed (nudged, σ)")
    S.panel_label(ax, "a", x=-0.2, y=1.05)

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
    bx.set_title("the downward rate separates the arms; the contrast does not", fontsize=6.5)
    S.panel_label(bx, "b", x=-0.2, y=1.08)
    S.save(fig, "fig2_selection")


if __name__ == "__main__":
    main()
