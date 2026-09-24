#!/usr/bin/env python3
"""
figED1_loeffel.py -- Extended Data Fig. 1: the lower-stratospheric precursor
couples to the surface just as strongly with no SSW.

Reads results/current/8_experiment/snapsi_loeffel_test.json (result N).
Recomputes nothing.

  a  per ensemble (centre x initialisation), within-member correlation of
     week-2 100 hPa polar-cap geopotential height with the days 15-25 surface
     response: SSW imposed (nudged) against the same ensemble with no SSW
     (control); short-lead initialisation s20190108 marked
  b  Fisher-pooled r per arm with 95% intervals, all matched ensembles and
     excluding s20190108, with the nudged-minus-control difference
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nature_style as S  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

SHORT = "s20190108"


def main():
    S.apply()
    d = json.loads((S.RESULTS / "8_experiment" / "snapsi_loeffel_test.json").read_text())
    p = d["primary"]
    nud = {(q["centre"], q["init"]): q["r"] for q in p["nudged"]["per_ensemble"]}
    con = {(q["centre"], q["init"]): q["r"] for q in p["control"]["per_ensemble"]}
    keys = sorted(set(nud) & set(con))
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(S.DOUBLE * 0.8, 62 * S.MM),
                                 gridspec_kw={"width_ratios": [1, 1], "wspace": 0.4})
    for k in keys:
        short = k[1] == SHORT
        ax.scatter(con[k], nud[k], s=10, marker="s" if short else "o",
                   facecolor="white" if short else S.C["vermillion"],
                   edgecolor=S.C["vermillion"], linewidths=0.6, zorder=3)
    lim = (-0.45, 0.55)
    ax.plot(lim, lim, color="k", lw=0.5, ls=(0, (2, 2)))
    ax.axhline(0, color="0.7", lw=0.4)
    ax.axvline(0, color="0.7", lw=0.4)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.set_xlabel("r within ensemble, no SSW (control)")
    ax.set_ylabel("r within ensemble, SSW imposed (nudged)")
    ax.scatter([], [], s=10, color=S.C["vermillion"], label="ensemble")
    ax.scatter([], [], s=10, marker="s", facecolor="white", edgecolor=S.C["vermillion"],
               label=f"{SHORT} (initialised after onset)")
    ax.legend(loc="lower right", fontsize=5.3)
    ax.set_title(f"{len(keys)} matched ensembles, 9 models", fontsize=6.5)
    S.panel_label(ax, "a", x=-0.22, y=1.04)

    rows = [("all", "arm_contrast", "pooled"),
            (f"excl. {SHORT}", "arm_contrast_excl_short_lead", "pooled_excl_short_lead")]
    for j, (lab, ck, pk) in enumerate(rows):
        for dx, arm in ((-0.12, "nudged"), (0.12, "control")):
            pl = p[arm][pk]
            y = j + dx
            bx.plot(pl["CI95"], [y, y], color=S.ARM[arm], lw=1.6, solid_capstyle="butt")
            bx.scatter([pl["r"]], [y], s=12, color=S.ARM[arm], zorder=3,
                       label=("SSW imposed" if arm == "nudged" else "no SSW") if j == 0 else None)
        ac = p[ck]
        bx.text(0.36, j, f"Δz {ac['delta_z_nudged_minus_control']:+.3f}\n"
                         f"[{ac['delta_z_CI95'][0]:+.3f}, {ac['delta_z_CI95'][1]:+.3f}]\n"
                         f"p = {ac['delta_p_two_sided']:.2f}, n = {ac['n_matched_ensembles']}",
                fontsize=5.3, va="center")
    bx.axvline(0, color="k", lw=0.4)
    bx.set_yticks([0, 1])
    bx.set_yticklabels([r[0] for r in rows])
    bx.set_ylim(-0.6, 1.6)
    bx.set_xlim(-0.05, 0.6)
    bx.invert_yaxis()
    bx.set_xlabel("pooled within-ensemble r (95% CI)")
    bx.legend(loc="lower right", fontsize=5.3)
    bx.set_title("Loeffel et al. 2026, across 18 events: r = 0.85", fontsize=6, color="0.35")
    S.panel_label(bx, "b", x=-0.3, y=1.04)
    S.save(fig, "figED1_loeffel")


if __name__ == "__main__":
    main()
