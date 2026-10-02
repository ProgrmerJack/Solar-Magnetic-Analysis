#!/usr/bin/env python3
"""
figED10_generality.py -- Extended Data Fig.: the class structure where the forcing
is strongest, in the Southern Hemisphere warming, and for regional temperature.

Reads (recomputes nothing): results/current/8_experiment/snapsi_contrast_generality.json

  a  back-shift test: contrast of SSW-imposed members displaced back by the
     imposed effect minus the control contrast, per centre-initialisation pair,
     against the imposed shift; vermilion, the 11 pairs whose SSW-imposed
     ensemble has < 3 NDW members; grey band, the tolerance fixed before the test
  b  the same differences against what a contrast proportional to the spread
     predicts from each pair's change in spread (post hoc); dashed, equality
  c  Southern Hemisphere minor warming: nudged contrast minus the contrast of
     control members shifted (shift only) or shifted and rescaled (location-scale)
     by the imposed change; NH pairs for comparison; 95% centre-bootstrap intervals
  d  northern-Eurasian temperature, days 8-24: DW - NDW contrast among SSW-imposed
     members against the contrast among control members displaced by the imposed
     NAM shift and classified identically, each keeping its own temperature
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nature_style as S  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

M = lambda s: s.replace("-", "−")                     # noqa: E731


def main():
    S.apply()
    d = json.loads((S.RESULTS / "8_experiment" / "snapsi_contrast_generality.json").read_text())
    p = pd.DataFrame(d["pairs_NH"]); rg = pd.DataFrame(d["pairs_regional"])
    tol = d["tolerance_sigma"]
    grey, verm = "0.55", S.C["vermillion"]
    fig = plt.figure(figsize=(S.DOUBLE, 112 * S.MM))
    gs = fig.add_gridspec(2, 2, wspace=0.34, hspace=0.62)

    ax = fig.add_subplot(gs[0, 0])
    ax.axhspan(-tol, tol, color="0.92", zorder=0)
    for st, col, lab in ((False, grey, "other 25 pairs"), (True, verm, "11 pairs without an NDW class")):
        q = p[p.strong == st]
        ax.scatter(q.s, q.Delta_b, s=12, color=col, zorder=3, label=lab)
    g1 = d["G1_strongly_forced"]
    sl = d["G2_dose_response"]["Delta_b_on_s_36"]
    xx = np.linspace(p.s.min(), p.s.max(), 10)
    ax.plot(xx, sl["intercept"] + sl["slope"] * xx, color="k", lw=0.6)
    m11 = g1["strong_11"]["Delta_b"]
    ax.axhline(0, color="k", lw=0.4)
    ax.set_xlabel("imposed shift of the NAM proxy, days 8-25 (σ of control)")
    ax.set_ylabel("back-shifted minus control contrast (σ)")
    ax.legend(loc="upper left", fontsize=4.8)
    ax.set_title(M(f"11 pairs: {m11['mean']:+.2f} [{m11['ci95'][0]:+.2f}, {m11['ci95'][1]:+.2f}]σ; "
                   f"slope {sl['slope']:+.2f} [{sl['ci95'][0]:+.2f}, {sl['ci95'][1]:+.2f}] per σ"), fontsize=5.8)
    S.panel_label(ax, "a", x=-0.2, y=1.04)

    bx = fig.add_subplot(gs[0, 1])
    pred = p.C_control * (p.sd_ratio - 1)
    lim = (min(pred.min(), p.Delta_b.min()) - 0.05, max(pred.max(), p.Delta_b.max()) + 0.05)
    bx.plot(lim, lim, color="k", lw=0.6, ls=(0, (3, 2)))
    bx.scatter(pred[~p.strong], p.Delta_b[~p.strong], s=12, color=grey, zorder=3)
    bx.scatter(pred[p.strong], p.Delta_b[p.strong], s=12, color=verm, zorder=3)
    gx = d["G1x_post_hoc_spread"]
    n11 = gx["strong_11"]["Delta_b_net_spread"]
    bx.text(0.97, 0.03, M(f"r = {gx['all_36']['r_Delta_b_vs_spread_pred']:.2f} (36 pairs)\n"
                          f"11 pairs net of spread: {n11['mean']:+.2f} [{n11['ci95'][0]:+.2f}, {n11['ci95'][1]:+.2f}]"),
            transform=bx.transAxes, va="bottom", ha="right", fontsize=5.4)
    bx.set_xlim(lim); bx.set_ylim(lim)
    bx.set_xlabel("predicted from the change in spread (σ)")
    bx.set_ylabel("back-shifted minus control contrast (σ)")
    bx.set_title("departures follow the spread (post hoc)", fontsize=6)
    S.panel_label(bx, "b", x=-0.24, y=1.04)

    cx = fig.add_subplot(gs[1, 0])
    g3 = d["G3_third_event_SH"]
    rows = [("SH, shift only", g3["SH"]["resid_shift"], verm), ("SH, location-scale", g3["SH"]["resid_locscale"], verm),
            ("NH 25 pairs, shift only", g3["NH_25"]["resid_shift"], grey),
            ("NH 25 pairs, location-scale", g3["NH_25"]["resid_locscale"], grey)]
    cx.axvspan(-tol, tol, color="0.92", zorder=0)
    for i, (lab, r, col) in enumerate(rows):
        cx.plot(r["ci95"], [i, i], color=col, lw=1.1); cx.scatter([r["mean"]], [i], s=14, color=col, zorder=3)
    cx.axvline(0, color="k", lw=0.5); cx.set_xlim(-0.32, 0.32)
    cx.set_yticks(range(len(rows))); cx.set_yticklabels([r[0] for r in rows], fontsize=5.4)
    cx.set_ylim(-0.6, len(rows) - 0.4); cx.invert_yaxis()
    cx.set_xlabel("nudged contrast minus null contrast (σ)")
    sh = g3["SH"]
    cx.set_title(M(f"SH warming: contrast {sh['C_nudged']['mean']:.2f} vs {sh['C_control']['mean']:.2f}σ, "
                   f"s.d. ratio {sh['sd_ratio']['mean']:.2f}"), fontsize=6)
    S.panel_label(cx, "c", x=-0.42, y=1.04)

    dx = fig.add_subplot(gs[1, 1])
    q = rg.dropna(subset=["CT_nudged", "CT_null"])
    lim = (min(q.CT_nudged.min(), q.CT_null.min()) - 0.1, max(q.CT_nudged.max(), q.CT_null.max()) + 0.1)
    dx.plot(lim, lim, color="k", lw=0.6, ls=(0, (3, 2)))
    dx.scatter(q.CT_null, q.CT_nudged, s=12, color=verm, zorder=3)
    g4 = d["G4_regional_temperature"]
    r, pl = g4["resid"], g4["resid_plant_0.5sigma"]
    dx.set_title(M(f"residual {r['mean']:+.2f}σ [{r['ci95'][0]:+.2f}, {r['ci95'][1]:+.2f}] "
                          f"(≈ {g4['resid_K_approx']:+.2f} K), {r['n_pairs']} pairs\n"
                          f"with a planted 0.5σ: {pl['mean']:+.2f} [{pl['ci95'][0]:+.2f}, {pl['ci95'][1]:+.2f}]"),
                 fontsize=5.6)
    dx.set_xlim(lim); dx.set_ylim(lim)
    dx.set_xlabel("matched shifted null from control members (σ)")
    dx.set_ylabel("N-Eurasian T, DW − NDW: SSW imposed (σ)")
    S.panel_label(dx, "d", x=-0.24, y=1.04)
    S.save(fig, "figED10_generality")


if __name__ == "__main__":
    main()
