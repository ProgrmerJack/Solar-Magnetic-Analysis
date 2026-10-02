#!/usr/bin/env python3
"""
fig2_selection.py -- Figure 2: in observations, event-free dates displaced by the SSW
shift reproduce the classes and predict held-out winters.

Reads (recomputes nothing):
  results/current/9_literature/era5_recompute_and_two_thirds.json   rate and contrast (H)
  results/current/9_literature/criterion_regional_contrast.json     regional contrast (new, revision 4)
  results/current/9_literature/criterion_transport_cv.json          held-out winters (Z2)

  a  share of the 39 observed SSWs meeting the published surface conditions, against
     event-free dates as they are and displaced by the measured shift (95% range of 200 sets)
  b  DW - NDW contrast of the 1000 hPa NAM (published criterion, days 8-52):
     observed against event-free dates (not shifted) classified identically
  c  DW - NDW contrast of 2 m temperature, days 8-24, in four regions -- a variable
     not used to define the classes: observed (95% event bootstrap for northern
     Eurasia, the registered target) against the matched null (mean, 95% range of
     4,000 replicates; shift estimated from other winters)
  d  whole-winter cross-validation: gain of the shifted null in the continuous ranked
     probability score of the days 8-52 NAM over unshifted climatology, and in the
     Brier score of the label over climatology and over the constant class rate
     (95% intervals from resampling winters); title, downward events observed and
     expected out of fold
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nature_style as S  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

M = lambda s: s.replace("-", "−")                     # noqa: E731  typographic minus


def main():
    S.apply()
    R = S.RESULTS / "9_literature"
    two = json.loads((R / "era5_recompute_and_two_thirds.json").read_text())
    rc = json.loads((R / "criterion_regional_contrast.json").read_text())
    tc = json.loads((R / "criterion_transport_cv.json").read_text())
    grey, verm = "0.6", S.C["vermillion"]
    fig = plt.figure(figsize=(S.DOUBLE, 108 * S.MM))
    gs = fig.add_gridspec(2, 2, wspace=0.55, hspace=0.8)

    # a -- the rate
    ax = fig.add_subplot(gs[0, 0])
    t = two["two_thirds_test"]
    vals = [t["null_pass_rate_c1c2"], t["shifted_pseudo_pass_rate_c1c2"], t["observed_pass_rate_c1c2"]]
    ax.bar(range(3), vals, color=["0.85", grey, verm], width=0.6)
    lo, hi = t["shifted_pseudo_CI95"]
    ax.errorbar(1, vals[1], yerr=[[vals[1] - lo], [hi - vals[1]]], color="k", lw=0.7, capsize=2)
    for i, v in enumerate(vals):
        ax.text(i - (0.2 if i == 1 else 0), v + 0.02, f"{v:.0%}", ha="center", fontsize=5.8)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["event-free\ndates", "event-free dates\n+ measured shift", f"observed\nSSWs ({two['n_events']})"], fontsize=5.6)
    ax.set_ylim(0, 1); ax.set_ylabel("share meeting the surface conditions")
    ax.set_title("the rate: 'about two thirds'", fontsize=6)
    S.panel_label(ax, "a", x=-0.18, y=1.04)

    # b -- the NAM contrast
    bx = fig.add_subplot(gs[0, 1])
    v = two["variants"]["Karpechko_1000hPa"]
    bx.bar([0, 1], [v["null_contrast"], v["contrast"]], color=[grey, verm], width=0.55)
    for i, x in enumerate([v["null_contrast"], v["contrast"]]):
        bx.text(i, x - 0.04, M(f"{x:.2f}σ"), ha="center", va="top", fontsize=5.8)
    bx.axhline(0, color="k", lw=0.5)
    bx.set_xticks([0, 1]); bx.set_xticklabels(["event-free dates\n(no shift)", "observed SSWs"], fontsize=5.6)
    bx.set_ylim(-0.95, 0.05); bx.set_ylabel("DW − NDW contrast, 1000 hPa NAM (σ)")
    bx.set_title(f"the contrast: {v['pct_reproduced_under_null']:.0f}% reproduced without SSWs", fontsize=6)
    S.panel_label(bx, "b", x=-0.22, y=1.04)

    # c -- a variable not used to define the classes
    cx = fig.add_subplot(gs[1, 0])
    blk = rc["primary_surface"]["cv"]
    regs = [("NEURASIA", "N. Eurasia\n(registered)"), ("HI_EUROPE", "high-lat.\nEurope"),
            ("MID_EASIA", "mid-lat.\nE. Asia"), ("MID_NAMER", "mid-lat.\nN. America")]
    for i, (k, lab) in enumerate(regs):
        r = blk[f"{k}_d8_24"]
        lo, hi = r["C_null_ci95_K"]
        cx.plot([i - 0.12] * 2, [lo, hi], color=grey, lw=2.2, solid_capstyle="butt")
        cx.scatter([i - 0.12], [r["C_null_mean_K"]], s=14, color="k", marker="_", zorder=3)
        cx.scatter([i + 0.12], [r["C_obs_K"]], s=12, color=verm, zorder=3)
        if "C_obs_ci95_K" in r:
            cx.plot([i + 0.12] * 2, r["C_obs_ci95_K"], color=verm, lw=0.8)
    cx.axhline(0, color="k", lw=0.5)
    cx.set_xticks(range(len(regs))); cx.set_xticklabels([r[1] for r in regs], fontsize=5.4)
    cx.set_ylabel("DW − NDW contrast, 2 m T, days 8-24 (K)")
    p = blk["NEURASIA_d8_24"]
    cx.set_title(M(f"regional temperature (not used to classify)\nnull {p['C_null_mean_K']:.2f} K, observed {p['C_obs_K']:.2f} K"),
                 fontsize=6)
    cx.text(0.98, 0.04, "grey: matched null (95%)\nvermilion: observed", transform=cx.transAxes, ha="right", fontsize=5.0)
    S.panel_label(cx, "c", x=-0.18, y=1.04)

    # d -- held-out winters
    dx = fig.add_subplot(gs[1, 1])
    pr, cr = tc["primary_published_surface"], tc["crps_window_mean_days8_52"]
    rows = [("CRPS, NAM days 8-52:\nshifted vs climatology", cr["climatology_minus_shifted"], cr["ci95"]),
            ("Brier, label:\nshifted vs climatology", pr["brier_climatology_minus_shifted"], pr["ci95_climatology_minus_shifted"]),
            ("Brier, label: shifted vs\nconstant class rate", pr["brier_rate_minus_shifted"], pr["ci95"])]
    for i, (lab, est, ci) in enumerate(rows):
        dx.plot(ci, [i, i], color=verm, lw=1.1); dx.scatter([est], [i], s=14, color=verm, zorder=3)
    dx.axvline(0, color="k", lw=0.5)
    dx.set_yticks(range(len(rows))); dx.set_yticklabels([r[0] for r in rows], fontsize=5.4)
    dx.set_ylim(-0.6, len(rows) - 0.4); dx.invert_yaxis()
    dx.set_xlabel("score gain of the shifted null (positive = better)")
    dx.set_title(f"held-out winters: {pr['observed_dw']} downward observed, {pr['expected_dw']:.1f} expected", fontsize=6)
    S.panel_label(dx, "d", x=-0.5, y=1.12)
    S.save(fig, "fig2_selection")


if __name__ == "__main__":
    main()
