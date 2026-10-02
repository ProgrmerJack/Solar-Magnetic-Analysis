#!/usr/bin/env python3
"""
fig1_summary.py -- Figure 1: an imposed SSW changes the outcome rate, not the class contrast.

Reads (recomputes nothing but a histogram and per-pair differences of stored values):
  results/current/8_experiment/snapsi_distribution_test.json   member outcomes (M)
  results/current/8_experiment/snapsi_selection_test.json       DW rates and contrasts (L)
  results/current/8_experiment/snapsi_causal_effect.json        causal shift per model (K)

  a  SNAPSI, all 36 NH ensembles: control members (grey) and SSW-imposed members,
     each re-centred on its own ensemble mean and placed at the mean shift
     (vermilion); the threshold of the downward criterion at zero
  b  share of members classified DW in each of the 36 centre x initialisation
     pairs, without and with the imposed SSW (lines join a pair); open markers, the
     11 pairs whose SSW-imposed ensemble has fewer than 3 NDW members, so that no
     class contrast can be formed there
  c  DW - NDW contrast of the 25 pairs in which both arms form one, no SSW against
     SSW imposed; dotted, the contrast a cut at zero gives on one Gaussian
  d  effect of imposing the SSW, in sigma of the control: on the mean polar-cap NAM
     (per model, and the mean over models) and on the DW - NDW contrast (per pair,
     and the paired mean with its centre-bootstrap interval)
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
    R = S.RESULTS
    dist = json.loads((R / "8_experiment" / "snapsi_distribution_test.json").read_text())
    sel = json.loads((R / "8_experiment" / "snapsi_selection_test.json").read_text())
    ce = json.loads((R / "8_experiment" / "snapsi_causal_effect.json").read_text())
    grey, verm = "0.6", S.C["vermillion"]
    fig = plt.figure(figsize=(S.DOUBLE, 112 * S.MM))
    gs = fig.add_gridspec(2, 2, wspace=0.34, hspace=0.62)

    # a -- one population, shifted, cut at zero
    ax = fig.add_subplot(gs[0, 0])
    MA = dist["member_A"]
    ctl = np.concatenate([np.array(v["control"]) - np.mean(v["control"]) for v in MA.values()])
    nud = np.concatenate([np.array(v["nudged"]) - np.mean(v["nudged"]) for v in MA.values()])
    shift = float(np.mean([np.mean(v["nudged"]) - np.mean(v["control"]) for v in MA.values()]))
    bins = np.linspace(-5, 4, 55)
    ax.hist(ctl, bins=bins, density=True, color=grey, alpha=0.55, label="no SSW (control)")
    ax.hist(nud + shift, bins=bins, density=True, histtype="step", color=verm, lw=1.2, label="SSW imposed (nudged)")
    ax.axvline(0, color="k", lw=0.8)
    vr = dist["pooled_all"]["variance_ratio"]; ci = dist["pooled_all"]["variance_ratio_CI95"]
    ax.text(1.4, 0.22, f"variance ratio\n{vr:.2f} [{ci[0]:.2f}, {ci[1]:.2f}]", fontsize=5.6, va="top")
    ax.text(0.08, 0.36, "DW | NDW", fontsize=5.4, ha="center")
    ax.set_xlabel("surface NAM proxy, days 8-25 (σ of control)"); ax.set_ylabel("density")
    ax.legend(fontsize=5.2, loc="upper right", bbox_to_anchor=(1.02, 1.0))
    ax.set_title("the imposed SSW moves the distribution across the class threshold", fontsize=6)
    S.panel_label(ax, "a", x=-0.14, y=1.04)

    # b -- the rate, all 36 pairs
    bx = fig.add_subplot(gs[0, 1])
    ens = {(e["centre"], e["init"], e["arm"]): e for e in sel["all_ensembles"] if e["hemisphere"] == "NH"}
    keys = sorted({(c, i) for c, i, _ in ens})
    n_open = 0
    for c, i in keys:
        n_, c_ = ens[(c, i, "nudged")], ens[(c, i, "control")]
        est = n_["contrast_estimable"] and c_["contrast_estimable"]
        bx.plot([0, 1], [c_["DW_rate"], n_["DW_rate"]], color="0.8", lw=0.4, zorder=1)
        bx.scatter([0, 1], [c_["DW_rate"], n_["DW_rate"]], s=9, zorder=2,
                   facecolor=([grey, verm] if est else "white"), edgecolor=[grey, verm], linewidths=0.6)
        n_open += not est
    sm = sel["summary"]
    for x, arm in ((0, "control"), (1, "nudged")):
        m = sm[arm]["mean_DW_rate_UNCONDITIONAL"]
        bx.plot([x - 0.22, x + 0.22], [m, m], color="k", lw=1.2)
        bx.text(x + (0.26 if x else -0.26), m, f"{m:.0%}", ha="left" if x else "right", va="center", fontsize=6)
    bx.text(0.5, 0.06, f"open: {n_open} of {len(keys)} pairs where the SSW-imposed ensemble has\n"
                       f"< 3 NDW members (DW share {sm['nudged']['dropped_DW_rate_range'][0]:.2f}–"
                       f"{sm['nudged']['dropped_DW_rate_range'][1]:.2f}), so no contrast exists",
            ha="center", fontsize=5.0, transform=bx.transAxes)
    bx.set_xticks([0, 1]); bx.set_xticklabels(["no SSW", "SSW imposed"])
    bx.set_xlim(-0.6, 1.6); bx.set_ylim(-0.2, 1.05)
    bx.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    bx.set_ylabel("share of members classified DW")
    bx.set_title(f"the DW rate responds in all {len(keys)} pairs", fontsize=6)
    S.panel_label(bx, "b", x=-0.16, y=1.04)

    # c -- the contrast, 25 estimable pairs
    cx = fig.add_subplot(gs[1, 0])
    rows = {}
    for r in sel["per_case"]:
        if r["hemisphere"] == "NH":
            rows.setdefault((r["centre"], r["init"]), {})[r["arm"]] = r["contrast_sigma"]
    pairs = {k: v for k, v in rows.items() if {"nudged", "control"} <= set(v)}
    xs = np.array([v["control"] for v in pairs.values()]); ys = np.array([v["nudged"] for v in pairs.values()])
    lim = (-2.45, -0.95)
    cx.plot(lim, lim, color="k", lw=0.6, ls=(0, (3, 2)))
    cx.axvline(REF, color=S.C["green"], lw=0.6, ls=":"); cx.axhline(REF, color=S.C["green"], lw=0.6, ls=":")
    cx.scatter(xs, ys, s=10, color=verm, edgecolor="white", linewidths=0.3, zorder=3)
    pr = sel["paired_NH"]
    cx.text(0.03, 0.97, f"{pr['n_pairs']} of {len(keys)} pairs, {pr['n_centres']} models", transform=cx.transAxes,
            va="top", fontsize=5.8)
    cx.text(lim[0] + 0.03, REF - 0.03, "dotted: one Gaussian\ncut at zero", color=S.C["green"], fontsize=5.2, va="top")
    cx.set_xlim(lim); cx.set_ylim(lim); cx.set_aspect("equal")
    cx.set_xlabel("DW − NDW contrast, no SSW (σ)"); cx.set_ylabel("DW − NDW contrast, SSW imposed (σ)")
    cx.set_title("the class contrast does not", fontsize=6)
    S.panel_label(cx, "c", x=-0.2, y=1.04)

    # d -- effect of the intervention on the mean and on the contrast, same units
    dx = fig.add_subplot(gs[1, 1])
    mod = ce["S_sigma_per_model"]
    rng = np.random.default_rng(0)                      # vertical jitter only
    mvals = -np.array(list(mod.values()))               # NAM sign: negative = downward
    dx.scatter(mvals, 1 + rng.uniform(-0.12, 0.12, len(mvals)), s=8, color=grey, zorder=2)
    ec = -mod["ECCC"]
    dx.text(ec, 0.78, "ECCC", ha="center", fontsize=5.2, color="0.4")
    gm = -ce["S_sigma_grand_mean"]
    dx.plot([gm, gm], [0.75, 1.25], color=verm, lw=1.6)
    dx.text(gm, 1.33, f"{gm:+.2f}σ".replace("-", "−"), ha="center", fontsize=5.8, color=verm)
    d_ = ys - xs
    dx.scatter(d_, 0 + rng.uniform(-0.12, 0.12, len(d_)), s=8, color=grey, zorder=2)
    lo, hi = pr["paired_difference_CI95_centre_bootstrap"]
    dx.plot([lo, hi], [0, 0], color=verm, lw=1.0, zorder=3)
    dx.plot([pr["paired_difference"]] * 2, [-0.25, 0.25], color=verm, lw=1.6, zorder=3)
    dx.text(pr["paired_difference"], -0.45, f"{pr['paired_difference']:+.2f}σ [{lo:+.2f}, {hi:+.2f}]".replace("-", "−"),
            ha="center", fontsize=5.8, color=verm)
    dx.axvline(0, color="k", lw=0.5)
    dx.set_yticks([0, 1]); dx.set_yticklabels(["DW − NDW\ncontrast\n(25 pairs)", "mean NAM\n(9 models)"])
    dx.set_ylim(-0.7, 1.55); dx.set_xlim(-3.7, 0.9)
    dx.set_xlabel("effect of imposing the SSW (σ of control)")
    dx.set_title("what the forcing changes, in one unit", fontsize=6)
    S.panel_label(dx, "d", x=-0.3, y=1.04)
    S.save(fig, "fig1_summary")


if __name__ == "__main__":
    main()
