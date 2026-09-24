#!/usr/bin/env python3
"""
fig4_archetypes.py -- Figure 4: the field's two archetypes are two draws.

Reads results/current/8_experiment/snapsi_archetype_test.json (result O) and
snapsi_causal_effect.json (result K). Recomputes nothing.

  a  every model's nudged members for Feb-2018 (s20180125) and Jan-2019
     (s20181213), the primary pair with comparable lead to onset, as the
     NAM proxy A (negative = downward-propagating sign); observed ERA5 marked
  b  percentile of the observed A_2018 - A_2019 among all member pairs
  c  Jan-2019 and Feb-2018 under four variants of the Karpechko criterion
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nature_style as S  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

PAIR = ("s20180125", "s20181213")
EV = {"s20180125": "Feb-2018", "s20181213": "Jan-2019"}


def main():
    S.apply()
    o = json.loads((S.RESULTS / "8_experiment" / "snapsi_archetype_test.json").read_text())
    k = json.loads((S.RESULTS / "8_experiment" / "snapsi_causal_effect.json").read_text())
    per = {(r["centre"], r["init"]): r for r in o["per_ensemble"]}
    centres = [c for c in o["centres"] if all((c, i) in per for i in PAIR)]

    fig = plt.figure(figsize=(S.DOUBLE, 92 * S.MM))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.55, 1], height_ratios=[1, 1],
                          wspace=0.32, hspace=0.55)
    ax = fig.add_subplot(gs[:, 0])
    rng = np.random.default_rng(0)            # jitter only; no statistic depends on it
    for yi, c in enumerate(centres[::-1]):
        for dj, init in ((0.18, PAIR[0]), (-0.18, PAIR[1])):
            a = np.array(o["member_A"][f"{c}|{init}"]["nudged"])
            y = yi + dj
            col = S.EVENT[EV[init]]
            ax.scatter(a, y + rng.uniform(-0.07, 0.07, len(a)), s=1.2, color=col,
                       alpha=0.45, linewidths=0, zorder=2)
            q = np.percentile(a, [25, 50, 75])
            ax.plot([q[0], q[2]], [y, y], color=col, lw=2.2, alpha=0.35,
                    solid_capstyle="butt", zorder=1)
            obs = per[(c, init)]["obs_A"]
            adj = per[(c, init)]["obs_A_adjusted"]
            ax.scatter([obs], [y], marker="D", s=14, facecolor=col,
                       edgecolor="k", linewidths=0.5, zorder=3)
            if abs(adj - obs) > 0.25:     # only where the model-ERA5 offset is visible
                ax.scatter([adj], [y], marker="D", s=14, facecolor="white",
                           edgecolor=col, linewidths=0.8, zorder=3)
    ax.axvline(0, color="k", lw=0.5, ls=(0, (2, 2)))
    ax.set_yticks(range(len(centres)))
    ax.set_yticklabels(centres[::-1])
    ax.set_xlabel("surface NAM proxy over days +8 to +25 (σ of control)")
    ax.set_xlim(-7.5, 3)
    ax.text(-0.25, len(centres) - 0.35, "DW  ←", ha="right", fontsize=6)
    ax.text(0.25, len(centres) - 0.35, "→  NDW", ha="left", fontsize=6)
    # Forced shift (K) and DW odds (O) AT THE PLOTTED INITIALISATIONS, so the
    # legend describes exactly what the panel shows.
    sk = k["S_sigma_per_init"]
    for init in PAIR:
        e, col = EV[init], S.EVENT[EV[init]]
        dw = np.mean([r["nudged_DW_rate"] for r in o["per_ensemble"] if r["init"] == init])
        ax.scatter([], [], s=12, marker="D", facecolor=col, edgecolor="k",
                   linewidths=0.5,
                   label=f"{e} ({init}): forced shift {sk[f'{e}|{init}']:+.2f}σ, "
                         f"DW odds {dw:.2f}")
    ax.scatter([], [], s=4, color="0.5", label="nudged members (bar: interquartile range)")
    ax.scatter([], [], s=12, marker="D", facecolor="0.5", edgecolor="k", linewidths=0.5,
               label="observed, ERA5")
    ax.scatter([], [], s=12, marker="D", facecolor="white", edgecolor="0.5",
               linewidths=0.8, label="observed, model-ERA5 initial offset removed")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, handletextpad=0.3,
              columnspacing=1.0, fontsize=5.5)
    S.panel_label(ax, "a", x=-0.2)

    # b: pair test
    bx = fig.add_subplot(gs[0, 1])
    pt = {r["centre"]: r for r in o["pair_test"]["primary"]["per_centre"]}
    ys = np.arange(len(centres))
    bx.axvspan(0.025, 0.975, color="0.92", zorder=0)
    bx.scatter([pt[c]["obs_diff_pct"] for c in centres[::-1]], ys, s=10,
               color="k", zorder=2)
    bx.set_yticks(ys)
    bx.set_yticklabels(centres[::-1], fontsize=5.5)
    bx.set_xlim(0, 1)
    bx.set_xlabel("percentile of observed $A_{2018}-A_{2019}$\namong all member pairs")
    bx.set_title(f"inside the central 95% in "
                 f"{len(centres) - o['pair_test']['primary']['centres_outside_central_95']}"
                 f" of {len(centres)} models", fontsize=6.5)
    S.panel_label(bx, "b", x=-0.33)

    # c: the label under four criterion variants
    cx = fig.add_subplot(gs[1, 1])
    labs = o["observed_labels"]
    rows = []
    for level in ("1000hPa", "850hPa"):
        for win in ([8, 52], [8, 25]):
            rows.append((level, win))
    for j, (level, win) in enumerate(rows[::-1]):
        for dy, ev in ((0.16, "Feb-2018"), (-0.16, "Jan-2019")):
            r = next(x for x in labs if x["event"] == ev and x["level"] == level
                     and x["window"] == win)
            cx.barh(j + dy, r["nam_mean"], height=0.28, color=S.EVENT[ev],
                    alpha=0.9 if r["class"] == "DW" else 0.35,
                    edgecolor="k" if r["class"] == "NDW" else "none", linewidth=0.6)
            if r["class"] == "NDW":
                cx.text(r["nam_mean"] + 0.04, j + dy, f"NDW ({r['nam_mean']:+.3f}σ)",
                        va="center", fontsize=5.5)
    cx.axvline(0, color="k", lw=0.5, ls=(0, (2, 2)))
    cx.set_yticks(range(len(rows)))
    cx.set_yticklabels([f"{lv[:-3]} hPa, d{w[0]}–{w[1]}" for lv, w in rows[::-1]],
                       fontsize=5.5)
    cx.set_xlabel("observed NAM mean (σ); criterion needs < 0")
    cx.set_xlim(-1.35, 0.6)
    cx.set_title("Jan-2019 is NDW in 1 of 4 variants", fontsize=6.5)
    S.panel_label(cx, "c", x=-0.33)
    S.save(fig, "fig4_archetypes")


if __name__ == "__main__":
    main()
