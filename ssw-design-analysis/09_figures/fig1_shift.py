#!/usr/bin/env python3
"""
fig1_shift.py -- Figure 1: an SSW shifts the surface distribution; it does not
split it.

Reads results/current/8_experiment/snapsi_causal_effect.json (result K) and
snapsi_distribution_test.json (result M). Recomputes nothing but a histogram.

  a  causal shift per model and initialisation, +-1.96 se, in K's HEADLINE
     units: each model's control spread pooled over its initialisations
     (S_sigma_pooled). K's per-init column S_sigma divides by a lead-dependent
     spread instead; the first draft plotted those points against model means
     in pooled units -- two yardsticks in one panel. Sign: NAM convention,
     negative = downward (-S, since S is a polar-cap pressure rise).
  b  pooled members, each ensemble re-centred on its own mean (as M requires):
     control, and nudged translated by the mean causal shift. Same shape,
     moved -- variance ratio and KS test from M.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nature_style as S  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

INIT_MARK = {"s20180125": "o", "s20180208": "s", "s20181213": "^", "s20190108": "v"}


def main():
    S.apply()
    k = json.loads((S.RESULTS / "8_experiment" / "snapsi_causal_effect.json").read_text())
    m = json.loads((S.RESULTS / "8_experiment" / "snapsi_distribution_test.json").read_text())
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(S.DOUBLE, 64 * S.MM),
                                 gridspec_kw={"width_ratios": [1.2, 1], "wspace": 0.28})

    pm = {c: -v for c, v in k["S_sigma_per_model"].items()}     # NAM sign
    models = sorted(pm, key=lambda c: -pm[c])
    ok = [pm[c] for c in models if c != "ECCC"]
    mu, sd = float(np.mean(ok)), float(np.std(ok, ddof=1))
    ax.axhspan(mu - sd, mu + sd, color="0.92", zorder=0)
    ax.axhline(mu, color="k", lw=0.7)
    ax.text(len(models) - 0.6, mu + sd + 0.05, f"mean excl. ECCC {mu:+.2f}σ (sd {sd:.2f})",
            ha="right", va="bottom", fontsize=5.5)
    for x, c in enumerate(models):
        for r in k["per_case"]:
            if r["centre"] != c:
                continue
            off = {"s20180125": -0.24, "s20180208": -0.08,
                   "s20181213": 0.08, "s20190108": 0.24}[r["init"]]
            col = S.EVENT["Feb-2018"] if r["event"] == "Feb-2018" else S.EVENT["Jan-2019"]
            v = -r["S_sigma_pooled"]
            e = 1.96 * r["se_Pa"] / r["sd_pooled_Pa"]
            ax.errorbar(x + off, v, yerr=e, fmt="none",
                        ecolor=col, elinewidth=0.5, alpha=0.7)
            ax.scatter(x + off, v, s=6, marker=INIT_MARK[r["init"]],
                       color=col, linewidths=0, zorder=3)
        ax.scatter(x, pm[c], s=40, marker="_", color="k", linewidths=1.4, zorder=4)
    ax.set_xticks(range(len(models)))
    ax.set_xticklabels(models, rotation=35, ha="right")
    ax.set_ylabel("causal surface NAM shift\n(σ of the model's control)")
    ax.axhline(0, color="k", lw=0.4)
    ax.annotate("Tukey outlier", xy=(models.index("ECCC"), pm["ECCC"]),
                xytext=(models.index("ECCC") - 2.2, pm["ECCC"] + 0.3), fontsize=5.5,
                arrowprops=dict(arrowstyle="-", lw=0.4))
    for init, mk in INIT_MARK.items():
        ev = "Feb-2018" if init.startswith("s2018") and init != "s20181213" else "Jan-2019"
        ax.scatter([], [], s=6, marker=mk, color=S.EVENT[ev], label=f"{ev}, {init}")
    ax.scatter([], [], s=40, marker="_", color="k", label="model mean")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, columnspacing=0.8,
              handletextpad=0.2, fontsize=5.3)
    S.panel_label(ax, "a", x=-0.12, y=1.12)

    # b: shift without split
    MA = m["member_A"]
    Nc = np.concatenate([np.array(v["nudged"]) - np.mean(v["nudged"]) for v in MA.values()])
    Cc = np.concatenate([np.array(v["control"]) - np.mean(v["control"]) for v in MA.values()])
    pa = m["pooled_all"]
    shift = pa["shift_sigma"]
    bins = np.linspace(-5, 3.5, 60)
    bx.hist(Cc, bins=bins, density=True, histtype="stepfilled", color="0.8",
            label=f"no SSW (control), n={pa['n_control']}")
    bx.hist(Nc + shift, bins=bins, density=True, histtype="step", lw=1.0,
            color=S.ARM["nudged"],
            label=f"SSW imposed (nudged), n={pa['n_nudged']}")
    bx.hist(Cc + shift, bins=bins, density=True, histtype="step", lw=0.7,
            color="k", ls=(0, (2, 1.5)), label="control, translated by the shift")
    bx.axvline(0, color="k", lw=0.4, ls=(0, (2, 2)))
    bx.annotate("", xy=(shift, 0.47), xytext=(0, 0.47),
                arrowprops=dict(arrowstyle="->", lw=0.7))
    bx.text(shift / 2, 0.48, f"shift {shift:+.2f}σ", ha="center", va="bottom", fontsize=5.5)
    lo, hi = pa["variance_ratio_CI95"]
    bx.text(0.98, 0.62, f"variance ratio {pa['variance_ratio']:.3f}\n[{lo:.3f}, {hi:.3f}]\n"
                        f"KS (shape) p = {pa['pure_translation_KS_p']:.2f}",
            transform=bx.transAxes, ha="right", va="top", fontsize=5.5)
    bx.set_xlabel("member surface NAM proxy, days +8 to +25\n(σ of that ensemble's control)")
    bx.set_ylabel("density")
    bx.set_ylim(0, 0.56)
    bx.legend(loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=5.3)
    S.panel_label(bx, "b", x=-0.14, y=1.12)
    S.save(fig, "fig1_shift")


if __name__ == "__main__":
    main()
