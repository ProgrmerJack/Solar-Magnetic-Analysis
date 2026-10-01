#!/usr/bin/env python3
"""
figED7_shiftrule.py -- Extended Data Fig. 7: the shift rule as a forecast.

Reads results/current/6_predictability/shift_rule_forecast.json (result V).
Recomputes nothing.

  a  observed frequency after the observed SSWs (Wilson 95% interval) against the
     probability taken from the models alone: cold northern-Eurasian fortnight
     below the event-free 5th, 10th, 20th percentile (SNAPSI, 39 events), and a
     negative annular mode, days 8-52 (CMIP6, CPC AO, 43 events); open symbols,
     climatology
  b  relative economic value of the rule against climatology, for users with
     cost/loss ratio alpha, cold-fortnight thresholds: positive between the
     climatological rate and the observed rate after SSWs, negative between that
     and the rule probability (the rule slightly over-states the risk), zero
     elsewhere (rule and climatology make the same decision)
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nature_style as S  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402


def main():
    S.apply()
    d = json.loads((S.RESULTS / "6_predictability" / "shift_rule_forecast.json").read_text())
    v1, ao = d["verification"]["V1"], d["verification"]["AO_negative"]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(S.DOUBLE * 0.72, 62 * S.MM),
                                 gridspec_kw={"wspace": 0.45})
    cols = {"0.05": S.C["blue"], "0.1": S.C["sky"], "0.2": S.C["green"]}
    rows = [(f"cold, below {int(float(q) * 100)}th pct.", v1[q], cols[q]) for q in ("0.05", "0.1", "0.2")]
    rows.append(("negative annular mode", ao, S.C["vermillion"]))
    for lab, r, col in rows:
        lo, hi = r["f_wilson95"]
        ax.plot([r["p_rule"]] * 2, [lo, hi], color=col, lw=0.8)
        ax.scatter([r["p_rule"]], [r["f"]], s=14, color=col, zorder=3,
                   label=f"{lab} ({r['count']}/{r['n']})")
        ax.scatter([r["q"]], [r["f"]], s=14, facecolor="white", edgecolor=col, zorder=3, linewidths=0.6)
        ax.plot([r["q"], r["p_rule"]], [r["f"]] * 2, color=col, lw=0.4, ls=(0, (1, 1)))
    ax.plot((0, 1), (0, 1), color="k", lw=0.5, ls=(0, (2, 2)))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    ax.set_xlabel("probability from the models alone\n(open: climatology)")
    ax.set_ylabel("observed frequency after SSWs")
    ax.legend(loc="upper left", fontsize=4.8)
    S.panel_label(ax, "a", x=-0.3, y=1.04)

    for q in ("0.05", "0.1", "0.2"):
        rc = v1[q]["rev_curve"]
        a = [float(k) for k in rc]
        # evaluated on a 0.02 grid; the curve is discontinuous at the climatological
        # rate and at the rule probability, so points, not a line
        bx.scatter(a, [rc[k] for k in rc], s=4, color=cols[q],
                   label=f"below {int(float(q) * 100)}th pct.")
    bx.axhline(0, color="k", lw=0.5)
    bx.set_xlim(0, 0.6); bx.set_ylim(-0.45, 1.0)
    bx.set_xlabel("user cost/loss ratio")
    bx.set_ylabel("relative economic value\n(1 = perfect forecast)")
    bx.legend(loc="upper right", fontsize=5.0)
    S.panel_label(bx, "b", x=-0.3, y=1.04)
    S.save(fig, "figED7_shiftrule")


if __name__ == "__main__":
    main()
