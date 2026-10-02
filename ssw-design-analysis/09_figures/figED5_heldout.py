#!/usr/bin/env python3
"""
figED5_heldout.py -- Extended Data Fig. 5: the held-out SSWs, 2023-2026.

Reads results/current/6_predictability/s2s_heldout_diagnosis.json (exploratory,
post hoc) and s2s_heldout_test.json. Recomputes nothing.

  a  every event (17 of 1998-2021, 3 held-out in 2023-24, and 4 March 2026 scored
     with real-time forecasts, s2s_heldout2026_test.json): multi-system mean rank of the
     observed polar-cap outcome against the observed response (negative = downward,
     negative NAM); the rank follows the realized outcome
  b  ECMWF model years 2022 and 2025 on the same ten 1998-2021 events
  c  mean rank on dates more than 30 days from every SSW, by winter, for the three
     held-out systems
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nature_style as S  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402


def main():
    S.apply()
    d = json.loads((S.RESULTS / "6_predictability" / "s2s_heldout_diagnosis.json").read_text())
    fig, (ax, bx, cx) = plt.subplots(1, 3, figsize=(S.DOUBLE, 62 * S.MM),
                                     gridspec_kw={"wspace": 0.45, "width_ratios": [1, 0.8, 1.2]})
    for r in d["D4"]["table"]:
        ho = r["set"] == "held-out"
        ax.scatter(r["A_obs"], r["mean_rank"], s=14 if ho else 9, zorder=3,
                   facecolor=S.C["vermillion"] if ho else S.C["grey"],
                   edgecolor=S.C["vermillion"] if ho else S.C["grey"])
        if ho:
            off = {"2023-02": (4, 4), "2024-01": (4, 4), "2024-03": (-30, -9)}[r["event"][:7]]
            ax.annotate(r["event"][:7], (r["A_obs"], r["mean_rank"]), fontsize=5,
                        xytext=off, textcoords="offset points")
    h26 = json.loads((S.RESULTS / "6_predictability" / "s2s_heldout2026_test.json").read_text())["HO2026_polar_cap"]
    ax.scatter(h26["A_obs"], h26["mean_rank"], s=22, marker="D", facecolor="white", edgecolor=S.C["vermillion"],
               linewidths=0.8, zorder=4)
    ax.annotate("2026-03", (h26["A_obs"], h26["mean_rank"]), fontsize=5, xytext=(4, -9), textcoords="offset points")
    ax.axhline(0.5, color="k", lw=0.4); ax.axvline(0, color="0.7", lw=0.4)
    ax.scatter([], [], s=9, color=S.C["grey"], label="1998–2021 (17)")
    ax.scatter([], [], s=14, color=S.C["vermillion"], label="held out 2023–24 (3)")
    ax.scatter([], [], s=22, marker="D", facecolor="white", edgecolor=S.C["vermillion"], label="2026, real time")
    ax.legend(loc="lower right", fontsize=5.2)
    ax.set_xlabel("observed response (Pa; negative = downward)")
    ax.set_ylabel("mean rank in the ensembles")
    ax.set_ylim(0, 1)
    S.panel_label(ax, "a", x=-0.28, y=1.04)

    same = d["D2"]["ecmwf_same_events"]["per_event"]
    for k, (r22, r25) in same.items():
        bx.plot([0, 1], [r22, r25], color="0.6", lw=0.5)
        bx.scatter([0, 1], [r22, r25], s=7, color=S.C["blue"], zorder=3)
    m = d["D2"]["ecmwf_same_events"]
    bx.plot([0, 1], [m["mean_rank_2022"], m["mean_rank_2025"]], color=S.C["vermillion"], lw=1.4)
    bx.axhline(0.5, color="k", lw=0.4)
    bx.set_xticks([0, 1]); bx.set_xticklabels(["ECMWF\nmodel year 2022", "ECMWF\nmodel year 2025"])
    bx.set_xlim(-0.3, 1.3); bx.set_ylim(0, 1)
    bx.set_ylabel("rank, same 1998–2021 events")
    S.panel_label(bx, "b", x=-0.35, y=1.04)

    cols = {"cnrm": S.C["green"], "ecmwf2025": S.C["blue"], "cma2025": S.C["orange"]}
    for k, col in cols.items():
        w = d["D3"][k]["dates_30d_from_onsets"]
        yrs = sorted(int(y) for y in w)
        v = [w[str(y)][0] for y in yrs]
        cx.plot(yrs, v, color=col, lw=0.8, marker="o", ms=2, label=k.replace("2025", " 2025").upper())
    cx.axhline(0.5, color="k", lw=0.4)
    cx.set_xlabel("winter (ending)"); cx.set_ylabel("mean rank, dates > 30 d from SSWs")
    cx.set_ylim(0, 1)
    cx.legend(loc="lower left", fontsize=5.0)
    S.panel_label(cx, "c", x=-0.22, y=1.04)
    S.save(fig, "figED5_heldout")


if __name__ == "__main__":
    main()
