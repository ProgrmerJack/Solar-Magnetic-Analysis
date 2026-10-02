#!/usr/bin/env python3
"""
figED9_forecasts.py -- Extended Data Fig. 9: what transfers to forecasts.

Reads results/current/6_predictability/state_forecast.json and shift_rule_forecast.json.
Recomputes nothing (skills and intervals are those stored by the producers).

  a  57 ERA5 SSWs, 1940-2026, whole-winter cross-validation: CRPS skill over
     climatology of the fixed SNAPSI rule, the constant CMIP6 x ERA5 rule and the
     state model (95% intervals, winter resamples); inset numbers, the state model
     with an SSW indicator and with the realised future wind (oracle), relative to it
  b  cold 17-day periods: observed count and the count each forecast expected
  c  17 SSWs of 1998-2021, forecasts issued 0-7 days after onset: per-event CRPS of
     the calibrated and raw dynamical ensembles, the state model and climatology
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
    d = json.loads((S.RESULTS / "6_predictability" / "state_forecast.json").read_text())
    sk = d["E1_crps_skill"]
    fig, (ax, bx, cx) = plt.subplots(1, 3, figsize=(S.DOUBLE, 62 * S.MM), gridspec_kw={"wspace": 0.5,
                                                                                        "width_ratios": [1, 0.9, 1.2]})
    rows = [("fixed SNAPSI rule", sk["F1_vs_F0"]["all"], S.C["vermillion"]),
            ("CMIP6 × ERA5 rule", sk["F2_vs_F0"]["all"], S.C["green"]),
            ("state model", sk["F3_vs_F0"]["all"], S.C["blue"])]
    for i, (lab, r, col) in enumerate(rows):
        ax.plot(r["ci95"], [i, i], color=col, lw=1.2)
        ax.scatter([r["skill"]], [i], s=16, color=col, zorder=3)
    ax.axvline(0, color="k", lw=0.5)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows])
    ax.set_xlabel("CRPS skill over climatology")
    f4, fo = sk["F4_vs_F3"]["all"], sk["F3_oracle_vs_F3"]["all"]
    ax.set_title(f"+SSW indicator {f4['skill']:+.3f} [{f4['ci95'][0]:+.3f}, {f4['ci95'][1]:+.3f}]\n"
                 f"+future wind {fo['skill']:+.3f} [{fo['ci95'][0]:+.3f}, {fo['ci95'][1]:+.3f}]", fontsize=5.4)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    S.panel_label(ax, "a", x=-0.62, y=1.04)

    cc = d["E1_cold_count"]
    labs = ["clim.", "SNAPSI\nrule", "CMIP6\nrule", "state", "state\n+SSW"]
    vals = [cc["expected"][k] for k in ("F0", "F1", "F2", "F3", "F4")]
    cols = [S.C["grey"], S.C["vermillion"], S.C["green"], S.C["blue"], S.C["sky"]]
    bx.bar(range(5), vals, color=cols, width=0.65)
    bx.axhline(cc["observed"], color="k", lw=0.8, ls=(0, (3, 2)))
    bx.text(-0.4, cc["observed"] + 0.5, f"observed: {cc['observed']} of {cc['n']}", ha="left", fontsize=5.2, bbox=dict(facecolor="white", edgecolor="none", pad=0.6), zorder=5)
    bx.set_xticks(range(5)); bx.set_xticklabels(labs, fontsize=4.8)
    bx.set_ylabel("cold 17-day periods expected")
    S.panel_label(bx, "b", x=-0.3, y=1.04)

    ev = d["E2"]["events"]
    order = sorted(range(len(ev)), key=lambda i: ev[i]["onset"])
    for key, col, mk, lab in (("cal", S.C["blue"], "o", "calibrated ensemble"), ("raw", S.C["sky"], "o", "raw ensemble"),
                              ("F3", S.C["orange"], "s", "state model"), ("F0", S.C["grey"], "^", "climatology")):
        y = [ev[i].get(key, np.nan) for i in order]
        cx.scatter(range(len(order)), y, s=8, marker=mk, color=col, label=f"{lab} ({d['E2']['mean_crps'][key]:.2f})", zorder=3)
    cx.set_xticks(range(len(order))); cx.set_xticklabels([ev[i]["onset"][:7] for i in order], rotation=90, fontsize=4.5)
    cx.set_ylabel("CRPS (K)")
    cx.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=4.6, ncol=2, scatterpoints=1)
    S.panel_label(cx, "c", x=-0.25, y=1.16)
    S.save(fig, "figED9_forecasts")


if __name__ == "__main__":
    main()
