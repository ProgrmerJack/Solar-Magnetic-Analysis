#!/usr/bin/env python3
"""
figED4_sweep.py -- Extended Data Fig. 4: every version of the criterion tracks one
shifted population.

Reads results/current/9_literature/criterion_sweep.json (result T-a).
Recomputes nothing.

  a  ERA5, 108 versions: downward rate after the 39 SSWs against the rate of
     event-free dates displaced by that version's measured shift (mean and
     central 95% of 400 sets); open symbols with the 150 hPa condition; star, the
     published criterion
  b  ERA5: share of the SSW class contrast reproduced by event-free dates, per
     version, sorted
  c  CMIP6, 27 versions (conditions 1-2), 1,517 events
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
    d = json.loads((S.RESULTS / "9_literature" / "criterion_sweep.json").read_text())
    e, c = d["era5"], d["cmip6"]
    fig, (ax, bx, cx) = plt.subplots(1, 3, figsize=(S.DOUBLE, 62 * S.MM),
                                     gridspec_kw={"wspace": 0.45})
    for r in e["specs"]:
        lo, hi = r["rate_shifted_q025_q975"]
        col = S.C["vermillion"] if not r["cond3"] else S.C["blue"]
        ax.plot([lo, hi], [r["rate_ssw"]] * 2, color=col, lw=0.4, alpha=0.5)
        ax.scatter([r["rate_shifted_mean"]], [r["rate_ssw"]], s=6, zorder=3,
                   facecolor=col if not r["cond3"] else "white", edgecolor=col, linewidths=0.5)
    pub = e["published_spec"]
    ax.scatter([pub["rate_shifted_mean"]], [pub["rate_ssw"]], s=40, marker="*", color="k", zorder=4,
               label="published criterion")
    ax.scatter([], [], s=6, color=S.C["vermillion"], label="surface conditions (54)")
    ax.scatter([], [], s=6, facecolor="white", edgecolor=S.C["blue"], label="+ 150 hPa condition (54)")
    lim = (0.0, 1.0)
    ax.plot(lim, lim, color="k", lw=0.5, ls=(0, (2, 2)))
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
    ax.set_xlabel("downward rate, shifted event-free dates")
    ax.set_ylabel("downward rate, 39 observed SSWs")
    s1, s3 = e["S1_S2_conditions_1_2"], e["S3_with_condition_3"]
    cal = e["calibration_posthoc"]
    ax.set_title(f"ERA5: calibrated p = {cal['c12']['p_max_abs_z_calibrated']:.2f} "
                 f"({cal['c3']['p_max_abs_z_calibrated']:.2f})", fontsize=6.3)
    ax.legend(loc="upper left", fontsize=5.0)
    S.panel_label(ax, "a", x=-0.28, y=1.04)

    for j, c3 in enumerate((False, True)):
        v = np.sort([r["share_reproduced_pct"] for r in e["specs"]
                     if r["cond3"] == c3 and r["share_reproduced_pct"] is not None])
        col = S.C["vermillion"] if not c3 else S.C["blue"]
        bx.scatter(np.arange(len(v)), v, s=5, facecolor=col if not c3 else "white", edgecolor=col,
                   linewidths=0.5)
    bx.axhline(100, color="k", lw=0.5, ls=(0, (2, 2)))
    bx.set_xlabel("version (sorted)")
    bx.set_ylabel("class contrast reproduced by\nevent-free dates (%)")
    bx.set_ylim(0, 170)
    bx.set_title("ERA5", fontsize=6.3)
    S.panel_label(bx, "b", x=-0.3, y=1.04)

    for r in c["specs"]:
        cx.scatter([r["rate_shifted_mean"]], [r["rate_ssw"]], s=7, color=S.C["green"], zorder=3)
    lim = (0.3, 0.8)
    cx.plot(lim, lim, color="k", lw=0.5, ls=(0, (2, 2)))
    cx.set_xlim(lim); cx.set_ylim(lim); cx.set_aspect("equal")
    cx.set_xlabel("downward rate, shifted event-free dates")
    cx.set_ylabel(f"downward rate, {c['n_events']:,} CMIP6 SSWs")
    cx.set_title(f"CMIP6: max|z| p = {c['S1_S2']['p_max_abs_z']:.2f}", fontsize=6.3)
    S.panel_label(cx, "c", x=-0.3, y=1.04)
    S.save(fig, "figED4_sweep")


if __name__ == "__main__":
    main()
