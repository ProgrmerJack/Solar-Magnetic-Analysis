#!/usr/bin/env python3
"""
fig5_regional.py -- Figure 5: regional cold risk after an imposed SSW is the
circulation shift acting through the ordinary relation.

Reads results/current/8_experiment/snapsi_regional_test.json (result R, SNAPSI)
and results/current/2_event_study/era5_regional_test.json (result R, ERA5).
Recomputes nothing.

  a  SNAPSI: the SSW's effect on regional temperature (days +8..+24, sigma of the
     control), split into the part the member's polar-cap NAM implies through the
     SSW-free relation and the residual R (point, 95% centre-bootstrap interval)
  b  SNAPSI: probability of a cold fortnight (below the control 10th percentile):
     control (0.10 by construction), with the SSW imposed, and as predicted from
     the circulation shift alone
  c  ERA5, 39 observed SSWs: the same split in K; whisker on the residual, the
     95% range of the residual on event-free dates of the same calendar period
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nature_style as S  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

REG = [("NEURASIA", "N Eurasia"), ("HI_EUROPE", "high-lat.\nEurope"),
       ("MID_EASIA", "mid-lat.\nE Asia"), ("MID_NAMER", "mid-lat.\nN America")]


def main():
    S.apply()
    sn = json.loads((S.RESULTS / "8_experiment" / "snapsi_regional_test.json").read_text())["results"]
    er = json.loads((S.RESULTS / "2_event_study" / "era5_regional_test.json").read_text())["regions"]
    fig = plt.figure(figsize=(S.DOUBLE, 62 * S.MM))
    gs = fig.add_gridspec(1, 3, wspace=0.42)
    x = np.arange(len(REG))
    grey, verm = "0.72", S.C["vermillion"]

    ax = fig.add_subplot(gs[0, 0])
    w = 0.26
    for i, (k, _) in enumerate(REG):
        st, r = sn[k]["S_T"], sn[k]["R"]
        ax.bar(i - w, st["mean"], w, color=S.C["blue"], label="SSW effect" if i == 0 else None)
        ax.errorbar(i - w, st["mean"], yerr=[[st["mean"] - st["ci95"][0]], [st["ci95"][1] - st["mean"]]],
                    color="k", lw=0.6, capsize=1.2)
        ax.bar(i, st["mean"] - r["mean"], w, color=grey, label="implied by NAM shift" if i == 0 else None)
        ax.bar(i + w, r["mean"], w, color=verm, label="residual" if i == 0 else None)
        ax.errorbar(i + w, r["mean"], yerr=[[r["mean"] - r["ci95"][0]], [r["ci95"][1] - r["mean"]]],
                    color="k", lw=0.6, capsize=1.2)
    ax.axhline(0, color="k", lw=0.4)
    ax.set_xticks(x); ax.set_xticklabels([l for _, l in REG], fontsize=5.3)
    ax.set_ylabel("temperature, days 8-24 (σ of control)")
    ax.legend(fontsize=5, loc="upper left", ncol=1, bbox_to_anchor=(0.0, 1.0))
    ax.set_ylim(-1.55, 1.25)
    ax.set_title("SNAPSI: imposed SSW", fontsize=6)
    S.panel_label(ax, "a", x=-0.28, y=1.04)

    bx = fig.add_subplot(gs[0, 1])
    for i, (k, _) in enumerate(REG):
        c, o, pr = sn[k]["cold_control"]["mean"], sn[k]["cold_obs"]["mean"], sn[k]["cold_pred"]["mean"]
        bx.bar(i - 0.25, c, 0.24, color="0.85", label="no SSW" if i == 0 else None)
        bx.bar(i, o, 0.24, color=S.C["blue"], label="SSW imposed" if i == 0 else None)
        bx.bar(i + 0.25, pr, 0.24, color=grey, hatch="///", edgecolor="0.4", lw=0.3,
               label="predicted from shift" if i == 0 else None)
    bx.set_xticks(x); bx.set_xticklabels([l for _, l in REG], fontsize=5.3)
    bx.set_ylabel("P(cold fortnight)")
    bx.legend(fontsize=5, loc="upper right")
    bx.set_title("SNAPSI: cold risk", fontsize=6)
    S.panel_label(bx, "b", x=-0.25, y=1.04)

    cx = fig.add_subplot(gs[0, 2])
    for i, (k, _) in enumerate(REG):
        e = er[k]
        cx.bar(i - w, e["S_T_mean"], w, color=S.C["blue"])
        cx.bar(i, e["S_T_mean"] - e["R_mean"], w, color=grey)
        cx.bar(i + w, e["R_mean"], w, color=verm)
        lo, hi = e["R_null_q025_q975"]
        cx.plot([i + w, i + w], [lo, hi], color="0.35", lw=0.6)
        cx.plot([i + w - 0.07, i + w + 0.07], [lo, lo], color="0.35", lw=0.6)
        cx.plot([i + w - 0.07, i + w + 0.07], [hi, hi], color="0.35", lw=0.6)
    cx.axhline(0, color="k", lw=0.4)
    cx.set_xticks(x); cx.set_xticklabels([l for _, l in REG], fontsize=5.3)
    cx.set_ylabel("temperature anomaly (K)")
    n = json.loads((S.RESULTS / "2_event_study" / "era5_regional_test.json").read_text())["n_events"]
    cx.set_title(f"ERA5: {n} observed SSWs", fontsize=6)
    S.panel_label(cx, "c", x=-0.25, y=1.04)
    S.save(fig, "fig5_regional")


if __name__ == "__main__":
    main()
