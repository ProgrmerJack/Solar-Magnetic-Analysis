#!/usr/bin/env python3
"""
figED6_continuity.py -- Extended Data Fig. 6: the SSW is itself a threshold on a
continuum of vortex weakenings.

Reads results/current/5_mechanism/vortex_threshold_continuity.json (result U).
Recomputes nothing: the panels show the stored per-episode values and the stored
estimates.

  a  114 NCEP vortex-weakening episodes 1958-2025: CPC Arctic Oscillation, days
     8-52, against the minimum u(10 hPa, 60N); filled, wind reversed (SSW)
  b  the same for northern-Eurasian 2 m temperature, days 8-52
  c  global piecewise-linear estimates: response per 10 m/s stronger minimum wind
     (dose; positive = a weaker vortex gives a lower AO / colder Eurasia) and step
     at reversal (jump; SSW side minus the continuation of the trend), with 95% winter-cluster bootstrap
     intervals, in units of each outcome's standard deviation across episodes;
     CMIP6 (5,726 episodes) in sigma of its own outcome
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
    d = json.loads((S.RESULTS / "5_mechanism" / "vortex_threshold_continuity.json").read_text())
    U = d["unit_table"]
    ob = d["observations"]
    fig, axs = plt.subplots(1, 3, figsize=(S.DOUBLE, 60 * S.MM), gridspec_kw={"wspace": 0.5})
    for ax, key, lab, pl in ((axs[0], "Y1_AO", "Arctic Oscillation, days 8–52", "a"),
                             (axs[1], "Y3_T_NEURASIA", "N-Eurasian T anomaly, days 8–52 (K)", "b")):
        x = np.array([u["X"] for u in U], float)
        y = np.array([np.nan if u[key] is None else u[key] for u in U], float)
        rev = x < 0
        ax.scatter(x[~rev], y[~rev], s=6, facecolor="white", edgecolor=S.C["grey"], linewidths=0.5,
                   label=f"no reversal ({np.sum(~rev & np.isfinite(y))})")
        ax.scatter(x[rev], y[rev], s=6, color=S.C["vermillion"], label=f"SSW ({np.sum(rev & np.isfinite(y))})")
        ax.axvline(0, color="k", lw=0.5, ls=(0, (2, 2)))
        ax.axhline(0, color="k", lw=0.3)
        e = ob[key]["E1_E2"]
        ax.set_title(f"dose {e['dose_slope_per_10ms']:+.2f} [{e['dose_ci95'][0]:+.2f}, {e['dose_ci95'][1]:+.2f}] per 10 m/s\n"
                     f"step at 0: {e['jump']:+.2f} [{e['jump_ci95'][0]:+.2f}, {e['jump_ci95'][1]:+.2f}]",
                     fontsize=5.6)
        ax.set_xlabel("minimum u(10 hPa, 60° N) (m s$^{-1}$)")
        ax.set_ylabel(lab)
        ax.legend(loc="lower right", fontsize=4.8)
        S.panel_label(ax, pl, x=-0.32, y=1.04)

    cx = axs[2]
    rows = []
    for key, lab in (("Y1_AO", "AO"), ("Y2_NAM1000", "NAM1000"), ("Y3_T_NEURASIA", "N-Eurasia T")):
        vals = np.array([u[key] for u in U if u[key] is not None], float)
        sd = float(np.std(vals, ddof=1))
        rows.append((lab, ob[key]["E1_E2"], sd))
    ypos = np.arange(len(rows) + 1)[::-1]
    for (lab, e, sd), yy in zip(rows, ypos[:-1]):
        cx.plot(np.array(e["dose_ci95"]) / sd, [yy + 0.12] * 2, color=S.C["blue"], lw=0.9)
        cx.scatter([e["dose_slope_per_10ms"] / sd], [yy + 0.12], s=10, color=S.C["blue"], zorder=3)
        cx.plot(np.array(e["jump_ci95"]) / sd, [yy - 0.12] * 2, color=S.C["vermillion"], lw=0.9)
        cx.scatter([e["jump"] / sd], [yy - 0.12], s=10, color=S.C["vermillion"], zorder=3)
    c6 = d["cmip6"]["E1_E2"]
    cx.plot(c6["dose_ci95"], [ypos[-1] + 0.12] * 2, color=S.C["blue"], lw=0.9)
    cx.scatter([c6["dose_slope_per_10ms"]], [ypos[-1] + 0.12], s=10, color=S.C["blue"], zorder=3,
               label="per 10 m/s stronger minimum wind")
    cx.plot(c6["jump_ci95"], [ypos[-1] - 0.12] * 2, color=S.C["vermillion"], lw=0.9)
    cx.scatter([c6["jump"]], [ypos[-1] - 0.12], s=10, color=S.C["vermillion"], zorder=3,
               label="step at reversal")
    cx.axvline(0, color="k", lw=0.5)
    cx.set_yticks(ypos)
    cx.set_yticklabels([r[0] for r in rows] + [f"CMIP6 ({d['cmip6']['n_units']:,})"])
    cx.set_xlabel("estimate (s.d. of outcome)")
    cx.set_ylim(-0.6, len(rows) + 0.5)
    cx.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), fontsize=4.8, ncol=2)
    S.panel_label(cx, "c", x=-0.45, y=1.04)
    S.save(fig, "figED6_continuity")


if __name__ == "__main__":
    main()
