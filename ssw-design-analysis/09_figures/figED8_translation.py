#!/usr/bin/env python3
"""
figED8_translation.py -- Extended Data Fig. 8: is the imposed SSW a translation?

Reads results/current/8_experiment/snapsi_residual_quantiles.json. Recomputes nothing.

  a  polar-cap NAM, days 8-25: residual quantile R(q) = Q_nudged - Q_control - shift,
     pooled over 36 ensembles (black, 95% interval) and per event (colours), with the
     tolerance declared before the test (shaded)
  b  the same for northern-Eurasian temperature, days 8-24
  c  tail-probability error of the translated control, E(q), both outcomes
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nature_style as S  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402


def panel(ax, block, tol, key, ylabel):
    ev = {"pooled": ("k", 0.0, "36 ensembles"), "Feb-2018": (S.C["vermillion"], -0.012, "February 2018"),
          "Jan-2019": (S.C["blue"], 0.012, "January 2019")}
    for name, (col, dx, lab) in ev.items():
        b = block[name][key]
        q = np.array([float(k) for k in b]); v = np.array([b[k] for k in b])
        ax.errorbar(q + dx, v[:, 0], yerr=[v[:, 0] - v[:, 1], v[:, 2] - v[:, 0]], fmt="o-" if name == "pooled" else "o",
                    ms=2.5 if name == "pooled" else 2, lw=0.9 if name == "pooled" else 0.6, color=col, capsize=0,
                    elinewidth=0.6, label=lab, zorder=3 if name == "pooled" else 2)
    ax.axhspan(-tol, tol, color="0.9", zorder=0)
    ax.axhline(0, color="k", lw=0.4)
    ax.set_xlabel("quantile")
    ax.set_ylabel(ylabel)


def main():
    S.apply()
    d = json.loads((S.RESULTS / "8_experiment" / "snapsi_residual_quantiles.json").read_text())
    tR, tE = d["tolerance"]["R_sigma"], d["tolerance"]["E_probability"]
    A, T = d["A_polar_cap_NAM_days8_25"], d["T_northern_Eurasia_days8_24"]
    fig, (ax, bx, cx) = plt.subplots(1, 3, figsize=(S.DOUBLE, 60 * S.MM), gridspec_kw={"wspace": 0.45})
    panel(ax, A, tR, "R", "residual quantile (s.d.)")
    ax.set_title(f"polar-cap NAM: {A['pooled']['reading']}", fontsize=6)
    ax.set_ylim(-0.6, 0.6); ax.legend(loc="lower left", fontsize=4.8)
    S.panel_label(ax, "a", x=-0.3, y=1.04)
    panel(bx, T, tR, "R", "residual quantile (s.d.)")
    bx.set_title(f"N-Eurasian temperature: {T['pooled']['reading']}", fontsize=6)
    bx.set_ylim(-0.6, 0.6)
    S.panel_label(bx, "b", x=-0.3, y=1.04)
    for blk, col, lab, dx in ((A, "k", "polar-cap NAM", -0.004), (T, S.C["green"], "N-Eurasian temperature", 0.004)):
        b = blk["pooled"]["E"]
        q = np.array([float(k) for k in b]); v = np.array([b[k] for k in b])
        cx.errorbar(q + dx, v[:, 0], yerr=[v[:, 0] - v[:, 1], v[:, 2] - v[:, 0]], fmt="o", ms=2.5, color=col,
                    elinewidth=0.7, capsize=0, label=lab)
    cx.axhspan(-tE, tE, color="0.9", zorder=0); cx.axhline(0, color="k", lw=0.4)
    cx.set_xlabel("tail probability q"); cx.set_ylabel("error of translated control")
    cx.set_ylim(-0.08, 0.08); cx.legend(loc="upper left", fontsize=4.8)
    S.panel_label(cx, "c", x=-0.32, y=1.04)
    S.save(fig, "figED8_translation")


if __name__ == "__main__":
    main()
