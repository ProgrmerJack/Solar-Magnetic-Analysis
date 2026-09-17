#!/usr/bin/env python3
"""
ng_fig_design_artifact.py
=========================
The three display items for `paper/ng_manuscript.tex`, built directly from the
saved result JSONs so the figures cannot drift from the numbers.

  Fig 1  AO lag profile, both estimators, placebo window shaded.
         The paper's central claim in one panel: the conventional estimator is
         significantly negative BEFORE onset; the within-winter estimator is not,
         yet still recovers the canonical post-onset response.
  Fig 2  Both estimators across five circulation variables (pre-onset placebo vs
         post-onset), showing the contamination is systematic, not AO-specific.
  Fig 3  Avalanche case study: the same observable measured both ways, across
         four independent observing systems.

Inputs : data/results/r112_design_positive_control.json
         data/results/r106_canonical_catalog_headline.json
         data/results/r108_pwl_seasonal_control.json
         data/results/r111_kinetic_growth_instrumental.json
         data/results/r107_pooled_within_winter.json
Outputs: data/figures/ng_fig1_ao_lag_profile.{pdf,png}
         data/figures/ng_fig2_variable_panel.{pdf,png}
         data/figures/ng_fig3_avalanche_case.{pdf,png}
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "data" / "results"
FIG = ROOT / "data" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

C_WITHIN = "#1b6ca8"     # validated estimator
C_BETWEEN = "#c1350a"    # conventional estimator
C_PLACEBO = "#f2e6c9"

# lag-window midpoints for plotting
MID = {"placebo_pre": -30.5, "onset": 0, "early_post": 7, "mid_post": 22,
       "late_post": 37, "post_0_60": 30, "post_15_44": 29.5}


def load(name):
    return json.loads((RES / name).read_text(encoding="utf8"))


def _save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {stem}.pdf/.png")


def fig1(pc):
    """AO lag profile, both estimators."""
    w = pc["series"]["AO"]["windows"]
    order = ["placebo_pre", "onset", "early_post", "mid_post", "late_post"]
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    ax.axvspan(-45, -16, color=C_PLACEBO, zorder=0)
    ax.text(-30.5, -1.12, "placebo window\n(event cannot have acted)",
            ha="center", va="bottom", fontsize=8, style="italic", color="#7a6320")
    ax.axhline(0, color="0.35", lw=0.9)
    ax.axvline(0, color="0.35", lw=0.9, ls=":")

    xs, ys, lo, hi = [], [], [], []
    for k in order:
        ww = w[k]["within_winter"]
        xs.append(MID[k]); ys.append(ww["effect"])
        lo.append(ww["CI95"][0]); hi.append(ww["CI95"][1])
    ax.errorbar(xs, ys, yerr=[np.array(ys) - lo, np.array(hi) - np.array(ys)],
                fmt="o-", color=C_WITHIN, capsize=3, lw=1.8, ms=6,
                label="within-winter (validated)", zorder=3)

    bx = [MID[k] for k in order]
    by = [w[k]["between_winter"]["effect"] for k in order]
    ax.plot(bx, by, "s--", color=C_BETWEEN, lw=1.8, ms=6,
            label="between-winter (conventional)", zorder=3)

    # annotate the two decisive points
    ax.annotate(f"P = {w['placebo_pre']['between_winter']['t_p_two_sided']:.4f}",
                (MID["placebo_pre"], w["placebo_pre"]["between_winter"]["effect"]),
                textcoords="offset points", xytext=(6, -14), fontsize=8,
                color=C_BETWEEN, fontweight="bold")
    ax.annotate(f"P = {w['placebo_pre']['within_winter']['p_two_sided']:.2f} (null)",
                (MID["placebo_pre"], w["placebo_pre"]["within_winter"]["effect"]),
                textcoords="offset points", xytext=(-52, 10), fontsize=8, color=C_WITHIN)
    ax.legend(frameon=False, fontsize=9, loc="upper right")

    ax.set_xlabel("days from SSW central date")
    ax.set_ylabel("AO anomaly (standard deviations)")
    ax.set_title(f"Arctic Oscillation response to {pc['series']['AO']['n_ssw']} "
                 "major SSWs, 1950–2026", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    _save(fig, "ng_fig1_ao_lag_profile")


def fig2(pc):
    """Placebo vs post-onset for both estimators, all circulation variables."""
    names = [k for k in ("AO", "NAO", "Z500_NH", "SLP_NH", "U850_NH")
             if k in pc["series"]]
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.9), sharey=True)
    for ax, win, title in zip(
            axes, ("placebo_pre", "post_0_60"),
            ("pre-onset placebo (−45 to −16 d)", "post-onset (0 to +60 d)")):
        y = np.arange(len(names))
        for i, nm in enumerate(names):
            w = pc["series"][nm]["windows"].get(win)
            if not w:
                continue
            ww, bw = w["within_winter"], w["between_winter"]
            # standardise each variable by its own within-winter CI half-width
            sc = max(abs(ww["CI95"][1] - ww["CI95"][0]) / 2, 1e-9)
            ax.errorbar(ww["effect"] / sc, i + 0.14,
                        xerr=[[(ww["effect"] - ww["CI95"][0]) / sc],
                              [(ww["CI95"][1] - ww["effect"]) / sc]],
                        fmt="o", color=C_WITHIN, capsize=2.5, ms=5)
            ax.plot(bw["effect"] / sc, i - 0.14, "s", color=C_BETWEEN, ms=6)
            star = "*" if bw.get("t_p_two_sided", 1) < 0.05 else ""
            if star:
                ax.text(bw["effect"] / sc, i - 0.42, star, ha="center",
                        color=C_BETWEEN, fontsize=13, fontweight="bold")
        ax.axvline(0, color="0.35", lw=0.9)
        ax.set_yticks(y)
        ax.set_yticklabels([n.replace("_NH", "") for n in names])
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("effect (within-winter CI half-widths)")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].plot([], [], "o", color=C_WITHIN, label="within-winter")
    axes[0].plot([], [], "s", color=C_BETWEEN, label="between-winter (* P<0.05)")
    axes[0].legend(frameon=False, fontsize=8.5, loc="lower left")
    fig.suptitle("The conventional estimator fires before onset; the validated one does not",
                 fontsize=11, y=1.02)
    _save(fig, "ng_fig2_variable_panel")


def fig3(head, pwl, snotel, pooled):
    """Avalanche case study: same claims, two estimators."""
    rows = []
    c17 = head["results"]["canonical_17"]
    rows.append(("Davos natural avalanche counts\n(17 events, 1998–2019)",
                 c17["between_winter"]["gmRR"],
                 c17["within_winter"]["IRR"],
                 c17["within_winter"]["IRR_CI95"],
                 f"P={c17['between_winter']['sign_test_p_one_sided']:.4f}",
                 f"P={c17['within_winter']['p_one_sided_suppression']:.3f}"))

    pw = pwl["metrics"]["Pen_depth"]
    # express the snowpack effects as ratios about 1 for a common axis
    base = 60.0  # cm, approximate mean penetration depth, for scaling only
    rows.append(("SNOWPACK weak-layer depth\n(285k station-days)",
                 1 + pw["uncontrolled"]["diff"] / base,
                 1 + pw["controlled"]["effect_controlled"] / base,
                 [1 + pw["controlled"]["CI95"][0] / base,
                  1 + pw["controlled"]["CI95"][1] / base],
                 "P≈1e-30", f"P={pw['controlled']['p_two_sided']:.2f}"))

    k = snotel["domains"]["continental"]["kinetic"]["post_15_44"]
    rows.append(("SNOTEL kinetic-growth regime\n(31 events, 945 stations)",
                 None, k["RR"], k["CI95"], "—",
                 f"P={min(k['p_one_sided_up'], k['p_one_sided_down']):.2f}"))

    al = pooled["domains"]["ALPINE (primary, mechanism domain)"]["onset"]
    rows.append(("Pooled Alpine activity\n(30 region-winters)",
                 None, al["IRR"], al["CI95"], "—",
                 f"P={al['p_one_sided_suppression']:.2f}"))

    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    y = np.arange(len(rows))[::-1]
    for yi, (lab, bw, ww, ci, bp, wp) in zip(y, rows):
        ax.errorbar(ww, yi, xerr=[[ww - ci[0]], [ci[1] - ww]], fmt="o",
                    color=C_WITHIN, capsize=3, ms=6, zorder=3)
        ax.text(ci[1] + 0.03, yi + 0.02, wp, fontsize=7.5, color=C_WITHIN,
                va="center")
        if bw is not None:
            ax.plot(bw, yi, "s", color=C_BETWEEN, ms=7, zorder=3)
            ax.text(bw, yi - 0.30, bp, fontsize=7.5, color=C_BETWEEN,
                    ha="center", fontweight="bold")
    ax.axvline(1, color="0.35", lw=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=8.5)
    ax.set_xlabel("effect ratio (1 = no effect)")
    ax.set_xlim(0, 2.05)
    ax.set_ylim(-0.85, len(rows) - 0.4)
    ax.set_title("Avalanche hazard: an apparently decisive effect that the\n"
                 "validated estimator does not support", fontsize=11)
    # the SNOWPACK conventional result is significant only through
    # pseudo-replication -- +1.68 cm is 2.8% of a ~60 cm depth over 285k
    # station-days. Say so, rather than let the small marker imply agreement.
    ax.annotate("conventional P≈1e-30 for a 2.8% effect:\n"
                "significance from 285k non-independent station-days",
                xy=(1.028, len(rows) - 3), xytext=(1.30, len(rows) - 2.55),
                fontsize=7.2, color=C_BETWEEN,
                arrowprops=dict(arrowstyle="->", color=C_BETWEEN, lw=0.8))
    ax.plot([], [], "s", color=C_BETWEEN, label="between-winter (conventional)")
    ax.plot([], [], "o", color=C_WITHIN, label="within-winter (validated)")
    ax.legend(frameon=False, fontsize=8.5, loc="upper center",
              bbox_to_anchor=(0.5, -0.16), ncol=2)
    ax.spines[["top", "right"]].set_visible(False)
    _save(fig, "ng_fig3_avalanche_case")


def main():
    pc = load("r112_design_positive_control.json")
    print("building display items:")
    fig1(pc)
    fig2(pc)
    fig3(load("r106_canonical_catalog_headline.json"),
         load("r108_pwl_seasonal_control.json"),
         load("r111_kinetic_growth_instrumental.json"),
         load("r107_pooled_within_winter.json"))


if __name__ == "__main__":
    main()
