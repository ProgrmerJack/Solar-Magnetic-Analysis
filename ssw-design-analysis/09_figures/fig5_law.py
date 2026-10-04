#!/usr/bin/env python3
"""
fig5_law.py -- Figure 5: one law for outcome-defined impact classes.

Reads (recomputes only the Gaussian curve and per-unit values stored in the files):
  results/current/8_experiment/snapsi_selection_test.json        SNAPSI rates and contrasts (L)
  results/current/8_experiment/snapsi_contrast_generality.json   G1x, G3, G5, G6
  results/current/8_experiment/davis2021_intervention_test.json  January 2021
  results/current/8_experiment/snapsi_free_natural_experiment.json  F and L (if present)
  results/current/8_experiment/damip_class_law.json              second domain

  a  back-shift contrast change against the change a contrast proportional to the
     spread predicts, C_ref (k - 1), unit by unit: SNAPSI zonal pairs and full
     nudging (post hoc), the free-arm natural experiment (registered); lines, one
     population (slope 1) and two types 2 s.d. apart (0.58)
  b  law residual (forced contrast minus the location-scale null of the reference)
     in every intervention set, with the shift-only residual for comparison;
     grey band, the tolerance fixed before the tests
  c  elasticity of the contrast with respect to the spread: one population
     predicts 1, forced mixing of two response types 0.58, 0.35, 0.14 for
     separations of 2, 3, 4 s.d. (threshold_law_theory.json)
  d  second domain (CanESM5 p1/p2, MIROC6, HadGEM3; historical vs natural forcing,
     1985-2014): observed law residual per model and region against its expectation
     when the law holds exactly with the model's own distribution (post hoc; grey,
     95% range); dashed, equality
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nature_style as S  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

M = lambda s: s.replace("-", "−")                     # noqa: E731


def J(name):
    p = S.RESULTS / "8_experiment" / name
    return json.loads(p.read_text()) if p.exists() else None


def g_abs(rate):
    t = norm.ppf(np.clip(rate, 1e-3, 1 - 1e-3))
    return norm.pdf(t) / (norm.cdf(t) * (1 - norm.cdf(t)))


def e_note(note):
    return bool(note)


def main():
    S.apply()
    sel, gen, dav, fre, dam = (J("snapsi_selection_test.json"), J("snapsi_contrast_generality.json"),
                               J("davis2021_intervention_test.json"), J("snapsi_free_natural_experiment.json"),
                               J("damip_class_law.json"))
    dmm = J("damip_class_law_multimodel.json")
    grey, verm, blue = "0.55", S.C["vermillion"], S.C["blue"]
    fig = plt.figure(figsize=(S.DOUBLE, 118 * S.MM))
    gs = fig.add_gridspec(2, 2, wspace=0.36, hspace=0.62)

    # a -- the contrast change follows the spread change
    ax = fig.add_subplot(gs[0, 0])
    pts = []
    for e in gen["pairs_NH"]:
        pts.append(("SNAPSI zonal (post hoc)", e["C_control"] * (e["sd_ratio"] - 1), e["Delta_b"]))
    if isinstance(gen.get("G5_nudged_full"), dict):
        for e in gen["G5_nudged_full"]["control"]["pairs"]:
            pts.append(("SNAPSI full 3-D (post hoc)", e["C_control"] * (e["sd_ratio"] - 1), e["Delta_b"]))
    if fre:
        for e in fre.get("ensembles", []):
            if e.get("mixed"):
                pts.append(("SNAPSI free, natural (registered)", e["law_pred_Delta_b"], e["Delta_b"]))
    u = pd.DataFrame(pts, columns=["set", "pred", "obs"]).dropna()
    lim = (min(u.pred.min(), u.obs.min()) - 0.05, max(u.pred.max(), u.obs.max()) + 0.05)
    xx = np.array(lim)
    ax.plot(xx, xx, color="k", lw=0.7, label="one population (slope 1)")
    e2 = J("threshold_law_theory.json")["mixture_elasticity"]["2.0"]
    ax.plot(xx, e2 * xx, color="0.5", lw=0.7, ls=":", label=f"two types, Δ = 2 (slope {e2:.2f})")
    cols = {"SNAPSI zonal (post hoc)": verm, "SNAPSI full 3-D (post hoc)": S.C["orange"],
            "SNAPSI free, natural (registered)": S.C["green"]}
    for k, g in u.groupby("set", sort=False):
        ax.scatter(g.pred, g.obs, s=8, color=cols[k], zorder=3, label=k)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("contrast change predicted by the spread change (σ)")
    ax.set_ylabel("observed contrast change (σ)")
    ax.legend(fontsize=4.2, loc="lower right", framealpha=0.9)
    ax.set_title("contrast changes are changes of spread", fontsize=6)
    S.panel_label(ax, "a", x=-0.2, y=1.04)

    # b -- law residual by set
    bx = fig.add_subplot(gs[0, 1])
    rows = []
    if fre and "L_pooled_law" in fre:
        Lp = fre["L_pooled_law"]
        names = {"SNAPSI_zonal_NH_25": "SNAPSI zonal-mean", "SNAPSI_SH": "SNAPSI SH warming",
                 "SNAPSI_full_3D": "SNAPSI full 3-D", "CESM2_Jan2021": "CESM2, January 2021",
                 "SNAPSI_free_natural": "SNAPSI free, natural"}
        for k, lab in names.items():
            ls_m, sh_m = Lp["location_scale"]["per_set"][k][0], Lp["shift_only"]["per_set"][k][0]
            rows.append((lab, {"mean": ls_m, "ci95": Lp["location_scale"]["per_set_ci95"].get(k)},
                         {"mean": sh_m, "ci95": Lp["shift_only"]["per_set_ci95"].get(k)}))
        rows.append(("pooled, five sets", Lp["location_scale"], Lp["shift_only"]))
    g6 = gen.get("G6_regional_classes")
    if g6:
        rows.append(("24 regional cold classes", g6["resid_law"], g6["resid_shift"]))
    if dmm:
        rows.append(("human forcing, 3 models", dmm["M1_CEUR_JJA_three_models"], None))
    bx.axvspan(-0.25, 0.25, color="0.92", zorder=0)
    for i, (lab, law, sh) in enumerate(rows):
        if sh is not None:
            if sh.get("ci95"):
                bx.plot(sh["ci95"], [i + 0.18] * 2, color=grey, lw=0.8)
            bx.scatter([sh["mean"]], [i + 0.18], s=10, facecolor="white", edgecolor=grey, zorder=3)
        if law.get("ci95"):
            bx.plot(law["ci95"], [i] * 2, color=verm, lw=1.1)
        bx.scatter([law["mean"]], [i], s=12, color=verm, zorder=3)
    bx.axvline(0, color="k", lw=0.5)
    bx.set_yticks(range(len(rows))); bx.set_yticklabels([r_[0] for r_ in rows], fontsize=5.0)
    bx.set_ylim(-0.6, len(rows) - 0.3); bx.invert_yaxis(); bx.set_xlim(-0.5, 0.65)
    bx.set_xlabel("forced contrast minus null contrast (σ)")
    bx.set_title("law (filled) vs shift only (open)", fontsize=6)
    S.panel_label(bx, "b", x=-0.62, y=1.10)

    # c -- elasticity
    cx = fig.add_subplot(gs[1, 0])
    ent = []
    e36 = gen["G1x_post_hoc_spread"]["elasticity_slope_36"]
    ent.append(("SNAPSI zonal (36)", e36["slope"], e36["ci95"], "post hoc"))
    if fre and fre.get("F3_Delta_b_vs_law"):
        f3 = fre["F3_Delta_b_vs_law"]
        ent.append(("SNAPSI free (registered)", f3["slope"], f3["slope_ci95"], ""))
    if dmm:
        e = dmm["M2_elasticity_all"]
        ent.append(("human forcing, 3 models (registered)", e["slope"], e["ci95"], ""))
        cal = dmm["post_hoc_calibration"]["expected_elasticity_under_law"]
        ent.append(("  expected under the exact law", cal["mean"], cal["range95"], "post hoc"))
    th = J("threshold_law_theory.json")["mixture_elasticity"]
    for v, lab in ((1.0, "one population"), (th["2.0"], "two types, Δ = 2"), (th["3.0"], "Δ = 3"), (th["4.0"], "Δ = 4")):
        cx.axvline(v, color="k" if v == 1 else "0.6", lw=0.6, ls="-" if v == 1 else ":")
        cx.text(v + 0.02, len(ent) - 0.5, lab, rotation=90, fontsize=4.6, ha="left", va="bottom", color="0.35")
    for i, (lab, sl, ci, note) in enumerate(ent):
        if ci:
            cx.plot(ci, [i] * 2, color=blue, lw=1.0)
        cx.scatter([sl], [i], s=12, color=grey if e_note(note) else blue, zorder=3)
    cx.set_yticks(range(len(ent))); cx.set_yticklabels([e[0] + (f" ({e[3]})" if e[3] else "") for e in ent], fontsize=5.0)
    cx.set_ylim(len(ent) - 0.5, -0.6); cx.set_xlim(0, 1.8)
    cx.set_xlabel("elasticity of the contrast with respect to the spread")
    cx.set_title("elasticity: one population 1; forced mixing of types ≥ 2 s.d. apart below 0.6", fontsize=6)
    S.panel_label(cx, "c", x=-0.62, y=1.04)

    # d -- second domain: observed law residual against its expectation under the exact law
    dx = fig.add_subplot(gs[1, 1])
    if dmm:
        cal = dmm["post_hoc_calibration_empirical_shape"]["per_model_target"]
        cols = {"CanESM5_p1": verm, "CanESM5_p2": S.C["orange"], "MIROC6": S.C["green"], "HadGEM3-GC31-LL": blue}
        mk = {"JJA": "o", "DJF": "s"}
        for k, v in cal.items():
            m, t = k.split("|")
            lo, hi = v["range95_under_law"]
            dx.plot([v["expected_under_law"]] * 2, [lo, hi], color="0.75", lw=0.6, zorder=1)
            dx.scatter(v["expected_under_law"], v["observed"], s=14, color=cols[m], marker=mk[t.split("_")[1]], zorder=3)
        lim = (-0.3, 0.3)
        dx.plot(lim, lim, color="k", lw=0.6, ls=(0, (3, 2)))
        for m, c in cols.items():
            dx.scatter([], [], s=12, color=c, label=m.replace("-GC31-LL", ""))
        dx.scatter([], [], s=12, color="0.4", marker="o", label="summer (hot class)")
        dx.scatter([], [], s=12, color="0.4", marker="s", label="winter (cold class)")
        dx.legend(fontsize=4.4, loc="upper left")
        dx.set_xlim(lim); dx.set_ylim(lim)
        n_in = sum(v["inside_range"] for v in cal.values())
        dx.set_title(f"human forcing: {n_in} of {len(cal)} model-region cases inside the law's 95% range", fontsize=5.8)
    dx.set_xlabel("law residual expected under the exact law (σ)")
    dx.set_ylabel("observed law residual (σ)")
    S.panel_label(dx, "d", x=-0.2, y=1.04)
    S.save(fig, "fig5_law")


if __name__ == "__main__":
    main()
