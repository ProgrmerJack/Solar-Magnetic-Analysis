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
  d  forced mixing of response types (class_law_general.json, T3, registered): curvature
     of the ensemble variance in the ensemble mean across natural or early-historical,
     historical and SSP5-8.5 forcing, per model-region-season case, against how far the
     forcing moves the arm mean; fixed types give -1, one population with linear spread
     change >= 0; filled, mixing excluded by the post-hoc calibration against matched
     fixed-type mixtures; open, not excluded
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
    gen2 = J("class_law_general.json")
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
    if gen2:
        rows.append(("temperature, 3 models", gen2["T1_by_variable"]["tas"]["pooled"], None))
        rows.append(("precipitation, 4 models", gen2["T1_by_variable"]["pr"]["pooled"], None))
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

    # d -- forced mixing of fixed response types: variance-mean curvature
    dx = fig.add_subplot(gs[1, 1])
    if gen2:
        cols = {"tas": verm, "pr": blue}
        for k, u in gen2["units"].items():
            if "T3" not in u:
                continue
            t = u["T3"]; c = u.get("T3_post_hoc_calibration", {}); var = k.split("|")[1]
            ex = bool(c.get("mixing_excluded_calibrated")); x = c.get("forced_change_sd", t["mean_range_sd"])
            dx.plot([x] * 2, t["ci95"], color=cols[var], lw=0.5, alpha=0.5, zorder=1)
            dx.scatter(x, t["kappa"], s=9, zorder=3, color=cols[var] if ex else "white", edgecolor=cols[var], linewidth=0.6)
        dx.axhline(-1, color="k", lw=0.7, ls=":"); dx.axhline(0, color="k", lw=0.5)
        dx.text(0.22, -0.95, "fixed response types (−1)", fontsize=4.8, va="bottom")
        dx.text(0.22, 0.05, "one population", fontsize=4.8, va="bottom")
        dx.scatter([], [], s=9, color=verm, label="temperature"); dx.scatter([], [], s=9, color=blue, label="precipitation")
        dx.scatter([], [], s=9, color="white", edgecolor="0.4", label="mixing not excluded (calibrated)")
        dx.legend(fontsize=4.4, loc="upper right")
        dx.set_xscale("log"); dx.set_ylim(-1.3, 2.5)
        dx.set_title("forced mixing of fixed types would make variance fall with mean²", fontsize=5.8)
    dx.set_xlabel("forced change of the mean, lowest to highest forcing (s.d. of the reference)")
    dx.set_ylabel("curvature of ensemble variance in mean")
    S.panel_label(dx, "d", x=-0.2, y=1.04)
    S.save(fig, "fig5_law")


if __name__ == "__main__":
    main()
