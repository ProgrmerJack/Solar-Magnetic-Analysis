#!/usr/bin/env python3
"""
fig5_law.py -- Figure 5: one law for outcome-defined impact classes.

Reads (recomputes only the Gaussian curve and per-unit values stored in the files):
  results/current/8_experiment/snapsi_selection_test.json        SNAPSI rates and contrasts (L)
  results/current/8_experiment/snapsi_contrast_generality.json   G1x, G3, G5, G6
  results/current/8_experiment/davis2021_intervention_test.json  January 2021
  results/current/8_experiment/snapsi_free_natural_experiment.json  F and L (if present)
  results/current/8_experiment/damip_class_law.json              second domain

  a  contrast change against the change a contrast proportional to the spread
     predicts, C_ref (k - 1), unit by unit: SNAPSI zonal pairs and full nudging
     (post hoc), the free-arm natural experiment and the human-forcing years
     (registered); lines, one population (slope 1) and two types 2 s.d. apart (0.6)
  b  law residual (forced contrast minus the location-scale null of the reference)
     in every intervention set, with the shift-only residual for comparison;
     grey band, the tolerance fixed before the tests
  c  elasticity of the contrast with respect to the spread: one population
     predicts 1, forced mixing of two response types 0.60, 0.37, 0.16 for
     separations of 2, 3, 4 s.d.
  d  second domain: human forcing and extreme-season classes (HadGEM3-GC31-LL,
     historical vs natural forcing, 1985-2014): class rate against law residual
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
    if dam:
        for y in dam["targets"]["CEUR_JJA_tercile"]["years"]:
            pts.append(("human forcing, C. Europe (registered)", y["law_pred_diff"], y["diff"]))
    u = pd.DataFrame(pts, columns=["set", "pred", "obs"]).dropna()
    lim = (min(u.pred.min(), u.obs.min()) - 0.05, max(u.pred.max(), u.obs.max()) + 0.05)
    xx = np.array(lim)
    ax.plot(xx, xx, color="k", lw=0.7, label="one population (slope 1)")
    ax.plot(xx, 0.6 * xx, color="0.5", lw=0.7, ls=":", label="two types, Δ = 2 (slope 0.6)")
    cols = {"SNAPSI zonal (post hoc)": verm, "SNAPSI full 3-D (post hoc)": S.C["orange"],
            "SNAPSI free, natural (registered)": S.C["green"], "human forcing, C. Europe (registered)": blue}
    for k, g in u.groupby("set", sort=False):
        ax.scatter(g.pred, g.obs, s=8, color=cols[k], zorder=3, label=k)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("contrast change predicted by the spread change (σ)")
    ax.set_ylabel("observed contrast change (σ)")
    ax.legend(fontsize=4.4, loc="upper left")
    ax.set_title("contrast changes are changes of spread", fontsize=6)
    S.panel_label(ax, "a", x=-0.2, y=1.04)

    # b -- law residual by set
    bx = fig.add_subplot(gs[0, 1])
    rows = []
    g3 = gen["G3_third_event_SH"]
    rows.append(("SNAPSI zonal, 25 pairs", g3["NH_25"]["resid_locscale"], g3["NH_25"]["resid_shift"]))
    rows.append(("SNAPSI SH warming", g3["SH"]["resid_locscale"], g3["SH"]["resid_shift"]))
    g6 = gen.get("G6_regional_classes")
    if g6:
        rows.append(("24 regional cold classes", g6["resid_law"], g6["resid_shift"]))
    if fre and "L_pooled_law" in fre:
        ps = fre["L_pooled_law"]["location_scale"]
        if fre.get("F2_resid_law") and fre["F2_resid_law"].get("ci95"):
            rows.append(("SNAPSI free, natural", fre["F2_resid_law"], None))
        rows.append(("pooled, all interventions", ps, fre["L_pooled_law"]["shift_only"]))
    if dam:
        for t, lab in (("CEUR_JJA_tercile", "human forcing, C. Europe JJA"), ("WRUS_JJA_tercile", "W. Russia JJA"),
                       ("CNAM_JJA_tercile", "C. N. America JJA"), ("CEAS_DJF_tercile", "C. Eurasia DJF")):
            v = dam["targets"][t]
            rows.append((lab, v["resid_law"], v["resid_shift"]))
    bx.axvspan(-0.25, 0.25, color="0.92", zorder=0)
    for i, (lab, law, sh) in enumerate(rows):
        if sh and sh.get("ci95"):
            bx.plot(sh["ci95"], [i + 0.18] * 2, color=grey, lw=0.8)
            bx.scatter([sh["mean"]], [i + 0.18], s=10, facecolor="white", edgecolor=grey, zorder=3)
        bx.plot(law["ci95"], [i] * 2, color=verm, lw=1.1)
        bx.scatter([law["mean"]], [i], s=12, color=verm, zorder=3)
    bx.axvline(0, color="k", lw=0.5)
    bx.set_yticks(range(len(rows))); bx.set_yticklabels([r_[0] for r_ in rows], fontsize=5.0)
    bx.set_ylim(-0.6, len(rows) - 0.3); bx.invert_yaxis(); bx.set_xlim(-0.45, 0.45)
    bx.set_xlabel("forced contrast minus null contrast (σ)")
    bx.set_title("law (filled) vs shift only (open)", fontsize=6)
    S.panel_label(bx, "b", x=-0.62, y=1.04)

    # c -- elasticity
    cx = fig.add_subplot(gs[1, 0])
    ent = []
    e36 = gen["G1x_post_hoc_spread"]["elasticity_slope_36"]
    ent.append(("SNAPSI zonal (36)", e36["slope"], e36["ci95"], "post hoc"))
    if fre and fre.get("F3_Delta_b_vs_law"):
        f3 = fre["F3_Delta_b_vs_law"]
        ent.append(("SNAPSI free (registered)", f3["slope"], f3["slope_ci95"], ""))
    if dam:
        for t, lab in (("CEUR_JJA_tercile", "C. Europe JJA"), ("WRUS_JJA_tercile", "W. Russia JJA"),
                       ("CNAM_JJA_tercile", "C. N. America JJA"), ("CEAS_DJF_tercile", "C. Eurasia DJF")):
            e = dam["targets"][t].get("elasticity_slope")
            if e:
                ent.append((f"human forcing, {lab}", e["slope"], e["ci95"], ""))
    for v, lab in ((1.0, "one population"), (0.60, "two types, Δ = 2"), (0.37, "Δ = 3"), (0.16, "Δ = 4")):
        cx.axvline(v, color="k" if v == 1 else "0.6", lw=0.6, ls="-" if v == 1 else ":")
        cx.text(v + 0.02, len(ent) - 0.5, lab, rotation=90, fontsize=4.6, ha="left", va="bottom", color="0.35")
    for i, (lab, sl, ci, note) in enumerate(ent):
        if ci:
            cx.plot(ci, [i] * 2, color=blue, lw=1.0)
        cx.scatter([sl], [i], s=12, color=grey if e_note(note) else blue, zorder=3)
    cx.set_yticks(range(len(ent))); cx.set_yticklabels([e[0] + (f" ({e[3]})" if e[3] else "") for e in ent], fontsize=5.0)
    cx.set_ylim(len(ent) - 0.5, -0.6); cx.set_xlim(0, 1.8)
    cx.set_xlabel("elasticity of the contrast with respect to the spread")
    cx.set_title("the contrast follows the spread, as one population requires", fontsize=6)
    S.panel_label(cx, "c", x=-0.62, y=1.04)

    # d -- second domain
    dx = fig.add_subplot(gs[1, 1])
    if dam:
        for t, col in (("CEUR_JJA_tercile", verm), ("WRUS_JJA_tercile", S.C["orange"]),
                       ("CNAM_JJA_tercile", S.C["green"]), ("CEAS_DJF_tercile", blue)):
            y = pd.DataFrame(dam["targets"][t]["years"])
            dx.scatter(y.rate_hist - y.rate_nat, y.resid_law, s=7, color=col, alpha=0.8, label=t.split("_tercile")[0].replace("_", " "))
        dx.axhspan(-0.25, 0.25, color="0.92", zorder=0)
        dx.axhline(0, color="k", lw=0.5); dx.axvline(0, color="k", lw=0.5)
        dx.legend(fontsize=4.8, loc="lower right")
        v = dam["targets"]["CEUR_JJA_tercile"]
        dx.set_title(M(f"human forcing: hot-summer rate {v['rate_nat']['mean']:.2f} → {v['rate_hist']['mean']:.2f}, "
                       f"law residual {v['resid_law']['mean']:+.2f}"), fontsize=5.8)
    dx.set_xlabel("change in class rate (historical − natural forcing), per year")
    dx.set_ylabel("law residual (σ)")
    S.panel_label(dx, "d", x=-0.2, y=1.04)
    S.save(fig, "fig5_law")


if __name__ == "__main__":
    main()
