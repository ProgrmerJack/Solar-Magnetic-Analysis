#!/usr/bin/env python3
"""
fig5_operational.py -- Extended Data Fig. (operational reforecasts, formerly Fig. 5): operational reforecasts after SSWs.

Reads results/current/6_predictability/s2s_forecast_test.json (result P, ECMWF)
and s2s_multimodel_test.json (result A, all systems). Recomputes nothing.

  a  ECMWF, starts 2-9 d before onset: ensemble-mean against observed surface
     response (A = minus polar-cap msl anomaly, days +8..+25) for each SSW
  b  discrimination r per system and multi-model mean, against the 95% range of
     the conditional calendar-window null (event-free dates, matched variance)
  c  mean PIT of the observed outcome per system, against the 95% range of the
     same null; 0.5 is a calibrated forecast
  d  regional 2 m temperature (s2s_regional_test.json): mean rank of the observed
     temperature in the ensembles, per region, nine-system mean (diamond) and each
     system (dots), against the central 95% of event-free dates
"""
import json
import sys
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nature_style as S  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

REG = [("NEURASIA", "N Eurasia"), ("HI_EUROPE", "high-lat.\nEurope"),
       ("MID_EASIA", "mid-lat.\nE Asia"), ("MID_NAMER", "mid-lat.\nN America")]
LABEL = {"ecmwf": "ECMWF", "eccc": "ECCC", "cma": "CMA", "hmcr": "HMCR", "kma": "KMA",
         "cnrm": "CNRM", "jma": "JMA", "cnr_isac": "CNR-ISAC", "ncep": "NCEP",
         "cptec": "CPTEC"}


def main():
    S.apply()
    p = json.loads((S.RESULTS / "6_predictability" / "s2s_forecast_test.json").read_text())
    m = json.loads((S.RESULTS / "6_predictability" / "s2s_multimodel_test.json").read_text())
    op = json.loads((S.RESULTS / "6_predictability" / "s2s_regional_test.json").read_text())["regions"]
    fig = plt.figure(figsize=(S.DOUBLE, 128 * S.MM))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1.1, 1.1], wspace=0.55, hspace=0.75)

    ax = fig.add_subplot(gs[0, 0])
    ev = p["bins"]["short"]["events"]
    a = np.array([e["A_ens"] for e in ev]); o = np.array([e["A_obs"] for e in ev])
    ax.scatter(a, o, s=12, color=S.C["vermillion"], zorder=3)
    lim = 1.08 * max(np.abs(a).max(), np.abs(o).max())
    ax.plot([-lim, lim], [-lim, lim], color="0.6", lw=0.5, ls=":")
    ax.axhline(0, color="k", lw=0.4); ax.axvline(0, color="k", lw=0.4)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("forecast ensemble mean (Pa)"); ax.set_ylabel("observed (Pa)")
    dsc = p["bins"]["short"]["discrimination"]
    ax.set_title(f"ECMWF, {len(ev)} SSWs", fontsize=6)
    ax.text(0.03, 0.97, f"r = {dsc['r_ensmean_vs_obs']:.2f}\nevent-free {dsc['r_null_conditional_mean']:.2f}"
            f"\np = {dsc['p_r_conditional']:.2f}", transform=ax.transAxes, va="top", fontsize=5.5)
    S.panel_label(ax, "a", x=-0.3, y=1.06)

    rows = [(LABEL[c], v) for c, v in m["centres"].items() if "D_r" in v]
    rows.append(("multi-model", m["multimodel"]["confirmatory"]))   # excludes ECMWF
    y = np.arange(len(rows))[::-1]
    for panel, (key, qkey, ref, xl, tt) in zip("bc", [
            ("D_r", "D_r_null_cond_q025_q975", None, "discrimination r", "event-to-event differences"),
            ("H1_mean_pit", "H1_null_pit_q025_q975", 0.5, "mean PIT of the outcome",
             "rank of the outcome in the ensemble")]):
        bx = fig.add_subplot(gs[0, 1 if panel == "b" else 2])
        for yy, (lab, v) in zip(y, rows):
            mm = lab == "multi-model"
            q = v.get(qkey)
            if q:
                bx.plot(q, [yy, yy], color="0.75", lw=3, solid_capstyle="butt")
            few = (not mm) and v["n_events"] < 6            # e.g. CPTEC, 4 events: shown open
            bx.scatter([v[key]], [yy], s=14 if mm else 9, marker="D" if mm else "o",
                       facecolors="none" if few else (S.C["black"] if mm else S.C["vermillion"]),
                       edgecolors=S.C["black"] if mm else S.C["vermillion"], zorder=3)
            bx.text(1.02, yy, f"{v['n_events']}", transform=bx.get_yaxis_transform(),
                    va="center", fontsize=5, color="0.4")
        if ref is not None:
            bx.axvline(ref, color="k", lw=0.4, ls=":")
        bx.set_yticks(y); bx.set_yticklabels([("multi-model\n(excl. ECMWF)" if lab == "multi-model" else lab)
                            for lab, _ in rows] if panel == "b" else [])
        bx.set_xlabel(xl)
        bx.set_title(tt, fontsize=6, pad=8)
        bx.text(1.02, 1.0, "n", transform=bx.transAxes, va="bottom", fontsize=5, color="0.4")
        S.panel_label(bx, panel, x=-0.45 if panel == "b" else -0.16, y=1.06)
    dx = fig.add_subplot(gs[1, 0:3])
    for i, (k, _) in enumerate(REG):
        mm = op[k]["multimodel"]["confirmatory"]
        lo, hi = mm["H1_null_pit_q025_q975"]
        dx.plot([lo, hi], [i, i], color="0.78", lw=4, solid_capstyle="butt",
                label="event-free dates (95%)" if i == 0 else None)
        pits = [v["H1_mean_pit"] for v in op[k]["centres"].values() if "H1_mean_pit" in v]
        dx.scatter(pits, [i + 0.22] * len(pits), s=5, color=S.C["vermillion"], alpha=0.7,
                   label="each system" if i == 0 else None)
        dx.scatter([mm["H1_mean_pit"]], [i], s=22, marker="D", color="k", zorder=3,
                   label="nine-system mean" if i == 0 else None)
        # round half up, as in the text (0.0135 -> 0.014; float formatting gave 0.013)
        pv = Decimal(str(mm["H1_p"])).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP)
        dx.text(0.72, i, f"p = {pv}", va="center", fontsize=5.5)
    dx.axvline(0.5, color="k", lw=0.4, ls=":")
    dx.set_yticks(range(len(REG))); dx.set_yticklabels([l.replace("\n", " ") for _, l in REG])
    dx.set_ylim(len(REG) - 0.4, -0.6)
    allp = [v["H1_mean_pit"] for k, _ in REG for v in op[k]["centres"].values() if "H1_mean_pit" in v]
    dx.set_xlim(min(0.3, min(allp) - 0.03), 0.8)
    dx.set_xlabel("mean rank of the observed temperature in the forecast ensembles (0.5 = calibrated)")
    dx.legend(fontsize=5, loc="upper center", ncol=3, bbox_to_anchor=(0.5, -0.32))
    dx.set_title("Ten operational systems, 17 SSWs: colder than forecast?", fontsize=6)
    S.panel_label(dx, "d", x=-0.13, y=1.04)
    S.save(fig, "fig5_operational")


if __name__ == "__main__":
    main()
