"""
15_nature_figures.py — Publication-Quality Figures for Nature Geoscience
========================================================================
Generates 6 figures for the manuscript:
  Fig 1: SSW-avalanche coupling (Switzerland + Norway)
  Fig 2: Lag-resolved response structure
  Fig 3: Solar cycle modulation (tercile stratification)
  Fig 4: Dose-response and avalanche type specificity
  Fig 5: SOC power-law distributions
  Fig 6: Mechanistic summary schematic
"""
import sys, json, logging, gc
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyArrowPatch

sys.path.insert(0, str(Path(__file__).parent))
from _analysis_utils import PROCESSED, RESULTS, FIGURES, LOG

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

# Style
plt.rcParams.update({
    "font.size": 9,
    "font.family": "sans-serif",
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "figure.dpi": 300,
})

COLORS = {
    "swiss": "#1b9e77",
    "norway": "#d95f02",
    "significant": "#e41a1c",
    "nonsig": "#999999",
    "solar_low": "#4575b4",
    "solar_mid": "#fee090",
    "solar_high": "#d73027",
    "dry": "#8c510a",
    "wet": "#01665e",
    "ssw": "#7570b3",
}


def fig1_ssw_coupling():
    """SSW-avalanche coupling: the primary finding."""
    LOG.info("Creating Fig 1: SSW-Avalanche Coupling...")

    nature = json.loads(open(RESULTS / "nature_tier_analysis.json", encoding="utf-8").read())
    ssw_data = nature.get("section_g_ssw", {})

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.5), sharey=True)

    # Panel A: Rate ratios by time window
    ax = axes[0]
    windows = ["0_15d", "15_30d", "30_60d"]
    window_labels = ["0-15 d", "15-30 d", "30-60 d"]
    y_pos = np.arange(len(windows))

    for i, region in enumerate(["swiss", "norway"]):
        prefix = "ssw_" if region == "swiss" else "norway_ssw_"
        rrs, ci_lo, ci_hi = [], [], []
        for w in windows:
            key = prefix + w
            if key in ssw_data:
                r = ssw_data[key]
                for pname, pval in r.get("params", {}).items():
                    if "ssw" in pname:
                        rrs.append(pval.get("rate_ratio", 1))
                        ci_lo.append(pval.get("rr_ci_low", 1))
                        ci_hi.append(pval.get("rr_ci_high", 1))
                        break
                else:
                    rrs.append(1); ci_lo.append(1); ci_hi.append(1)
            else:
                rrs.append(1); ci_lo.append(1); ci_hi.append(1)

        rrs = np.array(rrs)
        errors = np.array([[rrs - np.array(ci_lo)], [np.array(ci_hi) - rrs]]).reshape(2, -1)
        offset = -0.15 if region == "swiss" else 0.15
        color = COLORS["swiss"] if region == "swiss" else COLORS["norway"]
        label = "Switzerland (SLF)" if region == "swiss" else "Norway (NVE)"
        ax.barh(y_pos + offset, rrs, 0.28, xerr=errors, color=color,
                alpha=0.8, label=label, capsize=3, ecolor="gray")

    ax.axvline(1.0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(window_labels)
    ax.set_xlabel("Rate Ratio (avalanche activity)")
    ax.set_title("a  SSW effect on avalanches", fontsize=10, fontweight="bold", loc="left")
    ax.legend(fontsize=7, loc="lower right")
    ax.set_xlim(0, 1.5)

    # Panel B: Superposed Epoch Analysis
    ax = axes[1]
    sea_data = ssw_data.get("sea_by_lag", {})
    if sea_data:
        lags = sorted([int(k) for k in sea_data.keys()])
        means = [sea_data[str(l)]["mean"] for l in lags]
        ses = [sea_data[str(l)]["se"] for l in lags]

        ax.fill_between(lags, [m - 1.96*s for m, s in zip(means, ses)],
                        [m + 1.96*s for m, s in zip(means, ses)],
                        alpha=0.2, color=COLORS["ssw"])
        ax.plot(lags, means, color=COLORS["ssw"], linewidth=1.5)
        ax.axvline(0, color="red", linewidth=1, linestyle="--", alpha=0.7, label="SSW onset")
        ax.axhline(ssw_data.get("sea_baseline", 0), color="gray", linewidth=0.8,
                   linestyle=":", alpha=0.5, label="Baseline")
        ax.set_xlabel("Days relative to SSW onset")
        ax.set_ylabel("Mean avalanche activity (AAI)")
        ax.set_title("b  Superposed epoch analysis", fontsize=10, fontweight="bold", loc="left")
        ax.legend(fontsize=7)

    plt.tight_layout()
    fig.savefig(FIGURES / "fig1_ssw_avalanche.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES / "fig1_ssw_avalanche.pdf", bbox_inches="tight")
    plt.close(fig)
    LOG.info("  Saved fig1_ssw_avalanche.png/pdf")


def fig2_lag_structure():
    """Lag-resolved response structure."""
    LOG.info("Creating Fig 2: Lag Structure...")

    nature = json.loads(open(RESULTS / "nature_tier_analysis.json", encoding="utf-8").read())
    lag_data = nature.get("section_f_lags", {}).get("lag_sweep", [])

    fig, ax = plt.subplots(figsize=(7.2, 3.0))

    lags = [r["lag"] for r in lag_data]
    rrs = [r["rr"] for r in lag_data]
    ci_lo = [r["rr_ci_low"] for r in lag_data]
    ci_hi = [r["rr_ci_high"] for r in lag_data]
    sig_fdr = [r.get("significant_fdr", False) for r in lag_data]

    colors = [COLORS["significant"] if s else COLORS["nonsig"] for s in sig_fdr]
    ax.bar(lags, rrs, color=colors, alpha=0.7, width=0.8)
    ax.errorbar(lags, rrs, yerr=[np.array(rrs) - np.array(ci_lo),
                                  np.array(ci_hi) - np.array(rrs)],
                fmt="none", ecolor="gray", capsize=2, linewidth=0.5)

    ax.axhline(1.0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
    ax.set_xlabel("Lag (days after geomagnetic event)")
    ax.set_ylabel("Rate Ratio")
    ax.set_title("Day-by-day lag response (red = significant after FDR correction)",
                 fontsize=9)

    # Annotate pathway windows
    for lo, hi, label, yoff in [(1, 3, "Fast", 0.15), (5, 21, "Stratospheric", 0.10)]:
        mid = (lo + hi) / 2
        ax.annotate(label, xy=(mid, 0.2), fontsize=7, ha="center",
                    color="blue", alpha=0.7)

    plt.tight_layout()
    fig.savefig(FIGURES / "fig2_lag_structure.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES / "fig2_lag_structure.pdf", bbox_inches="tight")
    plt.close(fig)
    LOG.info("  Saved fig2_lag_structure.png/pdf")


def fig3_solar_modulation():
    """Solar cycle modulation of the effect."""
    LOG.info("Creating Fig 3: Solar Cycle Modulation...")

    nature = json.loads(open(RESULTS / "nature_tier_analysis.json", encoding="utf-8").read())
    f107 = nature.get("section_b_f107", {})

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))

    # Panel A: Tercile stratification
    ax = axes[0]
    terciles = ["solar_low", "solar_mid", "solar_high"]
    labels = ["Low F10.7\n(<33rd %ile)", "Mid F10.7\n(33-66th)", "High F10.7\n(>66th %ile)"]
    rrs, pvals = [], []
    for t in terciles:
        key = "stratified_" + t
        if key in f107:
            r = f107[key]
            pe = r.get("params", {}).get("post_event_1_3d", {})
            rrs.append(pe.get("rate_ratio", 1))
            pvals.append(pe.get("p", 1))
        else:
            rrs.append(1); pvals.append(1)

    bar_colors = [COLORS[t] for t in terciles]
    bars = ax.bar(range(3), rrs, color=bar_colors, alpha=0.8, edgecolor="black", linewidth=0.5)
    ax.axhline(1.0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
    ax.set_xticks(range(3))
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("Rate Ratio")
    ax.set_title("a  Effect by solar activity level", fontsize=10, fontweight="bold", loc="left")

    # Add significance stars
    for i, (rr, p) in enumerate(zip(rrs, pvals)):
        star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."
        ax.text(i, rr + 0.05, star, ha="center", fontsize=8)

    # Panel B: Primary model progression
    ax = axes[1]
    primary = nature.get("section_a_primary", {})
    models = ["model1_minimal", "model2_meteorological", "model3_solar_cycle",
              "model4_snow_cover", "model5_full"]
    model_labels = ["Minimal", "+Meteo", "+F10.7", "+Snow", "Full\n(18 vars)"]
    model_rrs, model_ps = [], []
    for m in models:
        if m in primary:
            pe = primary[m].get("params", {}).get("post_event_1_3d", {})
            model_rrs.append(pe.get("rate_ratio", 1))
            model_ps.append(pe.get("p", 1))
        else:
            model_rrs.append(1); model_ps.append(1)

    colors_m = [COLORS["significant"] if p < 0.05 else COLORS["nonsig"] for p in model_ps]
    ax.bar(range(len(models)), model_rrs, color=colors_m, alpha=0.7, edgecolor="black", linewidth=0.5)
    ax.axhline(1.0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
    ax.set_xticks(range(len(models)))
    ax.set_xticklabels(model_labels, fontsize=7)
    ax.set_ylabel("Rate Ratio")
    ax.set_title("b  Effect stability across models", fontsize=10, fontweight="bold", loc="left")

    for i, (rr, p) in enumerate(zip(model_rrs, model_ps)):
        star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."
        ax.text(i, rr + 0.03, star, ha="center", fontsize=7)

    plt.tight_layout()
    fig.savefig(FIGURES / "fig3_solar_modulation.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES / "fig3_solar_modulation.pdf", bbox_inches="tight")
    plt.close(fig)
    LOG.info("  Saved fig3_solar_modulation.png/pdf")


def fig4_dose_response():
    """Dose-response and avalanche type specificity."""
    LOG.info("Creating Fig 4: Dose-Response...")

    nature = json.loads(open(RESULTS / "nature_tier_analysis.json", encoding="utf-8").read())
    dose = nature.get("section_e_dose_response", {})
    multi = nature.get("section_d_multi_region", {})

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))

    # Panel A: Dose-response by Dst
    ax = axes[0]
    categories = ["moderate", "strong"]
    cat_labels = ["Moderate\n(Dst -100 to -50)", "Strong\n(Dst -200 to -100)"]
    cat_rrs = []
    cat_ps = []
    for c in categories:
        if c in dose:
            pe = dose[c].get("params", {}).get("dose_" + c, {})
            cat_rrs.append(pe.get("rate_ratio", 1))
            cat_ps.append(pe.get("p", 1))
        else:
            cat_rrs.append(1); cat_ps.append(1)

    colors = [COLORS["significant"] if p < 0.05 else COLORS["nonsig"] for p in cat_ps]
    bars = ax.bar(range(len(categories)), cat_rrs, color=colors, alpha=0.7,
                  edgecolor="black", linewidth=0.5)
    ax.axhline(1.0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels(cat_labels, fontsize=7)
    ax.set_ylabel("Rate Ratio")
    ax.set_title("a  Dose-response by storm intensity", fontsize=10, fontweight="bold", loc="left")
    for i, (rr, p) in enumerate(zip(cat_rrs, cat_ps)):
        star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."
        y_pos = max(rr + 0.03, 0.12)
        ax.text(i, y_pos, "RR=%.3f\n%s" % (rr, star), ha="center", fontsize=7)

    # Panel B: Avalanche type specificity
    ax = axes[1]
    types = ["switzerland_dry_natural", "switzerland_wet_natural", "switzerland_all_natural"]
    type_labels = ["Dry\nnatural", "Wet\nnatural", "All\nnatural"]
    type_colors = [COLORS["dry"], COLORS["wet"], "gray"]
    type_rrs, type_ps = [], []
    for t in types:
        if t in multi:
            pe = multi[t].get("params", {}).get("post_event_1_3d", {})
            type_rrs.append(pe.get("rate_ratio", 1))
            type_ps.append(pe.get("p", 1))
        else:
            type_rrs.append(1); type_ps.append(1)

    bars = ax.bar(range(len(types)), type_rrs, color=type_colors, alpha=0.7,
                  edgecolor="black", linewidth=0.5)
    ax.axhline(1.0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
    ax.set_xticks(range(len(types)))
    ax.set_xticklabels(type_labels, fontsize=8)
    ax.set_ylabel("Rate Ratio")
    ax.set_title("b  Effect by avalanche type", fontsize=10, fontweight="bold", loc="left")
    for i, (rr, p) in enumerate(zip(type_rrs, type_ps)):
        star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."
        ax.text(i, rr + 0.03, star, ha="center", fontsize=8)

    plt.tight_layout()
    fig.savefig(FIGURES / "fig4_dose_type.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES / "fig4_dose_type.pdf", bbox_inches="tight")
    plt.close(fig)
    LOG.info("  Saved fig4_dose_type.png/pdf")


def fig5_soc():
    """SOC power-law distributions."""
    LOG.info("Creating Fig 5: SOC Power Laws...")

    # Load flare data
    fl = pd.read_parquet(PROCESSED / "solar" / "flares.parquet")
    def class_to_flux(ct):
        if pd.isna(ct) or not isinstance(ct, str):
            return np.nan
        ct = ct.strip()
        cls = ct[0].upper()
        try:
            mag = float(ct[1:])
        except (ValueError, IndexError):
            mag = 1.0
        base = {"A": 1e-8, "B": 1e-7, "C": 1e-6, "M": 1e-5, "X": 1e-4}.get(cls, np.nan)
        return base * mag if not np.isnan(base) else np.nan
    fl["flux"] = fl["classType"].apply(class_to_flux)
    fluxes = fl["flux"].dropna().values
    fluxes = fluxes[fluxes > 0]
    del fl; gc.collect()

    # Load avalanche data
    act = pd.read_parquet(PROCESSED / "cryosphere" / "slf_activity.parquet")
    aai = act["aai_all_natural"].dropna().values
    aai = aai[aai > 0]
    del act; gc.collect()

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))

    # Panel A: Solar flare CCDF
    ax = axes[0]
    sorted_flux = np.sort(fluxes)[::-1]
    ccdf_y = np.arange(1, len(sorted_flux) + 1) / len(sorted_flux)
    ax.loglog(sorted_flux, ccdf_y, ".", markersize=1, alpha=0.3, color="darkorange")

    # Fit line
    nature = json.loads(open(RESULTS / "nature_tier_analysis.json", encoding="utf-8").read())
    soc = nature.get("section_h_soc", {})
    alpha_f = soc.get("flare_power_law", {}).get("alpha", 2.27)
    xmin_f = soc.get("flare_power_law", {}).get("xmin", 1e-5)
    x_fit = np.logspace(np.log10(xmin_f), np.log10(sorted_flux[0]), 100)
    y_fit = (x_fit / xmin_f) ** (-(alpha_f - 1))
    y_fit *= ccdf_y[np.searchsorted(-sorted_flux, -xmin_f)] if xmin_f < sorted_flux[0] else ccdf_y[0]
    ax.loglog(x_fit, y_fit, "r-", linewidth=1.5, label=r"$\alpha$ = %.2f" % alpha_f)

    ax.set_xlabel("Peak X-ray flux (W/m$^2$)")
    ax.set_ylabel("P(X > x)")
    ax.set_title(r"a  Solar flares ($\alpha$ = %.2f)" % alpha_f,
                 fontsize=10, fontweight="bold", loc="left")
    ax.legend(fontsize=8)

    # Panel B: Avalanche CCDF
    ax = axes[1]
    sorted_aai = np.sort(aai)[::-1]
    ccdf_aai = np.arange(1, len(sorted_aai) + 1) / len(sorted_aai)
    ax.loglog(sorted_aai, ccdf_aai, ".", markersize=1, alpha=0.3, color=COLORS["swiss"])

    alpha_a = soc.get("avalanche_power_law", {}).get("alpha", 2.08)
    xmin_a = soc.get("avalanche_power_law", {}).get("xmin", 11)
    x_fit_a = np.logspace(np.log10(xmin_a), np.log10(sorted_aai[0]), 100)
    y_fit_a = (x_fit_a / xmin_a) ** (-(alpha_a - 1))
    idx = np.searchsorted(-sorted_aai, -xmin_a)
    if idx < len(ccdf_aai):
        y_fit_a *= ccdf_aai[idx]
    ax.loglog(x_fit_a, y_fit_a, "r-", linewidth=1.5, label=r"$\alpha$ = %.2f" % alpha_a)

    ax.set_xlabel("Daily avalanche activity index")
    ax.set_ylabel("P(X > x)")
    ax.set_title(r"b  Avalanches ($\alpha$ = %.2f)" % alpha_a,
                 fontsize=10, fontweight="bold", loc="left")
    ax.legend(fontsize=8)

    plt.tight_layout()
    fig.savefig(FIGURES / "fig5_soc_power_laws.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES / "fig5_soc_power_laws.pdf", bbox_inches="tight")
    plt.close(fig)
    LOG.info("  Saved fig5_soc_power_laws.png/pdf")


def fig6_mechanism():
    """Mechanistic summary schematic."""
    LOG.info("Creating Fig 6: Mechanism Summary...")

    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7)
    ax.axis("off")

    # Boxes for the causal chain
    boxes = [
        (1.0, 6.0, "Solar Flare / CME", "darkorange"),
        (1.0, 4.5, "Geomagnetic Storm\n(Dst < -100 nT)", "firebrick"),
        (5.0, 4.5, "EPP into\nPolar Atmosphere", "mediumpurple"),
        (5.0, 3.0, "Stratospheric Chemistry\n(HNO3 increase, P<0.001)", "steelblue"),
        (5.0, 1.5, "SSW / Vortex\nDisruption", "navy"),
        (1.0, 1.5, "Tropospheric Response\n(NAO shift, blocking)", "teal"),
        (1.0, 0.0, "Avalanche Activity\nDECREASE (RR=0.44-0.70)", "darkgreen"),
    ]

    for x, y, text, color in boxes:
        rect = plt.Rectangle((x - 0.8, y - 0.35), 2.8, 0.7,
                              facecolor=color, alpha=0.15, edgecolor=color, linewidth=1.5)
        ax.add_patch(rect)
        ax.text(x + 0.6, y, text, ha="center", va="center", fontsize=7,
                fontweight="bold", color=color)

    # Arrows
    arrows = [
        ((2.4, 5.65), (2.4, 4.85), "1-3 d"),        # Flare -> Storm
        ((3.8, 4.5), (4.2, 4.5), ""),                # Storm -> EPP
        ((6.4, 4.15), (6.4, 3.35), "0-5 d"),         # EPP -> Chemistry
        ((6.4, 2.65), (6.4, 1.85), "5-21 d"),        # Chemistry -> SSW
        ((4.2, 1.5), (3.8, 1.5), "15-60 d"),         # SSW -> Troposphere
        ((2.4, 1.15), (2.4, 0.35), ""),               # Trop -> Avalanche
    ]

    for start, end, label in arrows:
        ax.annotate("", xy=end, xytext=start,
                    arrowprops=dict(arrowstyle="->", lw=1.5, color="gray"))
        if label:
            mid_x = (start[0] + end[0]) / 2
            mid_y = (start[1] + end[1]) / 2
            ax.text(mid_x + 0.3, mid_y, label, fontsize=6, color="gray",
                    fontstyle="italic")

    # Key statistics box
    stats_text = (
        "KEY STATISTICS\n"
        "---\n"
        "SSW -> Swiss aval:  RR=0.44 (P<0.001)\n"
        "SSW -> Norway aval: RR=0.30 (P<0.001)\n"
        "Strong storm dose:  RR=0.06 (P=0.004)\n"
        "Dry aval specific:  RR=0.70 (P=0.0001)\n"
        "Solar mid only:     RR=0.60 (P=0.0005)\n"
        "Full model (18var): RR=0.68 (P=0.013)"
    )
    ax.text(8.0, 6.0, stats_text, fontsize=6, fontfamily="monospace",
            va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow",
                      edgecolor="gray", alpha=0.9))

    # Title
    ax.text(5, 6.8, "Proposed Mechanistic Chain: Solar Activity → Avalanche Modulation",
            ha="center", va="center", fontsize=11, fontweight="bold")

    fig.savefig(FIGURES / "fig6_mechanism.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES / "fig6_mechanism.pdf", bbox_inches="tight")
    plt.close(fig)
    LOG.info("  Saved fig6_mechanism.png/pdf")


def main():
    fig1_ssw_coupling()
    gc.collect()
    fig2_lag_structure()
    gc.collect()
    fig3_solar_modulation()
    gc.collect()
    fig4_dose_response()
    gc.collect()
    fig5_soc()
    gc.collect()
    fig6_mechanism()
    gc.collect()

    print("\nAll 6 publication figures saved to data/figures/")
    print("  fig1_ssw_avalanche.png/pdf")
    print("  fig2_lag_structure.png/pdf")
    print("  fig3_solar_modulation.png/pdf")
    print("  fig4_dose_type.png/pdf")
    print("  fig5_soc_power_laws.png/pdf")
    print("  fig6_mechanism.png/pdf")


if __name__ == "__main__":
    main()
