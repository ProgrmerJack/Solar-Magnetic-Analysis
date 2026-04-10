"""
08_figures.py — Generate publication-quality figures
=====================================================
All figures for the Nature Geoscience manuscript.
"""
import sys
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyArrowPatch
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).parent))
from _analysis_utils import PROCESSED, RESULTS, FIGURES, load_panel

import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

# Nature style
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 8,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})


def figure1_event_study(panel, out):
    """Event-study plot: avalanche activity anomaly centered on event day 0."""
    events = panel[(panel["geo_event"] == 1) & (panel["is_winter"] == 1)]
    event_dates = events.index

    # Compute anomalies (subtract winter seasonal mean by day_of_season)
    panel_w = panel[panel["is_winter"] == 1].copy()
    seasonal = panel_w.groupby("day_of_season")["aai_all_natural"].mean()
    panel_w["aai_anom"] = panel_w["aai_all_natural"] - panel_w["day_of_season"].map(seasonal)

    lags = np.arange(-15, 31)
    epoch_matrix = []
    for ed in event_dates:
        row = []
        for lag in lags:
            target = ed + pd.Timedelta(days=int(lag))
            if target in panel_w.index and not np.isnan(panel_w.loc[target, "aai_anom"]):
                row.append(panel_w.loc[target, "aai_anom"])
            else:
                row.append(np.nan)
        epoch_matrix.append(row)

    epoch = np.array(epoch_matrix)
    mean_prof = np.nanmean(epoch, axis=0)
    sem_prof = np.nanstd(epoch, axis=0, ddof=1) / np.sqrt(np.sum(~np.isnan(epoch), axis=0))

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3), gridspec_kw={"width_ratios": [3, 1]})

    # Panel A: Epoch plot
    ax = axes[0]
    ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
    ax.axvline(0, color="red", linewidth=0.5, alpha=0.5, label="Event day")
    ax.axvspan(1, 3, alpha=0.15, color="blue", label="Fast pathway (d1–3)")
    ax.axvspan(5, 21, alpha=0.08, color="green", label="Strat. pathway (d5–21)")

    ax.plot(lags, mean_prof, color="black", linewidth=1)
    ax.fill_between(lags, mean_prof - 1.96 * sem_prof, mean_prof + 1.96 * sem_prof,
                    alpha=0.2, color="steelblue")

    ax.set_xlabel("Days relative to geomagnetic event")
    ax.set_ylabel("Avalanche activity anomaly")
    ax.set_title("a  Superposed epoch: avalanche activity", loc="left", fontweight="bold")
    ax.legend(fontsize=6, loc="upper right")

    # Panel B: Rate ratios
    ax2 = axes[1]
    results = json.loads((RESULTS / "primary_endpoint.json").read_text())
    models = [
        ("Fast 1-3d", "fast_pathway_1_3d"),
        ("Strat 5-21d", "primary_strat_5_21d"),
        ("Strong events", "dose_strong"),
    ]

    y_pos = []
    for i, (label, key) in enumerate(models):
        r = results.get(key, {})
        rr = r.get("rate_ratio", 1)
        lo = r.get("rr_ci_lower", rr)
        hi = r.get("rr_ci_upper", rr)
        p = r.get("p_value", 1)
        y = len(models) - i - 1
        y_pos.append(y)

        color = "red" if p < 0.05 else "gray"
        ax2.errorbar(rr, y, xerr=[[rr - lo], [hi - rr]], fmt="o",
                     color=color, markersize=4, capsize=3, linewidth=1)
        sig = "**" if p < 0.01 else "*" if p < 0.05 else ""
        ax2.text(hi + 0.02, y, f" p={p:.3f}{sig}", fontsize=6, va="center")

    ax2.axvline(1, color="gray", linewidth=0.5, linestyle="--")
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([m[0] for m in models], fontsize=7)
    ax2.set_xlabel("Rate ratio")
    ax2.set_title("b  Rate ratios", loc="left", fontweight="bold")
    ax2.set_xlim(0.5, 1.5)

    plt.tight_layout()
    fig.savefig(out)
    plt.close()
    print(f"  Saved: {out.name}")


def figure2_falsification(panel, out):
    """Falsification suite summary figure."""
    fals = json.loads((RESULTS / "falsification.json").read_text())
    results = json.loads((RESULTS / "primary_endpoint.json").read_text())

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5))

    # Panel A: Sham event distribution
    ax = axes[0, 0]
    sham = fals.get("sham_events", {})
    real_rr = sham.get("real_rate_ratio", 1)
    ax.axvline(real_rr, color="red", linewidth=1.5, label=f"Real (RR={real_rr:.3f})")
    ax.axvline(1, color="gray", linewidth=0.5, linestyle="--")
    ax.set_xlabel("Rate ratio (5-21d window)")
    ax.set_ylabel("Frequency")
    ax.set_title("a  Sham events (5-21d)", loc="left", fontweight="bold")
    ax.legend(fontsize=6)
    ax.text(0.05, 0.85, f"Rank: {sham.get('real_beats_pct_of_shams', 0):.1f}%ile",
            transform=ax.transAxes, fontsize=7)

    # Panel B: LOWO cross-validation
    ax = axes[0, 1]
    lowo = fals.get("lowo_cv", {})
    rr_by_w = lowo.get("rr_by_excluded_winter", {})
    if rr_by_w:
        winters = sorted(rr_by_w.keys())
        rrs = [rr_by_w[w] for w in winters]
        ax.barh(range(len(winters)), rrs, color="steelblue", edgecolor="white", linewidth=0.3)
        ax.axvline(lowo.get("mean_rr", 1), color="red", linewidth=1, linestyle="--",
                   label=f"Mean={lowo.get('mean_rr', 1):.3f}")
        ax.axvline(1, color="gray", linewidth=0.5, linestyle="--")
        ax.set_yticks(range(len(winters)))
        ax.set_yticklabels([w.replace("/", "\n") for w in winters], fontsize=5)
        ax.set_xlabel("Rate ratio (leave-one-out)")
        ax.legend(fontsize=6)
    ax.set_title("b  Leave-one-winter-out", loc="left", fontweight="bold")

    # Panel C: Outcome comparison (natural vs accidents)
    ax = axes[1, 0]
    categories = []
    rr_vals = []
    ci_lo = []
    ci_hi = []
    colors = []

    for label, key, color in [
        ("Natural (5-21d)", "primary_strat_5_21d", "steelblue"),
        ("Natural (1-3d)", "fast_pathway_1_3d", "navy"),
        ("Accidents", None, "orange"),
    ]:
        if key and key in results:
            r = results[key]
        elif key is None:
            r = fals.get("control_outcome", {})
        else:
            continue
        rr = r.get("rate_ratio", 1)
        lo = r.get("rr_ci_lower", rr - 0.1)
        hi = r.get("rr_ci_upper", rr + 0.1)
        categories.append(label)
        rr_vals.append(rr)
        ci_lo.append(rr - lo)
        ci_hi.append(hi - rr)
        colors.append(color)

    if categories:
        y_pos = range(len(categories))
        ax.errorbar(rr_vals, y_pos, xerr=[ci_lo, ci_hi], fmt="o",
                    color="black", markersize=4, capsize=3, linewidth=1)
        for i, c in enumerate(colors):
            ax.plot(rr_vals[i], i, "o", color=c, markersize=6, zorder=5)
        ax.axvline(1, color="gray", linewidth=0.5, linestyle="--")
        ax.set_yticks(list(y_pos))
        ax.set_yticklabels(categories, fontsize=7)
        ax.set_xlabel("Rate ratio")
    ax.set_title("c  Outcome comparison", loc="left", fontweight="bold")

    # Panel D: Regional comparison
    ax = axes[1, 1]
    regions = []
    for label, data, color in [
        ("Swiss Alps\n(natural)", results.get("fast_pathway_1_3d", {}), "steelblue"),
        ("Norway\n(danger level)", fals.get("control_region", {}), "forestgreen"),
    ]:
        rr = data.get("rate_ratio", 1)
        if not np.isnan(rr):
            lo = data.get("rr_ci_lower", rr - 0.1)
            hi = data.get("rr_ci_upper", rr + 0.1)
            p = data.get("p_value", 1)
            regions.append((label, rr, rr - lo, hi - rr, color, p))

    for i, (label, rr, lo, hi, color, p) in enumerate(regions):
        ax.errorbar(rr, i, xerr=[[lo], [hi]], fmt="o", color=color,
                    markersize=6, capsize=3, linewidth=1)
        sig = "**" if p < 0.01 else "*" if p < 0.05 else ""
        ax.text(rr + hi + 0.02, i, f"p={p:.3f}{sig}", fontsize=6, va="center")

    ax.axvline(1, color="gray", linewidth=0.5, linestyle="--")
    if regions:
        ax.set_yticks(range(len(regions)))
        ax.set_yticklabels([r[0] for r in regions], fontsize=7)
    ax.set_xlabel("Rate ratio (1-3d window)")
    ax.set_title("d  Regional replication", loc="left", fontweight="bold")

    plt.tight_layout()
    fig.savefig(out)
    plt.close()
    print(f"  Saved: {out.name}")


def figure3_mechanism(out):
    """Schematic of proposed mechanism with results annotated."""
    fig, ax = plt.subplots(1, 1, figsize=(7.2, 2.5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3)
    ax.axis("off")

    boxes = [
        (0.8, 1.5, "Solar\nDisturbance"),
        (3.0, 1.5, "Geomagnetic\nStorm\n(Kp≥5/Dst≤-50)"),
        (5.3, 2.2, "Stratospheric\nPerturbation\n(5-21d)"),
        (5.3, 0.8, "Fast Atmospheric\nResponse\n(1-3d)"),
        (8.0, 1.5, "Alpine\nAvalanche\nActivity"),
    ]

    for x, y, text in boxes:
        bbox = dict(boxstyle="round,pad=0.3", facecolor="lightblue", edgecolor="navy", linewidth=0.5)
        ax.text(x, y, text, ha="center", va="center", fontsize=7, bbox=bbox)

    # Arrows
    ax.annotate("", xy=(2.2, 1.5), xytext=(1.5, 1.5),
                arrowprops=dict(arrowstyle="->", color="black", lw=1))
    ax.annotate("", xy=(4.6, 2.2), xytext=(3.7, 1.8),
                arrowprops=dict(arrowstyle="->", color="gray", lw=0.8, linestyle="--"))
    ax.annotate("", xy=(4.6, 0.8), xytext=(3.7, 1.2),
                arrowprops=dict(arrowstyle="->", color="black", lw=1))
    ax.annotate("", xy=(7.2, 1.8), xytext=(6.1, 2.2),
                arrowprops=dict(arrowstyle="->", color="gray", lw=0.8, linestyle="--"))
    ax.annotate("", xy=(7.2, 1.2), xytext=(6.1, 0.8),
                arrowprops=dict(arrowstyle="->", color="red", lw=1.5))

    # Result annotations
    ax.text(5.3, 2.8, "NOT SIGNIFICANT\n(RR=0.96, p=0.63)", ha="center",
            fontsize=6, color="gray", style="italic")
    ax.text(5.3, 0.2, "SIGNIFICANT DECREASE\n(RR=0.77, p=0.008)", ha="center",
            fontsize=6, color="red", fontweight="bold")

    ax.set_title("Hypothesized pathways and empirical results", fontsize=9, fontweight="bold")
    fig.savefig(out)
    plt.close()
    print(f"  Saved: {out.name}")


def figure4_temporal(panel, out):
    """Temporal context: event distribution and activity time series."""
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 5), sharex=True)

    pw = panel[(panel["is_winter"] == 1) & (panel["aai_all_natural"].notna())].copy()

    # Panel A: Kp/Dst time series
    ax = axes[0]
    ax.plot(pw.index, pw["kp_max"], color="steelblue", linewidth=0.3, alpha=0.5)
    events = pw[pw["geo_event"] == 1]
    ax.scatter(events.index, events["kp_max"], color="red", s=8, zorder=5, label="Events")
    ax.set_ylabel("Daily Kp max")
    ax.axhline(5, color="red", linewidth=0.5, linestyle="--", alpha=0.5)
    ax.legend(fontsize=6)
    ax.set_title("a  Geomagnetic activity (winter days)", loc="left", fontweight="bold")

    # Panel B: Avalanche activity
    ax = axes[1]
    ax.plot(pw.index, pw["aai_all_natural"], color="saddlebrown", linewidth=0.3, alpha=0.5)
    rolling = pw["aai_all_natural"].rolling(30, center=True, min_periods=5).mean()
    ax.plot(pw.index, rolling, color="saddlebrown", linewidth=1, label="30d rolling mean")
    ax.set_ylabel("Natural AAI")
    ax.legend(fontsize=6)
    ax.set_title("b  Swiss natural avalanche activity", loc="left", fontweight="bold")

    # Panel C: Event markers + post-event windows
    ax = axes[2]
    ax.fill_between(pw.index, pw["post_event_1_3d"], color="blue", alpha=0.3, label="1-3d post-event")
    ax.fill_between(pw.index, pw["post_event_5_21d"] * 0.5, color="green", alpha=0.2, label="5-21d post-event")
    ax.set_ylabel("Post-event\nindicator")
    ax.set_xlabel("Date")
    ax.legend(fontsize=6)
    ax.set_title("c  Post-event exposure windows", loc="left", fontweight="bold")

    plt.tight_layout()
    fig.savefig(out)
    plt.close()
    print(f"  Saved: {out.name}")


def main():
    panel = load_panel(winter_only=False)
    FIGURES.mkdir(parents=True, exist_ok=True)

    print("Generating publication figures...")
    figure1_event_study(panel, FIGURES / "fig1_event_study.png")
    figure2_falsification(panel, FIGURES / "fig2_falsification.png")
    figure3_mechanism(FIGURES / "fig3_mechanism.png")
    figure4_temporal(panel, FIGURES / "fig4_temporal.png")
    print(f"\nAll figures saved to: {FIGURES}")


if __name__ == "__main__":
    main()
