#!/usr/bin/env python3
"""Generate Extended Data geographic figure showing EAWS centre-level
SSW-window danger-level anomalies.

Reads the EAWS analysis results and produces a schematic map of Alpine
centres with their directional SSW response.  Muted colours and
prominent uncertainty annotation (n_eff ≈ 2) communicate that this is
hypothesis-generating, not inferential.

Output: data/figures/extended_eaws_map.pdf
"""

import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]

# EAWS centre approximate coordinates (lon, lat) and SSW-window anomaly
# direction from the two-event analysis.  Values derived from
# scripts/review_rounds/r42_eaws_geographic.py results.
# direction: -1 = suppression, +1 = enhancement, 0 = mixed/neutral
CENTRES = [
    # Country, Region, lon, lat, direction
    ("CH", "Bern", 7.5, 46.9, -1),
    ("CH", "Graubünden", 9.8, 46.8, -1),
    ("CH", "Valais", 7.6, 46.2, -1),
    ("CH", "Central CH", 8.6, 46.9, -1),
    ("AT", "Tyrol", 11.4, 47.3, +1),
    ("AT", "Vorarlberg", 9.9, 47.2, 0),
    ("AT", "Salzburg", 13.1, 47.3, +1),
    ("AT", "Kärnten", 13.8, 46.7, +1),
    ("AT", "Steiermark", 15.0, 47.3, +1),
    ("AT", "Oberösterreich", 14.0, 47.8, 0),
    ("AT", "Niederösterreich", 15.8, 47.8, 0),
    ("IT", "South Tyrol", 11.4, 46.7, +1),
    ("IT", "Trentino", 11.1, 46.1, +1),
    ("IT", "Lombardia", 9.5, 46.2, -1),
    ("IT", "Veneto", 12.0, 46.4, 0),
    ("IT", "Friuli", 13.2, 46.4, +1),
    ("IT", "Piemonte", 7.7, 45.1, -1),
    ("IT", "Valle d'Aosta", 7.3, 45.7, -1),
    ("FR", "Northern Alps", 6.1, 45.5, -1),
    ("FR", "Southern Alps", 6.5, 44.5, -1),
    ("DE", "Bavaria", 11.5, 47.5, +1),
    ("SI", "Slovenia", 14.0, 46.2, +1),
]


def main():
    fig, ax = plt.subplots(figsize=(10, 7))

    # Background Alpine outline (simplified polygon)
    alpine_lon = [5.5, 7.0, 10.0, 13.5, 16.5, 16.0, 14.0, 10.0, 7.0, 5.5]
    alpine_lat = [44.0, 43.8, 45.5, 46.0, 47.5, 48.0, 47.5, 47.8, 47.5, 44.0]
    ax.fill(alpine_lon, alpine_lat, color="#f0f0f0", edgecolor="#cccccc",
            linewidth=1, zorder=0)

    # Plot centres with muted colours
    colours = {-1: "#6699cc", +1: "#cc8866", 0: "#999999"}
    labels = {-1: "Suppression", +1: "Enhancement", 0: "Mixed / neutral"}
    marker_size = 120

    for country, region, lon, lat, direction in CENTRES:
        ax.scatter(lon, lat, c=colours[direction], s=marker_size,
                   edgecolors="white", linewidths=0.8, zorder=2,
                   marker="o")
        ax.annotate(region, (lon, lat), fontsize=5.5,
                    ha="center", va="bottom", xytext=(0, 6),
                    textcoords="offset points", zorder=3)

    # Country borders (very simplified)
    for country_code, colour in [("CH", "#333"), ("AT", "#333"),
                                  ("FR", "#333"), ("IT", "#333"),
                                  ("DE", "#333")]:
        pass  # borders omitted for clarity; the centres carry the information

    # Legend
    handles = [mpatches.Patch(facecolor=colours[d], edgecolor="white",
                              label=labels[d]) for d in [-1, +1, 0]]
    ax.legend(handles=handles, loc="lower left", fontsize=8,
              framealpha=0.9, title="SSW-window anomaly direction",
              title_fontsize=8)

    # Prominent uncertainty box
    textstr = ("Based on 2 SSW events only ($n_{\\rm eff}$ ≈ 2)\n"
               "Hypothesis-generating; no inferential weight")
    props = dict(boxstyle="round,pad=0.5", facecolor="#fff3cd",
                 edgecolor="#856404", alpha=0.95)
    ax.text(0.98, 0.02, textstr, transform=ax.transAxes, fontsize=8,
            verticalalignment="bottom", horizontalalignment="right",
            bbox=props, style="italic")

    ax.set_xlabel("Longitude (°E)", fontsize=9)
    ax.set_ylabel("Latitude (°N)", fontsize=9)
    ax.set_title("EAWS centre-level SSW-window danger anomalies (2 events)",
                 fontsize=11, fontweight="bold")
    ax.set_xlim(4.5, 17)
    ax.set_ylim(43.5, 48.5)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3, linewidth=0.5)

    outpath = ROOT / "data" / "figures" / "extended_eaws_map.pdf"
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, bbox_inches="tight", dpi=300)
    print(f"Saved: {outpath}")
    plt.close()


if __name__ == "__main__":
    main()
