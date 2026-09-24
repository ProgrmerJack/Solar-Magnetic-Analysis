"""
nature_style.py
===============
Shared look for every display item, so the figures agree with each other and
with Nature Portfolio's artwork guide: sans-serif type at 5-7 pt, 89 mm single
and 183 mm double column widths, vector PDF plus a 300-dpi PNG preview, and the
Okabe-Ito palette, which stays distinguishable under the common colour-vision
deficiencies. Liberation Sans is metric-compatible with Arial.

Every figure script reads results from results/current/ (never recomputes a
statistic) and writes to 09_figures/out/.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "current"
OUT = Path(__file__).resolve().parent / "out"

MM = 1 / 25.4
SINGLE = 89 * MM
DOUBLE = 183 * MM

# Okabe-Ito
C = {"orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73",
     "yellow": "#F0E442", "blue": "#0072B2", "vermillion": "#D55E00",
     "purple": "#CC79A7", "black": "#000000", "grey": "#999999"}
EVENT = {"Feb-2018": C["vermillion"], "Jan-2019": C["blue"]}
ARM = {"nudged": C["vermillion"], "control": C["grey"]}


def apply():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 6.5, "axes.titlesize": 7, "axes.labelsize": 6.5,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
        "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "ps.fonttype": 42,     # embed TrueType, editable text
        "savefig.dpi": 300, "figure.dpi": 150,
        "legend.frameon": False,
    })


def panel_label(ax, s, x=-0.12, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=8, fontweight="bold",
            va="bottom", ha="left")


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight")
    print(f"-> {OUT.relative_to(ROOT)}/{name}.pdf / .png")
