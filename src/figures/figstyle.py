"""Shared figure style for the four main-text figures: colorblind-safe (Okabe-Ito)
palette, consistent fonts, and a saver that writes both vector (PDF) and raster (PNG at
Nature Methods resolution). Kept tiny and dependency-light so every figure script imports
the same look."""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIGURES = os.path.join(ROOT, "figures")

# Okabe-Ito colorblind-safe palette
OI = {
    "black": "#000000", "orange": "#E69F00", "skyblue": "#56B4E9", "green": "#009E73",
    "yellow": "#F0E442", "blue": "#0072B2", "vermillion": "#D55E00", "purple": "#CC79A7",
    "grey": "#999999",
}
ORACLE = OI["vermillion"]
NAIVE = OI["blue"]
METHOD = {"Baysor": OI["green"], "Proseg": OI["orange"], "ComSeg": OI["purple"]}


def apply_style():
    plt.rcParams.update({
        "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.8, "figure.dpi": 120, "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03, "pdf.fonttype": 42, "ps.fonttype": 42,
        "font.family": "DejaVu Sans",
    })


def panel_tag(ax, tag, dx=-0.02, dy=1.04):
    ax.text(dx, dy, tag, transform=ax.transAxes, fontsize=11, fontweight="bold",
            va="bottom", ha="right")


def save(fig, name):
    os.makedirs(FIGURES, exist_ok=True)
    png = os.path.join(FIGURES, name + ".png")
    pdf = os.path.join(FIGURES, name + ".pdf")
    fig.savefig(png, dpi=400)
    fig.savefig(pdf)
    print(f"wrote {png} and {pdf}")
