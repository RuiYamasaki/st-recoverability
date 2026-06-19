"""Figure 1: the answerability frontier and the oracle ceiling.

(a) generator schematic: a small synthetic field (Voronoi cells, nuclei, a few
    transcripts with true -> observed displacement arrows).
(b) oracle construction: one cell's indicator mask convolved with the displacement
    kernel (the per-cell density that the Bayes-optimal assignment maximises).
(c) the answerability frontier: oracle assignment accuracy over packing x displacement
    at the median density, with the full-grid range annotated.
(d) oracle minus naive at all 125 grid points (oracle >= naive everywhere).

Every value from results/sweep.csv (Gate 0). Panels a,b regenerate tiny synthetic fields
with recorded seeds for illustration only (no result value depends on them).
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.ndimage import gaussian_filter

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402
from generator import build_field, generate_transcripts  # noqa: E402
from figstyle import apply_style, save, panel_tag, OI, ORACLE, NAIVE  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS = os.path.join(ROOT, "results")
SEED = config.MASTER_SEED + 940000


def panel_a(ax):
    f = build_field(5700.0, 4.0, SEED + 1, n_target=36)
    tx = generate_transcripts(f, 30.0, SEED + 2)
    # colour each pixel by its cell, light pastel
    rng = np.random.default_rng(SEED)
    colors = rng.uniform(0.55, 0.95, size=(f.n_cells, 3))
    img = colors[f.label]
    ax.imshow(img, origin="lower", extent=[0, f.L_um, 0, f.L_um], interpolation="nearest")
    ax.plot(f.centers[:, 0], f.centers[:, 1], ".", color="0.15", ms=4)
    # a handful of transcripts with displacement arrows
    sel = np.linspace(0, tx.obs_xy.shape[0] - 1, 14).astype(int)
    for i in sel:
        x0, y0 = tx.true_xy[i]; x1, y1 = tx.obs_xy[i]
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="->", color=OI["vermillion"], lw=0.8))
    ax.set_xlim(0, f.L_um); ax.set_ylim(0, f.L_um)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("a  Generator: cells, nuclei, displaced transcripts", loc="left", fontsize=8.5)
    ax.set_xlabel("Voronoi cells (shaded) - nuclei (dots) - displacement (arrows)", fontsize=6.5)


def panel_b(ax):
    grid = 90
    mask = np.zeros((grid, grid))
    mask[30:60, 30:60] = 1.0  # one square cell footprint
    dens = gaussian_filter(mask, 6.0, mode="constant")
    im = ax.imshow(dens, origin="lower", cmap="magma", interpolation="bilinear")
    ax.contour(mask, levels=[0.5], colors="white", linewidths=0.8)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("b  Oracle: cell footprint blurred by\nthe displacement kernel", loc="left", fontsize=8.5)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("per-cell density", fontsize=6.5); cb.ax.tick_params(labelsize=6)
    ax.set_xlabel(r"argmax over cells of $p(g|t)\,(f_c * \phi_\sigma)$", fontsize=6.5)


def panel_c(ax, df):
    packs = np.sort(df.packing_cells_per_mm2.unique())
    sigs = np.sort(df.sigma_um.unique())
    dens = np.sort(df.mean_tx_per_cell.unique())
    dmid = dens[len(dens) // 2]
    sub = df[df.mean_tx_per_cell == dmid]
    M = np.full((len(packs), len(sigs)), np.nan)
    for _, r in sub.iterrows():
        i = np.where(packs == r.packing_cells_per_mm2)[0][0]
        j = np.where(sigs == r.sigma_um)[0][0]
        M[i, j] = r.oracle_acc
    im = ax.imshow(M, origin="lower", cmap="viridis", aspect="auto", vmin=0.23, vmax=0.96)
    ax.set_xticks(range(len(sigs))); ax.set_xticklabels([f"{s:g}" for s in sigs])
    ax.set_yticks(range(len(packs))); ax.set_yticklabels([f"{int(p):,}" for p in packs])
    ax.set_xlabel("displacement sigma (um)"); ax.set_ylabel("packing (cells/mm$^2$)")
    for i in range(len(packs)):
        for j in range(len(sigs)):
            if np.isfinite(M[i, j]):
                ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=5.5,
                        color="white" if M[i, j] < 0.6 else "black")
    ax.set_title(f"c  Oracle accuracy (density {dmid:g} tx/cell)\nfull-grid range 0.234 to 0.957",
                 loc="left", fontsize=8.5)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("oracle accuracy", fontsize=6.5); cb.ax.tick_params(labelsize=6)


def panel_d(ax, df):
    sc = ax.scatter(df.naive_acc, df.oracle_acc, c=df.sigma_um, cmap="cividis",
                    s=14, edgecolor="none")
    lim = [0.0, 1.0]
    ax.plot(lim, lim, "--", color="0.4", lw=1, label="oracle = naive")
    ax.set_xlim(0.0, 1.0); ax.set_ylim(0.0, 1.0)
    ax.set_aspect("equal")
    ax.set_xlabel("naive accuracy"); ax.set_ylabel("oracle accuracy")
    mn = (df.oracle_minus_naive).min()
    ax.set_title(f"d  Oracle >= naive at all 125 points\n(min margin +{mn:.3f})", loc="left", fontsize=8.5)
    ax.legend(loc="lower right", fontsize=6.5)
    cb = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("sigma (um)", fontsize=6.5); cb.ax.tick_params(labelsize=6)


def main():
    apply_style()
    df = pd.read_csv(os.path.join(RESULTS, "sweep.csv"))
    fig = plt.figure(figsize=(7.2, 6.6))
    gs = GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.34,
                  left=0.08, right=0.97, top=0.93, bottom=0.08)
    panel_a(fig.add_subplot(gs[0, 0]))
    panel_b(fig.add_subplot(gs[0, 1]))
    panel_c(fig.add_subplot(gs[1, 0]), df)
    panel_d(fig.add_subplot(gs[1, 1]), df)
    save(fig, "fig1_frontier")


if __name__ == "__main__":
    main()
