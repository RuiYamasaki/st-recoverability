"""Generate Gate-0 figures from the committed result CSVs.

  figures/oracle_accuracy_surface.png : the oracle assignment-accuracy frontier
      over (packing, displacement), with the real-data anchors marked.
  figures/identity_vs_typing.png      : the decoupling, exact cell-identity
      assignment (the frontier) vs cell-type profile recovery (a high oracle
      ceiling), with the naive baseline for contrast.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
FIGURES = os.path.join(ROOT, "figures")
os.makedirs(FIGURES, exist_ok=True)

REP_DENSITY = 160.0   # representative density for the accuracy surface (acc is ~density-independent)


def fig_surface():
    df = pd.read_csv(os.path.join(RESULTS, "sweep.csv"))
    sub = df[df.mean_tx_per_cell == REP_DENSITY]
    packings = sorted(sub.packing_cells_per_mm2.unique())
    sigmas = sorted(sub.sigma_um.unique())
    Z = np.zeros((len(sigmas), len(packings)))
    for i, s in enumerate(sigmas):
        for j, p in enumerate(packings):
            Z[i, j] = sub[(sub.sigma_um == s) & (sub.packing_cells_per_mm2 == p)].oracle_acc.iloc[0]

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(13, 5.2))

    im = ax0.imshow(Z, origin="lower", aspect="auto", cmap="viridis", vmin=0.2, vmax=1.0)
    ax0.set_xticks(range(len(packings)))
    ax0.set_xticklabels([f"{int(p)}" for p in packings])
    ax0.set_yticks(range(len(sigmas)))
    ax0.set_yticklabels([f"{s:g}" for s in sigmas])
    ax0.set_xlabel("cell packing (cells/mm$^2$)")
    ax0.set_ylabel("displacement sigma (um)")
    ax0.set_title(f"Oracle assignment accuracy (density {int(REP_DENSITY)} tx/cell)")
    for i in range(len(sigmas)):
        for j in range(len(packings)):
            ax0.text(j, i, f"{Z[i, j]:.2f}", ha="center", va="center",
                     color="white" if Z[i, j] < 0.6 else "black", fontsize=8)
    fig.colorbar(im, ax=ax0, label="oracle accuracy")
    # mark the real-data anchor packings (nearest grid columns)
    real = pd.read_csv(os.path.join(RESULTS, "realism.csv"))
    for ds, mark in [("xenium_breast", "Xenium\n~13.5k"), ("merfish_hypothal", "MERFISH\n~2.6k")]:
        pk = real[real.dataset == ds].packing_cells_per_mm2.mean()
        j = int(np.argmin([abs(pk - p) for p in packings]))
        ax0.text(j, len(sigmas) - 0.62, mark, ha="center", va="top",
                 fontsize=7.5, color="red", fontweight="bold",
                 bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="red", alpha=0.7))

    for p in packings:
        d = sub[sub.packing_cells_per_mm2 == p].sort_values("sigma_um")
        ax1.plot(d.sigma_um, d.oracle_acc, "-o", label=f"{int(p)} cells/mm$^2$", ms=4)
    ax1.axhspan(0.0, 0.0)  # no-op to keep axis stable
    ax1.axvspan(0.5, 2.0, color="0.85", zorder=0, label="realistic sigma band")
    ax1.set_xlabel("displacement sigma (um)")
    ax1.set_ylabel("oracle assignment accuracy")
    ax1.set_title("Frontier vs displacement, by packing")
    ax1.set_ylim(0.1, 1.0)
    ax1.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    out = os.path.join(FIGURES, "oracle_accuracy_surface.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out}")


def fig_decoupling():
    df = pd.read_csv(os.path.join(RESULTS, "sweep.csv"))
    pk = 9500.0
    sub = df[(df.packing_cells_per_mm2 == pk) & (df.mean_tx_per_cell == REP_DENSITY)].sort_values("sigma_um")
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(sub.sigma_um, sub.oracle_acc, "-o", color="C0", label="oracle: exact cell-identity accuracy")
    ax.plot(sub.sigma_um, sub.naive_acc, "--o", color="C0", alpha=0.6, label="naive: exact cell-identity accuracy")
    ax.plot(sub.sigma_um, sub.profile_ari_vs_truetype, "-s", color="C2", label="oracle: cell-type profile ARI")
    ax.plot(sub.sigma_um, sub.naive_profile_ari_vs_truetype, "--s", color="C3", label="naive: cell-type profile ARI")
    ax.set_xlabel("displacement sigma (um)")
    ax.set_ylabel("metric value")
    ax.set_ylim(0, 1.02)
    ax.set_title(f"Identity assignment vs cell-type recovery (packing {int(pk)} cells/mm$^2$)")
    ax.legend(fontsize=8, loc="center left")
    fig.tight_layout()
    out = os.path.join(FIGURES, "identity_vs_typing.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out}")


def main():
    fig_surface()
    fig_decoupling()


if __name__ == "__main__":
    main()
