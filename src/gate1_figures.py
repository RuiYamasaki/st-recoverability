"""Gate 1 figures from the committed result CSVs.

  figures/gate1_realistic_frontier.png : oracle accuracy vs displacement under
      realistic expression, with the data-pinned sigma and the two real-anchored
      operating points (MERFISH-like sparse, Xenium-like dense) marked.
  figures/gate1_structural_sensitivity.png : the frontier under non-Voronoi geometry
      and non-Gaussian displacement vs the baseline, at the dense packing.
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
REP_DENSITY = 160.0


def fig_frontier():
    sweep = pd.read_csv(os.path.join(RESULTS, "gate1_sweep.csv"))
    leak = pd.read_csv(os.path.join(RESULTS, "gate1_leakage.csv"))
    sub = sweep[sweep.mean_tx_per_cell == REP_DENSITY]
    fig, ax = plt.subplots(figsize=(8, 5.4))
    for p in sorted(sub.packing_cells_per_mm2.unique()):
        d = sub[sub.packing_cells_per_mm2 == p].sort_values("sigma_um")
        ax.plot(d.sigma_um, d.oracle_acc, "-o", ms=4, label=f"{int(p)} cells/mm$^2$")
    sig_pin = float(leak.sigma_pinned_um.iloc[0])
    lo = float(leak.sigma_bracket_lo_um.iloc[0]); hi = float(leak.sigma_bracket_hi_um.iloc[0])
    ax.axvspan(lo, hi, color="0.85", zorder=0, label=f"sigma 25-75 bracket [{lo:.1f},{hi:.1f}] um")
    ax.axvline(sig_pin, color="k", ls="--", lw=1.2, label=f"data-pinned sigma = {sig_pin:.2f} um")
    for _, r in leak.iterrows():
        mk = "D" if r.anchor_dataset == "xenium_breast" else "s"
        lab = ("Xenium-like (dense): %.2f" % r.oracle_acc_at_pinned) if r.anchor_dataset == "xenium_breast" \
            else ("MERFISH-like (sparse): %.2f" % r.oracle_acc_at_pinned)
        ax.plot([sig_pin], [r.oracle_acc_at_pinned], mk, color="red", ms=11,
                markeredgecolor="k", zorder=5, label=lab)
    ax.axhline(0.9, color="0.4", ls=":", lw=1)
    ax.set_xlabel("displacement sigma (um)")
    ax.set_ylabel("oracle assignment accuracy")
    ax.set_title("Gate 1: recovery frontier under realistic expression, with data-pinned displacement")
    ax.set_ylim(0.15, 1.0)
    ax.legend(fontsize=7.5, loc="lower left", ncol=2)
    fig.tight_layout()
    out = os.path.join(FIGURES, "gate1_realistic_frontier.png")
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"wrote {out}")


def fig_structural():
    st = pd.read_csv(os.path.join(RESULTS, "gate1_structural.csv"))
    pk = st.packing_cells_per_mm2.max()
    fig, ax = plt.subplots(figsize=(7, 5))
    styles = {"baseline": ("-o", "C0"), "aniso": ("--s", "C1"), "mixture": ("-.^", "C2")}
    for cond, (ls, col) in styles.items():
        d = st[(st.condition == cond) & (st.packing_cells_per_mm2 == pk)].sort_values("sigma_um")
        ax.plot(d.sigma_um, d.oracle_acc, ls, color=col, label=f"oracle: {cond}")
        ax.plot(d.sigma_um, d.naive_acc, ls, color=col, alpha=0.35, label=f"naive: {cond}")
    ax.set_xlabel("displacement sigma (um)")
    ax.set_ylabel("assignment accuracy")
    ax.set_title(f"Gate 1 structural sensitivity (dense packing {int(pk)} cells/mm$^2$)")
    ax.set_ylim(0, 1.0)
    ax.legend(fontsize=7.5, ncol=2)
    fig.tight_layout()
    out = os.path.join(FIGURES, "gate1_structural_sensitivity.png")
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"wrote {out}")


def main():
    fig_frontier()
    fig_structural()


if __name__ == "__main__":
    main()
