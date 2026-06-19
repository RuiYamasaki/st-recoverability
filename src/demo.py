"""Quick demo (no network, no real data, < 1 minute).

Builds two synthetic fields at a realistic data-pinned displacement (sigma = 2 um),
one sparse and one dense, and compares the Bayes-optimal (oracle) transcript-to-cell
assignment accuracy with the naive nearest-nucleus baseline. It reproduces, in miniature,
the headline of the study: in dense tissue the best-possible assignment accuracy drops
well below 1, i.e. transcript-to-cell assignment becomes unrecoverable, regardless of
method.

Run:  python src/demo.py
Output:  printed accuracies + figures/demo_oracle_vs_naive.png
"""
from __future__ import annotations

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402
from generator import build_field, generate_transcripts  # noqa: E402
from oracle import build_oracle_maps, oracle_assign, naive_assign  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGURES = os.path.join(ROOT, "figures")
SIGMA_UM = 2.0          # close to the data-pinned displacement (breast 1.99 um)
DENSITY = 160           # transcripts per cell
SEED = config.MASTER_SEED + 990000
REGIMES = [("sparse\n(2,575 cells/mm$^2$)", 2575.0), ("dense\n(13,625 cells/mm$^2$)", 13625.0)]


def accuracy(true_cell, assigned, interior):
    return float((assigned[interior] == true_cell[interior]).mean())


def main():
    os.makedirs(FIGURES, exist_ok=True)
    labels, oracle_accs, naive_accs = [], [], []
    print(f"Demo: oracle vs naive transcript-to-cell assignment at sigma = {SIGMA_UM} um\n")
    for i, (name, packing) in enumerate(REGIMES):
        f = build_field(packing, SIGMA_UM, SEED + 100 * i)
        tx = generate_transcripts(f, DENSITY, SEED + 100 * i + 1)
        Dmax, argcell = build_oracle_maps(f)
        oa = oracle_assign(f, Dmax, argcell, tx.obs_xy, tx.gene)
        na = naive_assign(f, tx.obs_xy)
        o = accuracy(tx.true_cell, oa, tx.interior)
        n = accuracy(tx.true_cell, na, tx.interior)
        labels.append(name.replace("\n", " ")); oracle_accs.append(o); naive_accs.append(n)
        print(f"  {name.splitlines()[0]:<8} packing {packing:>7.0f}: "
              f"oracle = {o:.3f}, naive = {n:.3f}, oracle-minus-naive = {o - n:+.3f}")

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    x = np.arange(len(REGIMES)); w = 0.36
    ax.bar(x - w / 2, oracle_accs, w, label="oracle (best possible)", color="#D55E00")
    ax.bar(x + w / 2, naive_accs, w, label="naive (nearest nucleus)", color="#0072B2")
    ax.axhline(0.9, ls="--", color="0.4", lw=1)
    ax.text(len(REGIMES) - 0.5, 0.905, "0.9", color="0.4", fontsize=8, va="bottom", ha="right")
    for xi, (o, n) in enumerate(zip(oracle_accs, naive_accs)):
        ax.text(xi - w / 2, o + 0.01, f"{o:.2f}", ha="center", fontsize=8)
        ax.text(xi + w / 2, n + 0.01, f"{n:.2f}", ha="center", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels([r[0] for r in REGIMES])
    ax.set_ylim(0, 1.0); ax.set_ylabel("assignment accuracy")
    ax.set_title(f"Recoverability ceiling at sigma = {SIGMA_UM} um")
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    out = os.path.join(FIGURES, "demo_oracle_vs_naive.png")
    fig.savefig(out, dpi=150)
    print(f"\nwrote {out}")
    print("Takeaway: in dense tissue even the oracle (the information-theoretic best) falls "
          "well below 1, so transcript-to-cell assignment there is unrecoverable by any method.")


if __name__ == "__main__":
    main()
