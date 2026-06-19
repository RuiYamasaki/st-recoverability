"""Method-fairness figure: best-config methods vs oracle/naive across the dense sigma CI,
under each of the three metrics. Shows whether the methods-vs-naive picture depends on the
metric (one-to-one penalises over-segmentation; many-to-one and ARI treat it differently).

Regen: micromamba run -n st python src/headroom_fair_figures.py
Source: results/headroom_fair_metrics.csv. Writes figures/headroom_fair_metrics_dense.png.
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
DENSE = [("dense_sigma1.43", 1.43), ("dense_sigma1.99", 1.99), ("dense_sigma2.63", 2.63)]
METRICS = [("acc_one_to_one", "one-to-one matched"),
           ("acc_many_to_one", "many-to-one homogeneity"),
           ("ari", "co-assignment ARI")]
SERIES = [("oracle", "C2", "-o"), ("naive", "0.5", "-o"),
          ("Baysor", "C0", "-s"), ("ComSeg", "C4", "-D"), ("Proseg", "C1", "-^")]


def main():
    mt = pd.read_csv(os.path.join(RESULTS, "headroom_fair_metrics.csv"))

    def val(method, cfg, metric):
        r = mt[(mt.config == cfg) & (mt.method == method)]
        if not len(r) or r.iloc[0].get("status", "ok") != "ok":
            return np.nan
        return float(r.iloc[0][metric])

    sig = [s for _, s in DENSE]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharex=True)
    for ax, (mcol, mname) in zip(axes, METRICS):
        for method, color, style in SERIES:
            ys = [val(method, c, mcol) for c, _ in DENSE]
            ax.plot(sig, ys, style, color=color, label=method)
        ax.set_title(mname)
        ax.set_xlabel("displacement sigma (um), dense")
        ax.set_xticks(sig)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("accuracy / ARI vs known truth")
    axes[0].legend(fontsize=8, loc="upper right")
    fig.suptitle("Method fairness: best documented config vs oracle and naive, by metric (dense regime)")
    fig.tight_layout()
    out = os.path.join(FIGURES, "headroom_fair_metrics_dense.png")
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
