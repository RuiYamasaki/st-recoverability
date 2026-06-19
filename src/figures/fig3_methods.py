"""Figure 3: real segmentation methods against the ceiling (dense tissue).

(a) one-to-one assignment accuracy vs displacement in dense tissue: oracle, naive, and
    Baysor / Proseg / ComSeg at their best documented configuration.
(b) best-method-to-oracle gap across the three metrics and three dense displacements
    (0.066 to 0.110 below the oracle).
(c) best-method-minus-naive across the three metrics, with the run-to-run noise band
    (~0.012): one-to-one 0.004-0.009 below, many-to-one up to +0.009, co-assignment
    (ARI) up to +0.019.
(d) free-segmentation (no nuclei prior, Baysor only) vs naive and oracle in dense tissue.

Values from results/headroom_fair_metrics.csv and headroom_fair_free.csv. No recomputation.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import apply_style, save, OI, ORACLE, NAIVE, METHOD  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS = os.path.join(ROOT, "results")
DENSE = ["dense_sigma1.43", "dense_sigma1.99", "dense_sigma2.63"]
SIGS = [1.43, 1.99, 2.63]
METRICS = [("acc_one_to_one", "one-to-one"), ("acc_many_to_one", "many-to-one"), ("ari", "co-assignment (ARI)")]
NOISE = 0.012


def get(df, method, config, col):
    r = df[(df.method == method) & (df.config == config)]
    return float(r[col].iloc[0]) if len(r) else np.nan


def panel_a(ax, df):
    for nm, col in [("oracle", ORACLE), ("naive", NAIVE)]:
        ys = [get(df, nm, c, "acc_one_to_one") for c in DENSE]
        ax.plot(SIGS, ys, "o-" if nm == "oracle" else "s--", color=col, lw=1.8, ms=6,
                label=nm, zorder=3 if nm == "oracle" else 2)
    for m in ["Baysor", "ComSeg", "Proseg"]:
        ys = [get(df, m, c, "acc_one_to_one") for c in DENSE]
        ax.plot(SIGS, ys, "^-", color=METHOD[m], lw=1.2, ms=5, label=m, alpha=0.9)
    ax.set_xlabel("displacement sigma (um)"); ax.set_ylabel("one-to-one accuracy")
    ax.set_xticks(SIGS); ax.set_title("a  Dense tissue: methods vs ceiling", loc="left", fontsize=8.5)
    ax.legend(fontsize=6, ncol=2, loc="upper right")


def panel_b(ax, df):
    # best-method-to-oracle gap per metric per dense sigma
    width = 0.25
    xs = np.arange(len(METRICS))
    gaps_all = []
    for k, c in enumerate(DENSE):
        gaps = []
        for col, _ in METRICS:
            orc = get(df, "oracle", c, col)
            best = np.nanmax([get(df, m, c, col) for m in ["Baysor", "Proseg", "ComSeg"]])
            gaps.append(orc - best)
        gaps_all += gaps
        ax.bar(xs + (k - 1) * width, gaps, width, color=plt.cm.Greys(0.4 + 0.2 * k),
               label=f"sigma {SIGS[k]:g}")
    ax.set_xticks(xs); ax.set_xticklabels([m[1] for m in METRICS], fontsize=6.5)
    ax.set_ylabel("oracle - best method")
    lo, hi = min(gaps_all), max(gaps_all)
    ax.set_title(f"b  Best-method-to-oracle gap\n({lo:.3f} to {hi:.3f})", loc="left", fontsize=8.5)
    ax.legend(fontsize=6, loc="upper left")


def panel_c(ax, df):
    # best-method-minus-naive per metric, across dense sigmas (range shown)
    xs = np.arange(len(METRICS))
    ax.axhspan(-NOISE, NOISE, color="0.85", alpha=0.7, zorder=0, label=f"noise floor +/-{NOISE}")
    ax.axhline(0, color="0.4", lw=0.8)
    for col, _ in METRICS:
        pass
    for j, (col, _) in enumerate(METRICS):
        deltas = []
        for c in DENSE:
            naive = get(df, "naive", c, col)
            best = np.nanmax([get(df, m, c, col) for m in ["Baysor", "Proseg", "ComSeg"]])
            deltas.append(best - naive)
        ax.plot([j, j], [min(deltas), max(deltas)], color=OI["purple"], lw=4, solid_capstyle="round", alpha=0.5)
        ax.scatter([j] * len(deltas), deltas, color=OI["purple"], s=22, zorder=4)
        ax.text(j + 0.12, max(deltas), f"{min(deltas):+.3f}\nto {max(deltas):+.3f}", fontsize=6, va="center")
    ax.set_xticks(xs); ax.set_xticklabels([m[1] for m in METRICS], fontsize=6.5)
    ax.set_xlim(-0.5, 2.8)
    ax.set_ylabel("best method - naive")
    ax.set_title("c  Methods vs nearest-nucleus baseline", loc="left", fontsize=8.5)
    ax.legend(fontsize=6, loc="lower left")


def panel_d(ax, df, free):
    # free segmentation (Baysor) vs naive and oracle, dense
    rows = [("oracle\n(prior)", get(df, "oracle", "dense_sigma1.99", "acc_one_to_one"), ORACLE),
            ("naive", get(df, "naive", "dense_sigma1.99", "acc_one_to_one"), NAIVE),
            ("Baysor\nprior", get(df, "Baysor", "dense_sigma1.99", "acc_one_to_one"), METHOD["Baysor"]),
            ("Baysor\nfree", get(free, "Baysor", "dense_sigma1.99", "acc_one_to_one"), OI["yellow"])]
    xs = np.arange(len(rows))
    ax.bar(xs, [r[1] for r in rows], color=[r[2] for r in rows], width=0.66, edgecolor="0.3", lw=0.5)
    for x, r in zip(xs, rows):
        ax.text(x, r[1] + 0.012, f"{r[1]:.3f}", ha="center", fontsize=6.5)
    ax.axhline(get(df, "naive", "dense_sigma1.99", "acc_one_to_one"), ls=":", color=NAIVE, lw=1)
    ax.set_xticks(xs); ax.set_xticklabels([r[0] for r in rows], fontsize=6.5)
    ax.set_ylim(0, 0.85); ax.set_ylabel("one-to-one accuracy (dense, sigma 1.99)")
    ax.set_title("d  Free segmentation falls below naive\n(Baysor; Proseg/ComSeg free unsupported)",
                 loc="left", fontsize=8.5)


def main():
    apply_style()
    df = pd.read_csv(os.path.join(RESULTS, "headroom_fair_metrics.csv"))
    free = pd.read_csv(os.path.join(RESULTS, "headroom_fair_free.csv"))
    fig = plt.figure(figsize=(7.2, 6.6))
    gs = GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.32,
                  left=0.09, right=0.97, top=0.93, bottom=0.09)
    panel_a(fig.add_subplot(gs[0, 0]), df)
    panel_b(fig.add_subplot(gs[0, 1]), df)
    panel_c(fig.add_subplot(gs[1, 0]), df)
    panel_d(fig.add_subplot(gs[1, 1]), df, free)
    save(fig, "fig3_methods")


if __name__ == "__main__":
    main()
