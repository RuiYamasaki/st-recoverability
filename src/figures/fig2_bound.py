"""Figure 2: the bound on real tissue.

(a) marker-validity permutation null: real adjacent-minus-distant displacement signal vs
    permutation p-value; admissible markers (p<0.05) separate from biology-diluted ones.
(b) data-pinned displacement sigma with bootstrap CI: breast 1.99 [1.43, 2.63] (28/38),
    lung 1.83 [1.62, 2.08] (31/32).
(c) dense-tissue oracle accuracy at ~13,600 cells/mm^2 across the sigma CI: breast and
    lung, both below 0.9.
(d) robustness: cell-typing stability (sigma 1.99-2.65, dense optimistic <=0.803), with
    realistic-overlapping-expression and negative-binomial results annotated.

Values from results/gate3_validity_breast.csv, gate3_pin_breast.csv, gate3_pin_lung.csv,
gate3_typing.csv, gate2_nbinom.csv. No recomputation.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import apply_style, save, OI  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS = os.path.join(ROOT, "results")
BREAST = OI["vermillion"]
LUNG = OI["blue"]


def panel_a(ax):
    v = pd.read_csv(os.path.join(RESULTS, "gate3_validity_breast.csv")).dropna(subset=["pval", "real_signal"])
    nlp = -np.log10(v.pval.clip(lower=1e-4))
    adm = v.admissible.values
    ax.scatter(v.real_signal[~adm], nlp[~adm], s=16, color=OI["grey"], edgecolor="none",
               label="rejected (biology-diluted)")
    ax.scatter(v.real_signal[adm], nlp[adm], s=16, color=OI["green"], edgecolor="none",
               label="admissible (p<0.05)")
    ax.axhline(-np.log10(0.05), ls="--", color="0.4", lw=1)
    ax.text(ax.get_xlim()[0], -np.log10(0.05) + 0.05, " p = 0.05", fontsize=6, color="0.3", va="bottom")
    ax.axvline(0, ls=":", color="0.6", lw=0.8)
    # label a couple of failing loose markers
    for nm in ["AHSP", "CRHBP", "CYP1A1"]:
        row = v[v.gene_name == nm]
        if len(row):
            r = row.iloc[0]
            ax.annotate(nm, (r.real_signal, -np.log10(max(r.pval, 1e-4))), fontsize=5.5,
                        color="0.3", xytext=(2, 2), textcoords="offset points")
    ax.set_xlabel("real displacement signal (adjacent - distant)")
    ax.set_ylabel(r"$-\log_{10}$ permutation p")
    n_adm = int(adm.sum())
    ax.set_title(f"a  Marker-validity null (1,000 perms)\n{n_adm}/{len(v)} admissible", loc="left", fontsize=8.5)
    ax.legend(loc="upper left", fontsize=6)


def panel_b(ax):
    b = pd.read_csv(os.path.join(RESULTS, "gate3_pin_breast.csv")).iloc[0]
    l = pd.read_csv(os.path.join(RESULTS, "gate3_pin_lung.csv")).iloc[0]
    rows = [("breast", b.sigma_point_um, b.sigma_ci_lo_um, b.sigma_ci_hi_um,
             int(b.n_admissible), int(b.n_candidate), BREAST),
            ("lung", l.sigma_point_um, l.sigma_ci_lo_um, l.sigma_ci_hi_um,
             int(l.n_admissible), int(l.n_candidate), LUNG)]
    for i, (nm, pt, lo, hi, na, nc, col) in enumerate(rows):
        ax.errorbar(pt, i, xerr=[[pt - lo], [hi - pt]], fmt="o", color=col, capsize=4, ms=7, lw=1.6)
        ax.text(hi + 0.08, i, f"{pt:.2f} um  [{lo:.2f}, {hi:.2f}]\n{na}/{nc} markers",
                va="center", fontsize=6.5)
    ax.set_yticks([0, 1]); ax.set_yticklabels(["breast", "lung"])
    ax.set_ylim(-0.5, 1.5); ax.set_xlim(0, 3.6)
    ax.set_xlabel("data-pinned displacement sigma (um)")
    ax.set_title("b  Pinned displacement, bootstrap CI", loc="left", fontsize=8.5)


def panel_c(ax):
    b = pd.read_csv(os.path.join(RESULTS, "gate3_pin_breast.csv")).iloc[0]
    l = pd.read_csv(os.path.join(RESULTS, "gate3_pin_lung.csv")).iloc[0]
    # dense oracle at ~13,600: optimistic = at sigma CI low, pessimistic = at sigma CI high
    data = [("breast", b.dense_oracle_at_sigma_ci_hi, b.dense_oracle_point, b.dense_oracle_at_sigma_ci_lo, BREAST),
            ("lung", l.dense_oracle_at_sigma_ci_hi, l.dense_oracle_point, l.dense_oracle_at_sigma_ci_lo, LUNG)]
    xs = [0, 1]
    for x, (nm, lo, pt, hi, col) in zip(xs, data):
        ax.errorbar(x, pt, yerr=[[pt - lo], [hi - pt]], fmt="o", color=col, capsize=5, ms=8, lw=1.8)
        ax.text(x + 0.08, pt, f"{lo:.3f}-{hi:.3f}", va="center", fontsize=6.5)
    ax.axhline(0.9, ls="--", color=OI["vermillion"], lw=1.2)
    ax.text(1.45, 0.905, "0.9", color=OI["vermillion"], fontsize=7, va="bottom", ha="right")
    ax.set_xticks(xs); ax.set_xticklabels(["breast", "lung"])
    ax.set_xlim(-0.4, 1.8); ax.set_ylim(0.5, 0.95)
    ax.set_ylabel("dense oracle accuracy")
    ax.set_title("c  Dense tissue (~13,600 cells/mm$^2$)\nbelow 0.9 across the sigma CI", loc="left", fontsize=8.5)


def panel_d(ax):
    t = pd.read_csv(os.path.join(RESULTS, "gate3_typing.csv"))
    nb = pd.read_csv(os.path.join(RESULTS, "gate2_nbinom.csv")).iloc[0]
    labels = {"kmeans_K10": "K=10", "kmeans_K14": "K=14", "kmeans_K14_altseed": "K=14'", "kmeans_K20": "K=20"}
    t = t.set_index("variant").loc[["kmeans_K10", "kmeans_K14", "kmeans_K14_altseed", "kmeans_K20"]].reset_index()
    xs = np.arange(len(t))
    # dense oracle optimistic end (sigma CI low) per typing variant
    ax.bar(xs, t.dense_oracle_at_sigma_ci_lo, color=OI["skyblue"], width=0.62,
           label="dense oracle, optimistic end")
    for x, val, s in zip(xs, t.dense_oracle_at_sigma_ci_lo, t.sigma_point_um):
        ax.text(x, val + 0.012, f"{val:.3f}\n({s:.2f} um)", ha="center", fontsize=5.8)
    ax.axhline(0.9, ls="--", color=OI["vermillion"], lw=1.2)
    ax.text(len(t) - 0.5, 0.905, "0.9", color=OI["vermillion"], fontsize=7, va="bottom", ha="right")
    ax.set_xticks(xs); ax.set_xticklabels([labels[v] for v in t.variant])
    ax.set_ylim(0.0, 1.0); ax.set_ylabel("dense oracle (optimistic)")
    ax.set_title("d  Robustness: cell-typing (sigma 1.99-2.65,\ndense optimistic <= 0.803)", loc="left", fontsize=8.5)
    note = (f"realistic overlapping expression: frontier mean |diff| 0.006\n"
            f"negative-binomial variance (dispersion {nb.nb_dispersion_median:.2f}): "
            f"dense shift {nb.dense_oracle_shift:+.3f}")
    ax.text(0.5, -0.34, note, transform=ax.transAxes, ha="center", va="top", fontsize=6,
            bbox=dict(boxstyle="round,pad=0.3", fc="0.95", ec="0.7", lw=0.6))


def main():
    apply_style()
    fig = plt.figure(figsize=(7.2, 6.8))
    gs = GridSpec(2, 2, figure=fig, hspace=0.5, wspace=0.34,
                  left=0.09, right=0.97, top=0.93, bottom=0.13)
    panel_a(fig.add_subplot(gs[0, 0]))
    panel_b(fig.add_subplot(gs[0, 1]))
    panel_c(fig.add_subplot(gs[1, 0]))
    panel_d(fig.add_subplot(gs[1, 1]))
    save(fig, "fig2_bound")


if __name__ == "__main__":
    main()
