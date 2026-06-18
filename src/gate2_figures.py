"""Gate 2 figures from the committed result CSVs (recomputing the dense oracle curve).

  figures/gate2_dense_vs_sigma.png : dense-regime oracle accuracy vs displacement, with
      the Xenium-pinned sigma, the statistical 95% CI, the combined (CI + design) band,
      the MERFISH-pinned sigma, and the 0.95 decision line marked.
  figures/gate2_design_sigma.png   : pinned sigma across the design-sensitivity grid,
      coloured by marker-selection cutoff (shows the threshold-0.7 family is stable).
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gate2_pin import (build_xenium, DENSE_PACKING, ORACLE_RES_CELL, ORACLE_GRID_MAX,  # noqa: E402
                       MERFISH_PINNED_UM)
from gate1_leakage_anchor import oracle_acc_at  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
FIGURES = os.path.join(ROOT, "figures")


def fig_dense_vs_sigma():
    unc = pd.read_csv(os.path.join(RESULTS, "gate2_sigma_uncertainty.csv")).iloc[0]
    model, real, _ = build_xenium()
    dens = real["density_median_tx_per_cell"]
    sigmas = np.array([0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.3, 1.61, 2.0, 2.53, 3.0, 4.0])
    oa = [oracle_acc_at(model, DENSE_PACKING, dens, s, 12345,
                        res_cell=ORACLE_RES_CELL, grid_max=ORACLE_GRID_MAX)[0] for s in sigmas]

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    ax.plot(sigmas, oa, "-o", color="C3", ms=4, label="dense oracle accuracy (Xenium packing)")
    ax.axhline(0.95, color="k", ls=":", lw=1.2, label="0.95 decision line")
    ax.axvspan(unc.sigma_stat_ci_lo_um, unc.sigma_stat_ci_hi_um, color="C0", alpha=0.18,
               label=f"statistical 95% CI [{unc.sigma_stat_ci_lo_um:.2f}, {unc.sigma_stat_ci_hi_um:.2f}] um")
    ax.axvspan(unc.sigma_combined_lo_um, unc.sigma_combined_hi_um, color="0.6", alpha=0.12,
               label=f"combined band [{unc.sigma_combined_lo_um:.2f}, {unc.sigma_combined_hi_um:.2f}] um")
    ax.axvline(unc.sigma_point_um, color="C3", ls="--", lw=1.2,
               label=f"Xenium-pinned sigma = {unc.sigma_point_um:.2f} um")
    ax.axvline(MERFISH_PINNED_UM, color="C2", ls="-.", lw=1, label=f"MERFISH-pinned = {MERFISH_PINNED_UM} um")
    ax.set_xlabel("displacement sigma (um)")
    ax.set_ylabel("oracle assignment accuracy (dense regime)")
    ax.set_title("Gate 2: dense-regime recovery limit vs the Xenium-direct sigma uncertainty")
    ax.set_ylim(0.55, 1.0)
    ax.legend(fontsize=7.5, loc="lower left")
    fig.tight_layout()
    out = os.path.join(FIGURES, "gate2_dense_vs_sigma.png")
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"wrote {out}")


def fig_design_sigma():
    d = pd.read_csv(os.path.join(RESULTS, "gate2_design_grid.csv"))
    d = d[np.isfinite(d.sigma_um)]
    fig, ax = plt.subplots(figsize=(7, 4.8))
    colors = {0.6: "C0", 0.7: "C1", 0.8: "C3"}
    for thr, sub in d.groupby("excl_threshold"):
        jitter = (np.arange(len(sub)) - len(sub) / 2) * 0.0
        ax.scatter(sub.sigma_um, [thr] * len(sub) + np.random.default_rng(0).uniform(-0.015, 0.015, len(sub)),
                   c=colors.get(thr, "k"), s=40, label=f"cutoff {thr} ({sub.n_markers.iloc[0]} markers)")
    ax.set_yticks(sorted(d.excl_threshold.unique()))
    ax.set_xlabel("pinned sigma (um)")
    ax.set_ylabel("marker exclusivity cutoff")
    ax.set_title("Gate 2 design-sensitivity: pinned sigma by marker cutoff and adjacency bins")
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    out = os.path.join(FIGURES, "gate2_design_sigma.png")
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"wrote {out}")


def main():
    fig_dense_vs_sigma()
    fig_design_sigma()


if __name__ == "__main__":
    main()
