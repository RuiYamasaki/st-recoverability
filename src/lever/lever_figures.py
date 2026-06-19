"""Lever analysis figures: the worked Janesick triple-positive lever.

Panel A: oracle answerability ceiling at the triple-positive ROI's own local density,
across the sigma CI, versus the 0.90 and 0.95 bars (Leg 1).
Panel B: synthetic spurious triple-positive rate (known truth, zero true triples) versus
displacement sigma, field-average and conditional on DCIS-adjacent-to-PR-source, with the
sigma=0 control, the data-pinned sigma and CI, and the observed real rate (Leg 2).

Reads results/lever/lever_roi.csv and lever_mechanism.csv. Writes figures/lever_*.png.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS = os.path.join(ROOT, "results", "lever")
FIGURES = os.path.join(ROOT, "figures")
THR = 2


def main():
    os.makedirs(FIGURES, exist_ok=True)
    roi = pd.read_csv(os.path.join(RESULTS, "lever_roi.csv"))
    mech = pd.read_csv(os.path.join(RESULTS, "lever_mechanism.csv")).sort_values("sigma_um")
    params = pd.read_csv(os.path.join(RESULTS, "lever_mechanism_params.csv")).iloc[0]
    obs = float(params[f"obs_triple_rate_among_dcis_t{THR}"]) * 100

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.3))

    # Panel A: ROI ceiling
    labels = {"tp_neighborhood": "ROI tissue (50um of a\ntriple-positive cell)",
              "tp_cells_own": "triple-positive cells'\nown local packing",
              "dcis_nest_p90_proxy": "densest DCIS nest\n(p90 proxy)"}
    xs = np.arange(len(roi))
    pt = roi.oracle_point.values
    lo = roi.oracle_opt_ci_lo.values   # optimistic (low sigma)
    hi = roi.oracle_pess_ci_hi.values  # pessimistic (high sigma)
    axA.errorbar(xs, pt, yerr=[pt - hi, lo - pt], fmt="o", color="#b2182b",
                 capsize=5, ms=8, lw=2, label="oracle ceiling (point, sigma CI)")
    axA.axhline(0.90, ls="--", color="0.4", lw=1.2, label="0.90 bar")
    axA.axhline(0.95, ls=":", color="0.55", lw=1.2, label="0.95 bar")
    axA.set_xticks(xs)
    axA.set_xticklabels([labels.get(r, r) for r in roi.roi_def], fontsize=8)
    axA.set_ylim(0.6, 1.0)
    axA.set_ylabel("oracle (best-possible) assignment accuracy")
    axA.set_title("A. Ceiling at the triple-positive ROI\n(below 0.95 across the whole sigma CI)", fontsize=10)
    axA.legend(fontsize=7, loc="upper right")

    # Panel B: mechanism
    s = mech.sigma_um.values
    fld = mech[f"spurious_triple_rate_among_dcis_t{THR}"].values * 100
    adj = mech[f"spurious_triple_rate_dcis_adj_pr_t{THR}"].values * 100
    axB.plot(s, adj, "o-", color="#b2182b", lw=2, ms=7,
             label="adjacent to a PR source\n(the described micro-environment)")
    axB.plot(s, fld, "s--", color="#2166ac", lw=2, ms=6, label="field average (all DCIS)")
    axB.axhline(obs, ls="-", color="0.2", lw=1.4, label=f"observed real rate ({obs:.2f}%)")
    # pinned sigma CI band
    pin = mech[mech.sigma_name == "pinned"].iloc[0]
    axB.axvspan(1.43, 2.50, color="0.85", alpha=0.5, zorder=0, label="pinned sigma CI [1.43, 2.50]")
    axB.axvline(2.10, ls=":", color="0.3", lw=1)
    axB.set_xlabel("displacement sigma (um)")
    axB.set_ylabel(f"spurious triple-positive rate among DCIS (thr>={THR}), %")
    axB.set_title("B. Synthetic null: zero true triples\n(sigma=0 control near zero)", fontsize=10)
    axB.legend(fontsize=7, loc="upper left")
    axB.set_xlim(-0.1, 2.7)

    fig.tight_layout()
    out = os.path.join(FIGURES, "lever_worked_example.png")
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
