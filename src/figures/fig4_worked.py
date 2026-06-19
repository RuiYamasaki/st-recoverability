"""Figure 4: the worked diagnostic on the Janesick triple-positive.

(a) Leg 1: the oracle ceiling at the triple-positive ROI's own local density (evaluated
    at the locked breast pin 1.99 um, CI [1.43, 2.63]); below 0.90 and 0.95 across the CI.
(b) Leg 2 local-to-local (apples-to-apples): the synthetic spurious triple-positive rate
    from pure displacement (zero true triples; sigma=0 control) versus the observed rate
    in the matched local condition (DCIS cells adjacent to a PR source). The synthetic
    under-produces by about 19x, so displacement is CONTRIBUTORY, not sufficient.
(c) Leg 3 plus the artifact-consistent vs corroborated contrast: the orthogonal evidence
    (scFFPE-seq, clinical PR-negative) and the ceiling positions of the corroborated
    colorectal and hypothalamus cases.

Integrity: no oracle or segmenter is run on the real transcripts. Leg 1 and the contrast
use the method-independent ceiling on synthetic fields at measured packing/density; Leg 2
uses the known-truth synthetic null; Leg 3 uses the published orthogonal data and the
released matrix (read, not re-assigned).

Values from results/lever/lever_roi_pin199.csv, lever_mechanism.csv, leg2_local.csv,
lever_ceiling.csv. No gate re-run.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import apply_style, save, OI, ORACLE  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS = os.path.join(ROOT, "results", "lever")
SECTION_MEDIAN_PACK = 6534
PGR_POS_FRAC = 2.9   # % of cells PGR+ (released matrix), consistent with PR-negative status


def panel_a(ax):
    r = pd.read_csv(os.path.join(RESULTS, "lever_roi_pin199.csv"))
    names = {"tp_neighborhood": "ROI tissue\n(7,103)", "tp_cells_own": "triple+ cells'\nown (7,985)",
             "dcis_nest_p90_proxy": "densest nest\n(9,329)"}
    xs = np.arange(len(r))
    pt = r.oracle_point.values
    lo = r.oracle_opt_ci_lo.values     # optimistic
    hi = r.oracle_pess_ci_hi.values    # pessimistic
    ax.errorbar(xs, pt, yerr=[pt - hi, lo - pt], fmt="o", color=ORACLE, capsize=5, ms=8, lw=1.8)
    for x, p in zip(xs, pt):
        ax.text(x + 0.12, p, f"{p:.3f}", va="center", fontsize=6.5)
    ax.axhline(0.90, ls="--", color="0.35", lw=1.1); ax.axhline(0.95, ls=":", color="0.55", lw=1.1)
    ax.text(len(r) - 0.5, 0.905, "0.90", fontsize=6.5, color="0.35", va="bottom", ha="right")
    ax.text(len(r) - 0.5, 0.955, "0.95", fontsize=6.5, color="0.55", va="bottom", ha="right")
    ax.set_xticks(xs); ax.set_xticklabels([names[v] for v in r.roi_def], fontsize=6.5)
    ax.set_ylim(0.6, 1.0); ax.set_ylabel("oracle assignment accuracy")
    ax.set_title("a  Leg 1: ROI above the ceiling\n(pin 1.99 um; section median packing 6,534)",
                 loc="left", fontsize=8.5)
    ax.text(0.02, 0.04, f"PGR+ in only {PGR_POS_FRAC}% of cells\n(consistent with PR-negative)",
            transform=ax.transAxes, fontsize=6, va="bottom",
            bbox=dict(boxstyle="round,pad=0.25", fc="0.96", ec="0.7", lw=0.5))


def panel_b(ax):
    l2 = pd.read_csv(os.path.join(RESULTS, "leg2_local.csv")).iloc[0]
    obs = l2.obs_rate_adj_to_prsource * 100
    syn_pin = l2.syn_adj_pinned * 100
    syn_hi = l2.syn_adj_ci_hi * 100
    cats = ["synthetic\nsigma=0 (control)", "synthetic\npinned 2.10 um", "synthetic\nupper CI 2.50 um",
            "observed real\n(DCIS adj. PR source)"]
    vals = [0.0, syn_pin, syn_hi, obs]
    cols = [OI["grey"], OI["skyblue"], OI["blue"], ORACLE]
    ys = np.arange(len(cats))[::-1]
    ax.barh(ys, vals, color=cols, height=0.6)
    for y, v in zip(ys, vals):
        ax.text(v + 0.5, y, f"{v:.2f}%", va="center", fontsize=6.8)
    ax.set_yticks(ys); ax.set_yticklabels(cats, fontsize=6.3)
    ax.set_xlim(0, 31); ax.set_xlabel("apparent triple-positive rate among DCIS (thr>=2), %")
    ax.set_title("b  Leg 2 local-to-local: displacement is\ncontributory, not sufficient", loc="left", fontsize=8.5)
    # ~19x gap annotation in the empty middle-right band (between the short synthetic bars and the long observed bar)
    ax.annotate("", xy=(obs, 0.5), xytext=(syn_hi, 0.5),
                arrowprops=dict(arrowstyle="<->", color="0.35", lw=1.0))
    ax.text((syn_hi + obs) / 2, 0.78, f"~{l2.obs_over_syn_pinned_ratio:.0f}x", fontsize=8,
            color="0.2", ha="center", va="bottom")
    ax.text(0.97, 0.97, "flag: under-produces by\n>1 order of magnitude\n(sufficient -> contributory)",
            transform=ax.transAxes, fontsize=6, ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="#fdecea", ec="#d98880", lw=0.6))


def panel_c(ax):
    cd = pd.read_csv(os.path.join(RESULTS, "lever_ceiling.csv"))
    crc = cd[cd.dataset == "crc"].iloc[0]
    roi = pd.read_csv(os.path.join(RESULTS, "lever_roi_pin199.csv"))
    jan = float(roi[roi.roi_def == "tp_cells_own"].oracle_point.iloc[0])
    # bars: ceiling position colored by orthogonal direction
    rows = [("Janesick breast\ntriple-positive (ROI)", jan, "DISAGREES", ORACLE),
            ("colorectal tumor\n(packing 17,424)", float(crc.oracle_own_median_point), "AGREES", OI["green"]),
            ("MERFISH hypothalamus\nhybrid neurons", 0.836, "AGREES", OI["green"])]
    xs = np.arange(len(rows))
    ax.bar(xs, [r[1] for r in rows], color=[r[3] for r in rows], width=0.6, edgecolor="0.3", lw=0.5)
    for x, r in zip(xs, rows):
        ax.text(x, r[1] + 0.012, f"{r[1]:.3f}\n{r[2]}", ha="center", fontsize=6.5,
                color=("#b2182b" if r[2] == "DISAGREES" else "#1b7837"), fontweight="bold")
    ax.axhline(0.90, ls="--", color="0.4", lw=1.0)
    ax.text(len(rows) - 0.5, 0.905, "0.90", fontsize=6.5, color="0.4", va="bottom", ha="right")
    ax.set_xticks(xs); ax.set_xticklabels([r[0] for r in rows], fontsize=6.3)
    ax.set_ylim(0.4, 1.0); ax.set_ylabel("oracle ceiling (assignment accuracy)")
    ax.set_title("c  Leg 3 and the artifact-consistent vs corroborated contrast", loc="left", fontsize=8.5)
    txt = ("Leg 3 orthogonal evidence (Janesick): dissociated scFFPE-seq finds ~30 PGR+ cells, none co-expressing "
           "ESR1 or ERBB2; ERBB2 detected in 87% of cells; tissue clinically PR-negative. The spatial assay "
           "over-reports the co-expression; both orthogonal sources under-report it (the direction misassignment "
           "predicts). CRC and MERFISH sit above/in the hard regime too, but their orthogonal data AGREES, so they "
           "are corroborated, not artifact-consistent: the ceiling alone does not flag a claim.")
    ax.text(0.5, -0.42, txt, transform=ax.transAxes, ha="center", va="top", fontsize=6.2, wrap=True,
            bbox=dict(boxstyle="round,pad=0.4", fc="0.96", ec="0.7", lw=0.6))


def main():
    apply_style()
    fig = plt.figure(figsize=(7.2, 7.4))
    gs = GridSpec(2, 2, figure=fig, hspace=0.5, wspace=0.32,
                  left=0.09, right=0.97, top=0.94, bottom=0.2, height_ratios=[1, 1])
    panel_a(fig.add_subplot(gs[0, 0]))
    panel_b(fig.add_subplot(gs[0, 1]))
    panel_c(fig.add_subplot(gs[1, :]))
    save(fig, "fig4_worked")


if __name__ == "__main__":
    main()
