"""Render the 10 Extended Data display items from committed results (decision D14).

No new analysis: every value is read from a committed results CSV (or a committed report,
for cell counts). Tables are rendered as vector+raster table figures. Imports figstyle
read-only (does not modify it). Outputs: figures/ED01_*..ED10_*.{png,pdf}.

ED item -> committed source:
  ED01 full frontier            results/sweep.csv
  ED02 generator calibration    results/realism_oracle.csv
  ED03 overlapping expression   results/gate1_sweep.csv vs results/sweep.csv
  ED04 model shape              results/gate1_structural.csv, gate2_nbinom_frontier.csv, gate2_nbinom.csv
  ED05 marker validity          results/gate3_validity_breast.csv, gate3_validity_lung.csv
  ED06 cell-typing stability    results/gate3_typing.csv
  ED07 selection sensitivity    results/gate2_sigma_uncertainty.csv, gate3_pin_breast.csv, gate3_pin_lung.csv
  ED08 method headroom metrics  results/headroom_fair_metrics.csv, headroom_fair_free.csv
  ED09 per-dataset summary      results/lever/lever_ceiling.csv, gate3_pin_lung.csv, realism.csv + committed reports
  ED10 worked-example detail    results/lever/lever_roi.csv, lever_mechanism.csv, leg2_local.csv, lever_roi_pin199.csv
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
LEVER = os.path.join(RESULTS, "lever")


def _table(ax, cell_text, col_labels, fontsize=6, highlight_rows=None, scale_y=1.3, colw=None):
    ax.axis("off")
    tbl = ax.table(cellText=cell_text, colLabels=col_labels, loc="center", cellLoc="center",
                   colWidths=colw)
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(fontsize)
    tbl.scale(1, scale_y)
    hi = set(highlight_rows or [])
    for (r, c), cell in tbl.get_celld().items():
        cell.set_linewidth(0.3)
        if r == 0:
            cell.set_facecolor("#d9d9d9"); cell.set_text_props(fontweight="bold")
        elif (r - 1) in hi:
            cell.set_facecolor("#fde0dc")
    return tbl


# ----------------------------------------------------------------------------- ED01
def ed01():
    df = pd.read_csv(os.path.join(RESULTS, "sweep.csv"))
    packs = np.sort(df.packing_cells_per_mm2.unique())
    sigs = np.sort(df.sigma_um.unique())
    dens = np.sort(df.mean_tx_per_cell.unique())
    fig = plt.figure(figsize=(7.2, 4.2))
    gs = GridSpec(1, len(dens), figure=fig, wspace=0.12, left=0.08, right=0.9, top=0.78, bottom=0.16)
    vmin, vmax = df.oracle_acc.min(), df.oracle_acc.max()
    for k, d in enumerate(dens):
        ax = fig.add_subplot(gs[0, k])
        sub = df[df.mean_tx_per_cell == d]
        M = np.full((len(packs), len(sigs)), np.nan)
        for _, r in sub.iterrows():
            M[np.where(packs == r.packing_cells_per_mm2)[0][0], np.where(sigs == r.sigma_um)[0][0]] = r.oracle_acc
        im = ax.imshow(M, origin="lower", cmap="viridis", aspect="auto", vmin=0.23, vmax=0.96)
        ax.set_title(f"{int(d)} tx/cell", fontsize=7)
        ax.set_xticks(range(len(sigs))); ax.set_xticklabels([f"{s:g}" for s in sigs], fontsize=5.5)
        if k == 0:
            ax.set_yticks(range(len(packs))); ax.set_yticklabels([f"{int(p/1000)}k" for p in packs], fontsize=5.5)
            ax.set_ylabel("packing (cells/mm$^2$)", fontsize=7)
        else:
            ax.set_yticks([])
        ax.set_xlabel("sigma (um)", fontsize=6.5)
    cax = fig.add_axes([0.92, 0.16, 0.015, 0.62])
    plt.colorbar(im, cax=cax).set_label("oracle accuracy", fontsize=7)
    fig.suptitle("Extended Data 1 | Full answerability frontier: oracle accuracy over packing x displacement,\n"
                 "per density slice (full-grid range 0.234 to 0.957)", fontsize=9, x=0.04, ha="left")
    save(fig, "ED01_full_frontier")


# ----------------------------------------------------------------------------- ED02
def ed02():
    df = pd.read_csv(os.path.join(RESULTS, "realism_oracle.csv"))
    cmap = {"xenium_breast": OI["vermillion"], "merfish_hypothal": OI["blue"]}
    nm = {"xenium_breast": "Xenium breast", "merfish_hypothal": "MERFISH hypothalamus"}
    fig, (axt, axp) = plt.subplots(1, 2, figsize=(7.2, 3.6))
    for ds, col in cmap.items():
        s = df[df.dataset == ds]
        axt.scatter(s.real_median_tx_per_cell, s.realized_median_tx_per_cell, c=col, s=22, label=nm[ds], edgecolor="none")
        axp.scatter(s.real_packing_cells_per_mm2, s.realized_packing_cells_per_mm2, c=col, s=22, label=nm[ds], edgecolor="none")
    for ax, lab in [(axt, "median transcripts/cell"), (axp, "packing (cells/mm$^2$)")]:
        lo = min(ax.get_xlim()[0], ax.get_ylim()[0]); hi = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([lo, hi], [lo, hi], "--", color="0.5", lw=1)
        ax.set_xlabel(f"real {lab}"); ax.set_ylabel(f"generator {lab}"); ax.legend(fontsize=6)
    txmax = df.tx_match_rel_err.max() * 100
    pkmax = df.packing_match_rel_err.max()
    axt.set_title(f"a  Transcripts/cell (max rel. err {txmax:.2f}%)", loc="left", fontsize=8.5)
    axp.set_title(f"b  Packing (max rel. err {pkmax:.1e})", loc="left", fontsize=8.5)
    fig.suptitle("Extended Data 2 | Generator calibration to real Xenium and MERFISH fields", fontsize=9, x=0.04, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save(fig, "ED02_calibration")


# ----------------------------------------------------------------------------- ED03
def ed03():
    g0 = pd.read_csv(os.path.join(RESULTS, "sweep.csv"))
    g1 = pd.read_csv(os.path.join(RESULTS, "gate1_sweep.csv"))
    key = ["packing_cells_per_mm2", "sigma_um", "mean_tx_per_cell"]
    m = g0[key + ["oracle_acc", "oracle_minus_naive"]].merge(
        g1[key + ["oracle_acc", "oracle_minus_naive"]], on=key, suffixes=("_g0", "_g1"))
    mad = float((m.oracle_acc_g1 - m.oracle_acc_g0).abs().mean())
    gap_g0 = float(g0.oracle_minus_naive.mean()); gap_g1 = float(g1.oracle_minus_naive.mean())
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.2, 3.6))
    axa.scatter(m.oracle_acc_g0, m.oracle_acc_g1, c=m.sigma_um, cmap="cividis", s=20, edgecolor="none")
    axa.plot([0.2, 1], [0.2, 1], "--", color="0.5", lw=1)
    axa.set_xlabel("Gate 0 oracle (disjoint markers)"); axa.set_ylabel("Gate 1 oracle (overlapping)")
    axa.set_title(f"a  Frontier overlay\nmean |diff| {mad:.3f}", loc="left", fontsize=8.5)
    axb.bar([0, 1], [gap_g0, gap_g1], color=[OI["grey"], OI["green"]], width=0.6)
    for x, v in zip([0, 1], [gap_g0, gap_g1]):
        axb.text(x, v + 0.002, f"{v:.3f}", ha="center", fontsize=7)
    axb.set_xticks([0, 1]); axb.set_xticklabels(["Gate 0", "Gate 1"])
    axb.set_ylabel("mean oracle - naive gap")
    axb.set_title(f"b  Oracle-minus-naive gap preserved\n({gap_g1:.3f} vs {gap_g0:.3f})", loc="left", fontsize=8.5)
    fig.suptitle("Extended Data 3 | Robustness to realistic overlapping expression", fontsize=9, x=0.04, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save(fig, "ED03_overlapping_expression")


# ----------------------------------------------------------------------------- ED04
def ed04():
    st = pd.read_csv(os.path.join(RESULTS, "gate1_structural.csv"))
    nb = pd.read_csv(os.path.join(RESULTS, "gate2_nbinom_frontier.csv"))
    nbm = pd.read_csv(os.path.join(RESULTS, "gate2_nbinom.csv")).iloc[0]
    packs = [2575.0, 6000.0, 13625.0]
    cond_style = {"baseline": ("-", OI["black"]), "aniso": ("--", OI["orange"]), "mixture": (":", OI["purple"])}
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.4), sharey=True)
    for ax, pk in zip(axes, packs):
        for cond, (ls, col) in cond_style.items():
            s = st[(st.condition == cond) & (st.packing_cells_per_mm2 == pk)].sort_values("sigma_um")
            ax.plot(s.sigma_um, s.oracle_acc, ls, color=col, lw=1.5, ms=4, marker="o", label=cond)
        ax.set_title(f"{int(pk):,} cells/mm$^2$", fontsize=7.5)
        ax.set_xlabel("sigma (um)")
    axes[0].set_ylabel("oracle accuracy")
    axes[0].legend(fontsize=6, loc="upper right")
    # NB overlay on the dense panel
    nbd = nb[nb.packing_cells_per_mm2 == 13625.0].sort_values("sigma_um")
    axes[2].plot(nbd.sigma_um, nbd.oracle_acc, "-", color=OI["skyblue"], lw=1.4, marker="s", ms=4,
                 label="neg-binom")
    axes[2].legend(fontsize=6, loc="upper right")
    fig.suptitle(f"Extended Data 4 | Robustness to model shape: non-Voronoi (aniso), non-Gaussian (mixture),\n"
                 f"and negative-binomial variance (dispersion {nbm.nb_dispersion_median:.2f}, dense shift "
                 f"{nbm.dense_oracle_shift:+.3f}); oracle >= naive throughout", fontsize=8.5, x=0.04, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    save(fig, "ED04_model_shape")


# ----------------------------------------------------------------------------- ED05
def ed05():
    fails = {"AHSP", "CRHBP", "CYP1A1", "EGFL7", "SLC4A1", "MZB1"}

    def prep(path):
        d = pd.read_csv(path).dropna(subset=["pval"]).copy()
        d = d.sort_values(["admissible", "pval"], ascending=[False, True])
        txt = [[r.gene_name, f"{r.real_signal:+.3f}", f"{r.pval:.3f}", "Y" if r.admissible else "n"]
               for _, r in d.iterrows()]
        hl = [i for i, (_, r) in enumerate(d.iterrows()) if (not r.admissible)]
        return txt, hl, int(d.admissible.sum()), len(d)

    bt, bhl, ba, bn = prep(os.path.join(RESULTS, "gate3_validity_breast.csv"))
    lt, lhl, la, ln = prep(os.path.join(RESULTS, "gate3_validity_lung.csv"))
    fig = plt.figure(figsize=(7.2, 9.6))
    gs = GridSpec(1, 2, figure=fig, wspace=0.18, left=0.04, right=0.98, top=0.92, bottom=0.03)
    axb = fig.add_subplot(gs[0, 0]); axl = fig.add_subplot(gs[0, 1])
    cols = ["gene", "signal", "p", "adm"]
    _table(axb, bt, cols, fontsize=5.0, highlight_rows=bhl, scale_y=1.05)
    axb.set_title(f"breast: {ba}/{bn} admissible (rejected shaded)", fontsize=8)
    _table(axl, lt, cols, fontsize=5.0, highlight_rows=lhl, scale_y=1.05)
    axl.set_title(f"lung: {la}/{ln} admissible", fontsize=8)
    fig.suptitle("Extended Data 5 | Marker-validity test: adjacent-minus-distant signal and permutation p\n"
                 "(admissible if p < 0.05 over 1,000 permutations)", fontsize=9, x=0.04, ha="left")
    save(fig, "ED05_marker_validity")


# ----------------------------------------------------------------------------- ED06
def ed06():
    t = pd.read_csv(os.path.join(RESULTS, "gate3_typing.csv"))
    lab = {"kmeans_K10": "K=10", "kmeans_K14": "K=14", "kmeans_K14_altseed": "K=14 (alt seed)", "kmeans_K20": "K=20"}
    t = t.set_index("variant").loc[list(lab)].reset_index()
    txt = [[lab[r.variant], f"{int(r.n_admissible)}/{int(r.n_candidate)}",
            f"{r.sigma_point_um:.2f}", f"[{r.sigma_ci_lo_um:.2f}, {r.sigma_ci_hi_um:.2f}]",
            f"{r.dense_oracle_point:.3f}", f"{r.dense_oracle_at_sigma_ci_lo:.3f}"]
           for _, r in t.iterrows()]
    cols = ["typing", "admissible", "sigma\n(um)", "sigma 95% CI", "dense\noracle", "dense\n(optimistic)"]
    fig, ax = plt.subplots(figsize=(7.2, 2.7))
    _table(ax, txt, cols, fontsize=7, scale_y=1.4, colw=[0.20, 0.16, 0.12, 0.20, 0.16, 0.16])
    ax.set_title("Extended Data 6 | Cell-typing stability: pinned displacement 1.99 to 2.65 um,\n"
                 "dense optimistic oracle accuracy at most 0.803 across four cell-typing choices",
                 loc="left", fontsize=8.5, y=1.02)
    save(fig, "ED06_typing_stability")


# ----------------------------------------------------------------------------- ED07
def ed07():
    u = pd.read_csv(os.path.join(RESULTS, "gate2_sigma_uncertainty.csv")).iloc[0]
    b = pd.read_csv(os.path.join(RESULTS, "gate3_pin_breast.csv")).iloc[0]
    l = pd.read_csv(os.path.join(RESULTS, "gate3_pin_lung.csv")).iloc[0]
    # committed dense-oracle-vs-sigma points (Gate 2 converged grid)
    sig = np.array([u.sigma_combined_lo_um, u.sigma_stat_ci_lo_um, u.sigma_point_um, u.sigma_stat_ci_hi_um])
    orc = np.array([u.oracle_dense_at_combined_lo, u.oracle_dense_at_stat_ci_lo,
                    u.oracle_dense_at_point, u.oracle_dense_at_stat_ci_hi])
    order = np.argsort(sig); sig, orc = sig[order], orc[order]
    cross = float(np.interp(0.95, orc[::-1], sig[::-1]))  # sigma where dense oracle = 0.95 (reading the curve)
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ax.plot(sig, orc, "o-", color=ORACLE, lw=1.8, ms=6, label="dense oracle (committed points)")
    ax.axhline(0.95, ls=":", color="0.4", lw=1); ax.axhline(0.90, ls="--", color="0.5", lw=1)
    ax.text(2.5, 0.952, "0.95", fontsize=6.5, color="0.4"); ax.text(2.5, 0.902, "0.90", fontsize=6.5, color="0.5")
    ax.axvspan(b.sigma_ci_lo_um, b.sigma_ci_hi_um, color="0.88", alpha=0.6, label="breast pin CI [1.43, 2.63]")
    ax.axvline(b.sigma_point_um, color=OI["vermillion"], lw=1.4)
    ax.text(b.sigma_point_um + 0.03, 0.60, f"breast pin {b.sigma_point_um:.2f} um", color=OI["vermillion"], fontsize=6.5, rotation=90, va="bottom")
    ax.axvline(l.sigma_point_um, color=OI["blue"], lw=1.2, ls="--")
    ax.text(l.sigma_point_um - 0.06, 0.60, f"lung pin {l.sigma_point_um:.2f} um ({int(l.n_admissible)}/{int(l.n_candidate)})",
            color=OI["blue"], fontsize=6.5, rotation=90, va="bottom", ha="right")
    ax.annotate("", xy=(b.sigma_point_um, 0.985), xytext=(cross, 0.985),
                arrowprops=dict(arrowstyle="<->", color="0.3", lw=1.0))
    ax.text((cross + b.sigma_point_um) / 2, 0.99, f"~{b.sigma_point_um / cross:.1f}x fall to exit",
            ha="center", va="bottom", fontsize=7, color="0.2")
    ax.plot([cross], [0.95], "kx", ms=7); ax.text(cross, 0.93, f"0.95 at ~{cross:.2f} um", fontsize=6.5, ha="center")
    ax.set_xlabel("displacement sigma (um)"); ax.set_ylabel("dense oracle accuracy")
    ax.set_xlim(0, 2.8); ax.set_ylim(0.55, 1.02); ax.legend(fontsize=6.5, loc="lower left")
    ax.set_title("Extended Data 7 | Selection-on-significance margin: displacement must fall ~5x to lift dense\n"
                 "accuracy out of the constrained regime; lung (near-zero selection, 31/32) gives the same pin",
                 loc="left", fontsize=8.5)
    save(fig, "ED07_selection_sensitivity")


# ----------------------------------------------------------------------------- ED08
def ed08():
    df = pd.read_csv(os.path.join(RESULTS, "headroom_fair_metrics.csv"))
    fr = pd.read_csv(os.path.join(RESULTS, "headroom_fair_free.csv"))
    order_cfg = ["dense_sigma1.43", "dense_sigma1.99", "dense_sigma2.63", "sparse_sigma1.99", "representative_sigma1.99"]
    order_m = ["oracle", "naive", "Baysor", "ComSeg", "Proseg"]
    df = df[df.status == "ok"].copy()
    df["ck"] = df.config.apply(lambda c: order_cfg.index(c) if c in order_cfg else 9)
    df["mk"] = df.method.apply(lambda m: order_m.index(m) if m in order_m else 9)
    df = df.sort_values(["ck", "mk"])
    short = {"dense_sigma1.43": "dense 1.43", "dense_sigma1.99": "dense 1.99", "dense_sigma2.63": "dense 2.63",
             "sparse_sigma1.99": "sparse 1.99", "representative_sigma1.99": "repr 1.99"}
    txt, hl = [], []
    for _, r in df.iterrows():
        txt.append([short[r.config], r.method, f"{r.acc_one_to_one:.3f}", f"{r.acc_many_to_one:.3f}",
                    f"{r.ari:.3f}", f"{r.frac_assigned:.2f}"])
    # free-segmentation Baysor rows
    frb = fr[(fr.method == "Baysor") & (fr.status == "ok")].copy()
    frb["ck"] = frb.config.apply(lambda c: order_cfg.index(c) if c in order_cfg else 9)
    for _, r in frb.sort_values("ck").iterrows():
        txt.append([short.get(r.config, r.config), "Baysor (free)", f"{r.acc_one_to_one:.3f}",
                    f"{r.acc_many_to_one:.3f}", f"{r.ari:.3f}", f"{r.frac_assigned:.2f}"])
        hl.append(len(txt) - 1)
    cols = ["regime", "method", "one-to-one", "many-to-one", "co-assign (ARI)", "frac assigned"]
    fig, ax = plt.subplots(figsize=(7.2, 8.2))
    _table(ax, txt, cols, fontsize=5.6, highlight_rows=hl, scale_y=1.08)
    ax.set_title("Extended Data 8 | Method-headroom metrics: Baysor, Proseg, ComSeg vs oracle and nearest-nucleus\n"
                 "(one-to-one, many-to-one homogeneity, co-assignment ARI); free-segmentation rows shaded",
                 loc="left", fontsize=8.5)
    save(fig, "ED08_method_metrics")


# ----------------------------------------------------------------------------- ED09
def ed09():
    cd = pd.read_csv(os.path.join(LEVER, "lever_ceiling.csv"))
    br = cd[cd.dataset == "breast"].iloc[0]; cr = cd[cd.dataset == "crc"].iloc[0]
    lp = pd.read_csv(os.path.join(RESULTS, "gate3_pin_lung.csv")).iloc[0]
    rm = pd.read_csv(os.path.join(RESULTS, "realism.csv"))
    mer = rm[rm.dataset == "merfish_hypothal"]
    # cell counts: breast/lung from committed reports; CRC/MERFISH flagged as dataset provenance
    rows = [
        ["Xenium breast", "167,780*", f"{br.pack_median_cells_per_mm2:,.0f}",
         f"{br.pack_p90_cells_per_mm2:,.0f}", f"{br.median_nn_um:.2f}", f"{br.density_median_tx_per_cell:.0f}"],
        ["Xenium lung", "150,365*", f"{lp.pin_packing:,.0f}", "n/a", "n/a", f"{lp.density:.0f}"],
        ["Xenium colorectal", "388,175†", f"{cr.pack_median_cells_per_mm2:,.0f}",
         f"{cr.pack_p90_cells_per_mm2:,.0f}", f"{cr.median_nn_um:.2f}", f"{cr.density_median_tx_per_cell:.0f}"],
        ["MERFISH hypothal.", "100-108/FOV‡", f"{mer.packing_cells_per_mm2.median():,.0f}",
         "n/a", f"{mer.median_nn_distance_um.median():.2f}", f"{mer.median_tx_per_cell.median():.0f}"],
    ]
    cols = ["dataset", "cells", "packing median\n(cells/mm$^2$)", "packing p90", "median NN\n(um)", "median tx/cell"]
    fig, ax = plt.subplots(figsize=(7.6, 2.8))
    _table(ax, rows, cols, fontsize=6.5, scale_y=1.5, colw=[0.20, 0.16, 0.18, 0.14, 0.14, 0.18])
    ax.set_title("Extended Data 9 | Per-dataset summary statistics", loc="left", fontsize=8.5)
    note = ("* cell count from committed report (GATE2/GATE3_REPORT). † colorectal count: 10x portal "
            "(dataset provenance, not a committed CSV). ‡ MERFISH: committed per-FOV counts (realism.csv); "
            "section total not in a committed file. Lung median NN was computed but not committed to CSV (shown n/a). "
            "Packing/NN/tx are committed-CSV values (lever_ceiling.csv, gate3_pin_lung.csv, realism.csv).")
    ax.text(0.0, -0.18, note, transform=ax.transAxes, fontsize=5.6, va="top", wrap=True,
            bbox=dict(boxstyle="round,pad=0.3", fc="0.96", ec="0.7", lw=0.5))
    save(fig, "ED09_dataset_summary")


# ----------------------------------------------------------------------------- ED10
def ed10():
    roi = pd.read_csv(os.path.join(LEVER, "lever_roi_pin199.csv"))
    mech = pd.read_csv(os.path.join(LEVER, "lever_mechanism.csv")).sort_values("sigma_um")
    l2 = pd.read_csv(os.path.join(LEVER, "leg2_local.csv")).iloc[0]
    params = pd.read_csv(os.path.join(LEVER, "lever_mechanism_params.csv")).iloc[0]
    fig = plt.figure(figsize=(7.2, 6.4))
    gs = GridSpec(2, 1, figure=fig, height_ratios=[1.0, 1.25], hspace=0.5, left=0.1, right=0.95, top=0.9, bottom=0.08)

    # (a) ROI localization: oracle ceiling at locked pin vs packing
    ax = fig.add_subplot(gs[0, 0])
    names = {"tp_neighborhood": "ROI tissue\n(7,103)", "tp_cells_own": "triple+ cells\n(7,985)",
             "dcis_nest_p90_proxy": "densest nest\n(9,329)"}
    xs = np.arange(len(roi))
    ax.errorbar(xs, roi.oracle_point, yerr=[roi.oracle_point - roi.oracle_pess_ci_hi,
                roi.oracle_opt_ci_lo - roi.oracle_point], fmt="o", color=ORACLE, capsize=5, ms=7, lw=1.6)
    for x, p in zip(xs, roi.oracle_point):
        ax.text(x + 0.12, p, f"{p:.3f}", va="center", fontsize=6.5)
    ax.axhline(0.90, ls="--", color="0.4", lw=1); ax.axhline(0.95, ls=":", color="0.55", lw=1)
    ax.set_xticks(xs); ax.set_xticklabels([names[v] for v in roi.roi_def], fontsize=6.5)
    ax.set_ylim(0.6, 1.0); ax.set_ylabel("oracle accuracy")
    ax.set_title(f"a  ROI localization and ceiling (pin 1.99 um; section median packing "
                 f"{params.roi_pack_median:.0f} -> ROI {params.roi_density_median:.0f} tx/cell; "
                 f"PGR+ in 2.9% of cells)", loc="left", fontsize=8)

    # (b) full Leg 2 mechanism table
    axt = fig.add_subplot(gs[1, 0])
    nm = {"sigma0": "0 (control)", "ci_lo": "1.43 (CI lo)", "pinned": "2.10 (pinned)", "ci_hi": "2.50 (CI hi)"}
    txt = []
    for _, r in mech.iterrows():
        txt.append([nm.get(r.sigma_name, f"{r.sigma_um:g}"),
                    f"{r.spurious_triple_rate_among_dcis_t2*100:.2f}",
                    f"{r.spurious_triple_rate_among_all_t2*100:.2f}",
                    f"{r.spurious_triple_rate_dcis_adj_pr_t2*100:.2f}",
                    f"{int(r.true_triple_count_t2)}"])
    txt.append(["observed real (local)", "-", f"{l2.obs_rate_section_wide*100:.2f}",
                f"{l2.obs_rate_adj_to_prsource*100:.2f}", "n/a"])
    cols = ["displacement sigma (um)", "field avg\n(among DCIS) %", "section/all\n%", "adjacent to\nPR source %", "true\ntriples"]
    _table(axt, txt, cols, fontsize=6.2, highlight_rows=[len(txt) - 1], scale_y=1.45)
    axt.set_title(f"b  Leg 2 mechanism (thr>=2): synthetic null (zero true triples) vs observed; observed adjacency "
                  f"{l2.obs_rate_adj_to_prsource*100:.1f}% vs synthetic 1.26-1.96% ({l2.obs_over_syn_pinned_ratio:.0f}x): "
                  f"displacement is contributory, not sufficient", loc="left", fontsize=7.6)
    fig.suptitle("Extended Data 10 | Worked-example detail: ROI localization and the full Leg 2 mechanism table",
                 fontsize=9, x=0.04, ha="left")
    save(fig, "ED10_worked_detail")


def main():
    apply_style()
    only = sys.argv[1:]
    items = {"1": ed01, "2": ed02, "3": ed03, "4": ed04, "5": ed05,
             "6": ed06, "7": ed07, "8": ed08, "9": ed09, "10": ed10}
    for k, fn in items.items():
        if not only or k in only:
            fn()


if __name__ == "__main__":
    main()
