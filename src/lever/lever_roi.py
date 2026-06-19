"""Lever analysis Task 1: the answerability ceiling at the triple-positive ROI's OWN
local density.

The scouting evaluated the ceiling at the dataset median and p90. Here we localise the
Janesick triple-positive DCIS region in the Xenium coordinates and evaluate the ceiling
at THAT region's own local packing and density.

Locating the ROI (reading the published segmented cell-by-gene matrix and centroids;
this is NOT a re-assignment of transcripts): the paper describes a small DCIS region
positive for ESR1, PGR and ERBB2. We take the apparent triple-positive cells in the
released matrix (ESR1>=t and ERBB2>=t and PGR>=t) and measure the local cell packing at
their own locations (k=10 nearest-neighbour local density, the same estimator the
diagnostic uses). Because the exact ROI polygon is not recoverable from the paper, we
also report a conservative proxy: the densest DCIS-nest packing (high-percentile local
packing among ESR1+/ERBB2+ DCIS cells). The oracle ceiling is evaluated at both.

Sigma is the data-pinned value already measured on this dataset by the lever diagnostic
(results/lever/lever_ceiling.csv, breast row; ~2.1 um, CI [1.43, 2.50]); we reuse it
rather than re-pinning, and cite it.

This does not run the oracle or any segmenter on the real transcripts. It only reads the
published matrix to locate the region, then evaluates the method-independent ceiling on
synthetic fields at the measured packing/density. Output: results/lever/lever_roi.csv.
Seeds: LEVER_SEED + offsets, recorded.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config  # noqa: E402
from expression import build_realistic_model_from_xenium, load_xenium_cells_genes  # noqa: E402
from gate2_pin import oracle_acc_at, ORACLE_RES_CELL, ORACLE_GRID_MAX  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS = os.path.join(ROOT, "results", "lever")
DATA = os.path.join(ROOT, "data")
LEVER_SEED = config.MASTER_SEED + 910000

H5 = os.path.join(DATA, "xenium_breast_rep1_cell_feature_matrix.h5")
PARQ = os.path.join(DATA, "xenium_breast_rep1_cells.parquet")
MARKERS = ["ESR1", "PGR", "ERBB2"]
POS_THRESH = 2          # primary positivity threshold (counts); a sweep is also reported
KNN = 10                # neighbours for local packing (matches the diagnostic estimator)


def local_packing(coords, k=KNN):
    """k-NN local cell packing (cells/mm^2) per cell: k / (pi * d_k^2), um -> mm^2."""
    d, _ = cKDTree(coords).query(coords, k=k + 1)
    return k / (np.pi * d[:, k] ** 2) * 1e6


def pinned_sigma():
    df = pd.read_csv(os.path.join(RESULTS, "lever_ceiling.csv"))
    r = df[df.dataset == "breast"].iloc[0]
    return float(r.sigma_point_um), float(r.sigma_ci_lo_um), float(r.sigma_ci_hi_um)


def main():
    os.makedirs(RESULTS, exist_ok=True)
    X, genes, coords, total = load_xenium_cells_genes(H5, PARQ)
    genes = list(genes)
    iE, iP, iH = (genes.index(m) for m in MARKERS)
    esr1, pgr, erbb2 = X[:, iE], X[:, iP], X[:, iH]
    lpack = local_packing(coords)
    N = X.shape[0]
    section_med_pack = float(np.median(lpack))
    section_med_dens = float(np.median(total))
    print(f"section: N={N} median_local_pack={section_med_pack:.0f} median_tx={section_med_dens:.0f}")

    # apparent triple-positive and DCIS double-positive cells in the released matrix
    t = POS_THRESH
    tp = (esr1 >= t) & (erbb2 >= t) & (pgr >= t)
    dcis = (esr1 >= t) & (erbb2 >= t)              # ER+/HER2+ DCIS double-positive
    print(f"thr>={t}: triple-positive cells={int(tp.sum())} ({tp.mean()*100:.4f}%); "
          f"DCIS ER+/HER2+ double={int(dcis.sum())} ({dcis.mean()*100:.3f}%)")

    # ROI-local packing/density: at the triple-positive cells' own locations, and a small
    # spatial neighbourhood around them (cells within R of any triple-positive cell)
    R_ROI = 50.0  # um neighbourhood radius defining the triple-positive ROI tissue
    tp_xy = coords[tp]
    nbr = cKDTree(coords).query_ball_point(tp_xy, R_ROI)
    roi_idx = np.unique(np.concatenate([np.array(ix, dtype=int) for ix in nbr])) if len(tp_xy) else np.array([], int)
    roi_pack_med = float(np.median(lpack[roi_idx])) if roi_idx.size else np.nan
    roi_dens_med = float(np.median(total[roi_idx])) if roi_idx.size else np.nan
    tp_pack_med = float(np.median(lpack[tp])) if tp.sum() else np.nan
    tp_dens_med = float(np.median(total[tp])) if tp.sum() else np.nan
    # conservative proxy: densest DCIS-nest packing (high percentile among DCIS cells)
    dcis_pack = lpack[dcis]
    dcis_p90 = float(np.percentile(dcis_pack, 90)) if dcis.sum() else np.nan
    dcis_p95 = float(np.percentile(dcis_pack, 95)) if dcis.sum() else np.nan
    dcis_med = float(np.median(dcis_pack)) if dcis.sum() else np.nan
    print(f"ROI(<= {R_ROI:.0f}um of a triple-positive cell): n={roi_idx.size} "
          f"median_local_pack={roi_pack_med:.0f} median_tx={roi_dens_med:.0f}")
    print(f"triple-positive cells' own: median_local_pack={tp_pack_med:.0f} median_tx={tp_dens_med:.0f}")
    print(f"DCIS ER+/HER2+ nest packing: median={dcis_med:.0f} p90={dcis_p90:.0f} p95={dcis_p95:.0f}")

    # build the realistic breast model (same as the diagnostic) for the oracle ceiling
    model, _ = build_realistic_model_from_xenium(n_types=15, seed=config.MASTER_SEED,
                                                 h5_path=H5, cells_parquet=PARQ,
                                                 name="realistic_xenium_breast")
    sig_pt, sig_lo, sig_hi = pinned_sigma()
    print(f"data-pinned sigma (reused from lever diagnostic): {sig_pt:.2f} um CI [{sig_lo:.2f}, {sig_hi:.2f}]")

    def oa(pack, dens, sig, sd):
        return oracle_acc_at(model, pack, dens, sig, sd,
                             res_cell=ORACLE_RES_CELL, grid_max=ORACLE_GRID_MAX)

    # evaluate the ceiling at three ROI definitions x three sigmas
    roi_defs = {
        "tp_neighborhood": (roi_pack_med, roi_dens_med),   # tissue around the triple-positive cells
        "tp_cells_own":    (tp_pack_med, tp_dens_med),      # the triple-positive cells' own local packing
        "dcis_nest_p90_proxy": (dcis_p90, roi_dens_med if np.isfinite(roi_dens_med) else section_med_dens),
    }
    rows = []
    for j, (name, (pack, dens)) in enumerate(roi_defs.items()):
        if not (np.isfinite(pack) and np.isfinite(dens)):
            continue
        o_pt, n_pt = oa(pack, dens, sig_pt, LEVER_SEED + 100 * j + 1)
        o_lo, _ = oa(pack, dens, sig_lo, LEVER_SEED + 100 * j + 3)   # low sigma -> optimistic
        o_hi, _ = oa(pack, dens, sig_hi, LEVER_SEED + 100 * j + 5)   # high sigma -> pessimistic
        above_090 = o_lo < 0.90    # optimistic end still below 0.90
        above_095 = o_lo < 0.95
        rows.append({
            "roi_def": name, "packing_cells_per_mm2": pack, "density_tx_per_cell": dens,
            "sigma_point_um": sig_pt, "sigma_ci_lo_um": sig_lo, "sigma_ci_hi_um": sig_hi,
            "oracle_point": o_pt, "oracle_opt_ci_lo": o_lo, "oracle_pess_ci_hi": o_hi,
            "naive_point": n_pt,
            "above_ceiling_0p90_whole_ci": bool(above_090),
            "above_ceiling_0p95_whole_ci": bool(above_095),
            "pos_threshold": t, "roi_radius_um": R_ROI,
            "section_median_packing": section_med_pack, "section_median_density": section_med_dens,
            "n_triple_positive_cells": int(tp.sum()), "n_dcis_cells": int(dcis.sum()),
            "dcis_nest_pack_median": dcis_med, "dcis_nest_pack_p90": dcis_p90, "dcis_nest_pack_p95": dcis_p95,
            "seed": LEVER_SEED + 100 * j,
        })
        print(f"[{name}] pack={pack:.0f} dens={dens:.0f} -> oracle {o_pt:.3f} "
              f"(opt {o_lo:.3f}, pess {o_hi:.3f}); naive {n_pt:.3f}; "
              f"above 0.90 across CI: {above_090}; above 0.95 across CI: {above_095}")

    df = pd.DataFrame(rows)
    out = os.path.join(RESULTS, "lever_roi.csv")
    df.to_csv(out, index=False)
    print(f"\nwrote {out} ({len(df)} ROI definitions)")
    return df


if __name__ == "__main__":
    main()
