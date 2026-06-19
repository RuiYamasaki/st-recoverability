"""Figure-prep / Task 1 (D13): the only new computation in the figures phase.

Two presentation refinements, both append-only, neither a gate re-run:

(A) Leg 2 local-to-local. The committed mechanism comparison contrasted the synthetic
    ADJACENCY-conditioned spurious rate (DCIS cells next to a PR source) with the
    section-wide observed rate, which is not apples-to-apples (the artifact is local).
    Here we compute the OBSERVED apparent-triple rate matched to the same local
    condition: among real ER+/HER2+ DCIS cells that are spatially adjacent to a real
    PR-source cell (PGR+ non-tumor), the apparent-triple-positive rate at threshold >=2.
    This reads the released cell-by-gene matrix and centroids only; it is NOT a re-
    assignment of transcripts and does not adjudicate the real claim. Compared against
    the committed synthetic adjacency rate (1.26% pinned, 1.96% upper CI) and the
    zero-displacement control (0%). Order-of-magnitude sufficiency flag per the spec.

(B) ROI-ceiling reconciliation to the locked breast pin. lever_roi.csv evaluated the ROI
    ceiling at the lever diagnostic's re-pin (2.10 um, CI [1.43, 2.50]). The manuscript
    locks the breast pin at the Gate 3 value 1.99 um, CI [1.43, 2.63]
    (results/gate3_pin_breast.csv). For figure consistency we re-evaluate the oracle
    ceiling at the SAME ROI packings/densities at the locked pin. This is a reconciliation
    of a region-run value to the manuscript pin, not a new gate.

Outputs: results/lever/leg2_local.csv, results/lever/lever_roi_pin199.csv.
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
LEVER_SEED = config.MASTER_SEED + 930000

H5 = os.path.join(DATA, "xenium_breast_rep1_cell_feature_matrix.h5")
PARQ = os.path.join(DATA, "xenium_breast_rep1_cells.parquet")
MARKERS = ["ESR1", "PGR", "ERBB2"]
T = 2
KNN = 10

# committed synthetic comparison values (results/lever/lever_mechanism.csv, adj-to-PR, thr>=2)
SYN_CONTROL = 0.0           # sigma=0
SYN_PINNED = 0.0126         # sigma 2.10
SYN_CI_HI = 0.0196          # sigma 2.50

# locked Gate 3 breast pin (results/gate3_pin_breast.csv)
PIN_POINT, PIN_LO, PIN_HI = 1.99, 1.43, 2.63


def local_packing(coords, k=KNN):
    d, _ = cKDTree(coords).query(coords, k=k + 1)
    return k / (np.pi * d[:, k] ** 2) * 1e6


def leg2_local():
    X, genes, coords, total = load_xenium_cells_genes(H5, PARQ)
    genes = list(genes)
    iE, iP, iH = (genes.index(m) for m in MARKERS)
    esr1, pgr, erbb2 = X[:, iE], X[:, iP], X[:, iH]
    dcis = (esr1 >= T) & (erbb2 >= T)                 # ER+/HER2+ cells (denominator)
    triple = dcis & (pgr >= T)                        # apparent triple-positive
    prsrc = (pgr >= T) & ~dcis                        # PR-source: PGR+ non-tumor cells
    tree = cKDTree(coords)
    d_nn, _ = tree.query(coords, k=2)
    mnn = float(np.median(d_nn[:, 1]))
    # DCIS cells adjacent (<= 1.5x median NN, the pin's near-band) to a PR-source cell
    if prsrc.sum():
        d_pr, _ = cKDTree(coords[prsrc]).query(coords, k=1)
        adj = d_pr <= 1.5 * mnn
    else:
        adj = np.zeros(X.shape[0], bool)
    dcis_adj = dcis & adj
    obs_adj_rate = float((triple & dcis_adj).sum()) / float(dcis_adj.sum()) if dcis_adj.sum() else np.nan
    obs_section_rate = float(triple.sum()) / float(dcis.sum())

    # region-local: among DCIS within 50um of a triple-positive cell (reported for context;
    # mildly circular since the region is defined by triple cells)
    nbr = tree.query_ball_point(coords[triple], 50.0)
    roi = np.zeros(X.shape[0], bool)
    if triple.sum():
        roi[np.unique(np.concatenate([np.array(ix, int) for ix in nbr]))] = True
    dcis_roi = dcis & roi
    obs_roi_rate = float((triple & dcis_roi).sum()) / float(dcis_roi.sum()) if dcis_roi.sum() else np.nan

    ratio = obs_adj_rate / SYN_PINNED if SYN_PINNED > 0 else np.inf
    within_oom = (0.1 <= ratio <= 10.0)
    row = {
        "threshold": T, "median_nn_um": mnn,
        "n_dcis": int(dcis.sum()), "n_triple": int(triple.sum()), "n_prsource": int(prsrc.sum()),
        "n_dcis_adj_prsource": int(dcis_adj.sum()),
        "obs_rate_adj_to_prsource": obs_adj_rate,          # apples-to-apples with synthetic adj
        "obs_rate_region_local_50um": obs_roi_rate,
        "obs_rate_section_wide": obs_section_rate,
        "syn_control_sigma0": SYN_CONTROL, "syn_adj_pinned": SYN_PINNED, "syn_adj_ci_hi": SYN_CI_HI,
        "obs_over_syn_pinned_ratio": ratio, "within_one_order_of_magnitude": bool(within_oom),
        "seed": LEVER_SEED,
    }
    print("Leg 2 local-to-local (threshold >=2):")
    print(f"  observed, DCIS adjacent to a PR-source cell : {obs_adj_rate*100:.3f}%  "
          f"({(triple & dcis_adj).sum()} / {dcis_adj.sum()})")
    print(f"  observed, region-local (50um of a triple)   : {obs_roi_rate*100:.3f}%")
    print(f"  observed, section-wide                       : {obs_section_rate*100:.3f}%")
    print(f"  synthetic adjacency: control {SYN_CONTROL*100:.2f}%, pinned {SYN_PINNED*100:.2f}%, "
          f"ci_hi {SYN_CI_HI*100:.2f}%")
    print(f"  observed/synthetic(pinned) ratio = {ratio:.2f}  within 1 order of magnitude: {within_oom}")
    if not within_oom:
        print("  *** FLAG: synthetic under/over-produces by more than one order of magnitude "
              "(Leg 2 would retreat from 'sufficient' to 'contributory') ***")
    return row


def roi_reconcile():
    model, _ = build_realistic_model_from_xenium(n_types=15, seed=config.MASTER_SEED,
                                                 h5_path=H5, cells_parquet=PARQ,
                                                 name="realistic_xenium_breast")
    # ROI definitions and densities from the committed lever_roi.csv
    roi = pd.read_csv(os.path.join(RESULTS, "lever_roi.csv"))
    rows = []
    for j, r in roi.iterrows():
        pack, dens = float(r.packing_cells_per_mm2), float(r.density_tx_per_cell)
        def oa(sig, off):
            return oracle_acc_at(model, pack, dens, sig, LEVER_SEED + 500 + 100 * j + off,
                                 res_cell=ORACLE_RES_CELL, grid_max=ORACLE_GRID_MAX)
        pt, npt = oa(PIN_POINT, 1)
        lo, _ = oa(PIN_LO, 3)     # optimistic
        hi, _ = oa(PIN_HI, 5)     # pessimistic
        rows.append({
            "roi_def": r.roi_def, "packing_cells_per_mm2": pack, "density_tx_per_cell": dens,
            "pin": "gate3_breast_1.99", "sigma_point_um": PIN_POINT,
            "sigma_ci_lo_um": PIN_LO, "sigma_ci_hi_um": PIN_HI,
            "oracle_point": pt, "oracle_opt_ci_lo": lo, "oracle_pess_ci_hi": hi, "naive_point": npt,
            "below_0p90_across_ci": bool(lo < 0.90), "below_0p95_across_ci": bool(lo < 0.95),
            # the 2.10-pin values for reference (from lever_roi.csv)
            "ref_oracle_point_pin210": float(r.oracle_point),
            "ref_oracle_opt_pin210": float(r.oracle_opt_ci_lo),
            "seed": LEVER_SEED + 500 + 100 * j,
        })
        print(f"[{r.roi_def}] pack={pack:.0f} dens={dens:.0f} @ pin 1.99[1.43,2.63]: "
              f"oracle {pt:.3f} (opt {lo:.3f}, pess {hi:.3f}); "
              f"[2.10-pin ref: pt {r.oracle_point:.3f}, opt {r.oracle_opt_ci_lo:.3f}]")
    return pd.DataFrame(rows)


def main():
    os.makedirs(RESULTS, exist_ok=True)
    print("=== (A) Leg 2 local-to-local ===")
    l2 = leg2_local()
    pd.DataFrame([l2]).to_csv(os.path.join(RESULTS, "leg2_local.csv"), index=False)
    print("\n=== (B) ROI-ceiling reconciliation to locked breast pin 1.99 um ===")
    rec = roi_reconcile()
    rec.to_csv(os.path.join(RESULTS, "lever_roi_pin199.csv"), index=False)
    print(f"\nwrote results/lever/leg2_local.csv and lever_roi_pin199.csv")


if __name__ == "__main__":
    main()
