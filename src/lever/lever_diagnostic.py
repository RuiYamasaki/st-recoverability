"""Lever-scout diagnostic: place a candidate Xenium dataset on the answerability frontier.

For each candidate dataset we MEASURE (never assume) three things with the existing
pinning machinery, exactly as the gates did, then read off the oracle ceiling:

  1. packing   - median and 90th-percentile (dense-tail) local cell packing, cells/mm^2,
                 from the segmented cell centroids (cells.parquet).
  2. density   - median transcripts per cell, from the cell-feature matrix.
  3. sigma     - the data-pinned displacement, via the Gate 3 principled non-circular
                 marker-validity test + admissible-set bootstrap pin (gate3_common).

The oracle (Bayes-optimal, method-independent) assignment accuracy is then evaluated at
the dataset's OWN median and OWN dense-tail (p90) packing/density, at the pinned sigma
point and at both ends of its bootstrap CI, on the converged grid. The 'ceiling position'
is the oracle accuracy in the dataset's dense regime: if it is below ~0.9, transcript-to-
cell assignment is provably unrecoverable there, and any conclusion that needs correct
assignment in that regime sits ABOVE the ceiling.

This script is append-only: it reuses Gate 0-3 code unchanged and writes only to
results/lever/. It does not import or rewrite any prior-gate result file.

Datasets (standard 10x Xenium output bundle: cell_feature_matrix.h5 + cells.parquet):
  - breast : Xenium FFPE Human Breast Cancer Rep1 (10x release 1.0.1) = the Janesick et
             al. 2023 Nat Commun dataset. Already cached under data/.
  - crc    : Xenium V1 Human Colorectal Cancer Add-on FFPE (10x Onboard Analysis 2.0.0),
             a 10x public dense CRC tumor dataset. Cached under data/lever/.

Output: results/lever/lever_ceiling.csv. Seeds: LEVER_SEED + offsets, recorded.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config  # noqa: E402
from expression import build_realistic_model_from_xenium  # noqa: E402
from gate3_common import (marker_validity, admissible_owner, candidate_markers,  # noqa: E402
                          pin_admissible, N_PERM, P_ADMIT)
from gate2_pin import (oracle_acc_at, ORACLE_RES_CELL, ORACLE_GRID_MAX,  # noqa: E402
                       DENSE_PACKING, SPARSE_PACKING)

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS = os.path.join(ROOT, "results", "lever")
DATA = os.path.join(ROOT, "data")
LEVER_SEED = config.MASTER_SEED + 900000

CEILING_BAR = 0.9   # below this in the dense regime => assignment unrecoverable => above the ceiling

DATASETS = {
    "breast": {
        "h5": os.path.join(DATA, "xenium_breast_rep1_cell_feature_matrix.h5"),
        "parquet": os.path.join(DATA, "xenium_breast_rep1_cells.parquet"),
        "label": "Xenium FFPE Human Breast Cancer Rep1 (10x release 1.0.1; Janesick 2023)",
    },
    "crc": {
        "h5": os.path.join(DATA, "lever", "xenium_crc_addon_cell_feature_matrix.h5"),
        "parquet": os.path.join(DATA, "lever", "xenium_crc_addon_cells.parquet"),
        "label": "Xenium V1 Human Colorectal Cancer Add-on FFPE (10x Onboard Analysis 2.0.0)",
    },
}


def diagnose(tag, h5, parquet, label, seed):
    model, real = build_realistic_model_from_xenium(
        n_types=15, seed=config.MASTER_SEED, h5_path=h5, cells_parquet=parquet,
        name=f"realistic_xenium_{tag}")
    abundance = model.mean_expr.sum(axis=1)
    abundance = abundance / abundance.mean()

    own_med = real["packing_median_cells_per_mm2"]
    own_p90 = real["packing_p90_cells_per_mm2"]
    dens = real["density_median_tx_per_cell"]
    print(f"\n[{tag}] {label}")
    print(f"[{tag}] K={model.n_types} G={model.n_genes} median_nn={real['median_nn_um']:.2f}um "
          f"pack_median={own_med:.0f} pack_p90={own_p90:.0f} dens_median={dens:.0f}")

    # --- principled non-circular marker-validity test + admissible-set sigma pin ---
    owner_cand, candidates = candidate_markers(model)
    vdf = marker_validity(real, owner_cand, candidates, seed)
    n_adm = int(vdf.admissible.sum())
    owner_adm = admissible_owner(model, vdf)
    print(f"[{tag}] candidate markers (excl>=0.5): {len(candidates)}; "
          f"admissible (p<{P_ADMIT}, {N_PERM} perms): {n_adm}")
    if n_adm < 2:
        print(f"[{tag}] WARNING: <2 admissible markers; sigma pin not reliable.")
    pin = pin_admissible(model, real, abundance, owner_adm, seed + 100)
    sig_pt = pin["sigma_point_um"]
    sig_lo = pin["sigma_ci_lo_um"]   # low sigma  -> optimistic (highest) oracle accuracy
    sig_hi = pin["sigma_ci_hi_um"]   # high sigma -> pessimistic (lowest) oracle accuracy
    print(f"[{tag}] data-pinned sigma = {sig_pt:.2f} um  CI [{sig_lo:.2f}, {sig_hi:.2f}] "
          f"(status {pin['sigma_fit_status']})")

    # --- oracle ceiling at the dataset's OWN packing/density (converged grid) ---
    def oa(pack, sig, sd):
        return oracle_acc_at(model, pack, dens, sig, sd,
                             res_cell=ORACLE_RES_CELL, grid_max=ORACLE_GRID_MAX)

    # own median packing
    om_pt, nm_pt = oa(own_med, sig_pt, seed + 201)
    om_lo, _ = oa(own_med, sig_lo, seed + 203)
    om_hi, _ = oa(own_med, sig_hi, seed + 205)
    # own dense-tail (p90) packing
    op_pt, np_pt = oa(own_p90, sig_pt, seed + 211)
    op_lo, _ = oa(own_p90, sig_lo, seed + 213)
    op_hi, _ = oa(own_p90, sig_hi, seed + 215)

    above_med = om_lo < CEILING_BAR   # optimistic end still below the bar at own median packing
    above_p90 = op_lo < CEILING_BAR   # optimistic end still below the bar at own dense tail
    print(f"[{tag}] oracle @ own MEDIAN pack ({own_med:.0f}): point={om_pt:.3f} "
          f"CI[{om_hi:.3f},{om_lo:.3f}] (naive {nm_pt:.3f})")
    print(f"[{tag}] oracle @ own P90    pack ({own_p90:.0f}): point={op_pt:.3f} "
          f"CI[{op_hi:.3f},{op_lo:.3f}] (naive {np_pt:.3f})  -> above {CEILING_BAR} ceiling: {above_p90}")
    print(f"[{tag}] (ref) fixed-regime dense({DENSE_PACKING:.0f}) point={pin['dense_oracle_point']:.3f}; "
          f"sparse({SPARSE_PACKING:.0f}) point={pin['sparse_oracle_point']:.3f}")

    return {
        "dataset": tag, "label": label,
        "n_types": model.n_types, "n_genes": model.n_genes,
        "median_nn_um": real["median_nn_um"],
        "pack_median_cells_per_mm2": own_med,
        "pack_p90_cells_per_mm2": own_p90,
        "density_median_tx_per_cell": dens,
        "n_candidate_markers": len(candidates), "n_admissible_markers": n_adm,
        "sigma_point_um": sig_pt, "sigma_ci_lo_um": sig_lo, "sigma_ci_hi_um": sig_hi,
        "sigma_fit_status": pin["sigma_fit_status"],
        "oracle_own_median_point": om_pt, "oracle_own_median_ci_lo": om_lo,
        "oracle_own_median_ci_hi": om_hi, "naive_own_median_point": nm_pt,
        "oracle_own_p90_point": op_pt, "oracle_own_p90_ci_lo": op_lo,
        "oracle_own_p90_ci_hi": op_hi, "naive_own_p90_point": np_pt,
        "ceiling_bar": CEILING_BAR,
        "above_ceiling_own_median": bool(above_med),
        "above_ceiling_own_p90": bool(above_p90),
        "ref_fixed_dense_oracle_point": pin["dense_oracle_point"],
        "ref_fixed_sparse_oracle_point": pin["sparse_oracle_point"],
        "ref_fixed_dense_packing": DENSE_PACKING, "ref_fixed_sparse_packing": SPARSE_PACKING,
        "seed": seed,
    }


def main():
    os.makedirs(RESULTS, exist_ok=True)
    only = sys.argv[1:] if len(sys.argv) > 1 else list(DATASETS.keys())
    rows = []
    for i, tag in enumerate(only):
        d = DATASETS[tag]
        rows.append(diagnose(tag, d["h5"], d["parquet"], d["label"], LEVER_SEED + 1000 * i))
    df = pd.DataFrame(rows)
    out = os.path.join(RESULTS, "lever_ceiling.csv")
    # append-safe: merge with any existing rows for datasets not re-run this call
    if os.path.exists(out):
        old = pd.read_csv(out)
        old = old[~old.dataset.isin(df.dataset)]
        df = pd.concat([old, df], ignore_index=True)
    df.to_csv(out, index=False)
    print(f"\nwrote {out} ({len(df)} datasets)")
    return df


if __name__ == "__main__":
    main()
