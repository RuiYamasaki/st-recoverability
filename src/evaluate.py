"""One evaluation of the oracle + naive baseline + metrics at a single grid point,
given a prebuilt field and oracle maps. Shared by sweep.py and anchor.py."""
from __future__ import annotations

import os
import sys

import numpy as np
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402
from generator import Field, generate_transcripts  # noqa: E402
from oracle import oracle_assign, naive_assign  # noqa: E402
from metrics import assignment_accuracy, profile_ari  # noqa: E402


def _same_type_error_frac(field: Field, assigned, tx) -> float:
    m = tx.interior & (assigned != tx.true_cell)
    if not m.any():
        return float("nan")
    atype = field.types[np.clip(assigned, 0, field.n_cells - 1)]
    return float((atype[m] == tx.true_type[m]).mean())


def eval_point(field: Field, Dmax, argcell, mean_tx_per_cell: float, tx_seed: int,
               displacement: str = "gaussian", disp_epsilon: float = 0.0) -> dict:
    tx = generate_transcripts(field, mean_tx_per_cell, tx_seed,
                              displacement=displacement, disp_epsilon=disp_epsilon)
    oa = oracle_assign(field, Dmax, argcell, tx.obs_xy, tx.gene, disp_epsilon=disp_epsilon)
    na = naive_assign(field, tx.obs_xy)
    acc_o, n_int = assignment_accuracy(oa, tx.true_cell, tx.interior)
    acc_n, _ = assignment_accuracy(na, tx.true_cell, tx.interior)
    # assignment-partition ARI: agreement of the recovered transcript->cell labelling
    # with the true one (this is the "ARI vs true cell labels" sense; it tracks
    # assignment quality, distinct from the expression-profile ARI below).
    if n_int > 0:
        m = tx.interior
        assign_ari_o = float(adjusted_rand_score(tx.true_cell[m], oa[m]))
        assign_ari_n = float(adjusted_rand_score(tx.true_cell[m], na[m]))
    else:
        assign_ari_o = assign_ari_n = float("nan")
    ari = profile_ari(field.model, field.n_cells, field.interior_cell, oa, tx.true_cell,
                      tx.gene, field.types, tx_seed + 101)
    # naive profile ARI, to contrast cell-typing recovery (oracle errors are
    # type-preserving; naive errors are type-random)
    ari_n = profile_ari(field.model, field.n_cells, field.interior_cell, na, tx.true_cell,
                        tx.gene, field.types, tx_seed + 101)
    # realized summary statistics (for the realism 15% match) over interior cells
    true_counts = np.bincount(tx.true_cell[tx.interior], minlength=field.n_cells)
    interior_idx = np.where(field.interior_cell)[0]
    med_tx = float(np.median(true_counts[interior_idx])) if interior_idx.size else float("nan")
    realized_packing = field.n_cells / ((field.L_um ** 2) * 1e-6)
    return {
        "n_cells": field.n_cells,
        "n_interior_cells": int(field.interior_cell.sum()),
        "grid": field.grid,
        "dx_um": field.dx_um,
        "r_mean_um": field.r_mean_um,
        "L_um": field.L_um,
        "n_interior_transcripts": n_int,
        "oracle_acc": acc_o,
        "naive_acc": acc_n,
        "oracle_minus_naive": acc_o - acc_n,
        "oracle_ge_naive": bool(acc_o >= acc_n - 1e-9),
        "assign_ari_vs_truecells_oracle": assign_ari_o,
        "assign_ari_vs_truecells_naive": assign_ari_n,
        "same_type_error_frac_oracle": _same_type_error_frac(field, oa, tx),
        "same_type_error_frac_naive": _same_type_error_frac(field, na, tx),
        "profile_ari_vs_truetype": ari["profile_ari_vs_truetype"],
        "profile_ari_vs_trueclust": ari["profile_ari_vs_trueclust"],
        "perfect_ari_vs_truetype": ari["perfect_ari_vs_truetype"],
        "naive_profile_ari_vs_truetype": ari_n["profile_ari_vs_truetype"],
        "n_cells_clustered": ari["n_cells_clustered"],
        "realized_median_tx_per_cell": med_tx,
        "realized_packing_cells_per_mm2": realized_packing,
    }
