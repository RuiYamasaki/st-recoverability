"""Gate-0 metrics.

  - assignment_accuracy: fraction of (interior) transcripts assigned to their true
    source cell. The headline oracle ceiling and the naive baseline both use this.
  - profile_ari: adjusted Rand index of recovered cell-level expression profiles vs
    the truth. Per interior cell we build a gene-count vector from the transcripts
    ASSIGNED to it, cluster cells into K groups (deterministic KMeans), and compare:
      * vs_truetype : ARI of the recovered-profile clustering against true cell type
                      (does downstream cell typing survive the assignment?).
      * vs_trueclust: ARI of the recovered-profile clustering against the clustering
                      of the TRUE profiles (separates assignment damage from the
                      intrinsic difficulty of clustering at this transcript depth).
    perfect_ari is the same clustering on perfectly-assigned (true) profiles vs true
    type, i.e. the best ARI achievable at this density regardless of assignment.
"""
from __future__ import annotations

import os
import sys

import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402


def assignment_accuracy(assigned: np.ndarray, true_cell: np.ndarray,
                        interior: np.ndarray) -> tuple[float, int]:
    m = interior
    n = int(m.sum())
    if n == 0:
        return float("nan"), 0
    return float((assigned[m] == true_cell[m]).mean()), n


def _cell_profiles(cell_ids: np.ndarray, gene: np.ndarray, n_cells: int, n_genes: int) -> np.ndarray:
    """Counts[n_cells, G] of genes grouped by cell id (ids outside range ignored)."""
    counts = np.zeros((n_cells, n_genes), dtype=np.float64)
    ok = (cell_ids >= 0) & (cell_ids < n_cells)
    np.add.at(counts, (cell_ids[ok], gene[ok]), 1.0)
    return counts


def _cluster(profiles: np.ndarray, seed: int, n_types: int) -> np.ndarray:
    """Row-normalise to composition and KMeans into K clusters (deterministic)."""
    tot = profiles.sum(axis=1, keepdims=True)
    comp = np.divide(profiles, tot, out=np.zeros_like(profiles), where=tot > 0)
    km = KMeans(n_clusters=n_types, n_init=10, random_state=seed)
    return km.fit_predict(comp)


def profile_ari(model, field_n_cells: int, interior_cell: np.ndarray,
                assigned: np.ndarray, true_cell: np.ndarray, gene: np.ndarray,
                true_types: np.ndarray, seed: int) -> dict:
    """Compute the oracle profile-recovery ARIs over interior cells."""
    K, G = model.n_types, model.n_genes
    nan_out = {"profile_ari_vs_truetype": float("nan"),
               "profile_ari_vs_trueclust": float("nan"),
               "perfect_ari_vs_truetype": float("nan"),
               "n_cells_clustered": 0}
    cells = np.where(interior_cell)[0]
    if cells.size < K + 1:
        nan_out["n_cells_clustered"] = int(cells.size)
        return nan_out

    rec_counts = _cell_profiles(assigned, gene, field_n_cells, G)[cells]
    true_counts = _cell_profiles(true_cell, gene, field_n_cells, G)[cells]
    # drop cells with no transcripts in either representation (undefined profile)
    keep = (rec_counts.sum(1) > 0) & (true_counts.sum(1) > 0)
    cells_k = cells[keep]
    rec_counts, true_counts = rec_counts[keep], true_counts[keep]
    if cells_k.size < K + 1:
        nan_out["n_cells_clustered"] = int(cells_k.size)
        return nan_out

    rec_lab = _cluster(rec_counts, seed, K)
    true_lab = _cluster(true_counts, seed, K)
    truth = true_types[cells_k]
    return {
        "profile_ari_vs_truetype": float(adjusted_rand_score(truth, rec_lab)),
        "profile_ari_vs_trueclust": float(adjusted_rand_score(true_lab, rec_lab)),
        "perfect_ari_vs_truetype": float(adjusted_rand_score(truth, true_lab)),
        "n_cells_clustered": int(cells_k.size),
    }
