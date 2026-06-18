"""Shared configuration: seeds, the fixed cell-type / gene model, physical scale.

Everything random in this project derives from MASTER_SEED. No un-seeded
randomness anywhere. Values here are *modelling choices* fixed before any oracle
output is seen; they are NOT tuned to clear a Gate-0 threshold.
"""
from __future__ import annotations

import numpy as np

# --------------------------------------------------------------------------
# Master seed. Every per-config seed is derived deterministically from this.
# --------------------------------------------------------------------------
MASTER_SEED = 20260618  # the project date, recorded so the whole study is reproducible

# --------------------------------------------------------------------------
# Fixed cell-type / gene model (shared by generator and oracle).
# K cell types, G genes. Each type has a normalised expression composition
# p[t, :] (sums to 1 over genes). The composition is a shared background plus a
# type-specific marker block, giving partially-overlapping profiles (realistic:
# types share housekeeping genes and differ on a handful of markers). This single
# model is held fixed across the entire sweep; it is not a swept axis.
# --------------------------------------------------------------------------
N_TYPES = 4                                   # K
N_GENES = 20                                  # G
TYPE_PROPORTIONS = np.array([0.40, 0.30, 0.20, 0.10])  # fixed, sums to 1
MARKERS_PER_TYPE = 4                          # genes boosted per type (non-overlapping blocks)
MARKER_MASS_FRACTION = 0.55                   # fraction of a type's expression mass on its markers
GENE_MODEL_SEED = MASTER_SEED + 7             # seed for any randomness in building the model


def build_type_gene_composition() -> np.ndarray:
    """Return p[K, G]: each row a normalised expression composition for a type.

    Construction (fully determined, no free randomness beyond GENE_MODEL_SEED for
    the small background jitter):
      - A shared background over all G genes (mildly non-uniform, so even the
        non-marker channel carries a little information).
      - Each type t owns a contiguous block of MARKERS_PER_TYPE marker genes
        (blocks are disjoint: type t -> genes [t*M : (t+1)*M)).
      - MARKER_MASS_FRACTION of the row mass is placed on the type's markers,
        the remainder follows the shared background.
    """
    rng = np.random.default_rng(GENE_MODEL_SEED)
    M = MARKERS_PER_TYPE
    assert N_TYPES * M <= N_GENES, "marker blocks must fit in the gene panel"

    background = rng.uniform(0.5, 1.5, size=N_GENES)
    background = background / background.sum()  # normalised shared background

    p = np.zeros((N_TYPES, N_GENES), dtype=float)
    for t in range(N_TYPES):
        marker_idx = np.arange(t * M, (t + 1) * M)
        marker = np.zeros(N_GENES)
        marker[marker_idx] = 1.0 / M                      # uniform over this type's markers
        p[t] = MARKER_MASS_FRACTION * marker + (1.0 - MARKER_MASS_FRACTION) * background
        p[t] = p[t] / p[t].sum()
    return p


# Precompute once and reuse everywhere (generator emits genes from this; oracle
# uses log of this as the gene-evidence term).
TYPE_GENE_COMPOSITION = build_type_gene_composition()  # shape (N_TYPES, N_GENES)


def config_seed(packing_idx: int, sigma_idx: int) -> int:
    """Deterministic per-(packing, sigma) seed for the spatial field + transcripts."""
    return MASTER_SEED + 1000 * (packing_idx + 1) + 17 * (sigma_idx + 1)
