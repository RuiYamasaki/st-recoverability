"""Expression models for the generator: the Gate 0 disjoint-marker model (default,
reproduces Gate 0 exactly) and a realistic overlapping model derived from a real
reference (MERFISH cell-class mean profiles).

An ExpressionModel bundles everything the generator/oracle/metrics need so the model
can be swapped without rewriting them: number of types/genes, type proportions, the
per-type normalised expression composition p[t, g] (rows sum to 1), and the set of
exclusive-marker genes per type (used by the marker-leakage statistic in Gate 1).
"""
from __future__ import annotations

import os
import re
import sys
from dataclasses import dataclass, field

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")


@dataclass
class ExpressionModel:
    name: str
    n_types: int
    n_genes: int
    type_names: list
    proportions: np.ndarray          # (K,) sums to 1
    composition: np.ndarray          # (K, G) rows sum to 1
    gene_names: list
    excl_threshold: float            # exclusivity cutoff used to pick marker genes
    excl_owner: np.ndarray           # (G,) owner type per gene, -1 if not exclusive
    mean_expr: np.ndarray = field(default=None)  # (K, G) raw per-type means (abundance)

    def exclusive_genes(self, t: int) -> np.ndarray:
        return np.where(self.excl_owner == t)[0]

    @property
    def n_exclusive(self) -> int:
        return int((self.excl_owner >= 0).sum())


def _exclusivity(composition: np.ndarray, threshold: float) -> np.ndarray:
    """Owner type per gene if one type holds >= threshold of the gene's cross-type
    intensity, else -1."""
    K, G = composition.shape
    owner = np.full(G, -1, dtype=int)
    col = composition.sum(axis=0)
    for g in range(G):
        if col[g] <= 0:
            continue
        t = int(np.argmax(composition[:, g]))
        if composition[t, g] / col[g] >= threshold:
            owner[g] = t
    return owner


def build_disjoint_model() -> ExpressionModel:
    """The Gate 0 model, wrapped. Used as the default so Gate 0 reproduces exactly."""
    comp = config.TYPE_GENE_COMPOSITION
    owner = _exclusivity(comp, threshold=0.7)
    return ExpressionModel(
        name="disjoint_marker_gate0",
        n_types=config.N_TYPES, n_genes=config.N_GENES,
        type_names=[f"type{t}" for t in range(config.N_TYPES)],
        proportions=config.TYPE_PROPORTIONS.copy(),
        composition=comp.copy(),
        gene_names=[f"g{g}" for g in range(config.N_GENES)],
        excl_threshold=0.7, excl_owner=owner,
    )


def _merge_class(name: str) -> str:
    """Merge sub-numbered cell classes (e.g. 'Endothelial 3' -> 'Endothelial')."""
    return re.sub(r"\s*\d+$", "", str(name)).strip()


def build_realistic_model_from_merfish(
    h5ad_path: str = None, excl_threshold: float = 0.7, min_cells_per_type: int = 60,
    merge: bool = False
) -> ExpressionModel:
    """Per-type mean expression profiles from MERFISH (Moffitt 2018) cell classes.

    Real cell types -> genuinely overlapping profiles over the 155-gene panel. By
    default the dataset's native cell_class labels are used (merge=False), which
    include confusable subtypes (e.g. Endothelial 1 vs 3); these make cell typing
    genuinely imperfect so the profile ARI de-saturates below 1.0. merge=True collapses
    sub-numbered classes into broad classes (more separable). Types with fewer than
    min_cells_per_type cells are dropped. Proportions are the real class frequencies.
    """
    import anndata as ad
    if h5ad_path is None:
        h5ad_path = os.path.join(DATA, "merfish_moffitt.h5ad")
    a = ad.read_h5ad(h5ad_path)
    X = np.asarray(a.X, dtype=np.float64)
    genes = list(map(str, a.var_names))
    raw = a.obs["cell_class"].astype(str)
    classes = np.array([_merge_class(c) for c in raw]) if merge else raw.to_numpy()

    names, counts = np.unique(classes, return_counts=True)
    keep = sorted([n for n, c in zip(names, counts) if c >= min_cells_per_type])
    means, props = [], []
    for n in keep:
        m = classes == n
        means.append(X[m].mean(axis=0))
        props.append(int(m.sum()))
    mean_expr = np.vstack(means)                  # (K, G) raw means
    proportions = np.array(props, dtype=float)
    proportions = proportions / proportions.sum()
    row = mean_expr.sum(axis=1, keepdims=True)
    composition = np.divide(mean_expr, row, out=np.zeros_like(mean_expr), where=row > 0)
    owner = _exclusivity(composition, excl_threshold)
    return ExpressionModel(
        name="realistic_merfish_cellclass" + ("_merged" if merge else "_native"),
        n_types=len(keep), n_genes=len(genes), type_names=keep,
        proportions=proportions, composition=composition, gene_names=genes,
        excl_threshold=excl_threshold, excl_owner=owner, mean_expr=mean_expr,
    )


def overlap_metrics(model: ExpressionModel) -> dict:
    """Quantify how overlapping the type profiles are."""
    C = model.composition
    K = model.n_types
    # mean pairwise cosine similarity between type composition vectors
    norm = C / (np.linalg.norm(C, axis=1, keepdims=True) + 1e-12)
    cos = norm @ norm.T
    iu = np.triu_indices(K, k=1)
    mean_cos = float(cos[iu].mean())
    max_cos = float(cos[iu].max())
    # fraction of expression mass (proportion-weighted) on exclusive vs shared genes
    excl_mask = model.excl_owner >= 0
    mass_excl = float(np.sum(model.proportions[:, None] * C[:, excl_mask]))
    total_mass = float(np.sum(model.proportions[:, None] * C))
    frac_excl = mass_excl / total_mass
    return {
        "mean_pairwise_cosine": mean_cos,
        "max_pairwise_cosine": max_cos,
        "n_exclusive_markers": model.n_exclusive,
        "frac_mass_exclusive_genes": frac_excl,
        "frac_mass_shared_genes": 1.0 - frac_excl,
    }


if __name__ == "__main__":
    print("=== disjoint (Gate 0) model ===")
    d = build_disjoint_model()
    print(f"K={d.n_types} G={d.n_genes} n_excl={d.n_exclusive}")
    print(overlap_metrics(d))
    print("\n=== realistic MERFISH model ===")
    for thr in (0.7, 0.8, 0.9):
        r = build_realistic_model_from_merfish(excl_threshold=thr)
        om = overlap_metrics(r)
        print(f"thr={thr}: K={r.n_types} G={r.n_genes} types={r.type_names}")
        print(f"   proportions={np.round(r.proportions,3)}")
        print(f"   {om}")
