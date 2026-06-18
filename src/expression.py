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
    dispersion: np.ndarray = field(default=None)  # (G,) per-gene NB dispersion (Gate 2 Task 3)

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


def load_xenium_cells_genes(h5_path: str = None, cells_parquet: str = None):
    """Load the Xenium cell-feature matrix (Gene Expression features only) as a dense
    cells-by-genes float32 array, aligned to the cells.parquet centroids."""
    import h5py
    import pandas as pd
    import scipy.sparse as sp
    if h5_path is None:
        h5_path = os.path.join(DATA, "xenium_breast_rep1_cell_feature_matrix.h5")
    if cells_parquet is None:
        cells_parquet = os.path.join(DATA, "xenium_breast_rep1_cells.parquet")
    with h5py.File(h5_path, "r") as f:
        m = f["matrix"]
        shape = tuple(m["shape"][:])
        M = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=shape)
        ftype = np.array([x.decode() for x in m["features"]["feature_type"][:]])
        names = np.array([x.decode() for x in m["features"]["name"][:]])
    gmask = ftype == "Gene Expression"
    X = M[gmask, :].T.toarray().astype(np.float32)        # cells x genes
    genes = list(names[gmask])
    cells = pd.read_parquet(cells_parquet)
    coords = cells[["x_centroid", "y_centroid"]].to_numpy(np.float64)
    total = cells["transcript_counts"].to_numpy(np.float64)
    return X, genes, coords, total


def build_realistic_model_from_xenium(n_types: int = 15, excl_threshold: float = 0.7,
                                      seed: int = None, min_cells: int = 50,
                                      h5_path: str = None, cells_parquet: str = None,
                                      name: str = "realistic_xenium_kmeans"):
    """Cluster the Xenium cell-by-gene matrix into types (MiniBatchKMeans on
    median-normalised log1p counts), build per-type mean expression profiles, exclusive
    markers, proportions, and a per-gene negative-binomial dispersion estimated from the
    within-type count variance. Clusters smaller than min_cells are dropped (their cells
    removed). h5_path/cells_parquet select the dataset (default: breast). Returns
    (model, real) where real carries the per-cell counts, centroids, types and section
    geometry for the leakage measurement."""
    from sklearn.cluster import MiniBatchKMeans
    if seed is None:
        seed = config.MASTER_SEED
    X, genes, coords, total = load_xenium_cells_genes(h5_path, cells_parquet)
    G = X.shape[1]
    # normalise to median total, log1p, cluster
    sf = total.copy(); sf[sf == 0] = 1.0
    Xn = np.log1p(X / sf[:, None] * np.median(total))
    km = MiniBatchKMeans(n_clusters=n_types, random_state=seed, n_init=10, batch_size=4096)
    raw_labels = km.fit_predict(Xn).astype(np.int32)

    # drop clusters below min_cells, remap labels to 0..K-1 (cells of dropped clusters removed)
    keep_clusters = [t for t in range(n_types) if (raw_labels == t).sum() >= min_cells]
    remap = {t: i for i, t in enumerate(keep_clusters)}
    cell_keep = np.array([lbl in remap for lbl in raw_labels])
    X, coords, total = X[cell_keep], coords[cell_keep], total[cell_keep]
    labels = np.array([remap[lbl] for lbl in raw_labels[cell_keep]], dtype=np.int32)
    n_types = len(keep_clusters)

    mean_expr = np.vstack([X[labels == t].mean(axis=0) for t in range(n_types)])
    proportions = np.array([(labels == t).mean() for t in range(n_types)], dtype=float)
    row = mean_expr.sum(axis=1, keepdims=True)
    composition = np.divide(mean_expr, row, out=np.zeros_like(mean_expr), where=row > 0)
    owner = _exclusivity(composition, excl_threshold)
    dispersion = _estimate_dispersion(X, labels, n_types)
    model = ExpressionModel(
        name=name, n_types=n_types, n_genes=G,
        type_names=[f"xen{t}" for t in range(n_types)],
        proportions=proportions, composition=composition, gene_names=genes,
        excl_threshold=excl_threshold, excl_owner=owner, mean_expr=mean_expr,
    )
    model.dispersion = dispersion
    from scipy.spatial import cKDTree
    d, _ = cKDTree(coords).query(coords, k=11)
    median_nn = float(np.median(d[:, 1]))
    local_pack = 10.0 / (np.pi * d[:, 10] ** 2) * 1e6
    real = {
        "X": X, "coords": coords, "types": labels, "total": total, "genes": genes,
        "median_nn_um": median_nn,
        "packing_median_cells_per_mm2": float(np.median(local_pack)),
        "packing_p90_cells_per_mm2": float(np.percentile(local_pack, 90)),
        "density_median_tx_per_cell": float(np.median(total)),
    }
    return model, real


def _estimate_dispersion(X, labels, n_types, clip=(1e-3, 50.0)):
    """Per-gene negative-binomial dispersion phi (Var = mu + phi*mu^2) by within-type
    method of moments, pooled as the median over types."""
    G = X.shape[1]
    phis = np.full((n_types, G), np.nan)
    for t in range(n_types):
        Xt = X[labels == t]
        if Xt.shape[0] < 10:
            continue
        m = Xt.mean(axis=0)
        v = Xt.var(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            phi = (v - m) / (m ** 2)
        phis[t] = np.where(m > 0, phi, np.nan)
    phi_g = np.nanmedian(phis, axis=0)
    phi_g = np.clip(np.nan_to_num(phi_g, nan=clip[0]), *clip)
    return phi_g


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
