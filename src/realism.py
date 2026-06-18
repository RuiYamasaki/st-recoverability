"""Realism anchor: stream small fields of view from public datasets and compute
the summary statistics that anchor the generator's parameter axes.

Per field of view (FOV) we compute:
  - median transcripts per cell        (-> generator 'density' axis: mean tx/cell)
  - cell-packing density (cells/mm^2)  (-> generator 'packing' axis)
  - median nearest-neighbour cell distance (um)  (cross-check on packing)

Datasets (exact IDs and access routes recorded in results/realism_meta.json):
  - Xenium: 10x Genomics "Xenium FFPE Human Breast Cancer, Replicate 1"
            (release 1.0.1), standalone cells.parquet. Direct HTTP, no auth.
  - MERFISH: Moffitt et al. 2018 mouse hypothalamic preoptic region, the
             squidpy-packaged AnnData (MERFISH_0.24.h5ad) hosted on figshare.
             Direct HTTP, no auth. (Vizgen's own showcase CSVs are gsutil-gated
             and return 403 on anonymous HTTP, so the figshare mirror is used and
             that access route is recorded honestly.)

Raw downloads are cached under data/ (gitignored). Only the small derived
results/realism.csv and results/realism_meta.json are committed.
"""
from __future__ import annotations

import json
import os
import urllib.request
from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")
os.makedirs(DATA, exist_ok=True)
os.makedirs(RESULTS, exist_ok=True)

FOV_SIDE_UM = 200.0          # side length of each square field of view (micrometres)
N_FOV_PER_DATASET = 3        # number of FOVs sampled per dataset
FOV_SELECT_SEED = 20260618   # deterministic FOV placement

XENIUM_CELLS_URL = (
    "https://cf.10xgenomics.com/samples/xenium/1.0.1/"
    "Xenium_FFPE_Human_Breast_Cancer_Rep1/"
    "Xenium_FFPE_Human_Breast_Cancer_Rep1_cells.parquet"
)
XENIUM_ID = "Xenium_FFPE_Human_Breast_Cancer_Rep1 (10x Genomics release 1.0.1)"

MERFISH_H5AD_URL = "https://ndownloader.figshare.com/files/40038538"  # MERFISH_0.24.h5ad
MERFISH_ID = "Moffitt2018_mouse_hypothalamic_preoptic (squidpy MERFISH_0.24.h5ad, figshare)"


def _download(url: str, dest: str) -> str:
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"  cached: {dest} ({os.path.getsize(dest):,} bytes)")
        return dest
    print(f"  downloading {url}\n    -> {dest}")
    req = urllib.request.Request(url, headers={"User-Agent": "gate0-realism/1.0"})
    with urllib.request.urlopen(req, timeout=300) as r, open(dest, "wb") as f:
        f.write(r.read())
    print(f"  done: {os.path.getsize(dest):,} bytes")
    return dest


@dataclass
class FovStat:
    dataset: str
    fov_id: str
    n_cells: int
    median_tx_per_cell: float
    mean_tx_per_cell: float
    packing_cells_per_mm2: float
    median_nn_distance_um: float
    fov_side_um: float
    x0_um: float
    y0_um: float


def _fov_windows(x: np.ndarray, y: np.ndarray, side: float, n: int, seed: int):
    """Deterministically pick n square windows of given side that contain cells.

    Tile the bounding box into a coarse grid of candidate windows, rank by cell
    count, and take the densest n (seed only breaks exact ties, deterministically).
    """
    rng = np.random.default_rng(seed)
    xmin, xmax, ymin, ymax = x.min(), x.max(), y.min(), y.max()
    xs = np.arange(xmin, xmax - side, side)
    ys = np.arange(ymin, ymax - side, side)
    cands = []
    for x0 in xs:
        for y0 in ys:
            m = (x >= x0) & (x < x0 + side) & (y >= y0) & (y < y0 + side)
            cands.append((int(m.sum()), float(x0), float(y0)))
    cands = [c for c in cands if c[0] >= 30]  # need a non-trivial number of cells
    # rank by count desc, tiny deterministic jitter to break ties
    order = sorted(range(len(cands)),
                   key=lambda i: (-cands[i][0], rng.random()))
    picked = [cands[i] for i in order[:n]]
    return picked


def _stats_for_windows(dataset, x, y, tx_per_cell, windows):
    """Compute per-FOV stats. NN distance uses the full section's KDTree so that
    cells at a window edge are not penalised by the crop."""
    tree = cKDTree(np.column_stack([x, y]))
    d, _ = tree.query(np.column_stack([x, y]), k=2)
    nn_all = d[:, 1]  # nearest other cell, micrometres
    out = []
    for k, (cnt, x0, y0) in enumerate(windows):
        m = (x >= x0) & (x < x0 + FOV_SIDE_UM) & (y >= y0) & (y < y0 + FOV_SIDE_UM)
        n_cells = int(m.sum())
        area_mm2 = (FOV_SIDE_UM ** 2) * 1e-6
        out.append(FovStat(
            dataset=dataset,
            fov_id=f"{dataset}_fov{k}",
            n_cells=n_cells,
            median_tx_per_cell=float(np.median(tx_per_cell[m])),
            mean_tx_per_cell=float(np.mean(tx_per_cell[m])),
            packing_cells_per_mm2=float(n_cells / area_mm2),
            median_nn_distance_um=float(np.median(nn_all[m])),
            fov_side_um=FOV_SIDE_UM,
            x0_um=x0, y0_um=y0,
        ))
    return out


def process_xenium():
    print("[Xenium]")
    dest = _download(XENIUM_CELLS_URL, os.path.join(DATA, "xenium_breast_rep1_cells.parquet"))
    df = pd.read_parquet(dest)
    print(f"  loaded {len(df):,} cells; columns: {list(df.columns)}")
    x = df["x_centroid"].to_numpy(float)
    y = df["y_centroid"].to_numpy(float)
    tx = df["transcript_counts"].to_numpy(float)  # panel transcripts assigned to the cell
    wins = _fov_windows(x, y, FOV_SIDE_UM, N_FOV_PER_DATASET, FOV_SELECT_SEED)
    return _stats_for_windows("xenium_breast", x, y, tx, wins)


def process_merfish():
    print("[MERFISH]")
    import anndata as ad
    dest = _download(MERFISH_H5AD_URL, os.path.join(DATA, "merfish_moffitt.h5ad"))
    adata = ad.read_h5ad(dest)
    obs = adata.obs
    print(f"  loaded {adata.n_obs:,} cells x {adata.n_vars} genes; obs columns: {list(obs.columns)}")
    # centroid columns: try obsm['spatial'] first, then common obs names
    if "spatial" in adata.obsm:
        coords = np.asarray(adata.obsm["spatial"], float)
        x, y = coords[:, 0], coords[:, 1]
    else:
        cx = next(c for c in obs.columns if c.lower() in ("center_x", "centroid_x", "x"))
        cy = next(c for c in obs.columns if c.lower() in ("center_y", "centroid_y", "y"))
        x = obs[cx].to_numpy(float); y = obs[cy].to_numpy(float)
    # restrict to a single coherent 2D section if a slice/bregma column exists
    slice_col = None
    for c in obs.columns:
        if c.lower() in ("bregma", "slice", "z", "z_slice"):
            slice_col = c
            break
    meta_note = "whole object (no slice column found)"
    if slice_col is not None:
        vals, counts = np.unique(obs[slice_col].to_numpy(), return_counts=True)
        keep = vals[np.argmax(counts)]  # the most populated slice
        sel = obs[slice_col].to_numpy() == keep
        x, y = x[sel], y[sel]
        X = adata.X[sel]
        meta_note = f"single slice {slice_col}={keep} ({int(sel.sum())} cells)"
        print(f"  restricted to {meta_note}")
    else:
        X = adata.X
    tx = np.asarray(X.sum(axis=1)).ravel().astype(float)  # transcripts per cell = row sum
    wins = _fov_windows(x, y, FOV_SIDE_UM, N_FOV_PER_DATASET, FOV_SELECT_SEED)
    stats = _stats_for_windows("merfish_hypothal", x, y, tx, wins)
    return stats, meta_note


def main():
    rows = []
    meta = {
        "fov_side_um": FOV_SIDE_UM,
        "n_fov_per_dataset": N_FOV_PER_DATASET,
        "fov_select_seed": FOV_SELECT_SEED,
        "datasets": {},
    }
    try:
        xen = process_xenium()
        rows += xen
        meta["datasets"]["xenium"] = {
            "id": XENIUM_ID, "url": XENIUM_CELLS_URL, "access": "direct HTTP, no auth",
            "tx_per_cell_field": "transcript_counts (panel transcripts assigned to cell)",
        }
    except Exception as e:
        print(f"  XENIUM FAILED: {e!r}")
        meta["datasets"]["xenium"] = {"error": repr(e)}
    try:
        mer, note = process_merfish()
        rows += mer
        meta["datasets"]["merfish"] = {
            "id": MERFISH_ID, "url": MERFISH_H5AD_URL, "access": "direct HTTP (figshare), no auth",
            "tx_per_cell_field": "row sum of cell-by-gene matrix", "section": note,
        }
    except Exception as e:
        print(f"  MERFISH FAILED: {e!r}")
        meta["datasets"]["merfish"] = {"error": repr(e)}

    df = pd.DataFrame([asdict(r) for r in rows])
    out_csv = os.path.join(RESULTS, "realism.csv")
    df.to_csv(out_csv, index=False)
    with open(os.path.join(RESULTS, "realism_meta.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(f"\nwrote {out_csv} ({len(df)} FOVs)")
    if len(df):
        with pd.option_context("display.width", 160, "display.max_columns", None):
            print(df.to_string(index=False))
    return df


if __name__ == "__main__":
    main()
