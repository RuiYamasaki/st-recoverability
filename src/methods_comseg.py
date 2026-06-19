"""Run ComSeg (Defard et al. 2024) on an exported headroom config, in nuclei-prior mode:
the prior is per-transcript nucleus membership (in_nucleus, from small disks around the
true cell centres) plus the true cell centroids as landmarks; ComSeg builds a transcript
co-expression graph and associates each RNA to a cell. Returns the per-transcript assigned
cell id (-1 for unassigned). Documented run_all defaults are used.
"""
from __future__ import annotations

import os
import sys
import tempfile
import warnings

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

NUCLEUS_RADIUS_UM = 3.0
MEAN_CELL_DIAMETER_UM = 10.0


def _in_nucleus(transcripts_df, cells_df):
    """Per-transcript nucleus id (1..N) if within NUCLEUS_RADIUS of its nearest centre,
    else 0 (ComSeg's in_nucleus prior convention: 0 = not in a nucleus)."""
    centres = cells_df[["x_centroid", "y_centroid"]].to_numpy()
    xy = transcripts_df[["x", "y"]].to_numpy()
    d, idx = cKDTree(centres).query(xy, k=1)
    prior = np.where(d <= NUCLEUS_RADIUS_UM, idx + 1, 0).astype(np.int64)
    return prior, idx


def run_comseg(transcripts_df, cells_df, model, with_prior=True,
               mean_cell_diameter=None, max_radius_factor=2.5):
    """Optional params (for the fairness sweep; None = documented default):
    mean_cell_diameter (um, default 10; ComSeg's most impactful parameter, sets the
    co-expression graph radius), max_radius_factor (max_cell_radius = mcd * factor)."""
    import random
    from comseg import dataset, dictionary
    mcd = MEAN_CELL_DIAMETER_UM if mean_cell_diameter is None else float(mean_cell_diameter)
    np.random.seed(0); random.seed(0)   # ComSeg samples internally; seed for reproducibility
    prior, nearest = _in_nucleus(transcripts_df, cells_df)
    work = tempfile.mkdtemp(prefix="comseg_")
    img_dir = os.path.join(work, "spots"); os.makedirs(img_dir)
    cen_dir = os.path.join(work, "centroids"); os.makedirs(cen_dir)

    n = len(transcripts_df)
    df = pd.DataFrame({
        "x": transcripts_df["x"].to_numpy(), "y": transcripts_df["y"].to_numpy(),
        "z": np.zeros(n), "gene": transcripts_df["gene"].to_numpy(),
        "in_nucleus": prior if with_prior else np.zeros(n, dtype=np.int64),
    })
    df.to_csv(os.path.join(img_dir, "fov0.csv"), index=False)
    # centroid CSV needs an `in_nucleus` column holding each centroid's cell id (1..N),
    # matching the per-spot in_nucleus prior convention (cell_id + 1).
    cen = pd.DataFrame({"x": cells_df["x_centroid"].to_numpy(),
                        "y": cells_df["y_centroid"].to_numpy(),
                        "z": np.zeros(len(cells_df)),
                        "in_nucleus": np.arange(len(cells_df)) + 1})
    cen.to_csv(os.path.join(cen_dir, "fov0.csv"), index=False)

    # Free mode: ComSeg's valid prior-free graph-partition method is "louvain" with no
    # per-RNA prior (prior_name=None). Note ComSeg's run_all still associates RNA communities
    # to the centroid landmarks (classify_centroid + associate_rna2landmark), so the true
    # centroids are retained even in this mode; it is therefore "no per-RNA nucleus prior",
    # not fully landmark-free.
    ds = dataset.ComSegDataset(
        path_dataset_folder=img_dir, prior_name=("in_nucleus" if with_prior else None),
        image_csv_files=["fov0.csv"], centroid_csv_files=["fov0.csv"],
        path_cell_centroid=cen_dir, dict_scale={"x": 1, "y": 1, "z": 1},
        mean_cell_diameter=mcd, gene_column="gene")
    ds.compute_edge_weight()
    cd = dictionary.ComSegDict(
        dataset=ds, mean_cell_diameter=mcd,
        community_detection="with_prior" if with_prior else "louvain", seed=0)
    cd.run_all(max_cell_radius=mcd * max_radius_factor)

    # extract per-spot cell_index_pred from the single FOV graph (node i == input row i)
    key = [k for k in cd.keys()][0]
    G = cd[key].G
    assigned = np.full(n, -1, dtype=np.int64)
    for i in range(n):
        if i in G.nodes:
            cip = G.nodes[i].get("cell_index_pred", 0)
            # cell_index_pred is the centroid's in_nucleus id (cell_id + 1); 0 = background
            if cip is not None and cip >= 1:
                assigned[i] = int(cip) - 1
    return assigned, nearest
