"""Run pciSeq (Qian et al., Nature Methods 2020) on an exported headroom config, in
nuclei-prior mode: the prior segmentation is a nucleus label image (small disks around the
true cell centres), the reference is the generator's true per-type mean expression, and
pciSeq reassigns each transcript to a cell by expression + proximity. Returns the
per-transcript assigned cell id (matched-accuracy is computed by the caller).

Documented defaults are used (only the non-algorithmic viewer/diagnostics UI is disabled).
"""
from __future__ import annotations

import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

PIXEL_UM = 0.5            # label-image pixel size (micrometres)
NUCLEUS_RADIUS_UM = 3.0   # nucleus disk radius for the prior segmentation


def _nucleus_label_image(cells_df, pixel_um=PIXEL_UM, nucleus_radius_um=NUCLEUS_RADIUS_UM):
    xs = cells_df["x_centroid"].to_numpy()
    ys = cells_df["y_centroid"].to_numpy()
    W = int(np.ceil(xs.max() / pixel_um)) + 2
    H = int(np.ceil(ys.max() / pixel_um)) + 2
    lab = np.zeros((H, W), dtype=np.int32)
    r = int(round(nucleus_radius_um / pixel_um))
    yy, xx = np.ogrid[-r:r + 1, -r:r + 1]
    disk = (xx ** 2 + yy ** 2) <= r ** 2
    for cid, (cx, cy) in enumerate(zip(xs, ys)):
        px, py = int(round(cx / pixel_um)), int(round(cy / pixel_um))
        y0, y1 = max(0, py - r), min(H, py + r + 1)
        x0, x1 = max(0, px - r), min(W, px + r + 1)
        sub = disk[(y0 - (py - r)):(y0 - (py - r)) + (y1 - y0),
                   (x0 - (px - r)):(x0 - (px - r)) + (x1 - x0)]
        # later cells overwrite on overlap; nuclei are small so overlaps are rare
        region = lab[y0:y1, x0:x1]
        region[sub] = cid + 1   # label = cell_id + 1 (0 = background)
    return lab


def run_pciseq(transcripts_df, cells_df, model):
    """Return a per-transcript assigned cell id (-1 for background/unassigned)."""
    import pciSeq
    lab = _nucleus_label_image(cells_df)
    coo = coo_matrix(lab)
    spots = pd.DataFrame({
        "Gene": transcripts_df["gene"].to_numpy(),
        "x": transcripts_df["x"].to_numpy() / PIXEL_UM,
        "y": transcripts_df["y"].to_numpy() / PIXEL_UM,
    })
    ref = pd.DataFrame(model.mean_expr.T, index=list(model.gene_names),
                       columns=[f"t{t}" for t in range(model.n_types)])
    opts = {"launch_viewer": False, "launch_diagnostics": False}
    _, geneData = pciSeq.fit(spots=spots, coo=coo, scRNAseq=ref, opts=opts)
    if len(geneData) != len(transcripts_df):
        # align by rounded pixel coords + gene if order/length differs
        key = lambda df, xs, ys, g: pd.Series(
            np.round(df[xs] * 1000).astype(np.int64).astype(str) + "_" +
            np.round(df[ys] * 1000).astype(np.int64).astype(str) + "_" + df[g].astype(str))
        gd = geneData.copy()
        gd["_k"] = key(gd, "x", "y", "Gene")
        sp = spots.copy(); sp["_k"] = key(sp, "x", "y", "Gene")
        nb = gd.set_index("_k")["neighbour"].reindex(sp["_k"]).to_numpy()
    else:
        nb = geneData["neighbour"].to_numpy()
    # neighbour is the label image value (cell_id+1); background = max label (n_cells+1) or 0
    assigned = np.full(len(transcripts_df), -1, dtype=np.int64)
    valid = (nb >= 1) & (nb <= cells_df.shape[0])
    assigned[valid] = (nb[valid] - 1).astype(np.int64)
    return assigned
