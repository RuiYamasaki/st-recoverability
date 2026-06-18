"""Segmenter-headroom shared machinery: export known-truth synthetic data at the most
realistic Gate 3 configuration, and score any method's transcript-to-cell assignment
against the truth with the same accuracy metric used for the oracle.

Configuration (most realistic Gate 3 config): Xenium-breast realistic overlapping
expression, negative-binomial within-type emission, data-pinned displacement at the Gate 3
breast value 1.99 um and its CI ends 1.43 and 2.63 um, dense packing 13,625 cells/mm^2
(plus sparse 2,575 and representative 6,534 for context), realistic density 164 tx/cell.

Scoring: methods output their own cell labels, so transcript-assignment accuracy is the
fraction of interior transcripts whose method-cell is matched to their true cell under the
optimal one-to-one matching of method-cells to true-cells (Hungarian on the contingency
table). For the oracle and the naive baseline (which already use the true cell ids) this
matched accuracy equals the direct id-match accuracy. Background/unassigned method labels
count as errors.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402
from generator import build_field, generate_transcripts  # noqa: E402
from oracle import build_oracle_maps, oracle_assign, naive_assign  # noqa: E402
from gate2_pin import build_xenium, ORACLE_RES_CELL, ORACLE_GRID_MAX  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
DATA = os.path.join(ROOT, "data")
EXPORT_DIR = os.path.join(DATA, "headroom")

HEADROOM_SEED = config.MASTER_SEED + 600000
DENSE, SPARSE, REP = 13625.0, 2575.0, 6534.0
SIGMA_POINT, SIGMA_LO, SIGMA_HI = 1.99, 1.43, 2.63   # Gate 3 breast pinned + CI ends
DENSITY = 164.0
N_TARGET = 250          # FOV cell count; oracle accuracy is field-size-independent, so the
                        # ceiling is unchanged. Kept modest to keep the external methods tractable.
EMISSION = "nbinom"

# (config name, packing, sigma)
CONFIGS = [
    ("dense_sigma1.43", DENSE, SIGMA_LO),
    ("dense_sigma1.99", DENSE, SIGMA_POINT),
    ("dense_sigma2.63", DENSE, SIGMA_HI),
    ("sparse_sigma1.99", SPARSE, SIGMA_POINT),
    ("representative_sigma1.99", REP, SIGMA_POINT),
]


def matched_accuracy(method_label: np.ndarray, true_cell: np.ndarray, interior: np.ndarray) -> float:
    """Fraction of interior transcripts correctly assigned under the optimal one-to-one
    matching of method-cells to true-cells."""
    m = interior
    n = int(m.sum())
    if n == 0:
        return float("nan")
    ml, tl = method_label[m], true_cell[m]
    umeth, mi = np.unique(ml, return_inverse=True)
    utrue, ti = np.unique(tl, return_inverse=True)
    C = np.zeros((len(umeth), len(utrue)), dtype=np.int64)
    np.add.at(C, (mi, ti), 1)
    r, c = linear_sum_assignment(-C)
    return float(C[r, c].sum() / n)


def build_config(model, abundance, packing, sigma, seed):
    """Build the field + transcripts (NB emission) and the oracle/naive assignments."""
    field = build_field(packing, sigma, seed, model=model,
                        res_cell=ORACLE_RES_CELL, grid_max=ORACLE_GRID_MAX, n_target=N_TARGET)
    tx = generate_transcripts(field, DENSITY, seed + 1, emission=EMISSION)
    Dmax, argcell = build_oracle_maps(field)
    oa = oracle_assign(field, Dmax, argcell, tx.obs_xy, tx.gene)
    na = naive_assign(field, tx.obs_xy)   # nearest centre == nuclei-prior assignment
    return field, tx, oa, na


def export_config(model, field, tx, oa, na, name):
    """Write the transcript table (with the true-centre nuclei prior) and the cell-centre
    table for the external methods. Raw tables go under data/headroom/ (gitignored)."""
    d = os.path.join(EXPORT_DIR, name)
    os.makedirs(d, exist_ok=True)
    gene_names = np.array(model.gene_names)
    tdf = pd.DataFrame({
        "transcript_id": np.arange(tx.obs_xy.shape[0]),
        "x": tx.obs_xy[:, 0], "y": tx.obs_xy[:, 1],
        "gene": gene_names[tx.gene],
        "true_cell": tx.true_cell,
        "prior_cell": na,                 # nearest-centre = nuclei prior cell id
        "interior": tx.interior.astype(int),
        "oracle_cell": oa, "naive_cell": na,
    })
    tdf.to_csv(os.path.join(d, "transcripts.csv"), index=False)
    cdf = pd.DataFrame({
        "cell_id": np.arange(field.n_cells),
        "x_centroid": field.centers[:, 0], "y_centroid": field.centers[:, 1],
        "type": field.types, "interior": field.interior_cell.astype(int),
    })
    cdf.to_csv(os.path.join(d, "cells.csv"), index=False)
    return d
