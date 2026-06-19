"""Shared machinery for the method-fairness gate: over-segmentation-robust metrics computed
from any method's per-transcript assignment against the known truth, over interior
transcripts.

Metrics (all over interior transcripts; background/unassigned = method label < 0):
  - acc_one_to_one : the existing Hungarian one-to-one matched accuracy
    (headroom_common.matched_accuracy). Penalises over-segmentation: each true cell can be
    matched to at most one predicted cell, so transcripts in a split-off predicted cell count
    as errors.
  - acc_many_to_one : assignment-homogeneity accuracy. Each predicted cell is mapped to the
    true cell it most overlaps (its majority true label); a transcript is correct if its
    predicted cell's majority-true-label equals its true cell. Does NOT penalise splitting one
    true cell into several predicted cells. Background transcripts count as errors (denominator
    is all interior transcripts).
  - ari : adjusted Rand index of the transcript co-assignment partition vs truth (sklearn
    adjusted_rand_score), with background treated as its own cluster. Symmetric: penalises both
    over- and under-segmentation; not gameable by splitting.
  - pred_true_ratio : (number of distinct predicted cells over interior transcripts) /
    (number of true interior cells). Quantifies over- (>1) or under- (<1) segmentation.
  - frac_assigned : fraction of interior transcripts given a real cell (not background).
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from headroom_common import matched_accuracy  # noqa: E402  (reused one-to-one scorer)


def many_to_one_homogeneity(assigned, true_cell, interior):
    """Each predicted cell -> its majority true cell; transcript correct if its predicted
    cell's majority equals its true cell. Background (assigned < 0) always incorrect.
    Denominator = all interior transcripts."""
    m = interior
    a = assigned[m]
    t = true_cell[m]
    n = int(m.sum())
    if n == 0:
        return float("nan")
    valid = a >= 0
    if not valid.any():
        return 0.0
    df = pd.DataFrame({"a": a[valid], "t": t[valid]})
    # majority true label per predicted cell
    maj = df.groupby("a")["t"].agg(lambda s: s.value_counts().idxmax())
    pred_maj = df["a"].map(maj).to_numpy()
    correct = int((pred_maj == df["t"].to_numpy()).sum())
    return correct / n


def all_metrics(assigned, true_cell, interior):
    """Return the full metric dict for one per-transcript assignment."""
    m = interior
    a = assigned[m]
    t = true_cell[m]
    n = int(m.sum())
    n_pred = int(len(np.unique(a[a >= 0])))
    n_true = int(len(np.unique(t)))
    return {
        "acc_one_to_one": matched_accuracy(assigned, true_cell, interior),
        "acc_many_to_one": many_to_one_homogeneity(assigned, true_cell, interior),
        "ari": float(adjusted_rand_score(t, a)),
        "n_pred_cells": n_pred,
        "n_true_cells": n_true,
        "pred_true_ratio": (n_pred / n_true) if n_true else float("nan"),
        "frac_assigned": float((a >= 0).mean()) if n else float("nan"),
        "n_interior_transcripts": n,
    }
