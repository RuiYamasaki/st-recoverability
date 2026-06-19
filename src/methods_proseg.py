"""Run Proseg (Newell lab, Nature Methods 2025) on an exported headroom config.

Install route: bioconda, `micromamba install -c bioconda rust-proseg` (package `rust-proseg`,
binary `proseg`); version recorded by the driver. Native linux-aarch64 build.

Proseg is a re-segmentation method: it takes an initial nucleus (or cell) segmentation and
jointly refines cell morphology and reassigns transcripts via a voxel sampler. It REQUIRES
an initial segmentation (`--cell-id-column`); it has no prior-free mode (running without the
prior aborts with "Missing required argument: --cell-id-column"). So only nuclei-prior is
supported here; the free-segmentation column is recorded as unsupported by the driver.

Nuclei-prior: each transcript whose observed position lies within NUCLEUS_RADIUS_UM of its
nearest true cell centre is tagged with that cell id (the nucleus prior, the same 3 um disk
prior used for pciSeq/ComSeg); the rest are unassigned (-1). Proseg grows cells from these
nucleus seeds. Documented defaults are used for every algorithm parameter; only the column
mapping, 2D flag, and output paths are set.

Proseg is an MCMC sampler with no seed flag, so its output is not bit-for-bit reproducible
(run-to-run variation is quantified in the report). Returns (assigned_cell_id, nearest):
per-transcript assigned cell id (-1 for background/unassigned).
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

NUCLEUS_RADIUS_UM = 3.0


def proseg_version():
    try:
        out = subprocess.run(["proseg", "--version"], capture_output=True, text=True, timeout=60)
        return out.stdout.strip() or out.stderr.strip()
    except Exception as e:  # noqa: BLE001
        return f"unknown ({e!r})"


def _nucleus_prior(transcripts_df, cells_df):
    centres = cells_df[["x_centroid", "y_centroid"]].to_numpy()
    xy = transcripts_df[["x", "y"]].to_numpy()
    d, idx = cKDTree(centres).query(xy, k=1)
    nuc = np.where(d <= NUCLEUS_RADIUS_UM, idx, -1).astype(np.int64)   # 0-based cell id or -1
    return nuc, idx


def run_proseg(transcripts_df, cells_df, model, with_prior=True, extra_args=None):
    """extra_args: optional list of documented Proseg CLI flags appended for the fairness
    sweep (for example ["--voxel-size", "0.5"], ["--prior-seg-reassignment-prob", "0.7"],
    ["--no-diffusion"]). None = documented defaults."""
    if not with_prior:
        raise RuntimeError(
            "Proseg has no prior-free mode: it requires an initial segmentation "
            "(--cell-id-column). Free-segmentation is unsupported for Proseg.")
    nuc, nearest = _nucleus_prior(transcripts_df, cells_df)
    n = len(transcripts_df)
    work = tempfile.mkdtemp(prefix="proseg_")
    in_csv = os.path.join(work, "tx.csv")
    out_csv = os.path.join(work, "tx_meta.csv.gz")
    pd.DataFrame({
        "transcript_id": np.arange(n),
        "x": transcripts_df["x"].to_numpy(),
        "y": transcripts_df["y"].to_numpy(),
        "z": np.zeros(n),
        "gene": transcripts_df["gene"].to_numpy(),
        "nucleus": nuc,
    }).to_csv(in_csv, index=False)

    cmd = [
        "proseg", in_csv,
        "-x", "x", "-y", "y", "-z", "z", "--gene-column", "gene", "--ignore-z-coord",
        "--cell-id-column", "nucleus", "--cell-id-unassigned", "-1",
        "--overwrite",
    ]
    if extra_args:
        cmd += list(extra_args)
    cmd += ["--output-transcript-metadata", out_csv]
    subprocess.run(cmd, cwd=work, check=True, capture_output=True, text=True)

    out = pd.read_csv(out_csv)
    # Proseg preserves input row order and count; map assignment -> cell id, background -> -1.
    assert len(out) == n, f"proseg returned {len(out)} rows for {n} transcripts"
    a = out["assignment"].to_numpy()
    bg = out["background"].to_numpy().astype(bool)
    assigned = np.where(np.isnan(a) | bg, -1, np.nan_to_num(a, nan=-1)).astype(np.int64)
    return assigned, nearest
