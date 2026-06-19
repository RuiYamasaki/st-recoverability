"""Run Baysor (Petukhov et al., Nature Biotechnology 2022) on an exported headroom config,
in nuclei-prior and free-segmentation modes.

Install route: prebuilt Linux binaries on the Baysor GitHub releases are x86_64-only and the
newest (cpp-0.8.x) releases ship no binaries; on this linux-aarch64 host Baysor was built
natively from the Julia package at tag v0.7.1 (Pkg.add from git into an isolated Julia 1.10
depot), invoked via the wrapper ~/baysor_run.sh (Baysor.run_cli). Version recorded by the
driver/report.

Baysor jointly segments and assigns: it clusters molecules into cells from their spatial +
expression neighbourhood. It takes an expected-cell-radius scale and, optionally, a prior
segmentation.
  - nuclei-prior: a per-molecule prior label (the molecule's nucleus cell id if it lies
    within NUCLEUS_RADIUS_UM of its nearest true centre, else 0 = no prior) is passed via
    Baysor's documented `:column` prior-segmentation mechanism, with the default
    prior-segmentation-confidence. Same 3 um nucleus disks as the other methods.
  - free-segmentation: no prior; Baysor segments de novo. This is Baysor's native mode.
In both modes the expected-cell-radius `scale` is set to the true mean cell radius (derived
from the centroid spacing); Baysor requires this geometry input, it is the correct value and
is not tuned toward any accuracy. All other parameters are documented defaults.

Baysor's optimiser is stochastic; v0.7 exposes no CLI seed, so output is not bit-for-bit
reproducible (run-to-run variation quantified in the report). Returns (assigned, nearest):
per-molecule assigned cell id (-1 for noise/unassigned).
"""
from __future__ import annotations

import glob
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

NUCLEUS_RADIUS_UM = 3.0
BAYSOR = os.path.expanduser("~/baysor_run.sh")
# Documented Baysor Xenium preset (configs/xenium.toml shipped with Baysor v0.7.1):
MIN_MOLECULES_PER_CELL = 50
PRIOR_SEG_CONFIDENCE = 0.5


def baysor_version():
    try:
        out = subprocess.run([BAYSOR, "--version"], capture_output=True, text=True, timeout=600)
        return (out.stdout + out.stderr).strip().splitlines()[-1] if (out.stdout or out.stderr) else "unknown"
    except Exception as e:  # noqa: BLE001
        return f"unknown ({e!r})"


def _scale_um(cells_df):
    """Expected cell radius = true mean cell radius, area_per_cell = (centroid bounding-box
    area) / n_cells, r_mean = sqrt(area_per_cell / pi). This recovers the generator's true
    r_mean = sqrt((1e6/packing)/pi) to within ~3% (the centroid bounding box closely tracks
    the field domain). The median nearest-neighbour spacing under-estimates it by ~40% under
    the generator's cell jitter, so it is not used."""
    cen = cells_df[["x_centroid", "y_centroid"]].to_numpy()
    n = len(cen)
    w = cen[:, 0].max() - cen[:, 0].min()
    h = cen[:, 1].max() - cen[:, 1].min()
    return float(np.sqrt((w * h / n) / np.pi))


def _nucleus_prior(transcripts_df, cells_df):
    cen = cells_df[["x_centroid", "y_centroid"]].to_numpy()
    xy = transcripts_df[["x", "y"]].to_numpy()
    d, idx = cKDTree(cen).query(xy, k=1)
    prior = np.where(d <= NUCLEUS_RADIUS_UM, idx + 1, 0).astype(np.int64)  # 0 = no prior
    return prior, idx


def run_baysor(transcripts_df, cells_df, model, with_prior=True):
    prior, nearest = _nucleus_prior(transcripts_df, cells_df)
    n = len(transcripts_df)
    scale = _scale_um(cells_df)
    work = tempfile.mkdtemp(prefix="baysor_")
    mol = os.path.join(work, "mol.csv")
    df = pd.DataFrame({
        "x": transcripts_df["x"].to_numpy(),
        "y": transcripts_df["y"].to_numpy(),
        "gene": transcripts_df["gene"].to_numpy(),
    })
    if with_prior:
        df["prior"] = prior
    df.to_csv(mol, index=False)

    out_dir = os.path.join(work, "out")
    os.makedirs(out_dir, exist_ok=True)
    cmd = [BAYSOR, "run", "-x", "x", "-y", "y", "-g", "gene",
           "-m", str(MIN_MOLECULES_PER_CELL), "-s", f"{scale:.4f}",
           "-o", out_dir + os.sep]
    if with_prior:
        cmd += ["--prior-segmentation-confidence", str(PRIOR_SEG_CONFIDENCE)]
    cmd += [mol]
    if with_prior:
        cmd += [":prior"]
    env = dict(os.environ, JULIA_NUM_THREADS=os.environ.get("JULIA_NUM_THREADS", "4"))
    subprocess.run(cmd, check=True, capture_output=True, text=True, env=env, timeout=3600)

    seg = _find_segmentation(work)
    out = pd.read_csv(seg)
    assert len(out) == n, f"baysor returned {len(out)} rows for {n} molecules"
    assigned = _cell_labels_to_ids(out)
    return assigned, nearest


def _find_segmentation(work):
    cands = (glob.glob(os.path.join(work, "out", "*segmentation.csv")) +
             glob.glob(os.path.join(work, "out", "seg*.csv")) +
             glob.glob(os.path.join(work, "**", "*segmentation.csv"), recursive=True))
    cands = [c for c in cands if "cell_stats" not in os.path.basename(c)]
    if not cands:
        raise FileNotFoundError(f"no Baysor segmentation output under {work}/out")
    return cands[0]


def _cell_labels_to_ids(out):
    """Map Baysor's per-molecule cell label to integer ids; noise/unassigned -> -1.
    Baysor marks background via is_noise and/or cell==0 (or empty)."""
    cell = out["cell"]
    noise = out["is_noise"].to_numpy().astype(bool) if "is_noise" in out.columns else np.zeros(len(out), bool)
    # normalise the cell column to a string, blank/"0"/"0.0" = background
    s = cell.astype(str).str.strip()
    bg = noise | s.isin(["0", "0.0", "", "nan", "NaN", "-1"])
    codes = pd.factorize(s.where(~bg, other="__bg__"))[0].astype(np.int64)
    # ensure background rows are -1
    assigned = np.where(bg, -1, codes)
    return assigned.astype(np.int64)
