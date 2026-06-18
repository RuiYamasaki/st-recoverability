"""Bayes-optimal (oracle) transcript-to-cell assignment, plus the naive
nearest-nucleus baseline.

Oracle derivation
-----------------
A transcript of gene g observed at position x was produced by some true cell c.
Under the known generative model the joint density of (source = c, gene = g,
observed position = x) factorises as

    P(c, g, x) = P(c) * P(g | type(c)) * P(x | c)

with, in this generator:
  - P(c) uniform: every cell has the same expected count (mean_tx_per_cell), so the
    source prior is flat and drops out of the argmax.
  - P(g | type(c)) = p[type(c), g]   (the fixed TYPE_GENE_COMPOSITION).
  - P(x | c) proportional to (f_c * phi_sigma)(x): the cell's uniform-in-Voronoi
    density f_c convolved with the isotropic Gaussian displacement kernel phi_sigma.

The Bayes-optimal (minimum-error) assignment is therefore

    c*(x, g) = argmax_c  p[type(c), g] * (f_c * phi_sigma)(x).

No estimator can beat this expected accuracy: it is the information-theoretic
ceiling for this generative model.

Computation
-----------
(f_c * phi_sigma) is evaluated on the field's pixel grid: each cell's indicator
mask (it owns the pixels nearest its centre) is Gaussian-blurred and normalised by
the cell's pixel area, giving D_c(pixel) ~ (f_c * phi_sigma). Because p[t, g] is
constant within a type, the argmax over cells reduces to, per type t and pixel,
the single best type-t cell:

    Dmax_t(pixel)  = max over type-t cells c of D_c(pixel)
    argcell_t(pixel) = the cell attaining that max.

Then for a transcript (x, g): pick the type maximising p[t, g] * Dmax_t(pixel(x)),
and assign to that type's argcell. This is exact in the pixel limit and avoids any
per-cell global computation (each cell is blurred only inside its bounding window).
"""
from __future__ import annotations

import os
import sys

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402
from generator import Field  # noqa: E402

NEG_INF = -1e30


def build_oracle_maps(field: Field):
    """Return (Dmax, argcell): each of shape (K, grid, grid).

    Dmax[t]  = per-pixel max blurred density over type-t cells.
    argcell[t] = per-pixel cell index attaining that max (-1 if no type-t cell near).
    """
    K, grid = config.N_TYPES, field.grid
    sigma_px = field.sigma_um / field.dx_um
    Dmax = np.zeros((K, grid, grid), dtype=np.float64)
    argcell = np.full((K, grid, grid), -1, dtype=np.int32)

    pad = int(np.ceil(4.0 * sigma_px)) + 1  # gaussian_filter truncates at 4 sigma
    for c in range(field.n_cells):
        pix = field.cell_pixels[c]
        if pix.size == 0:
            continue
        t = field.types[c]
        ry, rx = pix // grid, pix % grid
        y0, y1 = max(0, ry.min() - pad), min(grid, ry.max() + pad + 1)
        x0, x1 = max(0, rx.min() - pad), min(grid, rx.max() + pad + 1)
        win = np.zeros((y1 - y0, x1 - x0), dtype=np.float64)
        win[ry - y0, rx - x0] = 1.0
        if sigma_px > 1e-3:
            blurred = gaussian_filter(win, sigma_px, mode="constant", truncate=4.0)
        else:
            blurred = win
        dens = blurred / field.cell_area_px[c]   # ~ (f_c * phi_sigma) up to area norm
        sub_max = Dmax[t, y0:y1, x0:x1]
        sub_arg = argcell[t, y0:y1, x0:x1]
        better = dens > sub_max
        sub_max[better] = dens[better]
        sub_arg[better] = c
    return Dmax, argcell


def _pixels_of(field: Field, xy: np.ndarray):
    px = np.clip((xy[:, 0] / field.dx_um).astype(np.int64), 0, field.grid - 1)
    py = np.clip((xy[:, 1] / field.dx_um).astype(np.int64), 0, field.grid - 1)
    return py, px


def oracle_assign(field: Field, Dmax: np.ndarray, argcell: np.ndarray,
                  obs_xy: np.ndarray, gene: np.ndarray) -> np.ndarray:
    """Vectorised Bayes-optimal assignment. Returns assigned cell index per transcript."""
    K = config.N_TYPES
    py, px = _pixels_of(field, obs_xy)
    D = Dmax[:, py, px]              # (K, T) blurred densities at each transcript's pixel
    A = argcell[:, py, px]          # (K, T) candidate cell per type
    pg = config.TYPE_GENE_COMPOSITION[:, gene]  # (K, T) gene evidence p[t, g]
    logD = np.full_like(D, NEG_INF)
    pos = D > 0
    logD[pos] = np.log(D[pos])
    score = np.where(pos, logD + np.log(pg), NEG_INF)  # (K, T)
    best_t = np.argmax(score, axis=0)                          # (T,)
    T = obs_xy.shape[0]
    assigned = A[best_t, np.arange(T)]
    # fallback for the (rare) transcript whose pixel has no type with positive density
    bad = assigned < 0
    if bad.any():
        tree = cKDTree(field.centers)
        _, nn = tree.query(obs_xy[bad], k=1)
        assigned = assigned.copy()
        assigned[bad] = nn.astype(np.int32)
    return assigned.astype(np.int32)


def naive_assign(field: Field, obs_xy: np.ndarray) -> np.ndarray:
    """Naive nearest-nucleus baseline: assign each transcript to the nearest cell
    centre of its observed position (ignores the displacement model and gene id)."""
    tree = cKDTree(field.centers)
    _, nn = tree.query(obs_xy, k=1)
    return nn.astype(np.int32)
