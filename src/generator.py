"""Minimal 2D synthetic spatial-transcriptomics generator.

Pipeline (all randomness seeded):
  1. Cell centres: a perturbed (jittered) square lattice. Spacing s sets the
     packing; Gaussian jitter (JITTER_FRAC * s) makes the tessellation irregular
     without producing degenerate slivers. Packing = N / area is exact.
  2. True cell boundaries: the Voronoi tessellation of the centres. We never need
     the polygons explicitly: a pixel's Voronoi cell == its nearest centre, so the
     boundaries are represented by a nearest-centre label image.
  3. Cell types: each cell drawn from the fixed TYPE_PROPORTIONS.
  4. Expression: cell of type t emits, per gene g, Poisson(mean_tx * p[t, g])
     transcripts (config.TYPE_GENE_COMPOSITION). mean_tx is the 'density' axis.
  5. Placement: each transcript's true position is uniform inside its cell's
     Voronoi region (a random pixel of that cell + sub-pixel jitter). The observed
     position adds isotropic Gaussian displacement N(0, sigma^2 I) (the diffusion).
  6. The true transcript-to-cell assignment is recorded for every transcript.

Three exposed parameters:
  - mean_tx_per_cell  (density)
  - packing_cells_per_mm2  (cell packing; equivalently mean cell radius)
  - sigma_um  (displacement / diffusion)

Units are micrometres throughout, so the axes map directly onto real summary
statistics (median transcripts/cell, cells/mm^2).
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field as dc_field

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402
import expression  # noqa: E402

# --- generator structural constants (fixed modelling choices) ---
N_TARGET = 400            # cells per field, held constant so resolution/power are uniform
JITTER_FRAC = 0.33        # lattice jitter as a fraction of lattice spacing
RES_CELL = 6.0            # pixels per mean cell radius (grid resolution target)
GRID_MIN, GRID_MAX = 128, 768
INTERIOR_MARGIN_SIGMA = 3.0   # interior cells must sit >= r_mean + 3*sigma from the edge
ANISO_ELONGATION = (1.5, 2.5)  # Gate 1 non-Voronoi geometry: per-cell elongation ratio range

# The default expression model reproduces Gate 0 exactly. Gate 1 passes a realistic model.
DEFAULT_MODEL = expression.build_disjoint_model()


@dataclass
class Field:
    packing_cells_per_mm2: float
    sigma_um: float
    seed: int
    L_um: float
    n_cells: int
    r_mean_um: float          # mean cell radius = sqrt(area_per_cell / pi)
    dx_um: float
    grid: int                 # H == W
    centers: np.ndarray       # (N, 2) micrometres
    types: np.ndarray         # (N,) int in [0, K)
    label: np.ndarray         # (grid, grid) int: nearest-centre cell index per pixel
    cell_pixels: list         # cell_pixels[c] = flat pixel indices belonging to cell c
    cell_area_px: np.ndarray  # (N,) pixel count per cell
    interior_cell: np.ndarray # (N,) bool
    model: object = None       # ExpressionModel used for this field
    geometry: str = "voronoi"  # "voronoi" (Gate 0) or "aniso" (Gate 1 structural)


def _pixel_centers(grid: int, dx: float) -> np.ndarray:
    c = (np.arange(grid) + 0.5) * dx
    xx, yy = np.meshgrid(c, c, indexing="xy")  # x varies along columns
    return np.column_stack([xx.ravel(), yy.ravel()])


def _aniso_label(centers, pts, grid, rng):
    """Non-Voronoi label image: each pixel assigned to the cell minimising an
    anisotropic (elongated, randomly oriented) Mahalanobis distance. Cells become
    irregular and elongated rather than convex Voronoi polygons."""
    n = centers.shape[0]
    theta = rng.uniform(0.0, np.pi, size=n)
    elong = rng.uniform(*ANISO_ELONGATION, size=n)   # ratio a/b
    a = np.sqrt(elong); b = 1.0 / np.sqrt(elong)     # area-preserving (a*b == 1)
    P = pts.shape[0]
    best = np.full(P, np.inf); arg = np.full(P, -1, dtype=np.int32)
    for c in range(n):
        ct, st = np.cos(theta[c]), np.sin(theta[c])
        d = pts - centers[c]
        u = ct * d[:, 0] + st * d[:, 1]              # along major axis
        v = -st * d[:, 0] + ct * d[:, 1]             # along minor axis
        q = (u / a[c]) ** 2 + (v / b[c]) ** 2
        upd = q < best
        best[upd] = q[upd]; arg[upd] = c
    return arg.reshape(grid, grid).astype(np.int32)


def build_field(packing_cells_per_mm2: float, sigma_um: float, seed: int,
                n_target: int = N_TARGET, model=None, geometry: str = "voronoi") -> Field:
    model = model or DEFAULT_MODEL
    rng = np.random.default_rng(seed)
    # domain side so that N_TARGET cells give the requested packing exactly
    area_um2 = n_target / packing_cells_per_mm2 * 1e6
    L = float(np.sqrt(area_um2))
    r_mean = float(np.sqrt((area_um2 / n_target) / np.pi))

    # perturbed square lattice of centres
    side = int(round(np.sqrt(n_target)))
    s = L / side                      # lattice spacing
    gx = (np.arange(side) + 0.5) * s
    X, Y = np.meshgrid(gx, gx, indexing="xy")
    centers = np.column_stack([X.ravel(), Y.ravel()])
    centers = centers + rng.normal(0.0, JITTER_FRAC * s, size=centers.shape)
    centers = np.clip(centers, 1e-6, L - 1e-6)
    n_cells = centers.shape[0]

    # grid resolution: ~RES_CELL pixels per mean cell radius, clamped
    dx_target = r_mean / RES_CELL
    grid = int(np.clip(round(L / dx_target), GRID_MIN, GRID_MAX))
    dx = L / grid
    pts = _pixel_centers(grid, dx)

    if geometry == "voronoi":
        # nearest-centre label image == Voronoi membership of each pixel
        _, nn = cKDTree(centers).query(pts, k=1)
        label = nn.reshape(grid, grid).astype(np.int32)
    elif geometry == "aniso":
        label = _aniso_label(centers, pts, grid, rng)
    else:
        raise ValueError(f"unknown geometry {geometry!r}")

    # per-cell pixel lists and areas
    order = np.argsort(label.ravel(), kind="stable")
    sorted_lab = label.ravel()[order]
    boundaries = np.searchsorted(sorted_lab, np.arange(n_cells + 1))
    cell_pixels = [order[boundaries[c]:boundaries[c + 1]] for c in range(n_cells)]
    cell_area_px = np.diff(boundaries).astype(np.int64)

    # cell types from the model's proportions
    types = rng.choice(model.n_types, size=n_cells, p=model.proportions).astype(np.int32)

    # interior cells: far enough from the domain edge that the clipped Voronoi
    # region and displaced transcripts are not edge artifacts
    margin = r_mean + INTERIOR_MARGIN_SIGMA * sigma_um
    interior = (
        (centers[:, 0] >= margin) & (centers[:, 0] <= L - margin) &
        (centers[:, 1] >= margin) & (centers[:, 1] <= L - margin)
    )

    return Field(
        packing_cells_per_mm2=packing_cells_per_mm2, sigma_um=sigma_um, seed=seed,
        L_um=L, n_cells=n_cells, r_mean_um=r_mean, dx_um=dx, grid=grid,
        centers=centers, types=types, label=label, cell_pixels=cell_pixels,
        cell_area_px=cell_area_px, interior_cell=interior,
        model=model, geometry=geometry,
    )


@dataclass
class Transcripts:
    obs_xy: np.ndarray        # (T, 2) observed (displaced) positions, micrometres
    true_xy: np.ndarray       # (T, 2) pre-displacement positions
    gene: np.ndarray          # (T,) int gene id
    true_cell: np.ndarray     # (T,) int source cell index
    true_type: np.ndarray     # (T,) int source cell type
    interior: np.ndarray      # (T,) bool: source cell is interior


def generate_transcripts(field: Field, mean_tx_per_cell: float, seed: int,
                         displacement: str = "gaussian", disp_epsilon: float = 0.0,
                         abundance: np.ndarray = None) -> Transcripts:
    """Emit transcripts. displacement: 'gaussian' (Gate 0) or 'gauss_uniform' (Gate 1
    structural: a fraction disp_epsilon of transcripts land uniformly in the domain,
    a heavy-tailed contamination). abundance (K,), if given, scales each type's mean
    total count (used by leakage fitting to match real per-type abundance); None keeps
    the uniform mean_tx_per_cell across types that the oracle's uniform prior assumes."""
    rng = np.random.default_rng(seed)
    model = field.model or DEFAULT_MODEL
    p = model.composition  # (K, G)
    dx = field.dx_um

    true_xy_list, gene_list, cell_list = [], [], []
    for c in range(field.n_cells):
        t = field.types[c]
        scale = mean_tx_per_cell if abundance is None else mean_tx_per_cell * abundance[t]
        # per-gene Poisson counts -> total ~ Poisson(scale)
        counts = rng.poisson(scale * p[t])     # (G,)
        n_c = int(counts.sum())
        if n_c == 0:
            continue
        genes_c = np.repeat(np.arange(model.n_genes), counts)
        pix = field.cell_pixels[c]
        if pix.size == 0:
            continue
        chosen = pix[rng.integers(0, pix.size, size=n_c)]
        py = chosen // field.grid
        px = chosen % field.grid
        # pixel centre + uniform sub-pixel jitter -> continuous uniform-in-cell
        cx = (px + 0.5) * dx + rng.uniform(-dx / 2, dx / 2, size=n_c)
        cy = (py + 0.5) * dx + rng.uniform(-dx / 2, dx / 2, size=n_c)
        true_xy_list.append(np.column_stack([cx, cy]))
        gene_list.append(genes_c)
        cell_list.append(np.full(n_c, c, dtype=np.int32))

    true_xy = np.concatenate(true_xy_list, axis=0)
    gene = np.concatenate(gene_list).astype(np.int32)
    true_cell = np.concatenate(cell_list)
    disp = rng.normal(0.0, field.sigma_um, size=true_xy.shape)
    obs_xy = true_xy + disp
    if displacement == "gauss_uniform" and disp_epsilon > 0:
        bg = rng.random(true_xy.shape[0]) < disp_epsilon
        nb = int(bg.sum())
        if nb:
            obs_xy[bg] = rng.uniform(0.0, field.L_um, size=(nb, 2))
    elif displacement not in ("gaussian", "gauss_uniform"):
        raise ValueError(f"unknown displacement {displacement!r}")
    true_type = field.types[true_cell]
    interior = field.interior_cell[true_cell]
    return Transcripts(obs_xy=obs_xy, true_xy=true_xy, gene=gene,
                       true_cell=true_cell, true_type=true_type, interior=interior)
