# Gate 0 report: does a transferable recoverability frontier exist?

Branch: `gate0`. All numbers below trace to committed result files and columns.
Regenerate everything with `python src/run_all.py` (first run streams ~7 MB of
public data; `python src/run_all.py --skip-realism` reuses the cached real-data
summary). Re-derive the table and headline scalars with `python src/report_table.py`.

This report states raw numbers only. It does not assert a PASS or KILL conclusion;
the locked thresholds are judged by the orchestrator.

## 1. Generator parameterization and ranges swept

Units are micrometres throughout, so the axes map directly onto real summary
statistics. Code: `src/generator.py`, fixed model in `src/config.py`.

- Master seed: `20260618` (`config.MASTER_SEED`). Every field/transcript seed is
  derived deterministically (`config.config_seed(packing_idx, sigma_idx)`; the
  transcript draw adds a per-density offset). All seeds are recorded per row in
  `results/sweep.csv` (columns `field_seed`, `transcript_seed`).
- Cells: a perturbed square lattice (jitter 0.33 of lattice spacing), N = 400 cells
  per field held constant so resolution and statistical power are uniform. True
  boundaries are the Voronoi tessellation of the centres (a pixel's Voronoi cell is
  its nearest centre, so boundaries are exact on the label grid).
- Cell types: K = 4, fixed proportions [0.40, 0.30, 0.20, 0.10].
- Genes: G = 20. Each type owns a disjoint 4-gene marker block; 55% of a type's
  expression mass sits on its markers, the rest on a shared background
  (`config.TYPE_GENE_COMPOSITION`). This single model is fixed across the whole
  sweep and was set before any oracle output was seen.
- Transcripts: cell of type t emits Poisson(mean_tx * p[t, g]) per gene; each
  transcript is placed uniformly inside its cell's Voronoi region, then displaced by
  isotropic Gaussian N(0, sigma^2 I). True transcript-to-cell assignment is recorded.

Three swept axes, 5 levels each (5 x 5 x 5 = 125 grid points), chosen to bracket and
extend the real-data anchors measured in step 4:

| axis | parameter | levels |
|---|---|---|
| density | mean transcripts/cell | 50, 90, 160, 280, 500 |
| packing | cells per mm^2 | 2000, 3400, 5700, 9500, 16000 |
| displacement | sigma (um) | 0.5, 1.0, 2.0, 4.0, 8.0 |

The oracle assignment map depends only on (packing, sigma), so the field and oracle
maps are built once per (packing, sigma) and reused across densities. This is also
verified empirically: oracle assignment accuracy is density-independent (max spread
across the 5 densities at fixed packing and sigma is 0.015; mean 0.005;
`results/sweep.csv`, column `oracle_acc`).

## 2. Oracle computation and the oracle >= naive check

Code: `src/oracle.py`. Under the known generative model the joint density factorises
as P(c, g, x) = P(c) P(g | type(c)) P(x | c), with:
- P(c) uniform (every cell has the same expected count, so the source prior cancels);
- P(g | type(c)) = p[type(c), g], the fixed expression composition;
- P(x | c) proportional to (f_c * phi_sigma)(x): the cell's uniform-in-Voronoi
  density convolved with the Gaussian displacement kernel.

The Bayes-optimal (minimum-error) assignment is therefore
c*(x, g) = argmax_c p[type(c), g] (f_c * phi_sigma)(x). No estimator can beat its
expected accuracy: it is the information-theoretic ceiling for this model. It is
computed exactly in the pixel limit: each cell's indicator mask is Gaussian-blurred
and normalised by area to give (f_c * phi_sigma); because p[t, g] is constant within
a type, the per-cell argmax reduces to a per-type maximum plus the gene-evidence term.

Naive baseline (`naive_assign`): assign each transcript to the nearest cell centre of
its observed position (ignores the displacement model and the gene identity).

Confirmation: oracle assignment accuracy >= naive accuracy at every one of the 125
grid points (`results/sweep.csv`, column `oracle_ge_naive` all True; minimum of
`oracle_minus_naive` is +0.0061). The oracle also dominates at every realism anchor
row (`results/realism_oracle.csv`, column `oracle_ge_naive`).

Two distinct adjusted Rand indices are reported, because the task text and the locked
threshold refer to different objects:
- `assign_ari_vs_truecells_oracle`: ARI of the recovered transcript-to-cell partition
  vs the true partition ("ARI vs true cell labels"). It tracks assignment quality.
- `profile_ari_vs_truetype`: ARI of the recovered cell-level expression-profile
  clustering vs true cell types ("expression profiles vs true profiles").

## 3. Representative grid points

Source for every number: `results/sweep.csv`. Density is fixed at 160 tx/cell here
(assignment accuracy is density-independent, so this is representative). Column names:
oracle assignment accuracy = `oracle_acc`; oracle adjusted Rand (assignment partition
vs true cell labels) = `assign_ari_vs_truecells_oracle`; oracle adjusted Rand
(expression-profile clustering vs true types) = `profile_ari_vs_truetype`;
naive-baseline accuracy = `naive_acc`.

| regime | packing (cells/mm^2) | sigma (um) | density (tx/cell) | oracle assignment accuracy | oracle ARI (assignment vs true cell labels) | oracle ARI (profiles vs true types) | naive accuracy |
|---|---|---|---|---|---|---|---|
| low-difficulty extreme | 2000 | 0.5 | 160 | 0.956 | 0.917 | 1.000 | 0.941 |
| MERFISH-like packing, realistic sigma | 2000 | 2.0 | 160 | 0.888 | 0.797 | 1.000 | 0.851 |
| mid | 5700 | 4.0 | 160 | 0.651 | 0.459 | 1.000 | 0.554 |
| mid | 9500 | 2.0 | 160 | 0.770 | 0.613 | 1.000 | 0.696 |
| Xenium-like packing, realistic sigma | 16000 | 2.0 | 160 | 0.708 | 0.528 | 1.000 | 0.614 |
| high-difficulty extreme | 16000 | 8.0 | 160 | 0.239 | 0.111 | 1.000 | 0.128 |

Full-grid ranges (`results/sweep.csv`): oracle assignment accuracy [0.234, 0.957];
oracle assignment-ARI [0.110, 0.919]; oracle profile-ARI [1.000, 1.000]; naive
profile-ARI [0.102, 1.000]. Assignment accuracy is monotonically non-increasing in
sigma (0 violations across 25 curves) and in packing (0 violations across 25 curves).

Note on the two ARIs (a finding, not an artifact): the oracle's errors are
type-preserving. It routes each transcript to a cell whose type matches the
transcript's gene, so cell-level expression profiles stay type-pure even when exact
cell identity is lost. The naive baseline's errors are type-random (its profile-ARI
collapses to 0.10). The recoverability limit therefore bites hard on exact
cell-identity assignment but barely on cell-type recovery. See
`figures/identity_vs_typing.png`.

## 4. Real-data summary statistics, dataset IDs, access route, realistic band

Source: `results/realism.csv` (per-FOV stats) and `results/realism_meta.json`
(dataset IDs and access routes). Code: `src/realism.py`. FOV = 200 x 200 um window;
nearest-neighbour distance computed over the whole section then medianed within each
FOV; FOV placement is deterministic (seed 20260618).

| dataset | FOV | n cells | median tx/cell | packing (cells/mm^2) | median NN distance (um) |
|---|---|---|---|---|---|
| Xenium breast | fov0 | 565 | 86.0 | 14125 | 6.68 |
| Xenium breast | fov1 | 540 | 71.0 | 13500 | 6.37 |
| Xenium breast | fov2 | 530 | 66.0 | 13250 | 6.55 |
| MERFISH hypothalamus | fov0 | 108 | 266.8 | 2700 | 10.72 |
| MERFISH hypothalamus | fov1 | 101 | 311.9 | 2525 | 12.52 |
| MERFISH hypothalamus | fov2 | 100 | 280.3 | 2500 | 11.96 |

Dataset IDs and access:
- Xenium: "Xenium FFPE Human Breast Cancer, Replicate 1" (10x Genomics release 1.0.1).
  Direct HTTP, no auth: the standalone `..._cells.parquet`
  (`https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_cells.parquet`).
  transcripts/cell field: `transcript_counts`.
- MERFISH: Moffitt et al. 2018 mouse hypothalamic preoptic region, the
  squidpy-packaged AnnData `MERFISH_0.24.h5ad` on figshare
  (`https://ndownloader.figshare.com/files/40038538`). Direct HTTP, no auth.
  transcripts/cell field: row sum of the cell-by-gene matrix. (Vizgen's own showcase
  CSVs are gsutil-gated and return HTTP 403 to anonymous callers, so the figshare
  mirror was used; this is recorded in `results/realism_meta.json`.)

Realistic band on the generator axes: the region matching a real FOV within 15% on
both median tx/cell and packing, evaluated at the realistic displacement sub-band
sigma in [0.5, 2.0] um. Displacement is the one axis not observable from the public
per-cell tables, so the sigma band is an explicit assumption; the full sigma sweep is
also reported (`results/realism_oracle.csv`). The band has two clusters:
- Xenium-like: packing 13250 to 14125 cells/mm^2, median tx/cell 66 to 86.
- MERFISH-like: packing 2500 to 2700 cells/mm^2, median tx/cell 267 to 312.
On the swept grid the nearest cells are packing 16000 (Xenium) and packing 2000
(MERFISH), at sigma in {0.5, 1.0, 2.0}; the realism-anchored evaluation
(`results/realism_oracle.csv`) uses the exact measured packing and density of each
FOV rather than the nearest grid cell.

## 5. Oracle accuracy inside the realistic band

Source: `results/realism_oracle.csv`, rows with `in_realistic_sigma_band == True`
(sigma <= 2.0 um), column `oracle_acc`.
- Minimum oracle assignment accuracy in the realistic band: 0.719.
- Maximum oracle assignment accuracy in the realistic band: 0.953.
- Xenium-like (dense) FOVs in band: [0.719, 0.922].
- MERFISH-like (sparse) FOVs in band: [0.866, 0.953].
- Oracle assignment-ARI in band: [0.545, 0.910]. Oracle profile-ARI in band:
  [1.000, 1.000].

For reference, extending past the realistic sigma band (sigma = 4 and 8 um) drives the
oracle to 0.26 to 0.50 (Xenium-like) and 0.55 to 0.76 (MERFISH-like).

## 6. Can the generator match real median tx/cell and packing within 15%?

Yes, at every anchor (`results/realism_oracle.csv`, columns `tx_match_rel_err`,
`packing_match_rel_err`, `tx_match_within_15pct`, `packing_match_within_15pct`).
- Packing match: exact by construction (max relative error 1.4e-16). The generator's
  nominal packing equals the FOV packing; with N = 400 cells the domain side is set to
  sqrt(N / packing).
- Median tx/cell match: max relative error 0.0152 (1.52%) across all anchors. The
  generator's mean_tx_per_cell is set to the FOV's median tx/cell; the realised median
  over interior cells (`realized_median_tx_per_cell`) matches within 1.6%.
Matching parameter values are exactly the measured FOV statistics in the table above
(for example Xenium fov1: packing 13500 cells/mm^2, density 71 tx/cell;
MERFISH fov0: packing 2700 cells/mm^2, density 267 tx/cell).

## 7. Committed file paths and regen command

- Generator code: `src/generator.py` (fixed model `src/config.py`).
- Oracle and metrics code: `src/oracle.py`, `src/metrics.py`, `src/evaluate.py`.
- Sweep code and CSV: `src/sweep.py` -> `results/sweep.csv`.
- Realism code and CSVs: `src/realism.py` -> `results/realism.csv`,
  `results/realism_meta.json`; `src/anchor.py` -> `results/realism_oracle.csv`.
- Figures: `src/figures.py` -> `figures/oracle_accuracy_surface.png`,
  `figures/identity_vs_typing.png`.
- Regen everything: `python src/run_all.py` (or `--skip-realism` to reuse the cached
  real-data summary). Re-derive the table and scalars: `python src/report_table.py`.

## 8. Honesty notes and limitations

- The pixel oracle is exact in the continuum limit. In the high-accuracy corner (small
  sigma relative to large cells, the low-packing low-sigma region), the displacement
  is under-resolved by the grid, so oracle accuracy there is mildly optimistic. This
  corner is far from the decision-relevant mid and low band, where sigma is well
  resolved.
- The realistic sigma band [0.5, 2.0] um is an assumption (displacement is not
  observable from the public per-cell tables). The full sigma sweep is reported so the
  sensitivity is visible.
- The gene model is a fixed modelling choice (4 disjoint marker blocks, 55% marker
  mass). It was not tuned to clear any threshold. Its main consequence is that the
  oracle profile-ARI saturates at 1.0; the assignment-accuracy frontier is the robust,
  density-independent decision signal.
