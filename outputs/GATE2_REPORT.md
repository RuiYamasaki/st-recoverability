# Gate 2 report: pin the displacement sigma directly on Xenium, with error bars

Branch `gate2` (off `gate1`). All numbers trace to committed result files and columns.
Regenerate with `python src/gate2_run_all.py`; re-derive these numbers with
`python src/gate2_report.py`. Gate 0 and Gate 1 are untouched and still regenerate
byte-identically. This report states raw numbers only; it does not assert PASS or KILL.

Seeds: Task 1 base `MASTER_SEED + 400000`; Task 2 base `+ 410000`; Task 3 base
`+ 420000`; per-row seeds are in the CSVs. The Xenium clustering uses
`MASTER_SEED = 20260618`.

## Task 1: sigma pinned directly on Xenium (no cross-platform borrowing)

Dataset: 10x Genomics "Xenium FFPE Human Breast Cancer, Replicate 1" (release 1.0.1).
Files (direct HTTP, no auth): `..._cell_feature_matrix.h5` (12.1 MB, 313 gene-expression
features x 167,780 cells) and `..._cells.parquet` (centroids), both under
`cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/`. Code:
`src/expression.py` (`build_realistic_model_from_xenium`), `src/gate2_pin.py`.

Cell types: MiniBatchKMeans (seed = MASTER_SEED) on median-normalised log1p counts,
clusters below 50 cells dropped, giving K=14 types over the 313-gene panel. Exclusive
markers as in Gate 1 (exclusivity cutoff 0.7): 8 markers, owned by endothelial
(CLEC14A, RAMP2, HOXD9), epithelial (KRT23, PIGR) and plasma-B (MZB1, DERL3, TNFRSF17)
clusters. Section median nearest-neighbour distance 9.43 um; representative (median
local) packing 6,534 cells/mm^2; median 164 transcripts/cell.

Leakage and pin (`results/gate2_xenium_pin.csv`):
- real distance-stratified spatial marker leakage on the whole section: 0.0679
  (`real_spatial_leakage_mean`). The literal "fraction in other-type cells" diagnostic is
  0.408 (biology-dominated, as in Gate 1).
- Xenium-pinned displacement sigma: 1.61 um (`xenium_pinned_sigma_um`), versus the
  MERFISH-pinned 2.67 um. The smaller Xenium value is consistent with Xenium's finer
  reported localisation.
- oracle assignment accuracy at the Xenium-pinned sigma (converged grid): dense (packing
  13,625 cells/mm^2) 0.784, naive 0.703; sparse (packing 2,575) 0.906. Both below 0.95;
  dense below 0.9.

## Task 2: statistical CI and design-sensitivity on sigma

Source: `results/gate2_sigma_uncertainty.csv`, `results/gate2_design_grid.csv`. Code:
`src/gate2_uncertainty.py`. Dense-regime oracle accuracy is reported on a converged grid
(`res_cell` 24, `grid_max` 1500); a resolution check is in the honesty notes.

Statistical bootstrap (B = 2000 resamples over markers and over the near/far cell sets):
- leakage point 0.0679, 95% CI [0.0340, 0.1099].
- propagated to sigma: point 1.61 um, statistical 95% CI [0.79, 2.53] um.

Design-sensitivity grid (27 choices: adjacency R_near in {1.0, 1.5, 2.0}, far R_far in
{3, 4, 5} median-NN units, marker cutoff in {0.6, 0.7, 0.8}):
- pinned sigma range [0.20, 2.48] um. The spread is driven almost entirely by the marker
  cutoff (`results/gate2_design_grid.csv`):
  - cutoff 0.7 (8 markers, the pre-registered Gate 1 choice): sigma [1.49, 1.70] um, stable.
  - cutoff 0.6 (20 looser markers): sigma [0.20, 2.48] um; the 0.20 um value comes from the
    loosest markers with the tightest adjacency bin (R_near 1.0).
  - cutoff 0.8: 0 markers (too strict for the 313-gene panel), degenerate and excluded.

Combined band [min, max] = [0.20, 2.53] um.

Dense-regime oracle assignment accuracy (converged grid):
- at the point sigma 1.61 um: 0.787.
- at the statistical-CI ends: 0.891 at the low sigma (0.79 um), 0.675 at the high sigma
  (2.53 um). The statistical 95% CI keeps dense accuracy in [0.675, 0.891], entirely
  below 0.95.
- at the combined-band ends: 0.972 at the optimistic low sigma (0.20 um), 0.672 at the
  high sigma (2.53 um). The optimistic end of the combined band, reached only at the
  loose-marker design corner, exceeds 0.95.
Sparse-regime oracle accuracy: 0.905 at the point; combined ends [0.849, 0.985].

Dense oracle accuracy crosses 0.95 at about sigma 0.35 um (`figures/gate2_dense_vs_sigma.png`).

## Task 3: non-idealized negative-binomial emission

Source: `results/gate2_nbinom.csv`, `results/gate2_nbinom_frontier.csv`. Code:
`src/gate2_nbinom.py`, generator `emission='nbinom'`. Each cell's per-gene counts are
Poisson(mu * f_g) with a per-cell, per-gene gamma factor of mean 1 and variance phi_g,
phi_g the per-gene dispersion estimated by within-type method of moments from the real
Xenium counts (Var = mu + phi*mu^2). The real leakage is unchanged; only the synthetic
emission changes, so sigma is re-pinned and the oracle re-evaluated.
- NB dispersion phi median 4.191.
- re-pinned sigma: Poisson 1.57 um -> negative-binomial 1.67 um (shift +0.10 um).
- dense oracle accuracy: Poisson 0.795 -> negative-binomial 0.779 (shift -0.016); naive
  0.695. Sparse: 0.905 -> 0.904.
- frontier under NB emission (reduced grid): oracle range [0.267, 0.930];
  monotonicity-violation counts (sigma, packing) = (0, 0); oracle >= naive at every point
  (minimum gap +0.025). The frontier does not flatten and the dense-regime bite is
  preserved.

## Committed paths and regen

- Xenium model + dispersion: `src/expression.py`.
- Task 1: `src/gate2_pin.py` -> `results/gate2_xenium_pin.csv`.
- Task 2: `src/gate2_uncertainty.py` -> `results/gate2_sigma_uncertainty.csv`,
  `results/gate2_design_grid.csv`.
- Task 3: `src/gate2_nbinom.py` -> `results/gate2_nbinom.csv`,
  `results/gate2_nbinom_frontier.csv`.
- Generator extended in place (`emission`, `res_cell`/`grid_max` options); Gate 0/1
  defaults preserved byte-identically.
- Figures: `src/gate2_figures.py` -> `figures/gate2_dense_vs_sigma.png`,
  `figures/gate2_design_sigma.png`.
- Regen all: `python src/gate2_run_all.py`. Re-derive numbers: `python src/gate2_report.py`.

## Honesty notes and limitations

- Where the dense-regime accuracy band sits relative to the thresholds: across the
  statistical 95% CI the dense accuracy is [0.675, 0.891], all below 0.95 and mostly below
  0.9. The only way the dense accuracy reaches above 0.95 (0.972) is the combined-band
  optimistic end at sigma 0.20 um, which arises solely from the loosest marker cutoff
  (0.6) paired with the tightest adjacency bin; at the pre-registered cutoff (0.7) the
  pinned sigma is stable (1.49 to 1.70 um) and the dense accuracy is 0.78 to 0.79.
- Sigma stability: the combined sigma band spans about 13x (0.20 to 2.53 um), but the
  instability is isolated to the loosest marker cutoff; a defensible central value exists
  (1.61 um at cutoff 0.7, statistical CI [0.79, 2.53]).
- Resolution: at small sigma the default oracle grid under-resolves and UNDER-estimates
  accuracy. A resolution check (res_cell 6, 12, 24, 48) at packing 13,625 showed the
  sigma 0.20 um accuracy rising 0.952 -> 0.965 -> 0.969 -> 0.972 (converging near 0.972),
  and the sigma 1.61 um accuracy stable at 0.78. Dense-regime accuracies in this report
  use the converged grid (res_cell 24, grid_max 1500).
- Sigma is a physical micron quantity pinned at Xenium's representative (median local)
  packing and applied at the dense packing; the real cells span a range of local packings
  and the synthetic field is a regular jittered lattice, so the absolute micron value
  carries a modelling uncertainty that the design-sensitivity grid partly probes.
- The MERFISH and Xenium X matrices differ (volume-normalised vs raw integer counts); the
  Xenium NB dispersion is estimated from the raw integer counts.
