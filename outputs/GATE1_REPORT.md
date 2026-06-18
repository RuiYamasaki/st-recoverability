# Gate 1 report: does the frontier survive realistic assumptions?

Branch `gate1` (off `gate0`). All numbers trace to committed result files and columns.
Regenerate with `python src/gate1_run_all.py`; re-derive these numbers with
`python src/gate1_report.py`. Gate 0 results are untouched and still regenerate
byte-identically (`python src/run_all.py --skip-realism`). This report states raw
numbers only; it does not assert PASS or KILL.

Seeds: realistic-expression sweep base `config.MASTER_SEED + 100000`; leakage/anchor
base `+ 200000`; structural base `+ 300000`. Every per-row seed is in the CSVs.

## Task 1: realistic overlapping expression

Source of the model: per-type mean expression over the 155-gene MERFISH panel
(Moffitt 2018 mouse hypothalamic preoptic, squidpy `MERFISH_0.24.h5ad`), using the
dataset's native `cell_class` labels (10 types with >= 60 cells; not merged, so
confusable subtypes such as Endothelial 1 vs 3 are kept). Code: `src/expression.py`.
Overlap metrics (`results/gate1_expression_meta.json`):
- mean pairwise cosine similarity between type profiles: 0.210 (max 0.734, the
  Excitatory vs Inhibitory pair).
- fraction of expression mass on type-exclusive genes: 0.045; on shared genes: 0.955.
  (The Gate 0 disjoint model put 0.45 of mass on exclusive genes: the realistic model
  is dominated by shared expression, as real panels are.)
- 12 exclusive-marker genes (exclusivity threshold 0.7).

Oracle assignment-accuracy frontier under realistic expression
(`results/gate1_sweep.csv`, column `oracle_acc`; 5 x 5 x 5 grid, same axes as Gate 0):
- full-grid range: [0.251, 0.956]  (Gate 0 was [0.234, 0.957]).
- monotonicity-violation counts: 0 in sigma, 0 in packing.
- oracle >= naive at all 125 grid points: True (min `oracle_minus_naive` = +0.0115).
- vs Gate 0 at matched grid points: mean absolute difference in `oracle_acc` 0.0058,
  max 0.0319. The assignment frontier is essentially unchanged: it is not a
  disjoint-marker artifact.

Both adjusted Rand indices (`results/gate1_sweep.csv`):
- oracle expression-profile ARI vs true types (`profile_ari_vs_truetype`): range
  [0.667, 1.000], median 1.000. Gate 0 was [1.000, 1.000]. DE-SATURATED: the sanity
  check (profile ARI must fall below 1.0 under realistic expression) is satisfied.
  This metric is noisy (KMeans on 10 types including rare classes); the de-saturation
  is the relevant point, and it confirms the realistic model is genuinely overlapping.
- naive expression-profile ARI (`naive_profile_ari_vs_truetype`): range [-0.006, 1.000],
  median 0.743. A naive (segmentation-like) method loses cell typing as displacement
  grows; the Gate 0 saturated profile ARI of 1.0 was indeed a disjoint-marker artifact.

## Task 2: data-pinned displacement sigma

Marker leakage statistic (`src/gate1_leakage_anchor.py`). Two versions on the exclusive
markers (owner type tau), both computable on real and synthetic data without ground
truth:
- literal leakage: fraction of a marker's signal in non-owner-type cells. DIAGNOSTIC
  only: it is biology-dominated. Real (MERFISH section) = 0.402, which equals the
  synthetic co-expression baseline at sigma=0 (the synthetic emission uses the same
  per-type means), so it cannot identify displacement.
- spatial leakage (used to pin sigma): excess marker signal in non-owner cells ADJACENT
  to owner cells over the level in DISTANT non-owner cells, normalised by the owner
  level. Biological co-expression is spatially uniform and cancels in the near-minus-far
  contrast, so only displacement contributes; it rises monotonically with sigma.

Real spatial marker leakage on the MERFISH section (200 um sub-FOVs are too sparse for
the per-marker owner/adjacency split, so the whole imaged section is the field of view):
- mean 0.0768; per-marker 25th/75th percentiles 0.0177 / 0.0909
  (`results/gate1_leakage.csv`, columns `real_spatial_leakage_*`).
- synthetic baseline at sigma=0: 0.0167 (`synthetic_baseline_leakage_sigma0`).

Inverting the monotone synthetic leakage curve (built at MERFISH packing and density):
- data-pinned sigma = 2.67 um (`sigma_pinned_um`), status ok.
- 25-75 bracket = [0.10, 3.03] um (`sigma_bracket_lo_um`, `sigma_bracket_hi_um`).
- The measured leakage is reproducible within the plausible range [0, 15] um, so the
  Gate 1 threshold-2 KILL condition (leakage not reproducible by any plausible sigma) is
  not triggered.

Oracle assignment accuracy at the pinned physical sigma (point + bracket), evaluated at
each real anchor's packing and density (`results/gate1_leakage.csv`):
- MERFISH-like (sparse: packing 2575 cells/mm^2, density 286 tx/cell):
  oracle 0.836 at 2.67 um (naive 0.776); across the bracket [0.816, 0.989].
- Xenium-like (dense: packing 13625 cells/mm^2, density 74 tx/cell):
  oracle 0.643 at 2.67 um (naive 0.544); across the bracket [0.614, 0.976].
  Note: the displacement was measured on MERFISH (the Xenium cell-by-gene matrix was not
  downloaded), and the same physical sigma is applied at Xenium's geometry; this is
  recorded in the column `leakage_source`. Platform-specific sigma could differ.

Both real-anchored operating points sit below 0.9, with the oracle above the naive
baseline; the limit bites in the real-data-matched regime and bites hardest in the dense
Xenium-like regime. See `figures/gate1_realistic_frontier.png`.

## Task 3: structural-assumption sensitivity (a check, not the primary decision)

Reduced grid (packing in {2575, 6000, 13625} cells/mm^2; sigma in {1, 2, 4, 8} um;
density 160 tx/cell), realistic expression, three conditions
(`results/gate1_structural.csv`). Code: `src/gate1_structural.py`.
- baseline (Voronoi cells, Gaussian displacement): oracle range [0.290, 0.925];
  monotonicity violations (sigma, packing) = (0, 0); oracle >= naive everywhere (min gap
  +0.020); dense-regime oracle at sigma 2 / sigma 8 = 0.740 / 0.290.
- aniso (non-Voronoi: irregular elongated cells via anisotropic nucleus expansion,
  Gaussian displacement): range [0.281, 0.924]; violations (0, 0); oracle >= naive (min
  gap +0.121); dense 0.712 / 0.281.
- mixture (Voronoi cells, non-Gaussian heavy-tailed displacement: 10% of transcripts land
  uniformly in the domain, oracle made Bayes-optimal for the mixture): range
  [0.254, 0.833]; violations (0, 0); oracle >= naive (min gap +0.019); dense 0.660 / 0.254.

Under both perturbations the frontier stays monotone, the oracle stays at or above naive,
and the limit still bites in the dense regime (dense sigma-2 accuracy 0.66 to 0.74, all
below 0.9). The mixture caps the easy end at 0.833 because 10% of transcripts are
unassignable by construction. None of the structural-sensitivity KILL conditions
(monotonicity lost, or dense-regime bite flipped) is triggered. See
`figures/gate1_structural_sensitivity.png`.

## Committed paths and regen

- Expression model: `src/expression.py` (+ `results/gate1_expression_meta.json`).
- Generator / oracle / metrics extended in place (model + geometry + displacement
  options; Gate 0 defaults preserved byte-identically): `src/generator.py`,
  `src/oracle.py`, `src/metrics.py`, `src/evaluate.py`.
- Task 1: `src/gate1_sweep.py` -> `results/gate1_sweep.csv`.
- Task 2: `src/gate1_leakage_anchor.py` -> `results/gate1_leakage.csv`.
- Task 3: `src/gate1_structural.py` -> `results/gate1_structural.csv`.
- Figures: `src/gate1_figures.py` -> `figures/gate1_realistic_frontier.png`,
  `figures/gate1_structural_sensitivity.png`.
- Regen all: `python src/gate1_run_all.py`. Re-derive these numbers:
  `python src/gate1_report.py`.

## Honesty notes and limitations

- The displacement sigma is measured on MERFISH only. The Xenium leakage was not
  measured independently (its cell-by-gene matrix was not downloaded); the same physical
  sigma is applied at Xenium's geometry, and platform-specific sigma could differ.
- The spatial-leakage statistic uses adjacency thresholds (1.5 and 4.0 median
  nearest-neighbour distances) and the dataset's cell-type labels; the 25-75 bracket
  reflects the per-marker spread, not all sources of uncertainty.
- The oracle profile-ARI is a noisy metric here (KMeans on 10 real types with rare
  classes). The robust, density-independent decision metric is the oracle assignment
  accuracy; the profile ARI is reported to show de-saturation, not as a precise surface.
- The MERFISH X matrix is volume-normalised expression, so leakage is an expression-mass
  fraction rather than a raw transcript-count fraction.
