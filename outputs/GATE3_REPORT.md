# Gate 3 report: a principled marker-validity test and independent replication

Branch `gate3` (off `gate2`). All numbers trace to committed result files and columns.
Regenerate with `python src/gate3_run_all.py`; re-derive these numbers with
`python src/gate3_report.py`. Gates 0, 1 and 2 are untouched and still regenerate
byte-identically. This report states raw numbers only; it does not assert PASS or KILL.

Seeds: Task 1 base `MASTER_SEED + 500000`; Task 2 base `+ 510000`; Task 3 base
`+ 520000`; clustering seed `MASTER_SEED = 20260618`. Per-row seeds are in the CSVs.

## The validity test (non-circular by construction)

Code: `src/gate3_common.py` (`marker_validity`). For each candidate marker (owner type
tau), the real statistic is the adjacent-minus-distant displacement signal
(mean gene count in non-owner cells within 1.5 median-NN of an owner cell, minus the mean
in non-owner cells beyond 4 median-NN, divided by the owner-cell mean). The null shuffles
that gene's per-cell counts across the non-owner cells (a position permutation relative to
counts) and recomputes the statistic on the fixed real adjacency geometry, 1000 times. A
marker is admissible only if its real signal exceeds the null at p < 0.05. Marker
selection never uses the sigma the markers yield; the candidate pool is broad (exclusivity
>= 0.5), and admissibility is decided purely by displacement-signal significance.

## Task 1: breast dataset (source: results/gate3_validity_breast.csv, gate3_pin_breast.csv)

- 38 candidate markers; 28 admissible at p < 0.05 (1000 permutations).
- The test separates markers by real displacement signal, not by exclusivity: mean signal
  of admissible markers +0.094 versus rejected -0.055 (rejected markers have no positive
  near-owner gradient).
- Exclusive markers (exclusivity >= 0.7): 7 of 8 admissible (the failure is MZB1, signal
  0.002, p = 0.104). Loose markers (exclusivity in [0.6, 0.7), the set added when the
  cutoff was loosened to 0.6 in Gate 2): 7 of 12 admissible. The 5 rejected loose markers
  are exactly the sigma-destabilisers, and they fail on independent grounds (signal
  indistinguishable from the permutation null): AHSP p = 0.599, CRHBP p = 0.242,
  CYP1A1 p = 0.996, EGFL7 p = 0.642, SLC4A1 p = 0.454. The previously destabilising
  markers therefore fail the validity test.
- Re-pinned sigma on the admissible set: 1.99 um, bootstrap 95% CI [1.43, 2.63]
  (B = 2000 over markers and cells).
- Dense-regime oracle assignment accuracy (converged grid, packing 13,625 cells/mm^2):
  point 0.739 (naive 0.642, oracle >= naive); across the sigma CI [0.663, 0.805]. The
  whole CI is below 0.95 and below 0.9. Sparse-regime point 0.880. Dense accuracy is
  monotone decreasing across the sigma CI (3-point check, 0 violations).

This resolves the Gate 2 breach: the 0.20 um corner that pushed dense accuracy to 0.972
arose from biology-diluted markers that fail the principled validity test; on the
admissible set the dense limit holds across the bootstrap CI.

## Task 2: sigma robustness to cell-typing (source: results/gate3_typing.csv)

Cell types re-derived under several reasonable clustering choices (k and clustering seed),
admissible markers re-selected by the validity test each time, sigma re-pinned:

| variant | k | admissible | sigma (um) | sigma CI (um) | dense oracle (point) | dense at CI-optimistic |
|---|---|---|---|---|---|---|
| kmeans_K10 | 10 | 58 | 2.65 | [2.27, 3.28] | 0.654 | 0.694 |
| kmeans_K14 | 14 | 27 | 2.24 | [1.76, 2.63] | 0.705 | 0.762 |
| kmeans_K14_altseed | 14 | 43 | 1.99 | [1.60, 2.29] | 0.740 | 0.785 |
| kmeans_K20 | 20 | 21 | 2.13 | [1.53, 2.76] | 0.728 | 0.803 |

- sigma spread across typing choices: [1.99, 2.65] um (a factor of about 1.3, versus the
  old marker-cutoff lever that swung sigma over 0.20 to 2.48 um, a factor of about 12).
- dense-accuracy (point) spread: [0.654, 0.740]; maximum dense accuracy at any
  sigma-CI optimistic end: 0.803. All below 0.95 and below 0.9. oracle >= naive at every
  evaluated point; dense accuracy monotone decreasing across each sigma CI (0 violations).

Unlike the old marker-cutoff lever (which swung sigma over 0.20 to 2.48 um), sigma is
stable across cell-typing choices under the validity test.

## Task 3: independent Xenium replication, lung cancer (source: results/gate3_validity_lung.csv, gate3_pin_lung.csv)

Dataset: 10x Genomics "Human Lung Cancer Preview Data, Xenium Human Multi-Tissue and
Cancer Panel" (sample Xenium_V1_hLung_cancer_section, release 1.5.0; 150,365 cells; 541
features). Files (direct HTTP, no auth): cell_feature_matrix.h5 and cells.parquet under
`cf.10xgenomics.com/samples/xenium/1.5.0/Xenium_V1_hLung_cancer_section/`. Code:
`src/gate3_replicate.py`. Lung geometry: median NN 8.79 um, median local packing 7,285
cells/mm^2, dense tail (p90) 14,635 cells/mm^2 (so the dense packing 13,625 used for the
oracle is realistic for this tissue too), median 132 transcripts/cell.

- 32 candidate markers; 31 admissible. All 15 exclusive (>= 0.7) and all 5 loose
  [0.6, 0.7) markers pass (p <= 0.002): the lung candidate pool contains no biology-diluted
  markers.
- Re-pinned sigma on the admissible set: 1.83 um, bootstrap 95% CI [1.62, 2.08].
- Dense-regime oracle accuracy: point 0.767 (naive 0.672, oracle >= naive); across the
  sigma CI [0.746, 0.796]. The whole CI is below 0.95 and below 0.9. Sparse point 0.899.
  Monotone decreasing across the sigma CI (0 violations).

The dense limit holds on an independent tissue: dense oracle accuracy stays below 0.95
across the sigma CI on both breast and lung.

## Committed paths and regen

- Validity test + pinning: `src/gate3_common.py`.
- Task 1: `src/gate3_validity.py` -> `results/gate3_validity_breast.csv`,
  `results/gate3_pin_breast.csv`.
- Task 2: `src/gate3_typing.py` -> `results/gate3_typing.csv`.
- Task 3: `src/gate3_replicate.py` -> `results/gate3_validity_lung.csv`,
  `results/gate3_pin_lung.csv`.
- Xenium model builder extended in place (dataset path + name arguments); Gates 0 to 2
  defaults preserved byte-identically.
- Figures: `src/gate3_figures.py` -> `figures/gate3_validity_separation.png`,
  `figures/gate3_dense_limit_forest.png`.
- Regen all: `python src/gate3_run_all.py`. Re-derive numbers: `python src/gate3_report.py`.

## Honesty notes and limitations

- The validity test admits markers with a real spatial near-owner gradient. That gradient
  is consistent with transcript displacement but a position permutation alone cannot fully
  separate displacement from biological co-localisation of cell types; the test removes the
  biology-diluted (no-gradient) markers, which is what the sigma stability required.
- The dense packing (13,625 cells/mm^2) is the dense-tissue condition; it is at the dense
  tail of both breast and lung. The sparse and median regimes show a much weaker limit.
- Dense-regime oracle accuracies use the converged grid (res_cell 24), as in Gate 2.
- Cell types are clustering-derived; a matched scRNA-seq reference annotation was not
  pursued (the clustering-choice robustness in Task 2 is the typing-sensitivity evidence).
