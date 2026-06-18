# Segmenter-headroom report: real methods vs the oracle ceiling

Branch `headroom` (off `gate3`). All numbers trace to committed result files and columns.
Regenerate with `python src/headroom_run_all.py`; re-derive the table with
`python src/headroom_report.py`. Gates 0 to 3 are untouched and still regenerate
byte-identically. This report states raw numbers only; it does not assert whether methods
are near the wall or have headroom. The orchestrator interprets the gap.

Seeds: data export base `MASTER_SEED + 600000`; ComSeg internal sampling seeded
(np.random.seed(0), random.seed(0), ComSegDict seed 0); pciSeq variational inference is
deterministic given its input. Per-config seeds are in results/headroom_oracle_naive.csv.

## Metric

Matched transcript-assignment accuracy against the known truth: the fraction of interior
transcripts whose method-cell is matched to their true cell under the optimal one-to-one
matching of method-cells to true-cells (Hungarian on the contingency table); background or
unassigned transcripts count as errors. For the oracle and the naive baseline (which use
the true cell ids) this equals the direct id-match accuracy. Code: `src/headroom_common.py`
(`matched_accuracy`).

## Configuration

Most realistic Gate 3 config: Xenium-breast realistic overlapping expression, negative-
binomial within-type emission, data-pinned displacement at the Gate 3 breast value 1.99 um
and its CI ends 1.43 and 2.63 um, dense packing 13,625 cells/mm^2 (plus sparse 2,575 and
representative 6,534 for context), density 164 tx/cell, 250-cell fields (the oracle
accuracy is field-size-independent, so the ceiling is unchanged; the field is kept modest
so the external methods are tractable). Code: `src/headroom_export.py`.

## Which methods ran (versions, install outcomes, substitutions)

- pciSeq 0.0.61 (Qian et al., Nature Methods 2020): installed via pip and RAN on Python
  3.13. Nuclei-prior mode (prior = a nucleus label image of 3 um disks around the true
  centres; reference = the generator's true per-type mean expression).
- ComSeg 1.8.5 (Defard et al. 2024): installed via pip and RAN on Python 3.13. Nuclei-prior
  mode (prior = per-transcript in_nucleus membership from the same nuclei, plus the true
  centroids as landmarks). Free-segmentation mode (community_detection without prior) was
  attempted and FAILED with an AssertionError (recorded in results/headroom_methods.csv,
  column status/error).
- Proseg (Newell lab, Nature Methods 2025): NOT installable on this machine. It is a Rust
  CLI with no prebuilt binaries and no bioconda win-64 build; `cargo install` needs a Rust
  toolchain plus an MSVC or GNU C linker, none of which are present (no cl.exe, no gcc).
  Substituted by ComSeg + pciSeq (the closest installable published assignment methods).
- Baysor (Petukhov et al., Nature Biotechnology 2022): NOT installable on this machine. It
  is a Julia tool (Julia 1.10 plus a C compiler) with no prebuilt Windows binary; neither
  Julia nor a C compiler is present. Substituted by ComSeg + pciSeq.

No method was silently dropped: the two free-segmentation tools the prompt named (Proseg,
Baysor) are toolchain-blocked here and documented; the two running methods are real
published nuclei-prior assignment tools; the one attempted free-segmentation run (ComSeg
without prior) is recorded as failed.

## Table: dense operating point (and the sigma CI)

Source: oracle and naive from results/headroom_oracle_naive.csv (columns
oracle_acc_matched, naive_acc_matched); pciSeq and ComSeg from results/headroom_methods.csv
(column matched_accuracy, rows method in {pciSeq, ComSeg}, mode nuclei_prior). All values
are matched transcript-assignment accuracy vs the known truth.

| config (packing, sigma) | oracle | naive | pciSeq (nuclei) | ComSeg (nuclei) | ComSeg (free) |
|---|---|---|---|---|---|
| dense, sigma 1.43 um | 0.803 | 0.731 | 0.667 | 0.722 | failed |
| dense, sigma 1.99 um | 0.730 | 0.639 | 0.604 | 0.634 | failed |
| dense, sigma 2.63 um | 0.651 | 0.546 | 0.504 | 0.543 | failed |
| sparse, sigma 1.99 um | 0.882 | 0.833 | 0.230 | 0.762 | not run |
| representative, sigma 1.99 um | 0.808 | 0.741 | 0.475 | 0.730 | not run |

Notes (frac_assigned, results/headroom_methods.csv): ComSeg assigns ~0.99 of transcripts;
pciSeq leaves a substantial fraction as background (0.83 to 0.86 assigned in the dense
regime, only 0.24 assigned in the sparse regime where the small nuclei sit inside large
cells), which is why pciSeq's accuracy is lowest in the sparse regime.

## Best-real-method-to-oracle gap

- Best of the two sophisticated published methods (ComSeg, which exceeds pciSeq at every
  config): oracle minus best-sophisticated = 0.081, 0.096, 0.108 across the dense sigma CI
  (1.43, 1.99, 2.63 um); 0.120 sparse; 0.078 representative.
- Best of all real options including the naive nearest-nucleus baseline (naive exceeds both
  sophisticated methods at every config): oracle minus best-real = 0.072, 0.091, 0.105
  across the dense sigma CI; 0.048 sparse; 0.067 representative.
- Neither sophisticated method (pciSeq, ComSeg) exceeds the naive nearest-nucleus baseline
  on this known-truth synthetic; ComSeg is within a few points of naive in the dense
  regime, pciSeq is below it everywhere.

See figures/headroom_methods_vs_oracle.png and figures/headroom_dense_vs_sigma.png.

## Committed paths and regen

- Export + scorer: `src/headroom_common.py`, `src/headroom_export.py` ->
  results/headroom_oracle_naive.csv (raw transcript/cell tables under data/headroom/,
  gitignored).
- Method runners: `src/methods_pciseq.py`, `src/methods_comseg.py`; driver
  `src/headroom_methods.py` -> results/headroom_methods.csv.
- Report and figures: `src/headroom_report.py`, `src/headroom_figures.py` ->
  figures/headroom_*.png.
- Regen all: `python src/headroom_run_all.py`. Re-derive the table:
  `python src/headroom_report.py`.

## Honesty notes and limitations

- The oracle, naive, and the matched-accuracy scorer are fully deterministic. pciSeq is
  deterministic given its input. ComSeg samples internally; it is seeded (np.random,
  random, ComSegDict seed 0) but its community detection may not be bit-for-bit
  reproducible across runs, so results/headroom_methods.csv is committed from a fixed-seed
  run and is the recorded result. The byte-identical guarantee is verified for Gates 0 to 3.
- pciSeq and ComSeg were given the true nuclei as the prior to isolate assignment headroom
  from segmentation error, as specified. The nuclei prior (3 um disks) and ComSeg
  mean_cell_diameter (10 um) are the only configuration choices; method algorithm
  parameters are documented defaults.
- Free-segmentation context (no nuclei prior) is not represented by a successful run: the
  free-segmentation tools named (Proseg, Baysor) are toolchain-blocked here, and ComSeg's
  without-prior mode failed. The cost of inferring boundaries from scratch is therefore not
  quantified in this analysis.
- The fields are 250 cells for method tractability; the oracle ceiling matches the
  larger-field gate values (dense sigma 1.99 oracle 0.730 here vs 0.74 in Gate 3),
  confirming field-size independence.
