# Segmenter-headroom redo on Linux: the field-standard methods (Baysor, Proseg) vs the oracle ceiling

Branch `headroom-linux` (off `headroom`). This is the Linux/WSL redo motivated by Decision
D10: the Windows run could not install Baysor or Proseg (toolchain wall) and the substitutes
(pciSeq, ComSeg) were likely handicapped, so where the methods the field actually uses sit
relative to the oracle stayed open. This run installs and runs Baysor and Proseg, recomputes
the oracle and naive baselines inside this Linux environment, and places the methods against
that ceiling in both the fair nuclei-prior mode and the native free-segmentation mode. Raw
numbers only; no conclusion about "near the wall" versus "headroom" is asserted. The
orchestrator interprets the gaps.

No prior-gate or prior-headroom source file was modified. All new code is namespaced
`headroom_linux_*` / `methods_baysor.py` / `methods_proseg.py`; all new results are
append-only (`results/headroom_linux_*.csv`). The committed Windows results remain canonical;
cross-platform byte-identity was explicitly not attempted (D10).

## Environment (the reason the install story differs from the prompt's assumption)

The host is **aarch64 (ARM64) WSL2**, not x86_64, and started bare: system Python 3.14 only,
no pip, conda, cargo, Julia, or C/C++ compiler, and no passwordless sudo. A user-space
toolchain was bootstrapped: micromamba 2.8.1, then a conda env (`st`) with Python 3.11.15,
numpy 2.4.6, scipy 1.17.1, pandas 2.3.3, scikit-learn 1.9.0, h5py 3.16.0, anndata 0.12.17.
The ARM64 architecture, not the operating system, is what shaped every install outcome below.

## Methods: versions, exact install routes, what failed

- **Proseg 3.1.1** (Newell lab, Nature Methods 2025): installed natively from **bioconda**,
  `micromamba install -c bioconda rust-proseg` (package `rust-proseg`, binary `proseg`),
  linux-aarch64 build. Ran cleanly (about 13 to 29 s per config). `proseg --version` ->
  `proseg 3.1.1`.
- **Baysor 0.7.1** (Petukhov et al., Nature Biotechnology 2022): **no native ARM64 prebuilt
  exists.** Every Baysor GitHub release binary is x86_64-only (`baysor-x86_x64-linux-*`), the
  newest releases (cpp-0.8.x) ship no binaries at all, and the bioconda `baysor` package is
  linux-64 only. It was therefore built natively from the Julia package at git tag **v0.7.1**
  (`Pkg.add(url="https://github.com/kharchenkolab/Baysor.git", rev="v0.7.1")`) into an
  isolated Julia **1.10.11** (aarch64) depot, and invoked through the wrapper `~/baysor_run.sh`
  (`Baysor.command_main`, Comonicon CLI). Note: Baysor's `deps/build.jl` (which builds a
  standalone PackageCompiler binary) fails here for lack of a C compiler, but that step is
  optional; the package itself precompiles and runs. Ran about 50 to 60 s per config.
- **pciSeq**: NOT installable here. pip dependency resolution is impossible
  (`ResolutionImpossible`) against the modern numpy/scipy/scikit-learn stack on this Python.
- **ComSeg 1.8.5** (Defard et al. 2024): installs only after adding a conda toolchain
  (compilers + `scikit-misc`), and then **crashes at runtime on aarch64** with an LLVM code
  generator abort (`UNREACHABLE ... llvm/lib/CodeGen/TargetSchedule.cpp`, SIGABRT / exit 134)
  inside its numba-JIT path. It runs only with `NUMBA_DISABLE_JIT=1` (pure-Python, slower). A
  single validated dense run (`dense_sigma1.99`, nuclei-prior) under that setting gives matched
  accuracy 0.634, frac_assigned 0.998, which matches the Windows ComSeg dense number (0.634).
  A full within-Linux ComSeg sweep under `NUMBA_DISABLE_JIT=1` is being run for one consistent
  table; ComSeg rows will be appended to `results/headroom_linux_methods.csv` when it
  completes. The Windows pciSeq/ComSeg numbers remain the reference for those two methods.

So the two priority field-standard methods (Baysor, Proseg) both ran; the two Windows
substitutes (pciSeq, ComSeg) are blocked or degraded on ARM64 for toolchain reasons, exactly
the kind of handicap this redo was meant to remove for Baysor and Proseg.

## Metric

Matched transcript-assignment accuracy against the known truth: the fraction of interior
transcripts whose method-cell is matched to their true cell under the optimal one-to-one
matching of method-cells to true-cells (Hungarian on the contingency table); background or
unassigned transcripts count as errors. Identical scorer to the Windows run:
`headroom_common.matched_accuracy` (reused unchanged). For the oracle and naive baselines
(which use the true cell ids) matched accuracy equals direct id-match accuracy.

`frac_assigned` is reported alongside every method number: the fraction of interior
transcripts the method gave a real cell (not background), so any background-dumping is visible
and not silently scored as assignment error.

## Configuration (identical to the Windows headroom run)

Most realistic Gate 3 config: Xenium-breast realistic overlapping expression (K=14 types, 313
genes, from the 10x FFPE Human Breast Cancer Rep1 release 1.0.1), negative-binomial
within-type emission, data-pinned displacement at the Gate 3 breast sigma 1.99 um and its CI
ends 1.43 and 2.63 um, dense packing 13,625 cells/mm^2 (plus sparse 2,575 and representative
6,534), density 164 tx/cell, 250-cell fields (the oracle accuracy is field-size-independent).
Code: `src/headroom_common.py`, `src/headroom_linux_export.py`. Seeds: data export base
`MASTER_SEED + 600000`, per-config seeds recorded in `results/headroom_linux_oracle_naive.csv`
(column `seed`).

## Within-Linux oracle and naive baselines (every config)

Source: `results/headroom_linux_oracle_naive.csv`, columns `oracle_acc_matched`,
`naive_acc_matched`. Recomputed in THIS environment; bit-reproducible within Linux (identical
CSV across two runs). They match the canonical Windows values to three decimals, confirming
the generator and oracle are stable across platforms here.

| config (packing, sigma)        | oracle | naive | interior transcripts |
|--------------------------------|--------|-------|----------------------|
| dense, sigma 1.43 um           | 0.803  | 0.731 | 29,931 |
| dense, sigma 1.99 um           | 0.730  | 0.639 | 29,813 |
| dense, sigma 2.63 um           | 0.651  | 0.546 | 27,125 |
| sparse, sigma 1.99 um          | 0.882  | 0.833 | 33,384 |
| representative, sigma 1.99 um  | 0.808  | 0.741 | 31,348 |

## Main table: methods vs the ceiling, both modes

All values are matched transcript-assignment accuracy vs the known truth, with frac_assigned
in parentheses. Source: oracle/naive from `results/headroom_linux_oracle_naive.csv`; methods
from `results/headroom_linux_methods.csv` (columns `matched_accuracy`, `frac_assigned`,
rows by `method` and `mode`).

| config (packing, sigma)       | oracle | naive | Baysor nuclei-prior | Baysor free | Proseg nuclei-prior | Proseg free |
|-------------------------------|--------|-------|---------------------|-------------|---------------------|-------------|
| dense, sigma 1.43 um          | 0.803  | 0.731 | 0.701 (1.00)        | 0.633 (1.00)| 0.625 (0.84)        | not supported |
| dense, sigma 1.99 um          | 0.730  | 0.639 | 0.625 (1.00)        | 0.570 (1.00)| 0.564 (0.80)        | not supported |
| dense, sigma 2.63 um          | 0.651  | 0.546 | 0.541 (1.00)        | 0.489 (1.00)| 0.482 (0.80)        | not supported |
| sparse, sigma 1.99 um         | 0.882  | 0.833 | 0.758 (1.00)        | 0.749 (1.00)| 0.620 (0.75)        | not supported |
| representative, sigma 1.99 um | 0.808  | 0.741 | 0.710 (1.00)        | 0.674 (0.99)| 0.661 (0.85)        | not supported |

Proseg free-segmentation is recorded as "not supported": Proseg requires an initial
segmentation (`--cell-id-column`) and aborts without one ("Missing required argument:
--cell-id-column"), so it has no prior-free mode. Its nuclei-prior run is already its native,
field-standard joint segmentation-and-assignment from a nucleus seed (it grows and reshapes
cells; default `prior-seg-reassignment-prob` 0.5). Status rows in
`results/headroom_linux_methods.csv`, column `status` = `unsupported`.

## Gaps to the oracle and comparison to naive

Best method at every config is **Baysor (nuclei-prior)**. Source: derived from the two CSVs
above (`src/headroom_linux_report.py`).

Best-method-to-oracle gap (oracle minus Baysor nuclei-prior):
- dense sigma CI (1.43 / 1.99 / 2.63 um): **0.102 / 0.105 / 0.109**
- sparse: 0.123; representative: 0.098

Free-segmentation-to-oracle gap (oracle minus Baysor free, the only method with a prior-free
mode):
- dense sigma CI: **0.169 / 0.160 / 0.161**
- sparse: 0.133; representative: 0.134

Does each method beat the naive nearest-nucleus baseline?
- **No method exceeds naive at any config in either mode** (beats_naive = N everywhere).
- But Baysor (nuclei-prior) sits very close to naive in the dense regime: naive minus Baysor
  = 0.030 / 0.014 / 0.004 across the dense sigma CI (0.075 sparse, 0.031 representative).
- Baysor (free) below naive by 0.098 / 0.069 / 0.056 dense (0.084 sparse, 0.067
  representative).
- Proseg (nuclei-prior) below naive by 0.106 / 0.075 / 0.063 dense (0.213 sparse, 0.080
  representative).

## Background-dumping and over-segmentation (so the gaps are read correctly)

- **frac_assigned**: Baysor assigns essentially every interior transcript to a cell (1.00,
  except 0.99 at representative). Proseg leaves 15 to 25 percent of interior transcripts as
  background (frac 0.75 to 0.85), which counts as error under the metric and is part of why
  Proseg trails Baysor.
- **Over-segmentation**: under the one-to-one matched metric, splitting one true cell into two
  method-cells leaves one of them unmatched, so its transcripts score as errors. At the dense
  operating point (dense_sigma1.99, 182 true interior cells) the methods produce more cells
  than that: Baysor about 275 to 280 method-cells, Proseg about 242 to 247 (results below).
  Both over-segment by roughly 1.3 to 1.5 times the interior cell count (the field has 256
  cells total), which caps their matched accuracy. Source:
  `results/headroom_linux_repro.csv`, columns `n_method_cells`, `n_true_interior_cells`.

## Reproducibility

- Oracle, naive, and the matched-accuracy scorer are fully deterministic; the oracle/naive CSV
  is bit-identical across two within-Linux runs.
- **Baysor and Proseg expose no CLI seed**, so their outputs are not bit-for-bit reproducible.
  Run-to-run variation was quantified by repeating each three times at the dense operating
  point (dense_sigma1.99), `results/headroom_linux_repro.csv`:
  - Baysor nuclei-prior: 0.622 to 0.630 (spread 0.008), n_method_cells 275 to 280
  - Baysor free: 0.564 to 0.569 (spread 0.006), n_method_cells 260 to 266
  - Proseg nuclei-prior: 0.554 to 0.566 (spread 0.012), n_method_cells 242 to 247
  Run-to-run spread is at or below about 0.012 in matched accuracy, well under the gaps above.

## Parameter choices (documented defaults; not tuned to approach or miss the oracle)

- Nucleus prior, all methods: a transcript is tagged with its nearest true cell centre if
  within 3.0 um (the same 3 um nucleus-disk prior used for the Windows pciSeq/ComSeg runs).
- Baysor: Baysor's own shipped Xenium config (`configs/xenium.toml`):
  `min_molecules_per_cell = 50`, `prior_segmentation_confidence = 0.5` (nuclei-prior mode).
  The expected-cell-radius `scale` was set to the true mean cell radius, estimated from the
  centroid bounding box (area_per_cell = bbox_area / n_cells, r = sqrt(area_per_cell / pi)),
  which recovers the generator's true r_mean to within about 3 percent. This is a geometry
  input Baysor requires (the correct cell radius), not a tuned value. (An earlier run used a
  nearest-neighbour-spacing estimator that under-estimated the radius by about 40 percent
  because of the generator's cell jitter; that handicapped Baysor and was corrected. The
  corrected scale per config: 4.71 / 4.78 / 4.75 / 10.73 / 6.84 um.)
- Proseg: nucleus cell-id prior via `--cell-id-column` with `--cell-id-unassigned -1`, 2D
  (`--ignore-z-coord`), all algorithm parameters at documented defaults.

## Committed paths and regeneration

- Oracle/naive recompute + known-truth export: `src/headroom_linux_export.py` (reuses
  `headroom_common.build_config`, `export_config`, `matched_accuracy`) ->
  `results/headroom_linux_oracle_naive.csv`; raw tables under `data/headroom_linux/` (gitignored).
- Method runners: `src/methods_baysor.py` (Baysor via `~/baysor_run.sh`),
  `src/methods_proseg.py` (Proseg); driver `src/headroom_linux_methods.py` ->
  `results/headroom_linux_methods.csv`.
- Reproducibility + over-segmentation probe: `src/headroom_linux_repro.py` ->
  `results/headroom_linux_repro.csv`.
- Table and gaps: `src/headroom_linux_report.py`.
- Regenerate: `micromamba run -n st python src/headroom_linux_export.py`, then
  `micromamba run -n st python src/headroom_linux_methods.py` (Baysor + Proseg), then
  `micromamba run -n st python src/headroom_linux_report.py`. (For the within-Linux ComSeg
  rows, prefix the methods call with `NUMBA_DISABLE_JIT=1` and pass `comseg`.)

## Honesty notes and limitations

- Baysor and Proseg numbers are single recorded runs (no seed); the run-to-run spread above
  bounds the noise (about 0.01).
- The matched one-to-one metric penalises over-segmentation; both methods over-segment
  relative to the interior cell count, so part of each gap is segmentation granularity rather
  than transcript-level assignment error. The metric is held identical to the Windows run by
  design.
- Proseg has no prior-free mode, so the free-segmentation column for Proseg is genuinely
  absent (not a failure to run). Baysor is the only method here with both modes.
- pciSeq and ComSeg could not be re-run cleanly within this ARM64 environment for the reasons
  documented above; the Windows numbers stand for those two, and a JIT-disabled within-Linux
  ComSeg sweep is being appended for context.
- This report states raw numbers and gaps only. Whether the gaps mean the field-standard
  methods are near the recoverability wall or have material headroom is left to the
  orchestrator.
