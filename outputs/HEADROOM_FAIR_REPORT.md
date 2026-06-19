# Method-fairness hardening: best documented configuration, over-segmentation-robust metrics, free-segmentation

Branch `headroom-fair` (off `headroom-linux`). Decision D11. Two attacks on the Linux headroom
result are closed here: (1) configuration-sensitivity (a 40 percent Baysor scale handicap was
already found and fixed, proving sensitivity), and (2) the one-to-one matched metric penalises
over-segmentation. This gate gives each method its best honest configuration over a sweep of
that method's own documented parameters, scores every run under three metrics that treat
over-segmentation differently, and fixes or scopes free-segmentation. Raw numbers only; no
conclusion is asserted. The orchestrator judges against the locked fairness criteria.

No prior-gate or prior-headroom result file was modified. The three method runners were
extended with backward-compatible optional parameters (defaults unchanged, so prior results
regenerate identically); all new results are append-only and namespaced `headroom_fair_*`.
Oracle and naive were recomputed within this environment (the committed `headroom_linux`
oracle/naive CSV, used here for reference, matches the Windows values to 3 dp).

## Metrics (all over interior transcripts; background/unassigned counts as error)

- **one2one**: Hungarian one-to-one matched accuracy (`headroom_common.matched_accuracy`, the
  existing scorer). Penalises over-segmentation (one true cell matches at most one predicted
  cell).
- **many2one**: assignment-homogeneity accuracy. Each predicted cell is mapped to its
  majority true cell; a transcript is correct if its predicted cell's majority equals its true
  cell. Does not penalise splitting one true cell into several predicted cells.
- **ari**: adjusted Rand index of the transcript co-assignment partition vs truth (background
  is its own cluster). Symmetric: penalises both over- and under-segmentation; not gameable by
  splitting.
- **ratio**: predicted/true cell-count ratio over interior transcripts (over-segmentation if
  > 1). Note even the oracle shows ratio about 1.3 because interior transcripts spill into
  neighbouring boundary cells, so 1.3, not 1.0, is the reference.
- **frac**: fraction of interior transcripts assigned a real cell (not background).

Code: `src/headroom_fair_common.py` (`all_metrics`). Sources for every number are named below.

## Task 1: fair configuration sweep (nuclei-prior, at the dense operating point sigma 1.99 um)

Swept the parameters each method's documentation/paper flags as most impactful. ComSeg is slow
on this aarch64 host (NUMBA_DISABLE_JIT=1, about 110 to 230 s per run), so all three methods
were swept at the single data-pinned dense point (dense_sigma1.99); the best config is then
evaluated at every operating point (Task 2). Source: `results/headroom_fair_sweep.csv`
(columns `config_label`, `acc_one_to_one`, `acc_many_to_one`, `ari`, `pred_true_ratio`,
`frac_assigned`). Reference at this point: oracle 0.730/0.746/0.568, naive 0.639/0.659/0.452
(one2one/many2one/ari).

Baysor (swept: scale as a multiple of the true mean cell radius, min-molecules-per-cell,
prior-segmentation-confidence):

| config | one2one | many2one | ari | ratio |
|---|---|---|---|---|
| scale1.0x m50 c0.5  | 0.624 | 0.671 | 0.471 | 1.54 |
| scale1.0x m50 c0.8  | 0.634 | 0.674 | 0.472 | 1.54 |
| scale1.25x m50 c0.5 | 0.621 | 0.650 | 0.461 | 1.38 |
| **scale1.25x m50 c0.8 (BEST)** | **0.635** | 0.666 | 0.471 | 1.42 |
| scale1.5x m50 c0.5  | 0.585 | 0.612 | 0.435 | 1.31 |
| scale1.5x m50 c0.8  | 0.624 | 0.653 | 0.458 | 1.36 |
| scale2.0x m50 c0.5  | 0.497 | 0.517 | 0.360 | 1.13 |
| scale2.0x m50 c0.8  | 0.593 | 0.623 | 0.430 | 1.35 |
| scale1.25x m20 c0.5 | 0.621 | 0.656 | 0.464 | 1.43 |
| scale1.5x m20 c0.5  | 0.604 | 0.631 | 0.447 | 1.35 |

Proseg (swept: voxel-size, prior-seg-reassignment-prob, nuclear-reassignment-prob, diffusion):

| config | one2one | many2one | ari | ratio | frac |
|---|---|---|---|---|---|
| defaults             | 0.556 | 0.570 | 0.085 | 1.34 | 0.807 |
| voxel0.5             | 0.530 | 0.545 | 0.051 | 1.35 | 0.751 |
| priorreassign0.3     | 0.560 | 0.576 | 0.092 | 1.36 | 0.814 |
| **priorreassign0.7 (BEST)** | **0.563** | 0.577 | 0.093 | 1.35 | 0.815 |
| nodiffusion          | 0.481 | 0.493 | 0.027 | 1.30 | 0.683 |
| nuclearreassign0.4   | 0.552 | 0.569 | 0.088 | 1.36 | 0.810 |

ComSeg (swept: mean_cell_diameter, the co-expression graph radius):

| config | one2one | many2one | ari | ratio |
|---|---|---|---|---|
| mcd8  | 0.624 | 0.647 | 0.439 | 1.34 |
| mcd10 | 0.634 | 0.654 | 0.447 | 1.35 |
| **mcd12 (BEST)** | **0.635** | 0.655 | 0.448 | 1.34 |

Best config per method (by one2one): Baysor scale1.25x/min50/conf0.8; Proseg
prior-seg-reassignment-prob 0.7; ComSeg mean_cell_diameter 12.

## Task 2: best config per method at every operating point, under all three metrics

Source: `results/headroom_fair_metrics.csv` (best config from the sweep, nuclei-prior;
oracle/naive recomputed within Linux). Each cell is one2one / many2one / ari, with (ratio,
frac) after.

| config | oracle | naive | Baysor (best) | Proseg (best) | ComSeg (best) |
|---|---|---|---|---|---|
| dense s1.43 | 0.803/0.814/0.669 | 0.731/0.746/0.564 | 0.719/0.748/0.573 (1.37, 1.00) | 0.617/0.628/0.116 (1.31, 0.83) | 0.722/0.737/0.557 (1.28, 1.00) |
| dense s1.99 | 0.730/0.746/0.568 | 0.639/0.659/0.452 | 0.634/0.666/0.469 (1.43, 1.00) | 0.564/0.578/0.093 (1.35, 0.81) | 0.635/0.655/0.448 (1.34, 1.00) |
| dense s2.63 | 0.651/0.672/0.467 | 0.546/0.572/0.354 | 0.540/0.581/0.373 (1.51, 1.00) | 0.480/0.500/0.068 (1.43, 0.80) | 0.540/0.568/0.352 (1.39, 1.00) |
| sparse s1.99 | 0.882/0.888/0.791 | 0.833/0.843/0.711 | 0.724/0.749/0.633 (1.21, 1.00) | 0.635/0.638/0.062 (1.20, 0.74) | 0.768/0.781/0.630 (1.22, 0.98) |
| representative s1.99 | 0.808/0.819/0.677 | 0.741/0.756/0.578 | 0.707/0.735/0.570 (1.38, 1.00) | 0.659/0.671/0.155 (1.29, 0.85) | 0.724/0.741/0.565 (1.26, 1.00) |

Notes:
- Proseg leaves 15 to 26 percent of interior transcripts as background (frac 0.74 to 0.85);
  this collapses its ARI (background forms one large cluster), which is why Proseg's ARI is far
  below the others under a metric that scores background as error. Baysor and ComSeg assign
  about 100 percent.
- Over-segmentation is modest after fair configuration: best-config ratios are 1.2 to 1.5,
  near the oracle's own 1.2 to 1.4 boundary-spill reference.

## Best-fair-method margins in the dense regime, under each metric

Source: `results/headroom_fair_metrics.csv`. "best" = the best of Baysor/Proseg/ComSeg at that
config and metric.

| metric | config | oracle | naive | best (method) | best minus naive | oracle minus best |
|---|---|---|---|---|---|---|
| one2one | dense s1.43 | 0.803 | 0.731 | 0.722 (ComSeg) | -0.009 | 0.080 |
| one2one | dense s1.99 | 0.730 | 0.639 | 0.635 (ComSeg) | -0.004 | 0.095 |
| one2one | dense s2.63 | 0.651 | 0.546 | 0.540 (ComSeg) | -0.005 | 0.110 |
| many2one | dense s1.43 | 0.814 | 0.746 | 0.748 (Baysor) | +0.002 | 0.066 |
| many2one | dense s1.99 | 0.746 | 0.659 | 0.666 (Baysor) | +0.007 | 0.080 |
| many2one | dense s2.63 | 0.672 | 0.572 | 0.581 (Baysor) | +0.009 | 0.091 |
| ari | dense s1.43 | 0.669 | 0.564 | 0.573 (Baysor) | +0.009 | 0.095 |
| ari | dense s1.99 | 0.568 | 0.452 | 0.469 (Baysor) | +0.018 | 0.098 |
| ari | dense s2.63 | 0.467 | 0.354 | 0.373 (Baysor) | +0.019 | 0.094 |

So, in the dense regime with each method at its best documented configuration:
- Under one2one, the best method (ComSeg) is below naive by 0.004 to 0.009.
- Under many2one, the best method (Baysor) is above naive by 0.002 to 0.009.
- Under ari, the best method (Baysor) is above naive by 0.009 to 0.019.
- The best method stays 0.066 to 0.110 below the oracle under every metric.

## Task 3: free-segmentation (no per-cell nucleus prior)

Source: `results/headroom_fair_free.csv`. Each Baysor cell is one2one / many2one / ari (ratio,
frac).

| config | Baysor free | Proseg free | ComSeg free |
|---|---|---|---|
| dense s1.43 | 0.625/0.660/0.508 (1.22, 1.00) | unsupported | failed |
| dense s1.99 | 0.521/0.543/0.407 (1.12, 1.00) | unsupported | failed |
| dense s2.63 | 0.453/0.483/0.344 (1.14, 1.00) | unsupported | failed |
| sparse s1.99 | 0.712/0.736/0.619 (1.16, 1.00) | unsupported | (not run) |
| representative s1.99 | 0.603/0.616/0.492 (1.05, 0.99) | unsupported | (not run) |

- **Baysor**: native de-novo joint segmentation runs (best scale/min-molecules from the
  sweep). Free-segmentation is well below both nuclei-prior Baysor and naive, as expected when
  boundaries are inferred from scratch.
- **Proseg**: no prior-free mode. It requires an initial segmentation (`--cell-id-column` from
  a nucleus or cellpose mask) and aborts otherwise ("Missing required argument:
  --cell-id-column"). Recorded `unsupported` at every config.
- **ComSeg**: no working prior-free mode in v1.8.5. The previous "without_prior" failure was an
  invalid argument (the valid graph-partition options are `with_prior` and `louvain`). The
  `louvain` mode requires `prior_name=None`, but `ComSegDataset.__init__` then raises
  `KeyError(None)` at `df_centroid[self.prior_name]` (dataset.py:137): it structurally requires
  a nucleus prior column on the centroids. Independently, ComSeg's `run_all` associates RNA
  communities to the centroid landmarks (`classify_centroid` + `associate_rna2landmark`), so a
  nucleus/centroid prior is required regardless. ComSeg cannot run fully prior-free here.

Free-segmentation is therefore scoped to Baysor; Proseg and ComSeg are both re-segmentation
methods that require a nucleus (or cell) prior in their documented usage.

## Reproducibility, seeds, sources

- Oracle, naive, and all metric code are deterministic; the within-Linux oracle/naive CSV is
  bit-identical across runs.
- Baysor and Proseg expose no CLI seed; ComSeg is seeded (np.random/random/ComSegDict seed 0)
  but with NUMBA_DISABLE_JIT=1 here. Measured run-to-run spread is about 0.012 in matched
  accuracy (`results/headroom_linux_repro.csv`); the fair eval reproduces the sweep at
  dense_sigma1.99 within about 0.002 (Baysor 0.634 vs 0.635, Proseg 0.564 vs 0.563, ComSeg
  deterministic), consistent with that spread.
- Sources: `results/headroom_fair_sweep.csv` (Task 1), `results/headroom_fair_metrics.csv`
  (Task 2), `results/headroom_fair_free.csv` (Task 3), `results/headroom_linux_oracle_naive.csv`
  (oracle/naive reference). Code: `src/headroom_fair_common.py`, `headroom_fair_sweep.py`,
  `headroom_fair_eval.py`, `headroom_fair_free.py`, `headroom_fair_report.py`; runners
  `src/methods_baysor.py`, `methods_proseg.py`, `methods_comseg.py` (extended with
  backward-compatible optional parameters). Seeds: data export base `MASTER_SEED + 600000`,
  per-config seeds in `results/headroom_linux_oracle_naive.csv`.
- Confirmation: no prior-gate or prior-headroom result file was modified (`git diff` touches
  only `headroom_fair_*` results and the method source files; all `headroom_linux_*` and
  earlier-gate result CSVs are unchanged).

## Honesty notes and limitations

- The fair sweep was performed at the single data-pinned dense operating point
  (dense_sigma1.99) for cost; the best config is then evaluated at all five operating points.
  A per-operating-point sweep could shift a method's best config slightly, but the scale that
  wins for Baysor is expressed relative to the true cell radius, so it transfers.
- Proseg's free-segmentation column is genuinely absent (no prior-free mode), and ComSeg's is
  absent for the same structural reason; only Baysor has a working de-novo mode.
- This report states raw numbers, margins, and gaps under each metric and configuration. It
  does not judge whether any method beats naive or reaches the wall; that is the orchestrator's
  call against the locked fairness criteria.
