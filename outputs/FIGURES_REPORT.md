# Main-text figures: value provenance and refinement flags

**Branch:** `figures` (off `lever-analysis`). Append-only: no prior result or source file modified or deleted. Regenerate: `python src/figures/leg2_local.py` then `python src/figures/fig1_frontier.py fig2_bound.py fig3_methods.py fig4_worked.py` (each run separately). Outputs: `figures/fig{1,2,3,4}_*.pdf` (vector) and `.png` (raster, 400 dpi). Palette: Okabe-Ito colorblind-safe (`src/figures/figstyle.py`).

**Integrity linchpin held.** No oracle or segmenter is run on the real transcripts anywhere in the figures. Leg 1 and the contrast use the method-independent ceiling on synthetic fields at measured packing/density; Leg 2 uses the known-truth synthetic null; Leg 3 uses the published orthogonal data and the released matrix (read, not re-assigned).

## Two refinements this phase (the only new computation: `src/figures/leg2_local.py`)

### Flag 1 (decision-relevant): Leg 2 retreats from "sufficient" to "contributory"
The committed worked-analysis compared the synthetic adjacency-conditioned spurious rate against the section-wide observed rate, which is not apples-to-apples. Recomputed local-to-local (`results/lever/leg2_local.csv`), at threshold >=2:

| quantity | value | source |
|---|---|---|
| observed real, DCIS cells adjacent to a PR source | 24.19% (30/124) | leg2_local.csv `obs_rate_adj_to_prsource` |
| synthetic, pure displacement, sigma=0 control | 0.00% | lever_mechanism.csv (adj, sigma0) |
| synthetic, pure displacement, pinned 2.10 um | 1.26% | lever_mechanism.csv `spurious_triple_rate_dcis_adj_pr_t2` (pinned) |
| synthetic, pure displacement, upper CI 2.50 um | 1.96% | lever_mechanism.csv (ci_hi) |
| observed / synthetic(pinned) ratio | 19.2 | leg2_local.csv `obs_over_syn_pinned_ratio` |
| within one order of magnitude? | **No** | leg2_local.csv `within_one_order_of_magnitude` |

The pre-committed flag fires: the pure-displacement artifact under-produces the matched observed-local rate by about 19x (more than one order of magnitude). Per the spec, Leg 2 must retreat from "sufficient" to "contributory": displacement misassignment at the data-pinned sigma is a real, displacement-driven contributor to the spurious triple-positive signal (zero at sigma=0), but it does not by itself reproduce the observed-local magnitude. The manuscript prose for Leg 2 needs this change. Legs 1 and 3 are unaffected.

### Flag 2 (minor, presentation): ROI ceiling reconciled to the locked breast pin
Fig 4a evaluates the ROI ceiling at the manuscript's locked breast pin 1.99 um, CI [1.43, 2.63] (`results/gate3_pin_breast.csv`), not the lever diagnostic's re-pin (2.10 um, [1.43, 2.50]). The point values shift up slightly because 1.99 < 2.10 (`results/lever/lever_roi_pin199.csv`):

| ROI definition | oracle at 2.10 pin (manuscript) | oracle at locked 1.99 pin (Fig 4a) | optimistic CI end |
|---|---|---|---|
| triple-positive neighbourhood (7,103) | 0.798 | 0.809 | 0.862 |
| triple-positive cells' own (7,985) | 0.787 | 0.799 | 0.851 |
| densest DCIS nest (9,329) | 0.768 | 0.783 | 0.840 |

The qualitative claim is unchanged: all below 0.90 and 0.95 across the whole sigma CI. The manuscript's previously-quoted Fig 4a points ("0.79", "nest 0.768") should be updated to the reconciled locked-pin values ("0.80", "nest 0.783"), or annotated as the lever-pin values; optimistic ends (~0.85-0.86) are unchanged.

## Figure-by-figure value provenance (all from committed result files)

### Fig 1 (frontier and ceiling) - `results/sweep.csv`
- (c) oracle accuracy heatmap over packing x sigma at median density; full-grid range **0.234 to 0.957** (`oracle_acc` min/max).
- (d) oracle vs naive at all 125 points; **oracle >= naive everywhere, min margin +0.006** (`oracle_minus_naive` min = 0.00605).
- (a, b) schematic fields regenerated for illustration (seed 20260618+940000); no result value depends on them.

### Fig 2 (bound on real tissue)
- (a) validity null, breast: **28/38 admissible** (`gate3_validity_breast.csv` `admissible`, `pval`, `real_signal`); failing loose markers AHSP/CRHBP/CYP1A1 labelled.
- (b) pinned sigma: breast **1.99 [1.43, 2.63], 28/38** (`gate3_pin_breast.csv`); lung **1.83 [1.62, 2.08], 31/32** (`gate3_pin_lung.csv`).
- (c) dense oracle at ~13,600 cells/mm^2: breast **0.663-0.805**, lung **0.746-0.796** (`dense_oracle_at_sigma_ci_hi`..`_lo` in the two pin files); both below 0.9.
- (d) cell-typing stability: dense optimistic ends 0.694/0.762/0.785/0.803, sigma 1.99-2.65 (`gate3_typing.csv`); NB dispersion **4.19**, dense shift **-0.016** (`gate2_nbinom.csv`); realistic-expression frontier mean |diff| 0.006 (Gate 1, locked).

### Fig 3 (methods vs ceiling) - `results/headroom_fair_metrics.csv`, `headroom_fair_free.csv`
- (a) dense one-to-one accuracy vs sigma: oracle, naive, Baysor, Proseg, ComSeg.
- (b) best-method-to-oracle gap across three metrics and three dense sigmas: **0.066 to 0.110**.
- (c) best-method-minus-naive: one-to-one **-0.009 to -0.004**, many-to-one **+0.002 to +0.009**, co-assignment (ARI) **+0.009 to +0.019**; noise floor +/-0.012.
- (d) free segmentation (Baysor) dense sigma 1.99: free **0.521** < naive **0.639** < Baysor-prior 0.634-ish < oracle 0.730; Proseg/ComSeg free unsupported/failed (`headroom_fair_free.csv` status).

### Fig 4 (worked diagnostic) - `results/lever/*`
- (a) Leg 1: ROI oracle 0.783-0.809 (point), optimistic 0.840-0.862, at locked pin 1.99 (`lever_roi_pin199.csv`); section median packing 6,534; PGR+ in 2.9% of cells.
- (b) Leg 2 local-to-local: see Flag 1.
- (c) Leg 3 + contrast: Janesick ROI 0.799 (DISAGREES); colorectal 0.769 median (packing 17,424; AGREES) (`lever_ceiling.csv` crc); MERFISH hypothalamus 0.836 (AGREES). Leg 3 orthogonal: scFFPE-seq ~30 PGR+ cells none co-expressing ESR1/ERBB2; ERBB2 detected in 87% of cells; clinical PR-negative (Janesick et al., Nat Commun 14:8353, PMC10730913; prevalences from the released matrix).

## Provenance
- New code: `src/figures/figstyle.py`, `leg2_local.py`, `fig1_frontier.py`, `fig2_bound.py`, `fig3_methods.py`, `fig4_worked.py`.
- New results: `results/lever/leg2_local.csv`, `lever_roi_pin199.csv`. New figures: `figures/fig{1,2,3,4}_*.{png,pdf}`.
- Seeds: leg2_local/reconciliation LEVER_SEED 20260618+930000 = 21190618; Fig 1 schematic seed 20260618+940000 = 21200618.
- No prior result or source file modified (`git status --porcelain`: additions only, plus the Task-0 STATUS.md append).
- Discrepancies vs manuscript-locked values: only the two flags above (Leg 2 sufficient->contributory; Fig 4a points reconciled to the 1.99 pin). All other printed values match the committed results and the manuscript exactly.
