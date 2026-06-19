# LOCKED NUMBERS AND GUARDRAILS

This is the document trusted when chat history and memory cannot be. A number enters the locked table only after it has been verified to a committed source file.

## Guardrails (non-negotiable)
- **Verify then write.** No number, caption, table legend, or claim is authored before its source file is confirmed this session. If you catch yourself writing an unseen number, stop and ask for it.
- **Provenance is not consistency.** Every headline number traces to a committed result file and column, not to "it matches what we said before."
- **Gate discipline.** Pre-committed thresholds (below) are honored. They are not renegotiated after the data is seen.
- **Audit everything.** Deep-research output and executor (Claude Code) output are leads, not verdicts. Reproduce or cross-check before trusting.
- **No overclaiming.** The modest, exact, honestly-bounded claim is the strongest one. A clean negative or boundary result is a real contribution.
- **Writing rules (paper and every caption; internal docs exempt):** no em dashes (use commas, colons, or separate sentences; en dashes in numeric ranges are fine; grep after every edit). Readable by a non-specialist without losing precision; every number and confidence interval stated exactly.

## Locked numbers (project headline numbers only)
Gate 0 measured values below are verified to committed result files on branch gate0. The Gate 0 verdict itself (PASS/KILL) is pending orchestrator audit and is NOT recorded here as a finding.

| number (what it is) | value | source result file | column / field | date verified |
|---|---|---|---|---|
| oracle assignment accuracy, full sweep range | [0.234, 0.957] | results/sweep.csv | oracle_acc | 2026-06-18 |
| oracle assignment-partition ARI vs true cell labels, full range | [0.110, 0.919] | results/sweep.csv | assign_ari_vs_truecells_oracle | 2026-06-18 |
| oracle expression-profile ARI vs true types, full range | [1.000, 1.000] | results/sweep.csv | profile_ari_vs_truetype | 2026-06-18 |
| naive expression-profile ARI vs true types, full range | [0.102, 1.000] | results/sweep.csv | naive_profile_ari_vs_truetype | 2026-06-18 |
| oracle >= naive at every grid point | True (min margin +0.0061) | results/sweep.csv | oracle_ge_naive, oracle_minus_naive | 2026-06-18 |
| oracle assignment accuracy inside realistic band (sigma<=2um) | [0.719, 0.953] | results/realism_oracle.csv | oracle_acc where in_realistic_sigma_band | 2026-06-18 |
| generator vs real median tx/cell max relative error at anchors | 0.0152 | results/realism_oracle.csv | tx_match_rel_err | 2026-06-18 |
| generator vs real packing max relative error at anchors | 1.4e-16 | results/realism_oracle.csv | packing_match_rel_err | 2026-06-18 |
| real median transcripts/cell (Xenium breast / MERFISH hypothalamus) | 66-86 / 267-312 | results/realism.csv | median_tx_per_cell | 2026-06-18 |
| real packing (Xenium breast / MERFISH hypothalamus), cells/mm^2 | 13250-14125 / 2500-2700 | results/realism.csv | packing_cells_per_mm2 | 2026-06-18 |
| GATE 1: oracle assignment accuracy under realistic expression, full range | [0.251, 0.957] | results/gate1_sweep.csv | oracle_acc | 2026-06-18 |
| GATE 1: oracle frontier vs Gate 0, mean/max abs diff at matched points | 0.006 / 0.032 | results/gate1_sweep.csv vs results/sweep.csv | oracle_acc | 2026-06-18 |
| GATE 1: oracle profile ARI under realistic expression (de-saturated) | [0.667, 1.000] | results/gate1_sweep.csv | profile_ari_vs_truetype | 2026-06-18 |
| GATE 1: realistic-expression overlap (mean pairwise cosine; exclusive-mass frac) | 0.210 ; 0.045 | results/gate1_expression_meta.json | overlap_metrics | 2026-06-18 |
| GATE 1: data-pinned displacement sigma (point; 25-75 bracket), um | 2.67 ; [0.10, 3.03] | results/gate1_leakage.csv | sigma_pinned_um, sigma_bracket_lo/hi_um | 2026-06-18 |
| GATE 1: oracle acc at pinned sigma, MERFISH-like (sparse) / Xenium-like (dense) | 0.836 / 0.643 | results/gate1_leakage.csv | oracle_acc_at_pinned | 2026-06-18 |
| GATE 1: structural-sensitivity oracle ranges (baseline / aniso / mixture) | [0.290,0.925] / [0.281,0.924] / [0.254,0.833] | results/gate1_structural.csv | oracle_acc | 2026-06-18 |
| GATE 2: Xenium-direct pinned displacement sigma (point) | 1.61 um | results/gate2_xenium_pin.csv | xenium_pinned_sigma_um | 2026-06-18 |
| GATE 2: sigma statistical 95% CI (bootstrap B=2000) | [0.79, 2.53] um | results/gate2_sigma_uncertainty.csv | sigma_stat_ci_lo/hi_um | 2026-06-18 |
| GATE 2: sigma design-sensitivity range (27 design choices) | [0.20, 2.48] um | results/gate2_sigma_uncertainty.csv | sigma_design_lo/hi_um | 2026-06-18 |
| GATE 2: dense oracle accuracy at pinned sigma 1.61 um (converged grid) | 0.787 | results/gate2_sigma_uncertainty.csv | oracle_dense_at_point | 2026-06-18 |
| GATE 2: dense oracle accuracy at statistical-CI ends (low/high sigma) | 0.891 / 0.675 | results/gate2_sigma_uncertainty.csv | oracle_dense_at_stat_ci_lo/hi | 2026-06-18 |
| GATE 2: dense oracle accuracy at combined-band optimistic end (sigma 0.20 um) | 0.972 | results/gate2_sigma_uncertainty.csv | oracle_dense_at_combined_lo | 2026-06-18 |
| GATE 2: pinned sigma by marker cutoff (0.7 stable / 0.6 wide / 0.8 degenerate) | [1.49,1.70] / [0.20,2.48] / none | results/gate2_design_grid.csv | sigma_um by excl_threshold | 2026-06-18 |
| GATE 2: NB-emission sigma shift; dense accuracy shift | +0.10 um ; -0.016 | results/gate2_nbinom.csv | sigma_shift_um, dense_oracle_shift | 2026-06-18 |
| GATE 2: NB-emission frontier monotonicity / oracle>=naive | (0,0) viol ; True | results/gate2_nbinom.csv | nb_frontier_mono_viol_*, nb_frontier_oracle_ge_naive_all | 2026-06-18 |
| GATE 3: breast admissible markers (validity test, p<0.05, 1000 perms) | 28 of 38 | results/gate3_validity_breast.csv | admissible | 2026-06-18 |
| GATE 3: breast admissible vs rejected mean displacement signal | +0.094 / -0.055 | results/gate3_validity_breast.csv | real_signal by admissible | 2026-06-18 |
| GATE 3: breast destabilising loose markers permutation p-values | 0.24-0.996 (all fail) | results/gate3_validity_breast.csv | pval (AHSP,CRHBP,CYP1A1,EGFL7,SLC4A1) | 2026-06-18 |
| GATE 3: breast re-pinned sigma (admissible set) point; bootstrap CI | 1.99 ; [1.43, 2.63] um | results/gate3_pin_breast.csv | sigma_point_um, sigma_ci_lo/hi_um | 2026-06-18 |
| GATE 3: breast dense oracle accuracy across sigma CI | [0.663, 0.805] | results/gate3_pin_breast.csv | dense_oracle_at_sigma_ci_lo/hi | 2026-06-18 |
| GATE 3: cell-typing sigma spread (K=10/14/14alt/20) | [1.99, 2.65] um | results/gate3_typing.csv | sigma_point_um | 2026-06-18 |
| GATE 3: cell-typing max dense accuracy at any sigma-CI optimistic end | 0.803 | results/gate3_typing.csv | dense_oracle_at_sigma_ci_lo | 2026-06-18 |
| GATE 3: lung admissible markers; re-pinned sigma point; bootstrap CI | 31 of 32 ; 1.83 ; [1.62, 2.08] um | results/gate3_pin_lung.csv | n_admissible_markers, sigma_point_um, sigma_ci_lo/hi_um | 2026-06-18 |
| GATE 3: lung dense oracle accuracy across sigma CI | [0.746, 0.796] | results/gate3_pin_lung.csv | dense_oracle_at_sigma_ci_lo/hi | 2026-06-18 |

### Can-claim update (2026-06-18, post Gate 3)
- Under a principled, non-circular marker-validity test (admissible only if the adjacent-minus-distant displacement signal exceeds a position-permutation null at p<0.05, never selected by the resulting sigma), the markers that destabilised sigma at the loose cutoff fail the test, and the re-pinned sigma keeps dense-regime best-possible assignment accuracy below 0.95 (and below 0.9) across its bootstrap CI on the breast dataset, on an independent Xenium lung-cancer dataset, and across four cell-typing choices.
- Displacement sigma measured under the validity test is about 1.8-2.0 um on both tissues (typing spread 1.99-2.65 um), a factor of ~1.3 across reasonable choices, versus the ~12x swing across the old arbitrary marker cutoffs.

### Can-claim update (2026-06-18, post Gate 3, expanded)
- Under a principled non-circular validity test, on Xenium data calibrated to realistic expression overlap and a data-pinned displacement of about 1.8 to 2.0 um, the best-possible transcript-to-cell assignment accuracy in dense tissue (about 13,600 cells/mm2, realistic in both breast and lung) is bounded below 0.9 (about 0.66 to 0.81 across the sigma CI), replicated on two independent tissues, stable across cell-typing, and robust to within-type expression variance and to non-Voronoi geometry and non-Gaussian displacement.

### Cannot-claim update (2026-06-18, post Gate 3)
- That real published methods reach this ceiling. The bound is on the best possible; where Baysor, Proseg, and similar methods sit relative to it is measured in the headroom analysis, not yet known.
- Anything about Xenium in general; the bound is the dense-tissue regime.

### Headroom measured values (2026-06-18, branch headroom; matched transcript-assignment accuracy vs known truth)
| number (what it is) | value | source result file | column / field | date verified |
|---|---|---|---|---|
| oracle ceiling, dense sigma CI (1.43/1.99/2.63 um) | 0.803 / 0.730 / 0.651 | results/headroom_oracle_naive.csv | oracle_acc_matched | 2026-06-18 |
| naive nearest-nucleus, dense sigma CI | 0.731 / 0.639 / 0.546 | results/headroom_oracle_naive.csv | naive_acc_matched | 2026-06-18 |
| ComSeg (nuclei prior), dense sigma CI | 0.722 / 0.634 / 0.543 | results/headroom_methods.csv | matched_accuracy (ComSeg, nuclei_prior) | 2026-06-18 |
| pciSeq (nuclei prior), dense sigma CI | 0.667 / 0.604 / 0.504 | results/headroom_methods.csv | matched_accuracy (pciSeq, nuclei_prior) | 2026-06-18 |
| best-sophisticated-method (ComSeg) to oracle gap, dense sigma CI | 0.081 / 0.096 / 0.108 | results/headroom_oracle_naive.csv, headroom_methods.csv | derived | 2026-06-18 |
| methods that ran / toolchain-blocked | pciSeq 0.0.61, ComSeg 1.8.5 ran; Proseg, Baysor blocked | outputs/HEADROOM_REPORT.md | install outcomes | 2026-06-18 |

### Can-claim update (2026-06-18, post headroom)
- On known-truth synthetic data matched to the realistic Gate 3 configuration, two real published nuclei-prior assignment methods (pciSeq, ComSeg) both fall below the oracle ceiling at every regime, with the best of them (ComSeg) 0.08 to 0.12 below the dense-regime ceiling, and neither exceeds the naive nearest-nucleus baseline. (Proseg and Baysor could not be installed on the test machine; the free-segmentation cost is not quantified.)

### Cannot-claim update (2026-06-18, post headroom)
- That methods are "near the wall" or "have headroom": that interpretation of the measured oracle-minus-method gap is for the orchestrator, not asserted by the executor.
- A free-segmentation (no nuclei prior) headroom number; that mode did not produce a successful run on this machine.

### Can-claim update (2026-06-18, post Gate 2)
- The displacement sigma can be measured directly on Xenium from spatial marker leakage: point 1.61 um, statistical 95% CI [0.79, 2.53] um.
- Across the statistical 95% CI the dense-regime best-possible assignment accuracy is [0.675, 0.891], all below 0.95; the limit bites under directly-measured Xenium conditions across the statistical uncertainty, and survives negative-binomial within-type emission (monotone, oracle at or above naive).

### Cannot-claim update (2026-06-18, post Gate 2)
- That the dense-regime limit bites across the ENTIRE combined uncertainty band: at the loosest marker cutoff (0.6) with the tightest adjacency bin, sigma reaches 0.20 um where dense accuracy is 0.972 (above 0.95). The optimistic-end claim is sensitive to the marker-selection design choice.
- A sigma stable to within a tight factor across all design choices (it spans about 13x across the full grid, though it is stable at the pre-registered cutoff 0.7).

### Can-claim update (2026-06-18, post Gate 2)
- Displacement sigma measured directly on Xenium is about 1.61 um (bootstrap 95% CI [0.79, 2.53]); under the pre-registered exclusive-marker statistic the statistical CI keeps dense-regime best-possible assignment accuracy below 0.95, and the limit survives negative-binomial within-type variance.

### Cannot-claim update 2 (2026-06-18, post Gate 2)
- The dense-regime limit as a robust result. It depends on the marker-exclusivity criterion: loosening the cutoff to 0.6 pushes dense oracle accuracy to 0.972, above 0.95. The claim must not be made until the exclusivity criterion is principled (non-circular) and the limit is replicated on an independent Xenium dataset.
- Any "Xenium in general" framing; the limit is specific to the dense tissue tail.

### Gate 3 pre-committed thresholds (locked 2026-06-18; orchestrator judges, executor does not)
- PASS if, under a principled non-circular marker-validity test: the previously destabilizing markers are shown to fail the validity test on independent grounds; the re-pinned sigma is stable and its bootstrap CI keeps dense-regime oracle accuracy below 0.95 on the breast dataset AND on at least one independent Xenium dataset of a different tissue; and sigma and the dense limit are stable across reasonable cell-typing choices.
- KILL if, under the principled validity test, the dense-regime oracle accuracy exceeds 0.95 across the sigma CI on either dataset.
- KILL if no principled validity test cleanly separates the sigma-destabilizing markers from the stable ones, meaning the cutoff dependence is genuine rather than a marker-validity artifact.
- KILL if sigma or the dense limit swings across reasonable cell-typing choices the way it did across the old marker cutoffs.
- RETREAT (not full kill) if the dense limit holds on breast but not on the independent tissue: the claim narrows and the paper reframes around the method, the frontier, and the bounded result rather than a general dense-tissue limit.

### Can-claim update (2026-06-18, post Gate 1)
- The recovery frontier survives realistic overlapping expression (assignment frontier essentially unchanged, oracle-minus-naive gap preserved) and survives non-Voronoi geometry and non-Gaussian displacement (monotone, oracle at or above naive, dense-regime limit biting).
- At a sigma pinned to measured MERFISH spatial marker-leakage, best-possible transcript-assignment accuracy is about 0.84 sparse and 0.64 dense, both below 0.9.

### Cannot-claim update (2026-06-18, post Gate 1)
- The dense-regime limit for Xenium specifically (sigma was borrowed from MERFISH).
- A precise sigma (CI uncharacterized; per-marker spread reaches near zero, where the limit would not bite).
- Final accuracy numbers (per-type-mean emission makes them upper bounds; the real limit is at least this severe).

### Gate 2 pre-committed thresholds (locked 2026-06-18; orchestrator judges, executor does not)
- KILL if the optimistic (low-sigma) end of the combined Xenium-direct sigma uncertainty (statistical CI plus design-sensitivity) gives dense-regime oracle assignment accuracy above 0.95. The data cannot then rule out the absence of a limit.
- KILL if sigma cannot be reliably measured on Xenium: the leakage statistic is degenerate or noise-dominated, or sigma swings more than an order of magnitude across reasonable design choices with no defensible central value.
- KILL if, under non-idealized negative-binomial emission, the frontier flattens to oracle accuracy above 0.95 across the realistic regime.
- PASS if the entire combined Xenium-direct sigma uncertainty band keeps dense-regime oracle accuracy below 0.95 (the limit robustly bites), sigma is stable across design choices within its CI, and the frontier and the dense-regime bite survive non-idealized emission. Report the dense-regime accuracy as a band with honest error bars, and note where that band sits relative to 0.9.

## Gate 0 pre-committed thresholds (locked 2026-06-18)
Decision rests on the oracle frontier shape and the realism match. Running real segmenters (Baysor / Proseg) is an optional headroom signal, not a kill criterion, so install friction does not gate the decision.
- **KILL** if the oracle assignment accuracy (or oracle adjusted Rand index vs true cell labels) is > 0.95 across the entire realistic regime (the band where the generator matches real summary statistics). No fundamental limit in the relevant range means no frontier.
- **KILL** if oracle accuracy is uniformly below about 0.55 with no monotone dependence on density or packing across the grid. No structured frontier to characterize.
- **KILL** if the generator cannot match the real-data summary statistics (median transcripts per cell and cell-packing density) to within 15% anywhere in its parameter range. The bound would not transfer; fix the generator before proceeding.
- **PASS** if the oracle frontier is non-trivial (oracle accuracy varies smoothly from below 0.6 to above 0.9 across the realistic density/packing range) AND the real-data-matched regime sits in the informative band (oracle below about 0.9, i.e., a fundamental limit actually bites under real conditions).

## Can-claim (as of 2026-06-18)
- The literature facts verified by primary-source search (with citations): Proseg (Nat Methods 2025), segger (bioRxiv 2025), TopACT generator (Nature 2024), the field's open admission that it lacks segmentation ground truth.
- Nothing empirical about this project. No frontier, no oracle number, no recovery number exists yet.

## Cannot-claim (as of 2026-06-18)
- That an answerability frontier exists. (To be tested in Gate 0.)
- Any recovery, oracle, or accuracy number.
- A Nature-main ceiling. Honest ceiling is Nature Methods; Nature main requires a paired biological discovery.
- That the niche is "uncontested." It is contested; the moat is framing and speed.

## Disavowed strong claims (may appear in the paper only as the wrong framing being corrected, never as a finding)
- "This will be a Nature paper." A target is not a claim.
- "The recoverability-frontier niche is uncontested."

### Can-claim update (2026-06-18, post Gate 0)
- A non-trivial, computable, method-independent recovery ceiling exists on synthetic data matched to real summary statistics, monotone in packing and displacement, with oracle at or above a naive baseline everywhere.

### Cannot-claim update (2026-06-18, post Gate 0)
- The real-world severity or location of the frontier (displacement sigma is unanchored until Gate 1).
- Any final accuracy number (the disjoint-marker model makes Gate 0 accuracies optimistic upper bounds).
- That cell-typing survives identity loss (the saturated profile ARI of 1.0 is a disjoint-marker artifact, not a finding).

### Gate 1 pre-committed thresholds (locked 2026-06-18; orchestrator judges, executor does not)
- KILL if, under realistic overlapping expression and the data-pinned sigma, oracle assignment accuracy in the real-data-anchored regime is greater than 0.95 (the limit does not bite once the model is realistic; the Gate 0 frontier was a gene-model artifact).
- KILL if the measured real marker-leakage cannot be reproduced by any sigma in a plausible range (0 to 15 micrometres), indicating the displacement model is misspecified and the oracle ceiling will not transfer.
- KILL if the frontier loses monotonicity, or the "limit bites in the dense regime" conclusion flips, under either realistic expression or the structural-sensitivity perturbations.
- PASS if, under realistic expression with data-pinned sigma, the frontier stays non-trivial and monotone, oracle stays at or above naive, and the real-anchored operating point sits below about 0.9 with a structured frontier present (an easy recoverable end and a hard end), and this survives the structural-sensitivity perturbations. A very low real-anchored accuracy with structure intact is a strong frontier, not a failure.
- Sanity: profile ARI must fall below 1.0 under realistic expression; if it does not, the expression model is not genuinely overlapping.

### Can-claim update (2026-06-18, preliminary, post Linux headroom)
- On realistic synthetic data with known truth, given true nuclei, three published methods (Baysor, Proseg, ComSeg) sit 0.08 to 0.11 below the oracle and tie or trail nearest-nucleus in dense tissue under the one-to-one matched metric and default-or-corrected configurations.

### Cannot-claim update (2026-06-18, post Linux headroom)
- That sophisticated methods do not beat nearest-nucleus, as a headline. It is configuration-sensitive (a 40 percent Baysor handicap was found and fixed) and rests on the one-to-one matched metric, which penalises over-segmentation. Not claimable until each method is given its best configuration and the picture holds under over-segmentation-robust metrics.

### Fairness-gate pre-committed criteria (locked 2026-06-18; orchestrator judges, executor does not)
- The "sophisticated methods do not beat nearest-nucleus in dense tissue" observation is supported only if, with each method at its best documented configuration AND under all three metrics, no method exceeds the naive baseline by more than 0.01 (within the measured method run-to-run spread of about 0.012) in the dense regime.
- The "real headroom below the achievable ceiling" claim is supported if the best fairly-configured method stays at least 0.05 below the oracle in the dense regime across the metrics.
- If any fairly-configured method beats naive by more than about 0.02, or closes the oracle gap to under about 0.03, drop the "below naive" headline and frame the paper around the bound, the frontier, the diagnostic, and the headroom to the oracle.

### Can-claim update (2026-06-19, post lever-scout, pending worked analysis)
- On the Janesick Xenium breast dataset, the triple-positive DCIS co-expression sits in a regime where transcript-to-cell assignment is above the recoverability ceiling (best-possible accuracy below 0.9), so the spatial co-expression is not separable from a misassignment artifact by any method; and the authors' own orthogonal scFFPE-seq plus the clinical PR- annotation independently indicate the co-expression is not genuine.

### Cannot-claim update (2026-06-19)
- That we have proven the triple-positive false. We have not, and cannot, via our own re-assignment; the orthogonal data is the only evidence of spuriousness, and it carries the dropout caveat. The contribution is the method-independent ceiling and the mechanism demonstration, on top of a discordance the original authors already disclosed.
- A Nature-main biological discovery from this lever. There is none; the corrected value is the pre-existing clinical annotation, and the finding is incidental.
