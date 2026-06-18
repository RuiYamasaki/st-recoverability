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
