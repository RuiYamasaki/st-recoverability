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
None yet. No project measurement has been made. Literature facts verified by search live in the references with their citations, not here.

| number (what it is) | value | source result file | column / field | date verified |
|---|---|---|---|---|
| (none yet) | | | | |

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
