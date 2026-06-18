# STATUS — source of truth

**Project (working title):** A recoverability frontier for transcript-to-cell assignment in imaging-based spatial transcriptomics.

**One-line objective:** Map the information-theoretic answerability frontier (when transcript-to-cell assignment is recoverable *at all*, as a function of transcript density, cell packing, and displacement/diffusion) using a by-construction synthetic ground-truth generator, and ship it as a calibrated benchmark plus an installable diagnostic.

**Why this project (floor x ceiling):**
- Floor: a sound generator + frontier + diagnostic is an adoptable benchmark the field openly lacks. Realistic floor venue: Bioinformatics, Cell Reports Methods, or PLOS Computational Biology. There is no realistic "unpublishable" outcome if the benchmark is competently built.
- Ceiling: Nature Methods (top of the realistic ladder for a methods-and-benchmark paper). Nature main only as a thin tail, and only if the frontier converts into a biological discovery (overturns a class of published spatial conclusions AND is paired with a concrete biological re-finding), the way TopACT reached Nature.

## Decision log
- **D1 (2026-06-18):** Target = spatial-transcriptomics recoverability frontier. Chosen over the MICrONS digital-twin recoverability idea (confirmatory: the structure-to-function degeneracy principle is already published by Beiran & Litwin-Kumar, Nat Neurosci 2025; plus high competition and self-overlap with prior neuro work) and over the physics options (TESS / materials / continuous-gravitational-wave), which the deep research itself caps below Nature.
- **D2 (2026-06-18):** Optimization objective = maximize floor x ceiling (not "Nature tier", not "uncontested", both of which were shown incompatible with the others).
- **D3 (2026-06-18):** Deliberate clean break from the under-review IBL neuro paper, to avoid self-overlap. Same signature method style (synthetic ground truth + recoverability framework + open diagnostic), new field, different reviewers.
- **D4 (2026-06-18):** Novelty moat = the information-theoretic answerability-frontier framing + the oracle ceiling, NOT the synthetic generator. The Voronoi+Poisson generator approach is already published (TopACT, Nature 2024), so it is validated and citable, not novel. The neighborhood is active and rising, so speed matters.
- **D5 (2026-06-18):** First gate = Gate 0 (does a non-trivial, transferable frontier exist). Thresholds locked in LOCKED_NUMBERS_AND_GUARDRAILS.md and not to be renegotiated after seeing data.

## Open items
- Pick the exact target journal on the ladder before writing a single section, then pull that journal's current rules (word/figure caps, code/reproducibility policy, reporting summary) and build a compliance checklist.
- Verify the remaining cited tools (FastReseg, Cellist, MisTIC) when building the related-work section.
- Confirm whether the under-review IBL paper is public as a preprint (affects how it is cited and how much the prior recoverability framing is already attributable).
- Choose the specific public datasets for the realism anchor (candidates: a 10x Xenium release, a Vizgen MERFISH release, a Stereo-seq release).

## Current plan / next action
Run Gate 0 on a branch (see gate0_executor_prompt). Report back in the specified format. Orchestrator audits the report against the locked thresholds; PASS proceeds to the full frontier sweep + diagnostic, KILL pivots.

## Competitor watchlist (re-check before posting anything public)
- segger (Marioni / Huber, EMBL); Proseg (Newell lab, Fred Hutch); TopACT and the Cambridge topological-data-analysis group (own a published synthetic ST generator); TRACER authors; the authors of the June 2026 segmentation review (arXiv 2606.09675). The moat against all of them is the answerability-frontier framing and execution speed, since the generator tooling already exists in the field.
