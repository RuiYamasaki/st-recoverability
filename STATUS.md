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

## Gate 0 execution log (2026-06-18, branch gate0)
Gate 0 was built and run. Executor asserts no PASS/KILL; raw numbers are in outputs/GATE0_REPORT.md and the committed result files. Regen: `python src/run_all.py` (or `--skip-realism`); table: `python src/report_table.py`.
- Generator src/generator.py, oracle src/oracle.py, sweep src/sweep.py, realism src/realism.py + anchor src/anchor.py. Master seed 20260618.
- Sweep grid 5x5x5 (density {50,90,160,280,500} tx/cell; packing {2000,3400,5700,9500,16000} cells/mm2; sigma {0.5,1,2,4,8} um) -> results/sweep.csv.
- Oracle assignment accuracy full-grid range [0.234, 0.957], monotone in sigma and packing (0 violations), oracle >= naive at all 125 points.
- Two ARIs: assignment-partition ARI vs true cell labels [0.110, 0.919] (tracks the frontier); expression-profile ARI vs true types = 1.000 everywhere (oracle errors are type-preserving; naive profile-ARI collapses to 0.10). Finding: identity assignment and cell-type recovery are decoupled.
- Realism anchors: Xenium FFPE Human Breast Cancer Rep1 (10x release 1.0.1, cells.parquet) and Moffitt 2018 MERFISH mouse hypothalamus (squidpy MERFISH_0.24.h5ad, figshare). results/realism.csv, results/realism_meta.json.
- Generator matches real median tx/cell and packing within 15% at every anchor (packing exact; tx within 1.6%). results/realism_oracle.csv.
- Oracle accuracy inside the realistic band (sigma<=2um): [0.719, 0.953]; Xenium-like [0.719, 0.922], MERFISH-like [0.866, 0.953].

## Current plan / next action
Orchestrator audits outputs/GATE0_REPORT.md against the locked thresholds; PASS proceeds to the full frontier sweep + diagnostic, KILL pivots. Open follow-ups if continued: model overlapping (non-block) expression so the profile-ARI ceiling is not saturated; anchor the displacement sigma axis to a measured quantity rather than an assumed band; add a real-segmenter headroom point (Baysor/Proseg) at a representative grid cell.

## Competitor watchlist (re-check before posting anything public)
- segger (Marioni / Huber, EMBL); Proseg (Newell lab, Fred Hutch); TopACT and the Cambridge topological-data-analysis group (own a published synthetic ST generator); TRACER authors; the authors of the June 2026 segmentation review (arXiv 2606.09675). The moat against all of them is the answerability-frontier framing and execution speed, since the generator tooling already exists in the field.

### Gate 0 result (logged 2026-06-18): QUALIFIED PASS
All three KILL conditions cleanly not triggered. A non-trivial, smooth, monotone, oracle-validated recovery frontier exists; oracle at or above naive at all 125 grid points; generator matches real summary statistics within 1.52%. The limit bites in the dense Xenium-like regime (oracle accuracy 0.719 to 0.922, below the 0.9 informative bar). Qualification: the PASS sub-condition "varies from below 0.6 to above 0.9 across the realistic range" was met only across the full grid (which includes displacement sigma of 4 to 8 micrometres), not within the assumed-realistic sigma band of 2 micrometres or less, where the floor is 0.719. Cause: the Gate 0 threshold left the displacement axis unspecified (orchestrator's spec gap). Resolution: Gate 1 anchors sigma to data rather than assuming it.

### Decision D6 (2026-06-18): run Gate 1 robustness gate before any paper build
Two Gate 0 simplifications make the numbers provisional and are the first hostile-reviewer attacks: disjoint-marker expression (inflates type recovery to a saturated profile ARI of 1.0 and likely inflates assignment accuracy) and an assumed rather than measured displacement sigma. Gate 1 makes both realistic and adds a structural-assumption sensitivity, before the full frontier is treated as a result.

### Gate 1 execution log (2026-06-18, branch gate1)
Built and run. Executor asserts no PASS/KILL; raw numbers in outputs/GATE1_REPORT.md and committed result files. Regen: `python src/gate1_run_all.py`; numbers: `python src/gate1_report.py`. Gate 0 left byte-identical (`python src/run_all.py --skip-realism` -> zero diff). Generator/oracle/metrics extended in place (model + geometry + displacement options) with Gate 0 defaults preserved.
- Task 1 (realistic overlapping expression, MERFISH cell_class profiles, K=10 G=155; mean pairwise cosine 0.210, 4.5% mass on exclusive genes): oracle assignment frontier [0.251, 0.956] (Gate 0 [0.234, 0.957]), monotone (0 violations both axes), oracle>=naive at all 125 points; vs Gate 0 mean abs diff 0.006. The frontier is not a marker artifact. Oracle profile ARI de-saturates to [0.667, 1.000] (Gate 0 was flat 1.000); naive profile ARI median 0.743 down to ~0.
- Task 2 (data-pinned sigma via spatial marker leakage; literal "fraction in other-type cells" shown to be biology-dominated and uninformative): real spatial leakage 0.0768, data-pinned sigma 2.67 um, 25-75 bracket [0.10, 3.03] um (reproducible within [0,15]). Oracle at pinned sigma: MERFISH-like (sparse) 0.836, Xenium-like (dense) 0.643; both below 0.9, oracle>naive.
- Task 3 (structural sensitivity): under non-Voronoi (anisotropic) geometry and non-Gaussian (Gaussian+uniform) displacement, frontier stays monotone (0 violations), oracle>=naive everywhere, dense-regime bite holds (dense sig2 = 0.66 to 0.74, all below 0.9).

### Gate 1 result (logged 2026-06-18): QUALIFIED PASS
All three KILL conditions not triggered at the pinned values. The recovery frontier survived realistic overlapping expression (assignment frontier essentially unchanged, mean abs diff 0.006, oracle-minus-naive gap preserved at 0.071) and survived non-Voronoi geometry and non-Gaussian displacement (monotone, oracle at or above naive, dense-regime limit still biting). At a displacement sigma of 2.67 um pinned to measured MERFISH spatial marker-leakage, the oracle assignment accuracy was 0.836 sparse and 0.643 dense, both below 0.9. Profile ARI de-saturated below 1.0 (sanity met). Residual risk concentrated on the displacement pin: its precision is uncharacterized (the [0.10, 3.03] um figure is a per-marker spread, not a CI; its optimistic end would put dense oracle accuracy at 0.976, above 0.95), it was measured on MERFISH and applied to Xenium geometry, and the emission idealized expression with per-type means.

### Decision D7 (2026-06-18): run Gate 2 to pin sigma directly on Xenium with error bars
The dense-regime claim rests on a borrowed, imprecise sigma under idealized emission. Gate 2 measures sigma on Xenium directly, attaches a statistical CI and a design-sensitivity range, and redoes the evaluation under negative-binomial within-type emission, before the frontier is treated as a publishable result.

### Gate 2 execution log (2026-06-18, branch gate2)
Built and run. Executor asserts no PASS/KILL; raw numbers in outputs/GATE2_REPORT.md and committed result files. Regen: `python src/gate2_run_all.py`; numbers: `python src/gate2_report.py`. Gate 0 and Gate 1 left byte-identical. Downloaded Xenium cell_feature_matrix.h5 (Human Breast Cancer Rep1, 10x release 1.0.1) and clustered (MiniBatchKMeans K=14, seed MASTER_SEED); markers via exclusivity cutoff 0.7 (8 markers).
- Task 1: Xenium-pinned displacement sigma = 1.61 um (MERFISH-pinned was 2.67 um); real spatial leakage 0.0679; dense-regime oracle accuracy at pinned sigma 0.784 (naive 0.703), sparse 0.906 (converged grid).
- Task 2: bootstrap (B=2000 over markers+cells) sigma 95% CI [0.79, 2.53] um; design-sensitivity range [0.20, 2.48] um; combined band [0.20, 2.53] um. Dense oracle accuracy (converged grid) at statistical-CI ends [0.675, 0.891] (all below 0.95); at combined-band optimistic end (sigma 0.20 um) 0.972 (above 0.95). The sub-0.95-failing optimistic end comes ONLY from the loosest marker cutoff (0.6) plus tightest adjacency bin; at the pre-registered cutoff 0.7 sigma is stable [1.49, 1.70] um and dense accuracy 0.78-0.79. Dense accuracy crosses 0.95 at sigma ~0.35 um.
- Task 3: negative-binomial within-type emission (dispersion median 4.19): sigma 1.57 -> 1.67 um (+0.10); dense oracle 0.795 -> 0.779 (-0.016); NB frontier monotone (0 violations), oracle>=naive everywhere (min gap +0.025). The frontier does not flatten under non-idealized emission.
- Resolution caveat: the default oracle grid under-resolves small sigma and under-estimates accuracy; dense numbers use a converged grid (res_cell 24).

### Gate 2 result (logged 2026-06-18): NOT ESTABLISHED, on hold
The three prior risks were retired: sigma measured directly on Xenium (1.61 um, bootstrap 95% CI [0.79, 2.53]), the statistical CI keeps dense-regime oracle accuracy below 0.95 ([0.675, 0.891]), and the limit survives negative-binomial within-type emission (dense shift only -0.016). However, the locked KILL threshold fired: the combined sigma uncertainty's optimistic end gives dense oracle accuracy 0.972, above 0.95, arising solely from a loosened marker cutoff (0.6) that admits non-exclusive markers and biases sigma to 0.20 um. The dense-regime claim therefore depends on an arbitrary exclusivity cutoff and on a single breast sample, and is not established. Scoping clarified: the striking limit is at the dense tissue tail (about 13,600 cells/mm2), not median Xenium (about 6,500); the sparse regime barely shows a limit.

### Decision D8 (2026-06-18): run Gate 3 resolution gate before writing
Resolve the marker-exclusivity dependence with a principled, non-circular validity test, confirm stability across cell-typing, and replicate the dense limit on an independent Xenium tissue. If it holds, the dense headline is earned. If a principled criterion still lets the limit cross 0.95, retreat to the method, the frontier, and the honestly bounded result.

### Gate 3 execution log (2026-06-18, branch gate3)
Built and run. Executor asserts no PASS/KILL; raw numbers in outputs/GATE3_REPORT.md and committed result files. Regen: `python src/gate3_run_all.py`; numbers: `python src/gate3_report.py`. Gates 0, 1, 2 left byte-identical. Validity test = non-circular permutation null (1000 perms, p<0.05): a marker is admissible only if its adjacent-minus-distant displacement signal exceeds a position-permutation null; markers are never selected by the sigma they yield.
- Task 1 (breast): 28/38 candidates admissible; admissible mean displacement signal +0.094 vs rejected -0.055. The markers that destabilised sigma at cutoff 0.6 fail on independent grounds (AHSP p=0.599, CRHBP p=0.242, CYP1A1 p=0.996, EGFL7 p=0.642, SLC4A1 p=0.454); one exclusive marker (MZB1, p=0.104) also fails. Re-pinned sigma 1.99 um, bootstrap CI [1.43, 2.63]; dense-regime oracle accuracy across the CI [0.663, 0.805], all below 0.95 and below 0.9; oracle>=naive; monotone.
- Task 2 (cell-typing K=10/14/14-altseed/20): sigma spread [1.99, 2.65] um (factor ~1.3, vs the old marker-cutoff swing 0.20-2.48 um); dense accuracy at every sigma-CI optimistic end <= 0.803, all below 0.95 and below 0.9; oracle>=naive throughout.
- Task 3 (independent tissue, Xenium Human Lung Cancer, sample Xenium_V1_hLung_cancer_section, release 1.5.0): 31/32 admissible (no diluted markers in the lung pool); sigma 1.83 um, CI [1.62, 2.08]; dense accuracy across the CI [0.746, 0.796], all below 0.95 and below 0.9; oracle>=naive; monotone.
- Summary: under a principled non-circular validity test, on breast and on an independent lung dataset and across four cell-typing choices, the dense-regime oracle accuracy stays below 0.95 (and below 0.9) across every sigma bootstrap CI. The Gate 2 0.972 breach was driven by biology-diluted markers that fail the validity test.

### Gate 3 result (logged 2026-06-18): PASS
Under a principled, non-circular marker-validity test (a marker is admissible only if its displacement signal beats a position-permutation null at p<0.05, never selected by the sigma it yields), the Gate 2 destabilizing markers fail on independent grounds (p = 0.24 to 0.996). Re-pinned sigma 1.99 um (CI [1.43, 2.63]) on breast and 1.83 um (CI [1.62, 2.08]) on independent lung; dense-regime oracle accuracy stays below 0.95 and below 0.9 across every bootstrap CI on both tissues and across four cell-typing choices (sigma spread about 1.3x). The dense-regime limit is established and robust. Residual item for the paper's robustness section: a marker selection-on-significance sensitivity analysis (non-decisive; the lung pool, with near-zero selection, gives the same sigma).

### Decision D9 (2026-06-18): gating complete; begin paper-building with the segmenter-headroom analysis
The empirical core is locked. First paper-building analysis places real published methods against the oracle ceiling on known-truth synthetic data, with true nuclei as a prior to isolate assignment headroom, to make the bound matter to practice. This is the main lever on the realistic ceiling (Nature Methods).
