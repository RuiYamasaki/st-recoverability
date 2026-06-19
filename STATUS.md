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

### Segmenter-headroom execution log (2026-06-18, branch headroom)
Built and run. Executor asserts no conclusion (near-the-wall vs headroom is the orchestrator's call); raw numbers in outputs/HEADROOM_REPORT.md and committed result files. Regen: `python src/headroom_run_all.py`; table: `python src/headroom_report.py`. Gates 0-3 left byte-identical.
- Methods that ran (Python 3.13, nuclei-prior, documented defaults): pciSeq 0.0.61 (Nature Methods 2020) and ComSeg 1.8.5 (2024). Proseg (Rust) and Baysor (Julia) are toolchain-blocked on this Windows machine (no MSVC/GCC linker, no Julia, no prebuilt binaries) and were substituted by ComSeg + pciSeq. ComSeg free-segmentation mode (no prior) was attempted and failed (AssertionError); free-segmentation is therefore not represented by a successful run.
- Metric: matched transcript-assignment accuracy vs known truth (Hungarian one-to-one matching of method-cells to true-cells; equals direct id-match for oracle/naive).
- Dense regime across the sigma CI (1.43, 1.99, 2.63 um): oracle 0.803/0.730/0.651; naive 0.731/0.639/0.546; ComSeg 0.722/0.634/0.543; pciSeq 0.667/0.604/0.504. Sparse: oracle 0.882, naive 0.833, ComSeg 0.762, pciSeq 0.230. Representative: oracle 0.808, naive 0.741, ComSeg 0.730, pciSeq 0.475.
- best-sophisticated-method (ComSeg) to oracle gap: 0.081/0.096/0.108 across the dense sigma CI (0.120 sparse, 0.078 representative). Neither sophisticated method exceeds the naive nearest-nucleus baseline on this known-truth synthetic; naive-to-oracle gap 0.072/0.091/0.105 dense.

### Headroom analysis, Windows run (logged 2026-06-18): inconclusive for the field-standard methods
Baysor and Proseg were toolchain-blocked on the Windows/Git-Bash machine (no Rust, Julia, MSVC, or GCC; no Windows binaries). Substitutes pciSeq and ComSeg ran but are not the field-standard Xenium methods and appeared handicapped by the nuclei-only setup (both at or below the naive baseline; pciSeq dumped most transcripts to background in the sparse regime). Preliminary only: two published methods given true nuclei did not beat nearest-nucleus and sat about 0.08 to 0.12 below the oracle in dense tissue. Where Baysor and Proseg sit is still open.

### Decision D10 (2026-06-18): redo headroom on Linux (WSL) with Baysor and Proseg
Move to a Linux environment where the field-standard methods install cleanly, recompute oracle and naive within that environment, and place Baysor and Proseg against the ceiling in both nuclei-prior and native free-segmentation modes. This is the lever on the realistic Nature Methods ceiling and possibly the headline. Cross-platform byte-identity of prior gates is explicitly not required; the committed Windows results remain canonical.

### Headroom on Linux (logged 2026-06-18): favorable, pending fairness hardening
With Baysor, Proseg, and ComSeg installed and run on Linux (oracle and naive recomputed within Linux, matching Windows to 3 dp), all three given true nuclei tie or trail the nearest-nucleus baseline (margin -0.003 to -0.03) and sit 0.08 to 0.11 below the oracle in dense tissue; Baysor free-segmentation widens the oracle gap to 0.16 to 0.17. No method beats the oracle anywhere (no scorer bug). Baysor and ComSeg assign about 100 percent; Proseg dumps 15 to 25 percent to background (its scores are a conservative floor). A 40 percent Baysor cell-radius handicap was found and corrected before reporting, which proves the numbers are configuration-sensitive.

### Decision D11 (2026-06-18): run the method-fairness hardening gate before writing
Give each method its best documented configuration, add over-segmentation-robust metrics, and fix or scope free-segmentation for Proseg and ComSeg, so the "current methods do not beat nearest-nucleus" observation cannot be dismissed as unfair benchmarking by a method author at review. This is the final hardening before the paper.

### Lever-scout result (logged 2026-06-19): one Nature-Methods-grade lever, no Nature-main lever
Verified against primary sources. Janesick et al. (Nat Commun 14:8353, 2023, DOI 10.1038/s41467-023-43458-x) report a small DCIS region triple-positive for ESR1/PGR/ERBB2. The diagnostic on the exact dataset puts that dense regime above the ceiling (sigma about 2.1 um, oracle 0.74 to 0.80, all below 0.9). Orthogonal data contradicts the co-expression in the misassignment-predicted direction: the authors' own dissociated scFFPE-seq shows the few PGR+ cells do not co-express ESR1/ERBB2, and the block is clinically annotated PR-. The logical structure is Tier-A-shaped, but the tier is Nature Methods, not Nature main, for two reasons: the orthogonal contradiction was published by the original authors themselves (we explain a visible discordance, we do not discover a hidden false conclusion), and the finding is a minor incidental ROI, not a load-bearing biological claim with downstream impact. CRC (10x public Xenium colorectal, above ceiling) and MERFISH hypothalamus (sparser) are corroborated controls: above/in regime but orthogonal data agrees.

### Decision D12 (2026-06-19): run the worked lever analysis, target stays Nature Methods
Pursue the Janesick triple-positive as the worked Nature-Methods lever, not a Nature-main play. Establish the ROI-local ceiling, demonstrate the artifact mechanism on synthetic ground truth, document the orthogonal contradiction, and present the artifact-consistent vs corroborated comparison. The correction rests only on orthogonal data and the synthetic null, never on our own re-assignment of the real transcripts. If the worked analysis weakens any leg, drop the lever and ship the bound-plus-tool.

### Decision D13 (2026-06-19): build the four main figures; recompute Leg 2 local-to-local
Build the four main figures from committed results; recompute Leg 2 local-to-local (observed apparent-triple rate inside the actual triple-positive region vs the synthetic adjacency-conditioned rate) so the mechanism comparison is apples-to-apples. The integrity linchpin holds: no re-assignment of real transcripts is used anywhere in the figures.

### Decision D14 (2026-06-19): render the 10 Extended Data items from committed results; no new analysis
Render the 10 Extended Data display items from results that are already committed. No new analysis: everything needed is already computed and committed. The integrity linchpin holds: no re-assignment of real transcripts; the worked-example item uses the synthetic-field ceiling, the known-truth synthetic null, and the released matrix read only.
