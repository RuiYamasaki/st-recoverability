# Recoverability frontier for transcript-to-cell assignment in spatial transcriptomics

A by-construction synthetic ground-truth study mapping when transcript-to-cell assignment is recoverable at all, as a function of transcript density, cell packing, and displacement. Solo reanalysis and methods work on open data; commodity hardware only.

## Repository structure
```
.
├── README.md
├── STATUS.md                         # source of truth: decisions, plan, open items
├── LOCKED_NUMBERS_AND_GUARDRAILS.md  # trusted numbers, guardrails, can/cannot-claim
├── requirements.txt                  # starter deps (pin after first install, see below)
├── requirements.lock.txt             # exact pinned env (created on first install)
├── src/                              # generator, oracle, metrics, sweep, figures
├── data/                             # NOT committed: streamed public FOVs cached locally
├── results/                          # committed: small result CSVs / JSON only
├── figures/                          # committed: generated figures
└── outputs/                          # final deliverables only
```

## Environment (pin from commit one)
Starter `requirements.txt`:
```
numpy
scipy
pandas
matplotlib
anndata
scikit-learn
```
After the first install on the laptop, freeze the exact versions and commit the lock file:
```
pip install -r requirements.txt
pip freeze > requirements.lock.txt
git add requirements.lock.txt && git commit -m "Pin environment"
```
External segmenters (Baysor, Proseg) are optional for Gate 0 and are installed separately only if the headroom comparison is run.

## Reproducibility
- One command regenerates every result and figure: `python src/run_all.py` (or a `make all` target) re-derives `results/` and `figures/` from committed code and cached/streamed data.
- All randomness uses fixed, recorded seeds.
- `results/` holds only small derived files (summary stats, sweep tables). Raw public data is streamed and cached under `data/`, never committed.
- `outputs/` contains real deliverables only (manuscript, supplementary, the diagnostic package).

## Data
Public datasets are streamed from their portals (10x Xenium, Vizgen MERFISH, Stereo-seq). Dataset IDs and access routes are recorded in STATUS.md and the manuscript. Only small derived summary statistics are committed.
