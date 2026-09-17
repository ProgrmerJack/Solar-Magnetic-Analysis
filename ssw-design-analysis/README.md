# Study design reshapes estimates of sudden stratospheric warming impacts

Clean-room rebuild. **Nothing here inherits from the superseded project** — no script,
result or event catalogue is carried over except `ssw_canonical.csv`, which was rebuilt
from the NOAA CSL compendium and independently verified.

Target: Nature Geoscience **Analysis** (≤3,000 words, 200-word abstract, 4–6 display items).

## Layout

**Code is grouped by pipeline stage. Every result goes in one folder.**

**Results do not live in this directory.** They are in the repo-root `results/`
tree, shared with the superseded project — see `../results/README.md`.

| path | contents |
|---|---|
| **`../results/current/<theme>/`** | **every computed result, 53 files in 9 themed folders.** No result is written anywhere else. |
| `../results/_ALL_RESULTS.json` | all 53 payloads in a single file — the read path when you want the numbers, not the provenance |
| `../results/_PROVENANCE.json` | producer + sha256 baseline; drives the staleness check |
| `01_search_and_screen/` | preregistration, search queries, PRISMA screening, coding sheet |
| `02_event_catalogues/` | catalogue spec, freeze record, `build_catalogue.py` |
| `03_data_ingestion/` | raw acquisition with checksums |
| `04_reproduce_published/` | reproduction of published estimates |
| `05_corrected_estimators/` | the common event-study framework |
| `06_simulation_validation/` | estimator behaviour under known truth |
| `07_physical_decomposition/` | mechanism, predictability, SNAPSI |
| `08_literature_audit/` | recomputation of published criteria |
| `08_figures/`, `09_tables/` | display items |

Scripts stay in their stage folders because 147 import sites resolve modules by
stage path (`sys.path.insert(... "05_corrected_estimators")`) and the numbering
encodes pipeline order. Only the **outputs** were centralised; every script
writes to `RESULTS / "name.json"` where
`RESULTS = HERE.parents[1] / "results" / "current" / "<theme>"`.

## Catalogue

`python tools/catalog_repo.py` (from the repo root) rebuilds `CATALOG.md`,
`CATALOG.json` and `results/_ALL_RESULTS.json` from the files themselves. It
reports **stale** results (producing script changed since the result was written,
detected by sha256 — not mtime) and **orphans** (no script writes them). Run
`--record` only *after* re-running an analysis, to re-baseline.

## Reading order

1. `CONSOLIDATED_RESULTS.md` — the source of truth for what is established
2. `REVIEW_BRIEF.md` — results A–E as framed for review, and the errors caught
3. `GO_NO_GO.md` — the venue verdict and the 12-condition rubric
4. `02_event_catalogues/FREEZE_RECORD.md` — reconciles every event count

## Gate status

Go/no-go conditions are tracked in `GO_NO_GO.md`. **Condition 3 (pseudo-onset calibration
centred on zero) is the current blocker** and is being worked in `06_simulation_validation/`.

## Rules carried forward

- Never fit a single lag window; always mutually exclusive bins with an explicit baseline.
- Never use model standard errors; winter-block bootstrap, resampling whole winters.
- Report the whole lead–lag profile, never a favourable subset.
- Report the number of events **and** the number of unique winters.
- No claim of "causal", "inflated" or "upper bound" without evidence for the sign.
