# Results with no producing script

Moved out of `results/current/` on 2026-09-17. `tools/catalog_repo.py` resolves a
producer by grepping every script for the result's literal basename; these two
matched nothing, so no script in the repo can regenerate them. A number that
cannot be regenerated cannot be defended, so they do not belong in the live tree.

Neither is cited by any FINDING document or by `CONSOLIDATED_RESULTS.md` — checked
before moving. They appeared only in the generated `CATALOG.md` and
`results/README.md`.

- **`ensemble_precursor.json`** (written 2026-07-31 13:48) — a CanESM5 run from
  before `ensemble_precursor.py` was changed to write one file per model. The
  script now writes `ensemble_precursor_{model}.json`; the superseding
  `ensemble_precursor_CanESM5.json` (15:33) is in `results/current/7_ensemble/`.
  The two disagree numerically — this file gives −0.6558 σ at −30..−16 against the
  current −0.5908 σ — so it is a superseded run, not a duplicate.

- **`R_null_calibration.json`** (216 bytes) — nothing in the repo writes it and its
  numbers appear in no document. The live R diagnostics are `R_profile.json` and
  `R_gap_sensitivity.json`.
