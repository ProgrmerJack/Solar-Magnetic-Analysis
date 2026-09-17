# Quarantined scripts

233 scripts moved here 2026-07-28 by the audit in
`scripts/utilities/audit_scripts.py`. **Nothing was deleted** — tracked files were
moved with `git mv` so history follows them, and they can be restored with
`git mv archive/superseded_scripts/<path> scripts/<path>`.

They are quarantined rather than removed because they remain the provenance for
188 result files and for the claim→script map in `docs/REVIEWER_INDEX.md`. That
audit trail matters more than a tidy tree.

| reason | n | meaning |
|---|---|---|
| SUPERSEDED | 148 | uses a demonstrably wrong SSW event definition, or a design shown unsound |
| LIMITED | 48 | no valid SSW exposure; descriptive only |
| abandoned-hypothesis downloads | 30 | acquire POES/GOES/PSP/geomagnetic data for the discontinued solar–magnetic study |
| figures | 6 | render withdrawn manuscripts |
| EXPOSURE_SOURCE | 1 | `01_build_daily_panel.py` — builds `ssw_within_15d` from the 16-event list |

## The four wrong SSW event definitions

| source | scripts | defect |
|---|---|---|
| `ssw_catalog.parquet` | 100 | contains spurious 2012-01-11, omits MAR 2000, ends 2021 |
| `analysis_panel.parquet` (transitive) | 16 | carries `ssw_within_15d` from the 16-event list |
| `ssw_event_catalog.csv` | 15 | spurious 2012-01-11; omits MAR 2000 and MAR 2010 |
| `ssw_within_15d` named directly | 12 | same origin |
| local Charlton–Polvani detector | 3 | over-detects (40 vs 31 canonical) |

The only correct definition is `data/processed/atmospheric/ssw_canonical.csv`, built by
`scripts/download_extra/download_ssw_compendium_canonical.py`.

## Reusing anything from here

Most of these scripts are wrong only in **which events they treat as SSWs**. The analysis
logic is often sound. To revive one: repoint it at `ssw_canonical.csv`, then re-express it
as an event study with mutually exclusive lag bins (see `r113_event_study.py`), because the
single-window specification these use is separately unsound.
