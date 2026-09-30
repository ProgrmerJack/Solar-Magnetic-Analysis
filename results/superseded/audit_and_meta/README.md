# Quarantined results

183 result files moved here 2026-07-28 by `scripts/utilities/audit_results.py`.
**Nothing deleted** — tracked files moved with `git mv`; restorable with `git mv` back.

| reason | files | citable numbers |
|---|---|---|
| SUPERSEDED — producer used a wrong SSW event definition or an unsound design | 167 | 7,164 |
| ORPHAN — no producing script exists in either tree | 16 | 321 |

Retained in `data/results/`: 2 files from the corrected pipeline (160 numbers),
2 with no statistical claim, and the audit artefacts themselves.

## The ratio that matters

**7,164 citable numbers in this directory come from superseded code, against 160
from the corrected pipeline.** Any figure quoted from the manuscript history is,
on those odds, far more likely to come from here than from validated work.

## Orphans

These 16 carry 321 numbers and have **no producing script anywhere**, including the
archive — the code that made them was deleted in an earlier round. They are not
reproducible and must not be cited. Largest: `threshold_amplification_model.json`
(213 numbers).

## Reusing anything here

Don't. Re-derive it. Every file inherits at least one of: the spurious 2012-01-11
event, the omission of MAR 2000 / MAR 2010, the over-detecting local detector, or
the single-window specification shown unsound in `r113_event_study.py`.
