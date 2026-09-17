# Catalogue freeze record

Built by `build_catalogue.py` from `data/atmospheric/ssw_catalog/majorevents_raw.html`
(NOAA CSL compendium). Rule written in `CATALOGUE_SPECIFICATION.md` **before** the
catalogue was built or any lead–lag profile inspected.

**Re-frozen 2026-07-30** after the consensus rule was found to be era-dependent. See
§"Why the rule changed" below; the original 39-event freeze is superseded.

## What the frozen rule yields

| quantity | value |
|---|---|
| **N events (primary)** | **43** |
| **N unique winters** | **36** |
| event catalogue period | **1958-01-30 .. 2024-03-04** |
| winters containing >1 event | 7 — 1966, 1969, 1971, 1988, 1999, 2010, 2024 |
| rows in audit table (incl. excluded) | 48 |

`N events ≠ N unique winters` (43 vs 36) is reported everywhere, and the winter is the
bootstrap block, never the event.

## Why the rule changed: consensus is a FRACTION, not a count

The original rule was "detected by ≥4 of 6 reanalyses". Only four reanalyses cover
1958–1978 (NCEP-NCAR, ERA40, JRA-55, ERA5); ERA-Interim starts 1979 and MERRA-2 in 1980.
So "≥4 of 6" silently demanded **unanimity** before 1979 while asking only 67% in
1979–2002. The rule became stricter the further back it looked.

The damage was measurable on the previous build: mean consensus 3.3 (1958–78) against 5.4
(1979–2002), and **6 of 8 excluded "marginal" events were pre-1979, Fisher P = 0.011**.
"Marginal" had become a proxy for "pre-satellite", which would have put an era confound
into every catalogue-sensitivity comparison in the project.

The rule is now **consensus ≥ 2/3 of the reanalyses that actually cover the date**
(`CONSENSUS_FRACTION = 2/3`, with a `COVERAGE` table in `build_catalogue.py`). This is
era-independent by construction.

Effect: four events admitted — JAN 1963, JAN 1968, MAR 1969, JAN 1977, each 3/4 = 0.75 —
and the exclusions are now era-balanced at **2 pre-1979 / 2 post-1979**.

Remaining exclusions, all with the reason recorded in the audit table:

| event | detections | available | fraction |
|---|---|---|---|
| NOV 1958 | 1 | 4 | 0.250 |
| MAR 1965 | 1 | 4 | 0.250 |
| FEB 1981 | 1 | 6 | 0.167 |
| FEB 2002 | 3 | 6 | 0.500 |

## Frozen sets

| set | n | rule |
|---|---|---|
| `primary` | 43 | consensus ≥ 2/3 of covering products, median central date |
| `primary_compendium_only` | 41 | primary minus the two literature-sourced 2024 events |
| `consensus_strict` | 36 | detected by every product covering the date (fraction ≥ 0.99) |
| `consensus_half` | 42 | fraction ≥ 0.5 |
| `union` | 47 | ≥1 reanalysis (permissive) |
| `era5` | 42 | ERA5 detections |
| `jra_55` | 41 | JRA-55 detections |
| `ncep_ncar` | 39 | NCEP-NCAR detections |
| `era40` / `era_interim` / `merra2` | 29 / 26 / 27 | shorter records |

`consensus5` (the old ≥5/6 set) no longer exists — a raw count is not meaningful under a
fraction rule. It was replaced by `consensus_strict` and `consensus_half`.

## Invariance tests (spec §9): all passed

no duplicate ids or dates · minimum 20-day separation respected · no final warming in
primary · every central date within 1 Nov–31 Mar · sources all recognised · set lengths
consistent · subset relations hold.

The subset test was corrected at the same time: literature-sourced events carry no
consensus fraction, so `primary ⊆ consensus_half` cannot hold. The test now checks
`primary_compendium_only ⊆ consensus_half`, which is the meaningful relation.

March 2025 is carried in the audit table with `final_warming_flag = True` and an explicit
`exclusion_reason`, so the exclusion is visible rather than silent.

## Checksums

```
source_html_sha256        4147ac39fc0226c21c997aa9d68b7a0517e4bb1e754e81c4cb9d1ee456dbf028
event_catalogue_csv_sha256 76ad09e49bb9ac5888b1b970aaa74bb0a5be0d5c9f68265c93e251251ec06729
```

Recorded in `catalogue_checksum.txt` — sha256 of both the raw source HTML and the built
table, so any drift is detectable.

**The table checksum changed on 2026-09-17, and the catalogue did not.** The frozen file
was built on Windows, where `pandas.to_csv` defaults to CRLF; rebuilt on Linux it wrote
LF, so the file differed by exactly one byte per row (7,943 → 7,894) and the old hash
`11c2d4096919bffa13ead923ffa651193346bb228b9b60bfb2e5c31a431120f7` no longer matched.
Verified identical as data before the hash was re-recorded: same 48 rows × 19 columns,
`DataFrame.equals` True, same 43 primary events, all §9 invariance tests passing, and
`catalogue_sets.json` rebuilt byte-identical because it was already written with an
explicit `\n`. `to_csv` now pins `lineterminator="\n"`, so the checksum is
platform-independent and two consecutive rebuilds are byte-identical.

## Standing rule

`load_catalogue(which)` in `build_catalogue.py` is the **only** permitted way to obtain
events. Inline date lists are forbidden. Four mutually incompatible catalogues in the
superseded project came from exactly that practice.

## Downstream

Every catalogue-dependent result was recomputed on the 43-event primary after the
re-freeze: `canonical_event_study.py`, `multi_index_event_study.py`,
`design_sensitivity.py`. Consequences are recorded in
`../FINDING_design_robustness.md`.
