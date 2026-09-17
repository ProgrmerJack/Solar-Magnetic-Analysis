# Event catalogue specification

**Written before any catalogue was built or any lead–lag profile inspected.**
Once frozen, changes require a dated entry in `AMENDMENTS.md` stating whether the change
was made before or after seeing affected results.

---

## 1. Event definition

A **major mid-winter sudden stratospheric warming** is defined by the Charlton–Polvani
(2007) criterion as implemented in the NOAA CSL SSW Compendium (Butler et al. 2017,
*ESSD* 9, 63–76):

| element | rule |
|---|---|
| diagnostic | zonal-mean zonal wind |
| level | 10 hPa |
| latitude | 60° N |
| criterion | reversal to easterly (u < 0 m s⁻¹) |
| central date | first day of the reversal |
| season | 1 November – 31 March (mid-winter) |
| separation | consecutive central dates must be separated by ≥ 20 consecutive westerly days; otherwise treated as one event |
| final warmings | excluded — winds must return to westerly for ≥ 10 consecutive days before 30 April |
| event type | displacement and split events are **both included**; type is recorded, not used for inclusion |

Type classification is carried as metadata only. No analysis conditions on it without a
separate preregistered hypothesis.

## 2. Primary catalogue rule

> **Primary = the cross-reanalysis consensus catalogue: events detected as major
> mid-winter warmings by ≥ 2/3 of the compendium reanalyses that COVER the date
> (NCEP-NCAR 1948–, ERA-40 1957–2002, ERA-Interim 1979–2019, JRA-55 1958–,
> MERRA-2 1980–, ERA5 1940–). Central date = median of the qualifying
> reanalysis dates, rounded down to the earlier day when the median falls between two.**

**Amended 2026-07-30.** The rule was originally "≥ 4 of 6". That is era-dependent and had
to be replaced: only four reanalyses cover 1958–1978, so an absolute count of 4 demanded
unanimity before 1979 while asking 67% in 1979–2002. Measured on the resulting build, mean
consensus was 3.3 (1958–78) against 5.4 (1979–2002), and 6 of 8 excluded events were
pre-1979 (Fisher P = 0.011) — "marginal" had become a proxy for "pre-satellite", which
would have placed an era confound inside every catalogue comparison in the project. The
fraction form is era-independent by construction. Consequences are recorded in
`FREEZE_RECORD.md`.

Justification, fixed in advance and independent of any outcome:

1. A consensus rule does not privilege one reanalysis's detection idiosyncrasies.
2. It is reproducible from a single published table.
3. It weights an event detected consistently across products above one detected by a single
   product — unlike a union catalogue.
4. The median central date is robust to a single reanalysis dating an event early or late.
5. Expressing consensus as a fraction of *available* products keeps the threshold the same
   in every era, so catalogue membership does not encode observing-system history.

**Rejected as primary, with reasons:**
- *"ERA5 where available"* — not a rule; behaviour is undefined when ERA5 does not detect an
  event. This ambiguity is exactly what invalidated the superseded catalogue.
- *Union across reanalyses* — treats a one-reanalysis detection as equal to a six-reanalysis
  detection. Retained only as a deliberately permissive sensitivity set.
- *ERA5-only* — defensible, but discards 1958–1978 and makes the modern and historical
  periods non-comparable. Retained as a sensitivity set.

## 3. Sensitivity catalogues (all frozen together with the primary)

| id | rule |
|---|---|
| `primary` | consensus ≥ 2/3 of covering products, median central date |
| `primary_compendium_only` | primary minus the literature-sourced 2024 events |
| `era5` | ERA5 detections only |
| `jra_55` | JRA-55 detections only |
| `ncep_ncar` | NCEP-NCAR detections only (the reanalysis behind this repo's wind series) |
| `merra2` | MERRA-2 detections only |
| `consensus_strict` | detected by every product covering the date (fraction ≥ 0.99) |
| `consensus_half` | fraction ≥ 0.5 |
| `union` | detected by ≥ 1 reanalysis (permissive) |

`consensus5` (≥5/6) was removed in the 2026-07-30 amendment: a raw count is not meaningful
once consensus is a fraction. `consensus_strict` and `consensus_half` replace it.

Every analysis is run on `primary`; the seven-way comparison is reported for the canonical
AO/NAO event study and for the headline replication results.

## 4. Central-date disagreement

Primary uses the median qualifying date. Two additional checks:

1. **Per-reanalysis dates** — repeat the canonical event study using each reanalysis's own
   central date, giving a spread of profiles.
2. **Date jitter** — shift every central date by −3, −2, −1, +1, +2, +3 days and refit, to
   quantify sensitivity to modest displacement.

The catalogue is fixed once. It is **never** varied by outcome.

## 5. Multiple events in one winter

- Both events are retained. Nothing is dropped for sharing a winter.
- Days are assigned to the bin of the **nearest** central date, so bins remain mutually
  exclusive; no day contributes to two events.
- If two central dates are closer than 20 westerly days apart they are one event by §1.
- **Bootstrap resampling is by winter, never by event.** The winter is the independent
  block; two events in one winter are not two independent observations.
- Both counts are reported everywhere: **N events and N unique winters.**

## 6. Periods — three distinct quantities, never conflated

| quantity | meaning |
|---|---|
| **event catalogue period** | first to last central date in the catalogue |
| **outcome-data availability period** | span of the outcome series |
| **analysis overlap period** | intersection, and what actually enters a fit |

The catalogue period ends at the last included major mid-winter warming. Outcome data
extending beyond it does **not** extend the catalogue period, and the catalogue is never
described by the outcome span.

## 7. Post-compendium events

The compendium table ends at FEB 2023. Later events enter only if documented as major
mid-winter warmings in the peer-reviewed literature, flagged by `source`, and satisfying §1.

Currently: **2024-01-16** and **2024-03-04** (Lee et al. 2025, *Weather*,
doi:10.1002/wea.7656). The March 2025 event is **excluded** — documented as that winter's
early *final* warming, which §1 excludes.

Literature-sourced events have no per-reanalysis detection counts and therefore cannot meet
the ≥4/6 consensus rule mechanically. They are included in `primary` with
`consensus_count = NA` and `source = literature`, and a sensitivity catalogue
`primary_compendium_only` excludes them so their influence is always visible.

## 8. Output

One machine-readable table, `event_catalogue.csv`, with columns:

```
event_id, winter, primary_central_date, event_type, final_warming_flag,
primary_inclusion, ncep_ncar, era40, era_interim, jra55, merra2, era5,
consensus_count, source, exclusion_reason, sha256_of_source
```

The manuscript table is generated from this file. **No second, manually edited event list
may exist** — the failure that produced four incompatible catalogues in the superseded work.

## 9. Invariance tests (must all pass before freezing)

1. no duplicate `event_id`, no duplicate central date
2. minimum-separation rule satisfied between consecutive events
3. no event flagged as a final warming is included
4. every central date falls in 1 Nov – 31 Mar
5. every `primary` event traces to the compendium or a cited paper
6. counts reproduce exactly on a clean re-run from the raw source
7. every analysis imports the catalogue from one module — no inline date lists
8. file checksum recorded and archived
