# SSW surface-impact methodology

Global rules in `~/.claude/CLAUDE.md` apply. This file covers what is specific to this repo.

## What this is

Whether the literature's split of sudden stratospheric warmings into
"downward-propagating" and not carries information about the Earth system, or is a
threshold applied to one continuous population. Evidence comes from observations
(n=43 events), CMIP6 (n=1888), and the SNAPSI nudged-ensemble experiment.

**The repository name is a fossil.** It was a solar/geomagnetic–avalanche project
that was abandoned; see "Superseded" below. Nothing solar is live.

## Layout

```
ssw-design-analysis/     THE LIVE PROJECT. Numbered stages, executed in order.
  01_search_and_screen/    preregistration; why the literature census was dropped
  02_event_catalogues/     build_catalogue.py -- the ONLY event source
  03_data_ingestion/       acquire_*.py, plus the reduced data they write
  05_corrected_estimators/ event studies on the frozen catalogue
  06_simulation_validation/ nulls, calibration, selection-bias simulation
  07_physical_decomposition/ mechanism, predictability, SNAPSI analyses
  08_literature_audit/     recomputing published criteria on our data
  09_figures/ 10_tables/   display items; scripts READ results, never compute them
  manuscript/              main.md, the submission draft (written from CONSOLIDATED_RESULTS.md)
  environment/             the pinned environment and its check
  run_logs/<stage>/        stdout of every run, mirroring the stage layout
  tests/smoke_test.py      the standing checks
results/
  current/<theme>/         56 live result JSONs, grouped by question answered
  superseded/<family>/     retired results, kept for provenance
  _ALL_RESULTS.json        every live payload in one file -- the read path
  _PROVENANCE.json         producer path + sha256 baseline, drives staleness
tools/catalog_repo.py      rebuilds CATALOG.md/json; never edit those by hand
data/                      54 GB, almost all of it the superseded project
archive/                   everything dead, kept rather than deleted
```

**Stage 04 does not exist.** It was `04_reproduce_published`, struck when the
literature census was dropped (`01_search_and_screen/REMOVING_THE_HUMAN_DEPENDENCY.md`).
The gap is deliberate: ~147 sites import sibling stages by directory name, so
renumbering to close it would break the imports.

## Data

- **Events**: NOAA CSL compendium HTML at `data/atmospheric/ssw_catalog/majorevents_raw.html`,
  sha256 `4147ac39…`, frozen. Built into `02_event_catalogues/event_catalogue.csv`.
- **Reanalysis**: NCEP 1958–2024 via THREDDS NCSS; ERA5 via WeatherBench2 public GCS
  zarr (anonymous, no CDS account); CMIP6 zonal means, 20 usable members.
- **Experiment**: SNAPSI at `https://dap.ceda.ac.uk/badc/snap/data/post-cmip6/SNAPSI/`.
  Browsing is anonymous; downloads need a free CEDA account.
- Raw data under `data/` is **read-only**. Derived outputs go to `results/current/<theme>/`
  and to the reduced caches under `03_data_ingestion/` (all gitignored).

## Run

```bash
python3 -m venv .venv
.venv/bin/pip install -r ssw-design-analysis/environment/requirements.lock.txt

.venv/bin/python ssw-design-analysis/environment/check_environment.py   # must pass first
.venv/bin/python ssw-design-analysis/tests/smoke_test.py                # standing checks
.venv/bin/python tools/catalog_repo.py                                  # rebuild CATALOG.*
```

There is no single pipeline command; each stage script is run individually and writes
one result JSON. Run with the venv interpreter — the system interpreter has an xarray
with no NetCDF backend and cannot open any `.nc` in this repo.

## The pattern every new analysis follows

1. **Survey first.** Grep `results/current/` and the stage dirs for something that
   already answers the question. Improve it in place rather than adding a second script.
2. **One script, one question, one JSON.** Script goes in the stage that matches what it
   does. Its result goes to `results/current/<theme>/<script_name>.json`.
3. **Get events only from `load_catalogue(which)`.** Inline date lists are forbidden —
   four mutually incompatible catalogues in the superseded project came from that.
   `primary` = 43 events / 36 winters is the frozen default; anything else is a
   sensitivity, never a substitution.
4. **Record reproducibility metadata in the JSON**: seed, `n_boot`, input paths, the
   counts that define the sample. `tests/smoke_test.py` reads these.
5. **Tee stdout to `run_logs/<stage>/<script>.log`.** That log is the evidence the run
   happened.
6. **Re-run `tools/catalog_repo.py`**, then `--record` only after the result is known good.
7. **Write the finding up** as `FINDING_<topic>.md` in the same stage, and fold the
   headline into `ssw-design-analysis/CONSOLIDATED_RESULTS.md`, which is the single
   source of truth. Individual FINDING docs go stale; that file is maintained.
8. **If it failed in a way that can recur, append four lines to `FAILURES.md`.**

## Conventions specific to this repo

- `N_BOOT >= 1000` for any reported CI or p-value; several thousand near a threshold.
- Seed with `zlib.crc32(name.encode()) % N`, never `hash(str)`.
- **Compare results parsed, never byte-wise.** Every frozen artifact here was written
  on Windows and is CRLF (50 of 54 live JSONs); a Linux re-run writes LF and differs in
  bytes while being numerically identical. Pin `lineterminator="\n"` on `to_csv` and
  `newline="\n"` on `write_text` in any script you touch — but do not mass-edit all 60
  write sites, which would restamp every producer hash and mark all 53 results stale.
- A negative control is judged against its own FDR family, never pooled with positives.
- `CATALOG.md` / `CATALOG.json` are generated. Edit the tool, not the output.
- Staleness is by producer **content hash**, not mtime — a bulk edit rewrites every
  mtime and would mark everything stale.

## Known pitfalls

- **`CONSOLIDATED_RESULTS.md` is the source of truth; `paper/ng_manuscript.tex` is not.**
  The manuscript predates every current result and still leads with the abandoned
  avalanche application.
- **Six FINDING docs are stale**, written against the pre-2026-07-30 39-event catalogue.
  `05_corrected_estimators/FINDING_canonical.md` is worse than stale — its headline is
  reversed by its own JSON. Listed in `CONSOLIDATED_RESULTS.md` §5b.
- **Two Gate 3 artifacts disagree.** `gate3_clean_null.json` is correct (0.930/0.070,
  passes). `gate3_final.json` is output of the flawed `finalise_gate3.py`.
- Full log of traps: `FAILURES.md`.

## Superseded

`paper/`, `notebooks/`, `scripts/`, `docs/`, `ZENODO_DEPOSITS.json`, `requirements.txt`
and nearly all of `data/` belong to the abandoned solar/geomagnetic–avalanche project.
Four files under `data/processed/` are still read by live code
(`ao_daily_cpc.txt`, `ncep_stratosphere.parquet`, `ssw_canonical.csv`,
`cryosphere/snotel_daily.parquet`); the rest is inert. `archive/` holds the retired
scripts, figures and documents.
