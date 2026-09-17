# SSW surface-impact methodology

Does the literature's split of sudden stratospheric warmings into
"downward-propagating" and not carry information about the Earth system, or is it a
threshold applied to one continuous population?

Evidence comes from three systems: observations (43 events, 1958–2024), CMIP6
(1,888 events, 20 members), and SNAPSI — a nudged-ensemble experiment in which the
stratospheric evolution is identical across members by construction, so the
classification can be tested against a driver that is fixed by design rather than
against a modelled null.

## Where things are

| | |
|---|---|
| **Conventions, layout, how to add an analysis** | [`CLAUDE.md`](CLAUDE.md) |
| **Results — the single source of truth** | [`ssw-design-analysis/CONSOLIDATED_RESULTS.md`](ssw-design-analysis/CONSOLIDATED_RESULTS.md) |
| **Traps, and the approaches that failed** | [`FAILURES.md`](FAILURES.md) |
| **Every script and result, generated** | [`CATALOG.md`](CATALOG.md) |
| **Venue assessment and the work remaining** | [`ssw-design-analysis/NATURE_GEO_PLAN.md`](ssw-design-analysis/NATURE_GEO_PLAN.md) |

Individual `FINDING_*.md` documents go stale; `CONSOLIDATED_RESULTS.md` is the one
that is maintained. Where a claim was withdrawn, §7 records it rather than dropping
it — nine headline claims died during this project, and the tests that killed them
are the reason the survivors are trustworthy.

## Setup

```bash
python3 -m venv .venv
.venv/bin/pip install -r ssw-design-analysis/environment/requirements.lock.txt

.venv/bin/python ssw-design-analysis/environment/check_environment.py
.venv/bin/python ssw-design-analysis/tests/smoke_test.py
.venv/bin/python tools/catalog_repo.py
```

Run everything with the venv interpreter. A system Python with an xarray that has no
NetCDF backend fails on every `.nc` file here with a message about missing
dependencies rather than a missing file, which is hard to diagnose and was a real
blocker — `check_environment.py` exists to catch exactly that.

## Data

Events come from the NOAA CSL compendium, frozen and checksummed; reanalysis from
NCEP (THREDDS) and ERA5 (WeatherBench2 public GCS zarr, anonymous); models from CMIP6
zonal means; the experiment from SNAPSI on CEDA, which is browsable anonymously and
needs a free account only to download. Everything under `data/` is read-only.

## The repository name is a fossil

This began as a solar/geomagnetic study of avalanche activity. That project was
abandoned: three headline claims collapsed under stricter specifications and both of
its manuscripts are marked do-not-submit. `paper/`, `scripts/`, `docs/`, most of
`data/` and everything in `archive/` belong to it. Four files under `data/processed/`
are still read by live code; the rest is inert. Nothing solar is live.
