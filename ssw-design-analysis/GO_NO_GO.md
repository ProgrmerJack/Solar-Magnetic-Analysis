# Go/no-go conditions for Nature Geoscience

**Verdict, 2026-07-30: NO GO for Nature Geoscience.** Not because conditions remain
unmet, but because the thesis those conditions were built to support has been tested and
found false. Design choices do not reshape SSW impact estimates — estimator choice moves
the AO estimate 0.002 σ and catalogue choice 0.115 σ on a common analysis window, against
a response of ≈1.0 σ. See `FINDING_design_robustness.md`.

The work is sound and publishable as a **robustness result**. Target *WCD*, *GRL*,
*JGR-Atmospheres* or *J. Climate*.

| # | condition | status |
|---|---|---|
| 1 | Literature selection preregistered and systematic | **STRUCK** — requires a second coder (unachievable solo; LLM coding does not substitute, see `REMOVING_THE_HUMAN_DEPENDENCY.md`). |
| 2 | ≥10 representative published results re-analysed | **STRUCK** — was to be narrowed to studies with public code+data. Moot: the field-wide claim it supported is withdrawn. |
| 3 | **Pseudo-onset calibration centred near zero** | **CLOSED — all three regimes.** index: 3 harmonics, +0.033 σ, cov 0.947 / FPR 0.053. aggregated: +0.026 σ, cov 0.983. station-level: per-station LOWO climatology, +0.009 σ, cov 0.950. Re-confirmed per-bin 2026-07-30: null mean −0.017 to +0.075 across all 8 AO bins. See `GATE3_FINAL.md`. |
| 4 | Conclusions survive alternative seasonal models | **MET** — requirement is non-monotone in flexibility: 3 harmonics best for indices, per-station climatology required at station level. |
| 5 | Results survive primary and alternative event catalogues | **MET.** Catalogue re-frozen 2026-07-30 at **43 events / 36 winters** under an era-independent rule (≥2/3 of *covering* products; the old "≥4 of 6" demanded unanimity pre-1979 and made "marginal" a proxy for "pre-satellite", Fisher P=0.011). Post-onset response −0.94 to −1.14 across all 9 sets, every CI excluding zero; spread 0.115 σ on the common 1980–2019 window. |
| 6 | Result not driven by AO/NAO alone | **NOT MET.** Across 4 indices: AO and NAO post-onset bins survive FDR, **PNA 0/8** — the response is annular/Atlantic, not hemisphere-wide. Non-circulation outcomes were blocked behind the struck census. |
| 7 | Direction of changes reported honestly, incl. increases and sign reversals | **MET** — PNA increases are reported alongside the AO/NAO decreases, and the AAO negative-control investigation is documented in full including the false alarm and its resolution. |
| 8 | Physical decomposition explains meaningful heterogeneity | **NOT MET, and now understood.** R²=0.217, permutation **P=0.28**, 0 of 5 predictors significant. The reason is not weak predictors: mean design sensitivity is **−0.002 σ**, so there is essentially nothing to explain. The earlier "marginal events" route to this condition was **RETRACTED** (confounded with pre-satellite era; monotonicity failed; u10 mechanism rested on 3 events) — see `FINDING_marginal_RETRACTION.md`. |
| 9 | No headline depends on a single region or few events | **MET for the surviving headline.** The post-onset response holds on 43 events, on 29 isolated events, on 29 satellite-era events, and on the 21 that are both. |
| 10 | Independent atmospheric scientist + statistician audit | **NOT MET** — genuinely requires people; no substitute identified. |
| 11 | Every result from the frozen clean pipeline | **MET.** `load_catalogue()` is the sole event source; catalogue rebuilds byte-for-byte; `canonical_event_study.json` verified byte-identical across processes after the seed fix. |
| 12 | Conclusion interesting without "causal"/"inflated"/"upper bound" | **MET, but the conclusion inverted.** The defensible statement is that the response is *robust* to the design choices commonly invoked against it — not that design inflates it. |

## What changed on 2026-07-30

Three defects were found and fixed, in the order they surfaced:

1. **Era-dependent catalogue rule** — consensus as an absolute count of 6 silently demanded
   unanimity before 1979. Replaced with a fraction of covering products; 43 events.
2. **Non-reproducible bootstrap seeds** — `hash(<str>)` is randomised per process, so no
   p-value from `canonical_event_study.py`, `multi_index_event_study.py` or
   `marginal_events.py` was reproducible. Replaced with `zlib.crc32`; reproduction now
   verified byte-identical. The `consensus_strict` pre-onset bin had been flipping
   significance on 8 of 40 seeds.
3. **Negative-control FDR family** — the AAO was pooled into a 32-test family with the
   positive outcomes, whose small p-values raise the BH threshold and make the *control*
   easier to flag. Judged on its own 8 bins it passes (min q = 0.24).

Two false alarms were chased and correctly dismissed, both recorded so they are not
repeated: the AAO "negative control failure" (multiple testing, not the estimator — ENSO
was tested as a common cause and attenuates it 0%), and an apparent 20–27% bootstrap
under-coverage that was entirely a sample-size confound (`ao_null_perbin.py`, marked
VOID; the like-for-like test in `boot_calibration.py` returns SE ratio 0.96–1.00).

## The one live scientific question

Whether the pre-onset AO anomaly is real. It is −0.82 (p=0.015) on the full record and
loses significance under every restriction, but **every restricted interval still contains
−0.82**, and neither mechanism test is significant (era p=0.115, clustering p=0.234,
independent at Fisher P=0.49). 43 events cannot resolve it.

SNAPSI can: `nudged` (stratosphere nudged to observed) minus `control` (nudged to
climatology) is the causal stratospheric contribution *by experimental design*, with 50
members per centre per case instead of one realisation of history.

**Correction, 2026-07-31: SNAPSI was never blocked by a CEDA outage.** The earlier
"unreachable, 404" conclusion came from probing `/badc/snapsi`, which is simply the wrong
path. The real root is `/badc/snap/data/post-cmip6/SNAPSI` and it is fully browsable
**anonymously**. Only the file download requires a free CEDA account.

`acquire_snapsi.py` now enumerates the archive without credentials and has written
`snapsi_manifest.csv`: **21,486 files, 11 centres, 5 experiments** (control, control-full,
free, nudged, nudged-full), **6 cases** (s20180125, s20180208, s20181213, s20190108,
s20190829, s20191001), 50 members each. Pattern-constructed URLs are sample-verified;
21,386 of 21,486 verified, the 2 failing nodes (CNR-ISAC/GLOBO/nudged/s20190108) flagged
rather than trusted.

Download tiers for `nudged`+`control`, `psl` at 6hrPt — grid resolution varies enormously
across centres (NRL 2.0 MB/file to UKMO 67.2 MB/file):

| scope | files | size |
|---|---|---|
| CCCma pilot — validates the pipeline end to end | 385 | **2.1 GB** |
| 4 smallest-grid centres (CCCma, KMA, NRL, SNU) — a real multi-model answer | 2,547 | **16.8 GB** |
| all 11 centres | 6,121 | 147 GB |

Blocked only on a free account at https://services.ceda.ac.uk/ — then set
`CEDA_USERNAME` / `CEDA_PASSWORD` and re-run.
