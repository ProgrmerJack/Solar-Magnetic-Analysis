# Gate 3, station level — the aggregated result does NOT transfer

Script: `calibrate_station_level.py` · Output: `station_calibration.json`
Metadata acquired by `03_data_ingestion/acquire_snotel_stations.py`
(887/945 stations resolved from the NRCS AWDB REST API, checksummed).

Sample: **103 stations stratified across elevation deciles, 510,269 station-days,
46 winters, 28 events**; elevation 9–3,560 m, latitude 36.1–68.1°.
Real-event influence removed first (189,350 rows, 37%), leaving 320,919 clean rows.

## The heterogeneity is real, and it is not elevation-ordered

| quantity | value |
|---|---|
| per-station seasonal amplitude | **5.94 – 39.33 °C** (median 11.77) — a 6.6× range |
| corr(elevation, seasonal amplitude) | −0.28 |
| corr(elevation, mean level) | −0.003 |

Stations differ enormously in the *amplitude* of their seasonal cycle, and elevation explains
only a modest part of it (r = −0.28). Mean level is uncorrelated with elevation, which is
unsurprising given that station × winter fixed effects absorb level anyway — the seasonal
*shape* is what the seasonal term has to carry.

## Result: at station level, shared seasonality costs a factor of four

| seasonal specification | null mean | in sd units |
|---|---|---|
| shared 3 harmonics | +0.2564 °C | +0.037 |
| per-elevation-band 3 harmonics | +0.2558 °C | +0.037 |
| **per-station LOWO climatology** | **+0.0621 °C** | **+0.009** |

Per-station climatology reduces the residual bias by **roughly four-fold** (0.256 → 0.062 °C).
Elevation banding achieves nothing at all (+0.2558 vs +0.2564) — consistent with the
diagnostic above: the heterogeneity is station-specific, not elevation-ordered, so grouping
stations by height does not capture it.

## This corrects the earlier reading

The state-aggregated test found shared and unit-specific seasonality indistinguishable, and
I reported that as contradicting the plan's assumption that a shared curve is inadmissible
for spatial data. **That reading was right for aggregated outcomes and wrong as a general
statement.** The two results together give the actual rule:

* **spatially aggregated outcomes** — shared seasonality is adequate; averaging within the
  unit has already removed the heterogeneity;
* **station-level outcomes** — shared seasonality leaves ~4× the bias, and only
  per-station climatology removes it. The plan's concern was correct here.

The residual +0.009 sd under per-station climatology is the smallest bias measured anywhere
in this project, smaller than the index case (+0.033 sd).

## Standing rule added

Any station-resolved or fine-gridded outcome must use **per-unit leave-one-winter-out
day-of-year climatology**, not a shared curve and not a banded curve. Coverage and
false-positive rate for that specification are being computed at 1,000 replicates; Gate 3
for station-level outcomes stays open until they report.

---

## Coverage: passes

Per-station LOWO climatology, 40 draws × 1,000 winter-block bootstrap replicates:

| metric | value | 95% CI | criterion | verdict |
|---|---|---|---|---|
| null mean | +0.009 sd | — | \|mean\| ≤ 0.05 | pass |
| coverage | **0.950** | (0.831, 0.994) | contains 0.95 | pass |
| false-positive rate | **0.050** | — | ≈ 0.05 | pass |

The point estimates are exactly nominal (38/40 covered, 2/40 rejected). The interval is
wide because n = 40 draws; each draw costs 1,000 fits over 320,919 rows, and more was not
affordable. The result should be read as *consistent with nominal calibration*, not as a
precise demonstration of it.

Second limitation: **103 stations, not 887.** The sample is stratified across elevation
deciles to preserve the heterogeneity under test, and the full seasonal-amplitude range
(5.94–39.33 °C) is present in it, but it is a sample.
