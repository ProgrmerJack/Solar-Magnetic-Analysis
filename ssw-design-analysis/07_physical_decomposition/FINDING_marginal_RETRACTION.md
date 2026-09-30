# Retraction and correction: the marginal-event finding does not hold up

The three strengthening tests were run. **All three weakened it, and the third
found a bug in my own catalogue rule.** This supersedes `FINDING_marginal_events.md`,
which must not be used.

## 1. Monotonicity — not supported

If marginal detections were misclassified events, the pre-onset anomaly should weaken
steadily as more reanalyses agree. Per-event anomalies against each event's own winter
baseline, all 45 union events, Spearman with 10,000-permutation test:

| index | Spearman(detection count, pre-onset) | perm P | post-onset | perm P |
|---|---|---|---|---|
| AO | +0.221 | 0.147 | −0.055 | 0.72 |
| NAO | +0.240 | 0.115 | −0.134 | 0.39 |

The sign is consistent across both indices and in the predicted direction, but **neither is
significant**, and the pattern is not monotone — events detected by 3/6 have a *larger*
pre-onset anomaly (−1.71) than those detected by 1/6 (−1.07).

## 2. NAO replication — same direction, same non-significance

NAO reproduces AO's shape (3/6 most negative pre-onset) and a nearly identical rank
correlation. That is mildly reassuring about consistency and does nothing for significance.

## 3. The vortex-wind mechanism — driven by three events, and not monotone

Mean 10 hPa wind at onset by detection count:

| detected by | 1/6 | 3/6 | 4/6 | 5/6 | 6/6 |
|---|---|---|---|---|---|
| u10 (m/s) | **+14.79** | +1.39 | −1.66 | −1.03 | +0.01 |

Spearman ρ = −0.017, **P = 0.93**. The striking "+8.09 m/s, vortex never reversed" contrast
I reported comes **entirely from the three events detected by a single reanalysis**. Events
detected by 3/6 already sit at +1.39, close to the strong events. The mechanism was inferred
from a group mean that a subgroup of three dominated.

## 4. The bug: the consensus rule is era-dependent

Not all six reanalyses cover all dates:

| era | reanalyses available | what ">=4 of 6" actually demands |
|---|---|---|
| 1958–1978 | 4 (NCEP, ERA40, JRA-55, ERA5) | **unanimity** |
| 1979–2002 | 6 | 67% |
| 2003–2019 | 5 | 80% |
| 2020–2024 | 4 | **unanimity** |

Mean achieved consensus is 3.3 for 1958–1978 against 5.4 for 1979–2002. The primary
catalogue rule is therefore **stricter in the early record than in the satellite era**, purely
because fewer products exist.

**Consequence:** "marginal" is substantially a proxy for "early".

| | not marginal | marginal |
|---|---|---|
| 1979 or later | 28 | 2 |
| pre-1979 | 9 | 6 |

**Fisher exact P = 0.011**; **6 of the 8 marginal events are pre-1979.**

The pre-onset anomaly attributed to marginal detections is therefore confounded with
pre-satellite-era data quality — in both the event detection *and* the AO index itself. The
two explanations cannot be separated with these data, and the classification explanation was
the one I offered.

## What must be fixed

1. **Redefine consensus as a fraction of reanalyses covering that date**, not an absolute
   count out of six. `CATALOGUE_SPECIFICATION.md` §2 and `build_catalogue.py` both need
   changing, and the catalogue must be rebuilt and re-frozen.
2. **Re-run every catalogue-dependent result** afterwards — the canonical event study, the
   multi-index study, and the design-sensitivity comparison all use the current primary set.
3. **Re-test the marginal-event contrast with era controlled**, e.g. restricted to 1979+
   where coverage is uniform. With only 2 marginal events after 1979 that test will be
   badly underpowered, which is itself the honest answer.

## Standing

The double dissociation in `FINDING_marginal_events.md` is real *in the data as split*, and
survived FDR. But the split is confounded with era, the proposed mechanism rests on three
events, and monotonicity fails. **It cannot support a publication claim**, and the
Nature Geoscience case I sketched on the back of it is withdrawn.
