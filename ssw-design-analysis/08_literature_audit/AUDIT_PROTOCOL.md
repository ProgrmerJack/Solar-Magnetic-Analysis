# Literature audit protocol — converting Result D into a paper with consequences

**Written 2026-08-01, BEFORE extracting any published number.** This ordering is
the whole point: the project's failure mode is choosing a target and then finding
a result to fit it. Nothing below is conditioned on what the audit will show.

## Why this exists

Result D currently proves a bias exists (outcome-based SSW classification
reproduces 30–44% of the reported surface anomaly under a true null) and derives
its correction (bias = beta x bias_selection, R2=0.999, generalises to ENSO at
R2=0.985). It has **no body count**. No published number has moved.

The audit's only job: recompute published SSW surface-composite effect sizes
under the correction and report **how many change conclusion**. If the field's
central estimate drops materially, D is a paper about a live literature. If it
does not, D is a methods note and goes to WCD as one. Both outcomes are
publishable; only one is a major paper. **This document does not assume which.**

## Prespecified inclusion rule

A study enters the audit if ALL hold:
1. reports a QUANTITATIVE surface effect size (NAO/AO/NAM index, SLP, T2m) for SSWs;
2. the exposure is an SSW event set, not a continuous stratospheric index;
3. enough is reported to reconstruct the estimator: window, baseline, and whether
   events were classified using post-onset surface information;
4. published 2015 or later, OR earlier if cited >100 times (the propagation path
   matters more than recency).

Excluded, and why, recorded per study — never dropped silently.

## The two strata, fixed in advance

- **S1 OUTCOME-CLASSIFIED**: events selected or split using post-onset surface
  information (e.g. Karpechko-type "downward propagating" criteria). D predicts
  these carry bias = beta x bias_selection.
- **S2 NOT OUTCOME-CLASSIFIED**: all events, or events split on stratospheric or
  pre-onset predictors only. D predicts **no** bias here.

S2 is the control arm. If corrected and published values differ in S2, the
correction is wrong, not the literature. **This is the falsification test and it
is stated before any extraction.**

## Recomputation

For each study, recompute its reported quantity on the frozen 43-event catalogue
using its own stated window and baseline, then apply the additive correction with
beta estimated on event-free days. Report published, reproduced, corrected, and
whether the study's stated conclusion survives.

"Conclusion change" is prespecified as ANY of:
- a significant effect becoming non-significant at the study's own threshold;
- effect size changing by >33%;
- a sign change.

## Declared in advance

- Framed as a general property of outcome-based classification, illustrated by a
  widely used scheme. Never "X et al. are wrong."
- Karpechko et al. (2017) will be read in full to establish what was actually
  CLAIMED before anything is built on it. If the classification was used
  descriptively and the composite never presented as an unbiased effect size, D's
  critique softens and that will be stated in the abstract, not buried.
- Expected yield is unknown. If <3 studies change conclusion, the honest headline
  is "the bias is real, quantified and correctable, but the published literature
  is largely robust to it" — which is a WCD paper and will be submitted as one.

## Candidate pool (identified, not yet screened)

| study | why it is a candidate | stratum (provisional) |
|---|---|---|
| Hall et al. 2022, JGR-Atmos, CMIP6 SSW surface impacts | multi-model composites with the exact estimator; 11 of these models already processed here through an independent pipeline | S2 likely |
| Karpechko et al. 2017 | origin of the outcome-based classification | S1 |
| WCD 4, 213, 2023 (large ensembles, SSW precursors and NAO) | quotable frequencies (65% negative-NAO) and magnitudes | S2 likely |
| IOP ERL 2023, NW-Europe surface hazards after SSW | applied impact estimates on surface-coupled variables | S1 candidate |
| ACP 24, 1389, 2024 (air quality impacts of downward-propagating SSWs) | uses dSSW status as exposure for a surface-coupled outcome | **S2 — reassigned after screening** |
| ACP 26, 3723, 2026 (sorting SSWs by downward influence) | recent application of the classification | S1 |

Screening against the inclusion rule comes next; this list is a starting pool,
not a result.

## Independent re-implementation

Five substantive bugs have been found in this project (non-reproducible hash
seeds, a mis-specified FDR family, a sample-size confound nearly reported as a
finding, CMIP6 _FillValue contamination, a noleap calendar collapsing 165 years
into one winter). Five found implies more remain. The D pipeline will be
re-implemented from this specification by an independent route and the numbers
checked to agree before submission. With no co-author, the re-implementation is
the substitute for one.
