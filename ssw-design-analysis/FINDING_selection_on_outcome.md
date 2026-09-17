# FINDING: the design choice that actually matters is selection on the outcome

**Date:** 2026-07-31
**Status:** positive result, internally validated. This is the research direction
that survives after the design-sensitivity thesis failed
(`FINDING_design_robustness.md`).

---

## 1. The setup

Two design axes were tested on the frozen 43-event catalogue and returned null:

| axis | effect on the AO estimate |
|---|---|
| estimator (conventional composite vs within-winter) | **0.002 σ** |
| event catalogue (9 sets, common 1980–2019 window) | **0.115 σ** |

Against a response of ≈1.0 σ, neither matters. The design choices the field
debates are not where the problem is.

This tests a third axis, embedded in standard practice rather than debated:
**selecting events on the outcome.**

## 2. The practice, as published

Karpechko et al. (2017) classify an SSW as "downward-propagating" (dSSW) if

1. mean NAM at 1000/850 hPa over days **+8..+52** is negative
2. fraction of those days with negative surface NAM **> 0.5**
3. fraction with negative NAM at 150/100 hPa > 0.7

and the classification remains in active use (ACP **24**, 1389, 2024; ACP **26**,
3723, 2026).

Criteria 1 and 2 *are* the surface response.

**This is not a criticism of the classification.** It is sound for what it was
built for — precursors, predictability, and what distinguishes coupling from
non-coupling events. A case-control design may legitimately select on outcome
and look back at predictors. The problem is narrower and specific: when the
**selected group's mean surface anomaly is reported as an impact**, part of that
number is manufactured by the selection, and the size has not been quantified.

## 3. The measurement

Pseudo-onsets are drawn on clean data — every day within −60..+75 of a real
onset deleted first — so the true effect is **exactly zero** by construction.
The published classification is then applied to them. Whatever the selected
group shows is pure selection.

Two corrections were needed before the numbers meant anything, both recorded
because the first version of this analysis got them wrong:

- **Contaminated reference.** The day-of-year climatology was built from the full
  record, which contains depressed real-event days, so clean pseudo-onsets scored
  +0.279 σ against it instead of zero. Predicted offset from the clean-minus-full
  climatology difference: +0.2296. Rebuilding the climatology from clean days
  only brings the null to **+0.003**.
- **Unmatched selection rate.** The criterion selects 70% of real events but only
  38% of pseudo-events, and selecting a *smaller* fraction is a *stronger*
  selection. The uncorrected comparison overstated the bias (it gave 65–91%).
  The null selection is now rate-matched to the real classification rate.

## 4. Result

Selection bias under a true null, rate-matched, 1,500 draws:

| outcome | real dSSW reported | selection bias (true effect = 0) | reproduced by selection |
|---|---|---|---|
| AO — the selection variable itself | −1.433 | −0.529 | **37%** |
| NAO — different index, correlated | −0.464 | −0.140 | **30%** |
| PNA — different index, weakly correlated | +0.186 | +0.082 | **44%** |

Roughly a third to a half of the reported "downward-propagating SSW impact" is
reproducible from selection alone, with no real effect whatsoever — **and it
propagates to outcomes that are not the selection variable.**

## 5. Internal validation: the bias is additive and recoverable

If the bias is measured correctly it should subtract cleanly, returning the
all-events estimate:

| outcome | dSSW reported | − bias | = implied | all-events | residual |
|---|---|---|---|---|---|
| AO | −1.433 | −0.529 | −0.904 | −0.895 | **−0.010** |
| NAO | −0.464 | −0.140 | −0.325 | −0.300 | **−0.024** |
| PNA | +0.186 | +0.082 | +0.104 | +0.119 | **−0.015** |

Residuals of 0.010–0.024 σ across three independent outcomes. The decomposition
is not a rhetorical device; it closes numerically.

## 6. What IS uncontaminated: the classification rate

The selected group's *magnitude* is compromised. The *rate* is not:

| | real events | pseudo-events | |
|---|---|---|---|
| classified downward-propagating | **70%** | **38%** | P < 0.0001 |

Real SSWs are genuinely far more likely to satisfy the downward-propagation
criteria than chance allows. That is a real, strong, uncontaminated signal — and
it is the quantity the field should be reporting.

## 6b. The bias obeys a projection law — sign and size are predictable

Extending the test to the 945-station SNOTEL network (`selection_snotel.py`)
produced a surprise: for western-US temperature and snow depth the selection
bias has the **opposite sign** — it attenuates rather than inflates. So
"selecting on the outcome inflates impacts" is false as a general claim.

The correct statement is a projection law, and it is derivable rather than
empirical. For jointly normal window means (Y_w, A_w), selecting on A_w gives

    E[Y_w | selected] − E[Y_w] = ρ · (σ_Y/σ_A) · (E[A_w | selected] − E[A_w])

The bracket is exactly the bias measured for the selection variable itself, and
ρ·σ_Y/σ_A is the **ordinary regression slope β of Y_w on A_w**. So

> **bias(outcome) = β(outcome window on selection-variable window) × bias(selection variable)**

Measured against 34,800 clean pseudo-event windows, with β estimated on
event-free data so the real events never inform the prediction of their own bias:

| outcome | ρ | β | predicted bias | measured bias | error |
|---|---|---|---|---|---|
| AO (anchor, β≡1) | +1.000 | +1.000 | −0.660 | −0.660 | −0.000 |
| NAO | +0.625 | +0.200 | −0.132 | −0.133 | **−0.001** |
| PNA | −0.460 | −0.179 | +0.118 | +0.126 | +0.008 |
| SNOTEL mean temperature | −0.109 | −0.035 | +0.023 | +0.009 | −0.014 |
| SNOTEL snow depth | +0.075 | +0.064 | −0.042 | −0.056 | −0.014 |

**R² = 0.999**, fitted slope −0.657 against the derived −0.660, intercept −0.005.

Two specification errors were made before reaching this, both recorded because
each looked like evidence against the law:

1. **Daily instead of window correlation.** The selection acts on the +8..+52
   window mean. Daily correlation gave R²=0.868, NAO off by 66%.
2. **Correlation instead of regression slope** — the σ_Y/σ_A factor was dropped.
   Window means of persistent variables (AO) spread much more than those of less
   persistent ones (NAO), so ρ alone over-predicts every non-selection outcome by
   roughly a factor of three. R²=0.823, NAO predicted −0.412 against −0.133.

Both were mine. The law was never in doubt; the predictor was mis-specified.

Honest limits: five outcomes, of which AO is definitional (β≡1 by construction),
leaving four free points. Absolute accuracy is ±0.014 σ across β from −0.18 to
+1.0. For outcomes whose bias is near zero the *relative* error is larger
(SNOTEL temperature: 0.023 predicted vs 0.009 measured), which is what near-zero
quantities always do.

## 7. The recommendation

- **Do not report the dSSW-selected group's mean surface anomaly as an impact.**
  A third to a half of it is selection for AO/NAO/PNA, and the contamination
  carries into any outcome that co-varies with the selection variable — snowpack,
  temperature, air quality, energy demand.
- **Never assume the bias is conservative.** Its sign follows β, not intuition:
  for western-US snowpack it runs *opposite* to the effect, so "selection can
  only inflate, so our result is a lower bound" is unsafe.
- **Do report the classification rate against a matched null.** 70% vs 38% at
  P < 0.0001 is a strong result that costs nothing in rigour.
- **If the selected-group magnitude is wanted, correct it.** The bias is additive
  (§5, residuals 0.010–0.024 σ) and predictable from one regression slope on
  event-free data (§6b, R²=0.999) — no per-outcome null rerun required.

## 8. Why this is the paper

The honest arc: the two design choices the field argues about (estimator,
catalogue) do not move SSW impact estimates. The one nobody flags — selection on
the outcome, built into a standard and widely used classification — accounts for
30–44% of a headline number. The bias is quantified, validated, and correctable,
and an uncontaminated alternative is available and stronger.

Extension in progress (`selection_snotel.py`): the same test on the 945-station
SNOTEL network, so the demonstration lands on a genuine applied impact variable
rather than a circulation index.
