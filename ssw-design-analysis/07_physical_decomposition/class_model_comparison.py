#!/usr/bin/env python3
"""
class_model_comparison.py
=========================
ONE SHIFTED POPULATION, A CONTINUOUSLY STATE-DEPENDENT ONE, OR TWO REGIMES?

Plan approved 2026-09-30. Written before any of these models was fitted.

QUESTION
  Which statistical description of the surface response after an SSW predicts
  held-out events best: one shifted distribution, a distribution whose centre and
  spread vary continuously with the event's state, or two regimes ("downward" and
  "not") whose membership depends on the event's state? The class reading of the
  literature is the third; our reading is the first (with small continuous
  departures). Proper scores on held-out data decide, on the class model's own
  terms, rather than a threshold applied after the outcome.

A. CMIP6 (1,517 Charlton-Polvani events, 20 members, 10 models; the data, outcome
   and predictors of forecast_value_test.py: annular mode days +8..+52 in member
   sigma units, not demeaned; predictors of tier P2, information up to onset,
   standardised within member, plus day-of-year sine and cosine)
   M0 SHIFT       Gaussian; mean = calendar regression; spread = residual s.d.
                  (forecast_value_test F0).
   M1 CONTINUOUS  Gaussian; mean linear in predictors + calendar (ridge, as F1);
                  log spread linear in the same predictors (fitted by maximum
                  likelihood on the ridge residuals, L2 penalty chosen by inner
                  grouped cross-validation).
   M2 REGIMES     two Gaussian components with their own constant means (plus the
                  calendar terms, shared) and variances; membership probability
                  logistic in the predictors (EM, 20 random starts, best likelihood).
   M2b TWO POPULATIONS, state-blind: M2 with constant membership probability.
   Cross-validation: leave one MODEL out (as forecast_value_test). Scores per
   event: CRPS (closed form for Gaussians and Gaussian mixtures) and the ignorance
   (negative log predictive density). Skill of each model against M0; the class
   model's gain over the continuous one, G = score(M1) - score(M2), with a
   model-cluster bootstrap interval (2,000). The identical pipeline on 200
   size-matched sets of event-free pseudo-onsets (forecast_value_test draws).
   Primary readout: G_SSW - mean(G_pseudo); p = P(G_pseudo >= G_SSW), CRPS
   primary, ignorance secondary; and M2b against M0 at SSWs relative to pseudo.
   Positive control: synthetic outcomes with a planted two-regime structure
   (two thirds / one third, regime means 1.0 sigma apart, membership logistic in
   the first predictor with slope 1) with the real predictors; detection rate of
   G > 0 at p < 0.05 over 100 datasets. Negative control: the same with no regimes.

B. SNAPSI (36 Northern Hemisphere ensembles per arm, nine models; members' days
   +8..+25 polar-cap NAM proxy standardised by the control, as in the
   distribution test). In each ensemble the forcing is shared, so two kinds of
   response would appear as a mixture within the nudged ensemble. Five-fold
   cross-validated mean log predictive density of a two-component Gaussian mixture
   minus that of one Gaussian, per ensemble; nudged minus control on paired
   ensembles, with a 10,000-resample centre bootstrap. Positive control: control
   members with two thirds of them displaced by 1.0 sigma; detection rate over
   100 resamples.

Reading, fixed now: the class reading is supported if M2 beats M1 after SSWs by
more than on ordinary days (p < 0.05) or if the SNAPSI mixture gain is larger with
the SSW imposed (interval excluding zero). If neither, and the positive controls
detect the planted structure, the two-regime description adds nothing detectable.

Output: results/current/6_predictability/class_model_comparison.json
"""
raise SystemExit("design registered 2026-09-30; implementation follows")
