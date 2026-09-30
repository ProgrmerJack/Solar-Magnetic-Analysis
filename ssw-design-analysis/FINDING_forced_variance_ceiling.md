# Finding: SSWs do not differ from one another in their forced surface response

**Status: established in a 1888-event CMIP6 ensemble with a placebo control that
removes 91% of the naive signal. Observations (n=43) are consistent but have
almost no power and do NOT establish it.**

Script: `07_physical_decomposition/forced_variance_ceiling.py`
Data: `forced_variance_ceiling.json`

## The question the field skipped

The literature asks *which* SSWs propagate downward and builds classifiers to
find out. That presupposes an answer to a prior question nobody asked: **do SSWs
differ from one another in their forced surface response at all?**

If every SSW forces the same expected surface anomaly, and the spread between
events is internal tropospheric noise, then there is nothing to classify. Every
classifier — Karpechko's, a stratospheric one, a machine-learned one, one not yet
invented — is sorting noise.

## The decomposition

Write the per-event surface anomaly as

```
Y_i = mu + S + f_i + eps_i     at real SSWs
Y_j = mu     +       eps_j     at event-free pseudo-onsets
```

`S` is the mean forced shift, common to every event — the real, large SSW effect,
which is not in dispute. `f_i` is the **event-specific** part of the forced
response, with `Var(f) = s_f^2`. `eps` is internal noise. With `f` independent of
`eps`,

```
s_f^2 = Var(Y | SSW) - Var(Y | pseudo)
```

`s_f` is the only quantity a classifier can recover. Everything else is noise no
predictor can see, in principle.

## The result, and the control that changed it

| | CMIP6 (n=1888) | observations (n=43) |
|---|---|---|
| mean shift S | −0.509 sigma | −0.904 sigma |
| Var(Y \| SSW) | 0.5189 | 1.1167 |
| Var(Y \| pseudo) | 0.4480 | 1.1788 |
| naive s_f^2 | **+0.0709** [+0.0354, +0.1074] | −0.0620 [−0.4332, +0.3136] |
| **placebo s_f^2** (pre-onset, true value 0) | **+0.0642** [+0.0222, +0.1075] | +0.2832 [−0.4412, +1.1262] |
| **corrected s_f^2** | **+0.0066** [−0.0460, +0.0570] | −0.3453 [−1.1665, +0.4057] |
| corrected s_f | 0.081 (upper 95% **0.239**) | 0 (upper 95% 0.637) |

**The placebo is the whole story in CMIP6.** The naive estimate looks decisively
positive. But running the identical decomposition on the pre-onset window
(days −52..−8), where `s_f = 0` by construction because the event has not happened
yet, returns +0.0642 — **91% of the apparent signal is already there before onset.**

The cause is that SSWs are not randomly timed. They cluster in disturbed winters
whose tropospheric variance is already elevated (pre-onset variance ratio 1.13),
and pseudo-onsets matched only on day-of-year do not reproduce that. Reporting the
naive number as "event-specific forced response" would have been a straight false
positive, and this project has now caught three of those with controls.

After the paired correction, **`s_f^2` is indistinguishable from zero.**

## The ceiling on every classifier

For a perfect classifier — one that knows `f_i` exactly — splitting at the
q-quantile of `f` gives, for Gaussian `f`,

```
contrast(q) = s_f * phi(z_q) * (1/q + 1/(1-q))
```

**A first draft of this used 2*sqrt(2/pi) = 1.596 and called it the maximum. That
was wrong, and wrong in the direction that flattered the argument.** The
expression is *minimised* at q = 0.5 and grows without bound for more extreme
splits (1.596 at q=0.5, 1.655 at q=0.7, 1.950 at q=0.9). A classifier free to pick
its split point has no finite ceiling at all. So the ceiling is evaluated at **the
split fraction each study itself uses**, which for this literature is 54–70% DW,
putting the factor in 1.60–1.66.

Using the **uncorrected** `s_f` upper bound — the conservative arm, since it still
contains the confound:

| reported contrast | q | \|value\| | ceiling (95%) | ratio | matched |
|---|---|---|---|---|---|
| CMIP6 DW−NDW (this project) | 0.50 | 1.143 | 0.523 | **2.2x** | yes |
| Karpechko-criterion AO, obs | 0.70 | 1.782 | 0.927 | **1.9x** | yes |
| ACP 26,3723 (2026) published NAO | 0.59 | 0.850 | 0.900 | 0.9x | yes |
| ERA5 1000 hPa NAM, Karpechko | 0.54 | 0.684 | 0.895 | 0.8x | yes |
| ERA5 850 hPa NAM, ACP criterion | 0.59 | 0.639 | 0.900 | 0.7x | yes |

With the placebo-corrected CMIP6 bound the ceiling falls to **0.381 sigma**, and
the CMIP6 contrast of 1.143 is **3.0x** above it.

## What this does and does not establish

**Establishes, in models:** in a 1888-event ensemble the reported DW−NDW contrast
is 2.2–3.0x larger than anything a perfect classifier could produce. Most of the
contrast therefore cannot be classification; it must be selection.

**Does NOT establish, in observations:** at n=43 the ceiling is 0.89–0.93 sigma,
and three of the four published observational contrasts (0.64–0.85) fall *below*
it. **This test does not exclude them.** Only the AO contrast (1.782) exceeds the
observational ceiling. The observational arm is reported as consistency, not
evidence — exactly as in `FINDING_no_class.md`.

The known-truth calibration (C2) shows why: at n=43 the estimator returns
0.197 ± 0.253 when the true `s_f` is 0. It is upward-biased and noisy at that
sample size, so only the upper bound is meaningful.

## Independent convergence with the bias law

`FINDING_stratifier_law.md` found by a completely different route — the OLS
projection `beta x dS` — that the CMIP6 AM_post contrast is **96% selection**,
leaving a residual of 0.046. This decomposition says a perfect classifier could
produce at most 0.523. The residual sits far below the ceiling, as it must.

**Two methods with different assumptions — regression projection versus variance
decomposition — give the same answer.** That is the strongest internal check the
project has.

## The definitive test, which needs one thing

`s_f` is measured here indirectly, by differencing variances against pseudo-events
and correcting with a placebo. **An initialised ensemble measures it directly:**
with many members per event, the ensemble-mean response *is* `mu + S + f_i`, and
`Var(f)` is read straight off the between-event spread of ensemble means, with no
pseudo-events and no placebo correction needed.

SNAPSI is exactly that dataset — 6 SSW cases, 11 centres, 50 members each, and the
archive path is already mapped (`03_data_ingestion/acquire_snapsi.py`,
21,486-row manifest). Browsing needs no credentials; **downloading needs a free
CEDA account** from https://services.ceda.ac.uk/. That is the only blocker.

This also settles the tension with Loeffel et al. (2025), who report from ICON
ensemble re-forecasts that events *do* differ in their forced response
("individual SSW events differ significantly in their likelihood to induce a
canonical tropospheric response"). Their design can see `f_i` directly and mine
cannot; if `s_f` is genuinely small but non-zero, both results are correct and the
reconciliation is that **the forced differences are real but far too small to be
recovered from a single realisation** — which is all the observational literature
ever has.

## Limits

- Established in MODELS. Ten CMIP6 models, with known spread in how well they
  reproduce the SSW surface composite.
- Assumes `f` independent of `eps`. If a large forced response also suppressed
  internal variance, the decomposition misattributes. Untested.
- The placebo correction may be over-conservative: if preconditioning genuinely
  makes some events force a larger response, part of the pre-onset excess is
  signal, and subtracting it removes real effect. That is why the ceiling is
  quoted from the *uncorrected* bound.
- Gaussian `f` assumed in the split-factor algebra.
