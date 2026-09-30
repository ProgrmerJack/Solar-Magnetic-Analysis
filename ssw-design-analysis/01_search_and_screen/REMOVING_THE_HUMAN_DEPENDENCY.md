# Can the second-coder dependency be removed? Research as of 2026 Q2

## The dependency is attached to one claim, not to the paper

A second coder is required because inter-rater reliability validates **subjective
classification of what published papers did**. That is needed only for a claim *about the
literature* — "N% of studies use design X". It is not needed for a claim about **estimators**.

So the dependency can be removed, but only by changing what is claimed. It cannot be
removed by automating the coding, for reasons below.

---

## Option A — LLM-assisted coding: a reduction, not an elimination

State of the art as of 2026 Q2:

* **PRISMA-trAIce (2025)**, a checklist extension for AI-assisted evidence synthesis,
  requires reporting the tool, version, prompts **and how human verification was performed**
  ([JMIR AI 2025](https://ai.jmir.org/2025/1/e80247)). Human verification is a reporting
  requirement, not an optional step.
* Dual-model LLM ensembles that mimic the two-reviewer process (GPT-4-turbo + Claude-3-Opus
  with cross-critique) are published and reach ~95% extraction accuracy with Cohen's κ
  reported ([Frontiers Digit. Health 2026](https://www.frontiersin.org/journals/digital-health/articles/10.3389/fdgth.2026.1799623/full);
  [PubMed 39836495](https://pubmed.ncbi.nlm.nih.gov/39836495/)).
* But the same literature reports that LLM reliability **"degrades on tasks requiring
  interpretive judgment"**.

Classifying a paper's baseline construction, its treatment of pre-onset periods, or whether
its language is causal or associative **is** interpretive judgment — precisely the regime
where the method is weakest. LLM assistance would reduce the human workload but would still
require documented human verification, and would be weakest on the fields that matter most.

**Verdict: does not remove the dependency.**

---

## Option B — SNAPSI: experimental ground truth, no coding at all

**This is the strongest option and it removes the dependency entirely.**

SNAPSI (Stratospheric Nudging And Predictable Surface Impacts,
[GMD 15, 5073, 2022](https://gmd.copernicus.org/articles/15/5073/2022/)) is a community
protocol in which **11 operational centres** ran ~50-member sub-seasonal ensembles
initialised near three events (NH 2018, NH 2019, SH 2020), in three configurations:

* free-running;
* **stratosphere nudged to observed conditions**;
* **stratosphere nudged to climatology**.

The difference between the last two **is the causal stratospheric contribution, by
experimental construction.** That is a ground truth obtained from a designed experiment, not
from a simulation whose assumptions I would be choosing myself.

This makes the central question directly answerable without touching the literature:

> Apply each candidate estimator — conventional composite, single-window, within-winter
> event study — to the free-running ensembles, and compare each estimate against the
> nudged-minus-climatology benchmark. The discrepancy is the estimator's bias, measured
> against a physical standard.

It also satisfies the plan's own item-6 requirement for "a large model ensemble or
nudging/scrambling experiment where the incremental stratospheric contribution can be
independently estimated".

**Access:** CEDA catalogue, record
[0a5a1ce22fb047749e040879efa8e9b5](https://catalogue.ceda.ac.uk/uuid/0a5a1ce22fb047749e040879efa8e9b5/).
Direct BADC paths 404 without authentication; CEDA requires a **free account**. That is a
registration step, not a collaborator — a different order of dependency entirely.

**Limitation, stated plainly:** three events, two of them the same hemisphere. It gives
strong evidence about estimator behaviour, not about how often the field is affected.

---

## Option C — reproduction instead of classification

Code-and-data sharing in the relevant literature rose **from 11% in 2014 to 64% in 2024**
(sixfold), while papers sharing neither fell from 29% to 4%
([EGUsphere 2025](https://egusphere.copernicus.org/preprints/2025/egusphere-2025-5210/egusphere-2025-5210.pdf)).

Restricting the replication sample to studies with public code and data replaces subjective
coding with an objective operation: the code runs and yields a number, or it does not. No
inter-rater reliability is required because no rater judgment is involved.

**Cost:** the sample is no longer representative — it is the subset that shares code, which
is younger and probably better-practised than the field average. That biases *against*
finding problems, which is at least the conservative direction, and must be stated.

---

## What is lost, and what is not

**Lost:** any claim about prevalence. Without a survey there is no defensible way to say
"most SSW-impact studies use this design". Go/no-go conditions 1 and 2, as written, become
unachievable by a solo author and should be struck rather than fudged.

**Not lost:** the methodological result itself. "These estimators differ by X against
experimental ground truth, and the difference has this dynamical structure" is a complete,
falsifiable, publishable claim that needs no census.

## Novelty check

The closest recent work, [ACP 26, 3723, 2026](https://acp.copernicus.org/articles/26/3723/2026/),
classifies SSWs by their downstream tropospheric impact (Eurasia / North America / both /
non-downward). Verified by direct reading: it applies **one** classification framework
consistently and does **not** compare compositing or event-selection methods, nor quantify
how method choice changes the estimated response. **No conflict.**

---

## Recommendation

1. **Acquire SNAPSI** (free CEDA registration). It converts the estimator question from
   "what do papers do?" to "what is each estimator's bias against a designed experiment?"
2. **Drop the census** and strike go/no-go conditions 1 and 2, rather than satisfying them
   weakly with single-coder or LLM-only coding.
3. **Keep a small, citable exemplar set** of designs — sufficient to show the estimators are
   in real use, without claiming prevalence.
4. **Retarget** accordingly: without a field-wide survey this is not a Nature Geoscience
   *Analysis* about the literature. It is a strong methods paper for *Weather and Climate
   Dynamics*, *GRL*, *JGR-Atmospheres* or *QJRMS* — which is where the internal review
   already judged the honest deliverable to sit.

The human dependency is removable. The Nature Geoscience framing largely is not, because
that framing rested on the field-wide claim the census was there to support.
