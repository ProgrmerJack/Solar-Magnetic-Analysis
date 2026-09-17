export const meta = {
  name: 'ssw-audit-round3',
  description: 'Remaining literature audit: forecast circularity, prior art closer than 2019, S1 body count, stratospheric-split usage',
  phases: [
    { title: 'Audit' },
    { title: 'Verify' },
  ],
}

const RULES = `
CONTEXT. Single-author observational SSW study. A methods result is being written up:
classifying SSWs into "downward-propagating" (DW) vs not, using the Karpechko et al.
(2017) criterion (VERIFIED from source: over days +8..+52 after onset, (1) mean NAM at
1000 hPa negative, (2) fraction of those 45 days with negative 1000 hPa NAM > 0.5,
(3) fraction with negative 150 hPa NAM > 0.7), produces a DW-minus-NDW surface contrast
that is 93-97% reproducible by applying the identical criterion to random winter dates
with no SSW at all.

PRIOR ART ALREADY FOUND AND VERIFIED FROM SOURCE. Do NOT re-report these. Find things
CLOSER, NEWER, or in adjacent literatures:
1. White, Garfinkel, Gerber, Jucker, Hitchcock et al., J. Climate 32, 85 (2019), section
   3b: differences "qualitatively similar to that found in the DW - NDW differences (but
   with differing magnitudes)"; "the differences at positive lags in the troposphere are
   entirely there by construction." Not quantified, no correction, not extended to
   stratospheric classification.
2. Coughlin & Gray (2009), JAS 66, 531, "A Continuum of Sudden Stratospheric Warmings" --
   continuum on the SSW-DEFINITION axis (major vs minor).
3. Maury, Claud, Manzini, Hauchecorne, Keckhut (2016), JGR-Atmos, 10.1002/2015JD024226 --
   "the idea of a 'warming continuum'"; "there is no statistical difference between SWEs
   with regard to their feedbacks on planetary waves and hence their potential influence
   into the troposphere."
4. Baldwin et al. (2021), Rev. Geophys. 59, e2020RG000708, section 7.2: "about two thirds
   ... of SSW events are characterized as having a visible downward impact", no
   uncertainty range, no null comparison. Also records that split-vs-displacement already
   failed replication.

HARD RULES - violating these is worse than returning nothing.
1. Quote VERBATIM from sources you actually retrieved. Exact sentences only.
2. NEVER write from memory. Not retrieved => say so. A truthful "not retrieved" is worth
   far more than a confident paraphrase. This project has lost four headlines to
   unverified assumptions.
3. Report every URL you actually opened.
4. Load web tools FIRST, one call:
   ToolSearch query "select:WebSearch,WebFetch" max_results 5
5. Paywall ladder: doi.org -> publisher -> Semantic Scholar API
   (api.semanticscholar.org/graph/v1/paper/DOI:<doi>?fields=title,abstract,openAccessPdf)
   -> OpenAlex (api.openalex.org/works/doi:<doi>) -> CORE -> ResearchGate -> institutional
   repositories -> PhD theses (they restate methods in full and are nearly always open)
   -> papers that QUOTE the passage verbatim.
6. Distinguish what authors CLAIMED from what their method mechanically DOES from what
   later papers say they showed.
`

const SWEEP = {
  type: 'object',
  properties: {
    findings: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          citation: { type: 'string' },
          doi: { type: 'string' },
          url: { type: 'string' },
          what_it_does: { type: 'string' },
          verbatim: { type: 'string' },
          quantitative_claim: { type: 'string' },
          relevance: { type: 'string' },
        },
        required: ['citation', 'url', 'what_it_does', 'verbatim'],
      },
    },
    summary: { type: 'string' },
    could_not_access: { type: 'string' },
  },
  required: ['findings', 'summary', 'could_not_access'],
}

phase('Audit')

const TASKS = [
  {
    label: 'forecast:circularity',
    prompt: `THE MOST IMPORTANT QUESTION IN THIS ROUND. Determine whether SUBSEASONAL FORECAST
SKILL claims rest on outcome-based SSW classification, and whether that makes them circular.

Be careful and fair: predicting a well-defined target is NOT circular merely because the
target is defined by an outcome. Circularity bites when a composite is presented as the
causal surface effect of the event, or when a reported skill gain is computed only over
events selected on the very outcome being predicted. Classify every claim you find into
one of those two categories and JUSTIFY the assignment explicitly.

Find and quote with numbers:
- "forecast skill is enhanced following SSWs" claims (ACC, RPSS, CRPSS, % improvement,
  lead time gained). Include Sigmond, Scinocca, Kharin, Shepherd (2013) Nature Geoscience
  6, 98, doi 10.1038/ngeo1698 -- how are events selected there? Quote the selection rule.
- Karpechko, Charlton-Perez, Balmaseda, Tyrrell, Vitart (2018) GRL on predicting the
  surface after SSWs -- doi 10.1029/2018GL079859. Quote how skill is conditioned.
- Domeisen, Butler, Charlton-Perez et al. (2020) JGR-Atmos, the two-part S2S review of
  SSW predictability -- what does it say about conditioning skill on downward propagation?
- claims about what fraction of SSWs are followed by negative NAO/NAM
- whether any operational centre (ECMWF, Met Office, NOAA CPC) uses the DW/NDW
  distinction, or "downward propagation", in seasonal/monthly outlooks or in published
  operational documentation

THE KEY QUESTION FOR THE PAPER, answer it directly and unambiguously in summary:
does the selection-bias correction threaten any OPERATIONAL or FORECAST-SKILL claim, or
does it only affect retrospective composites? If it only affects retrospective composites,
say so plainly -- a truthful negative here is exactly as valuable as a positive.`,
  },
  {
    label: 'prior-art:closer-than-2019',
    prompt: `Find prior art CLOSER than the White et al. (2019) passage. I need to know if anyone has:
(a) QUANTIFIED the selection-induced share of a DW/NDW surface contrast (any number, any
    percentage, any effect-size decomposition);
(b) proposed a CORRECTION for it;
(c) noted that STRATOSPHERIC classification (Hitchcock PJO, reversal depth, vortex
    strength) carries the same bias -- this is the part I believe nobody has flagged;
(d) applied collider bias / endogenous selection / post-treatment conditioning formally
    to atmospheric composite analysis.

Search terms to actually run: "composite" bias "event selection" stratosphere;
"selection effect" annular mode composite; "by construction" SSW composite; circularity
stratosphere troposphere coupling; "regression to the mean" SSW; conditioning on the
response atmospheric composites; Elwert Winship collider climate; "post-treatment"
conditioning geoscience.

Also check COMMENT/REPLY exchanges and errata in J. Climate, JAS, QJRMS, GRL, JGR-Atmos,
WCD, ACP 2015-2026. And check EGUsphere / EarthArXiv / ESSOAr preprints from 2024-2026,
where such a critique would most likely appear first.

Report the CLOSEST hit even if it only partially overlaps, and say plainly what remains
unclaimed.`,
  },
  {
    label: 'bodycount:S1-studies',
    prompt: `Build the body count. Find published studies that (a) classify SSW events using
POST-ONSET SURFACE information, AND (b) then report a surface composite or between-group
contrast as a quantitative result.

Already confirmed, do not re-derive: Karpechko et al. 2017 QJRMS (origin, 229 cites);
J. Climate 32, 85, 2019 (DW anomalies "around twice that of the total"); ACP 26, 3723,
2026 (NAO -0.762 DW vs +0.088 NDW, follows Karpechko at 850 hPa, and sub-classifies
using 2 m temperature over +40 days).

FIND MORE. Use the OpenAlex cited-by endpoint to enumerate recent citers of Karpechko:
api.openalex.org/works?filter=cites:<openalex_id_for_10.1002/qj.3017>,from_publication_date:2020-01-01
and read meta.count plus the result titles. Then screen the promising ones.

Domains to cover: surface climate/extremes, cold-air outbreaks, air quality/ozone, wind
and solar energy, electricity demand, sea ice, precipitation, snow, health.

For each hit report the actual NUMBER claimed, with units, and the verbatim sentence.
In summary give a COUNT of confirmed studies and the range of effect sizes they report.`,
  },
  {
    label: 'stratospheric:pjo-usage',
    prompt: `Map how widely the STRATOSPHERIC classifications are used, since I have evidence they
carry the same bias (a purely stratospheric split manufactures +0.889 sigma of surface
contrast on random dates with no SSW).

Targets:
- Hitchcock, Shepherd, Manney (2013), J. Climate 26, 2096, doi 10.1175/JCLI-D-12-00202.1,
  the PJO definition. Quote the definition verbatim: variable, level(s), latitude band,
  time window. Does ANY surface/near-surface field enter the DEFINITION? Answer explicitly
  from quoted text.
- How many papers use PJO/non-PJO to compare SURFACE composites? Get citation counts from
  OpenAlex/Semantic Scholar and screen recent citers.
- Other stratosphere-only splits used to predict/compare surface response: vortex
  strength at onset, reversal duration, descent depth, vortex geometry (split vs
  displacement), absorbing vs reflecting events (Kodera, Perlwitz-Harnik).
- Specifically: Maycock & Hitchcock (2015) GRL on split vs displacement annular mode
  signatures, and any paper comparing surface composites between stratospherically
  defined SSW subtypes.

summary: how large is the stratospherically-classified literature compared with the
surface-classified one? Which is bigger?`,
  },
]

const audit = (await parallel(TASKS.map(t => () =>
  agent(RULES + '\n\n' + t.prompt, { label: t.label, phase: 'Audit', schema: SWEEP })
))).filter(Boolean)

log('Audit returned ' + audit.length + '/' + TASKS.length)

const VERDICT = {
  type: 'object',
  properties: {
    refuted: { type: 'boolean' },
    confidence: { type: 'string', enum: ['high', 'medium', 'low'] },
    sources_checked: { type: 'array', items: { type: 'string' } },
    what_i_found: { type: 'string' },
    correction: { type: 'string' },
  },
  required: ['refuted', 'confidence', 'sources_checked', 'what_i_found', 'correction'],
}

phase('Verify')

const verified = await parallel([
  () => agent(RULES + `
YOU ARE AN ADVERSARIAL VERIFIER. Default to REFUTING. If you cannot independently
confirm, that is refutation-by-default.

CLAIM UNDER TEST: "Nobody has quantified the selection-induced share of the DW/NDW
surface contrast, nor proposed a correction, nor noted that stratospheric classification
carries the same bias."

Evidence gathered:
` + JSON.stringify(audit[1], null, 1).slice(0, 20000) + `

Find the scoop. Angles the previous agent likely missed: the ozone/chemistry-climate
community; ENSO and MJO composite methodology, where the identical problem arises and may
be solved under another name; statistical-meteorology journals (J. Climate's statistics
papers, Monthly Weather Review, Environmetrics, Advances in Statistical Climatology);
textbooks on compositing; the detection-and-attribution literature; and any paper using
the phrase "selection effect" or "sampling bias" about atmospheric event composites.

refuted=true means substantive prior art exists. State precisely what remains novel.`,
    { label: 'verify:quantification-scoop', phase: 'Verify', schema: VERDICT, effort: 'high' }),

  () => agent(RULES + `
YOU ARE AN ADVERSARIAL VERIFIER. Default to REFUTING.

CLAIM UNDER TEST: "The selection-bias correction threatens only retrospective composites,
not operational forecast-skill claims."

Evidence gathered:
` + JSON.stringify(audit[0], null, 1).slice(0, 20000) + `

Try to REFUTE it, i.e. find a real forecast-skill or operational claim that IS
contaminated by outcome-based event selection. Angles: verification studies that report
skill only over "downward-propagating" events; papers computing NAO skill composited on
post-onset surface state; ECMWF/Met Office technical memoranda and newsletters; the S2S
project verification literature; any seasonal-outlook methodology document conditioning
on stratospheric-to-surface coupling.

refuted=true means an operational or forecast-skill claim IS affected. Name it and quote
the selection rule verbatim. If after real searching you find none, say so -- that is a
publishable scope limit, not a failure.`,
    { label: 'verify:forecast-scope', phase: 'Verify', schema: VERDICT, effort: 'high' }),
])

return { audit, verified: verified.filter(Boolean) }
