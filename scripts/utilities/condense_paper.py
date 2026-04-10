#!/usr/bin/env python3
"""Condense main.tex Results and Discussion for Nature Geoscience format."""
import re

with open('paper/main.tex', 'r', encoding='utf-8') as f:
    content = f.read()

# Find boundaries
results_start = content.index('\\section*{Results}')
methods_start = content.index('% ============================================================\n%  METHODS')

before_results = content[:results_start]
after_discussion = content[methods_start:]

# Condensed Results + Discussion
condensed = r"""\section*{Results}

\subsection*{SSW events suppress dry slab avalanche activity across three mountain regions}

Of 16~SSW events during the Swiss study period (1999--2019), 14 show reduced dry slab avalanche counts relative to matched seasonal controls (sign test $P = 0.002$; $\mathrm{BF}_{10} = 64$; median $\Delta = -0.91$~avalanches~day$^{-1}$; Fig.~\ref{fig:ssw_matched}). The two exceptions---February 2018 (split-vortex ``Beast from the East'') and January 2019---are physically explained by SSW-type-dependent surface cooling (see below). Leave-one-out analysis confirms stability: 15/16~jackknife sign tests remain significant at $P < 0.01$ (Extended Data Table~\ref{tab:loo}).

The same matched-comparison design applied to independent databases yields concordant results. In Utah (4~SSW events, 2012--2025), all four show decreased dry slab activity (mean $\Delta = -0.98$; $P = 0.041$), including \textit{decreases} for the 2018 and 2019 events that produced increases in Switzerland---demonstrating region-specific heterogeneity rather than SSW-level inconsistency. In Norway (NVE danger level forecasts; 5~inland mountain regions, 62--70$^\circ$N; $n = 3{,}255$~region-days), SSW-window danger levels are significantly lower than controls ($\Delta = -0.27$; Mann--Whitney $P < 0.0001$; 4/4~events concordant), using a measurement system entirely independent of occurrence counts. Pan-European EAWS danger levels (2011--2015; 141{,}339~region-days) reveal geographically structured responses: Switzerland and Italy show concordant decrease (3/3~events each) while Austria and Germany show increase---consistent with the known SSW-type dependence of European surface impacts\cite{Mitchell2013}. Random-effects meta-analysis of Swiss and Utah data yields a pooled $\Delta = -0.75$~day$^{-1}$ ($P = 0.012$; $I^2 = 32\%$). Across all 24~SSW--region pairs, 22 show decrease ($P < 0.0001$; $\mathrm{BF}_{10} > 1{,}000$).

\subsection*{Temporal structure reveals a common-cause pathway}

A critical diagnostic for causal interpretation is the temporal profile relative to SSW onset. The avalanche anomaly is present \textit{before} SSW onset: 5-day rolling rate ratios reveal peak suppression at lag~$-13$~days (Switzerland) and lag~$-12$~days (Norway), with progressive weakening through onset and return to baseline by day~$+15$. Phase-level aggregation to five non-overlapping SSW lifecycle periods yields near-perfect cross-country concordance (Spearman $r = 0.975$, $P = 0.005$; $n = 5$~phases), demonstrating identical temporal fingerprints across countries separated by 15$^\circ$~latitude using independent measurement systems.

This timing is predicted by the common-cause hypothesis. SSW events are preceded by 1--3~weeks of enhanced tropospheric planetary wave activity\cite{Polvani2004}; under the common-cause framework, this same wave forcing simultaneously reorganises surface weather. The avalanche anomaly should therefore emerge \textit{before} the stratospheric wind reversal---precisely as observed. Phase-resolved ERA5 analysis confirms: pre-SSW warming ($+0.45$~K, $P = 0.15$), peak post-SSW warming ($+0.92$~K, $P = 0.001$), followed by cold-air outbreak after day~$+16$ ($-0.84$~K, $P = 0.006$). Formal mediation analysis rules out annular mode pathways: neither NAO nor AO mediate the SSW--avalanche association (proportion mediated $< 7\%$; bootstrap CIs include zero; Supplementary Information).

\subsection*{SNOWPACK simulations reveal a loading paradox resolved by slab bridging}

SNOWPACK model output across 129~IMIS stations\cite{Bartelt2002,Lehning2002} (2001--2020; 29{,}296~station-day records) reveals a physical paradox: all eight stability indices at persistent weak layers \textit{decrease} after SSW onset (natural stability $\mathrm{sn38\_pwl}$: $-18\%$, $P = 0.004$, Cohen's $d = -0.94$; critical crack length: $-17\%$, $P = 0.008$), while snow depth \textit{increases} ($+38\%$, 15/15~events, $P < 0.0001$) and SWE increases ($+44\%$, $P < 0.0001$). More snow, more weak layers, yet fewer avalanches.

This loading paradox is resolved by slab bridging: the ratio of snow depth to critical crack length ($\mathrm{HS_{mod}/ccl_{pwl}}$) increases in 12/13~events ($P = 0.0005$, $d = 0.98$). Deeper burial distributes stress more uniformly across weak-layer interfaces, raising the natural-release threshold despite reduced point stability\cite{Schweizer2003}. Independent Rutschblock field tests ($n = 4{,}433$) confirm increased bulk slab stability during SSW windows ($P = 0.003$), while ECT tests (probing weak-layer fracture propagation) show no signal ($P = 0.88$)---indicating the effect operates on slab-scale properties rather than weak-layer fracture initiation, consistent with the bridging mechanism.

A process-based sintering model (Arrhenius bond growth, Clausius--Clapeyron vapour transport, Mellor densification) forced with ERA5 temperature quantifies the slab-strengthening pathway. Full-sample results are non-significant ($P = 0.17$), but stratifying by SSW morphological type\cite{Charlton2007,Mitchell2013} resolves this: displacement-type SSWs ($n = 10$), which produce Alpine surface warming, yield $+12.5\%$ slab strength enhancement ($P = 0.0006$; 9/10~positive). Split-vortex SSWs ($n = 6$), which produce surface cooling, yield uniformly negative sintering. The two Swiss outlier events fall precisely in the cooling category, transforming them from unexplained exceptions into physically predicted outcomes.

\subsection*{Planetary wave forcing operates independently of SSW events}

A wave activity index (WAI) constructed from 7-day 10~hPa wind deceleration predicts dry slab counts at 10-day lag even on \textit{non-SSW days} ($r_s = 0.064$, $P = 0.0007$; $n = 2{,}822$), with physically consistent lag structure: no concurrent effect ($P = 0.27$), emergence at 7--10~days, dissipation by lag~14~d. Daily-level atmospheric chain analysis ($n > 3{,}000$~winter days) identifies each mechanistic link: stratospheric warming $\to$ surface warming at lag 1--7~d ($r_s = 0.076$, $P < 0.0001$); surface warming $\to$ reduced dry slabs at lag 5--10~d ($r_s = -0.064$, $P = 0.001$). Each link explains $< 1\%$ of daily variance---consistent with published stratosphere--surface coupling effect sizes\cite{Baldwin2001,Sigmond2013,Kidston2015}. Model comparison establishes a non-linear, regime-dependent association: quintile models strongly outperform both binary SSW and continuous vortex specifications ($\Delta$AIC~$> 149$), while SSW intensity metrics show no dose--response (all $|r| < 0.1$). Together, these results establish that the planetary wave forcing $\to$ avalanche pathway operates continuously, with SSW events as high-amplitude markers of an ongoing wave-driven modulation of surface hazard.

Trigger-type specificity further constrains the mechanism: natural-release avalanches decrease during SSW windows (RR~$= 0.91$) while human-triggered rates are unchanged (RR~$= 1.01$; $P = 0.024$), ruling out confounding by recreational exposure patterns (Supplementary Information).

\subsection*{Robustness}

Event-level permutation inference (1,000~shuffled date sets) yields $P = 0.083$ for the mean anomaly magnitude---higher than the sign test $P = 0.002$ because two large positive outliers inflate the mean. For heavy-tailed distributions (kurtosis~$= 4.2$), the sign test is the appropriate primary inference\cite{Wilcox2017}; the Bayesian sign test ($\mathrm{BF}_{10} = 64$) confirms very strong evidence. Window sensitivity analysis shows robustness across 10--30~day windows; the 20-day weakening is explained by the post-day-$+16$ cold-air outbreak reversing sintering-favourable conditions. LOOCV across 21~winters yields RR~$< 1$ in all folds. Specification curve analysis across 16~model variants shows 69\% with decreased activity. Full robustness details, including NAO-adjusted regression, ZIP model validation, and statistical power calculations, are in the Supplementary Information.

% ============================================================
%  DISCUSSION
% ============================================================
\section*{Discussion}

We report that SSW events mark a discrete reduction in dry slab avalanche hazard---a previously unrecognised connection between stratospheric variability and surface geophysical hazard. The finding is robust across three mountain regions (22/24~SSW--region pairs decrease; $P < 0.0001$; $\mathrm{BF}_{10} > 1{,}000$), multiple measurement systems (occurrence counts, expert danger forecasts, field stability tests, physics-based simulations), and complementary statistical frameworks (non-parametric tests, Bayesian inference, meta-analysis, zero-inflated regression). Five converging lines of evidence---the pre-SSW temporal anomaly, cross-country phase concordance ($r = 0.975$, $P = 0.005$), absent SSW-intensity dose--response, rejected annular mode mediation, and continuous planetary wave--avalanche coupling on non-SSW days ($P = 0.0007$)---collectively support a common-cause interpretation: the planetary wave forcing that triggers SSW simultaneously reorganises mid-latitude surface weather in ways that stabilise the snowpack.

The mechanism operates through a multi-scale snowpack response. At the weak-layer interface, SSW-associated snowfall increases point-instability (SNOWPACK $\mathrm{sn38\_pwl}$: $-18\%$, $d = -0.94$); at the slab scale, deeper burial enhances load-spreading (bridge ratio increase: 12/13~events, $P = 0.0005$, $d = 0.98$), and temperature-driven sintering strengthens slab cohesion (displacement events: $+12.5\%$, $P = 0.0006$). Rutschblock field tests confirm the slab-scale effect ($P = 0.003$) while ECT tests (weak-layer specific) show no signal. This multi-scale reconciliation---increased point-instability at weak layers but increased bulk stability at the slab scale---provides a physically complete resolution of the loading paradox. The SSW-type stratification (displacement vs.\ split\cite{Charlton2007,Mitchell2013}) transforms the full-sample heterogeneity from a limitation into mechanistic evidence, as the two Swiss outlier events are precisely those for which the sintering model predicts reduced slab strength.

Several limitations constrain interpretation. (1)~The Swiss analysis rests on $n = 16$~events from a single 21-year record; the Norwegian NVE data ($n = 4$~events) uses predicted danger rather than observed occurrences. Formal replication with Austrian (ZAMG), French (ANENA), and Japanese (NIED) occurrence-count databases remains a priority. (2)~ERA5 event-level analysis ($n = 8$) is underpowered; confirmation with station-level MeteoSwiss data is desirable. (3)~The 15-day pre-SSW anomaly precludes definitive separation of top-down forcing from common-cause confounding with observational data alone; stratospheric nudging experiments\cite{Hitchcock2013} could provide a causal test. (4)~Our sintering model is simplified (1D, forced with grid-point temperature); full SNOWPACK/Crocus simulations with isolated SSW-conditional forcing would quantify the relative contributions of sintering, bridging, and wind redistribution. (5)~The daily-level chain explains $< 1\%$ of variance per link, reflecting the intrinsic signal-to-noise ratio of multi-step geophysical chains.

Our results carry two practical implications. First, the SSW-type classification alone predicts Swiss dry slab decrease with 90\% accuracy for displacement events (9/10), and SSW morphological type can be determined within 1--2~days of onset, providing a potential 2-week extended-range hazard signal. Second, the planetary wave independence result suggests that continuous monitoring of wave activity indices could flag reduced avalanche hazard $\sim$10~days before the surface effect materialises, regardless of whether a formal SSW event is triggered. More broadly, the non-linear, threshold-like character of this association suggests that similar discrete regime-dependent signals may exist in other surface hazard categories but would be missed by analyses seeking only continuous linear relationships.

"""

# Reconstruct the file
new_content = before_results + condensed + "\n" + after_discussion

with open('paper/main.tex', 'w', encoding='utf-8') as f:
    f.write(new_content)

# Word count check
words = new_content.split()
print(f"New total words: {len(words)}")

sections = re.split(r'\\section\*\{', new_content)
for s in sections[1:]:
    name = s.split('}')[0]
    wc = len(s.split())
    print(f"  {name}: ~{wc} words")
