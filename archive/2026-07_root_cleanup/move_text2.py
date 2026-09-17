import re

with open('paper/main.tex', 'r', encoding='utf-8') as f:
    text = f.read()

abstract_match = re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}', text, re.DOTALL)
old_abstract = abstract_match.group(1).strip()
new_abstract = r'''Stratospheric sudden warming (SSW) events produce persistent surface weather anomalies, yet their consequences for geophysical hazards remain unexplored.
Here we identify a previously unknown connection: SSW episodes drive a large, reproducible reduction in natural dry slab avalanche activity across four countries and two continents.
In Switzerland (21~winters, 16~SSW events), 14 of 16 events show reduced counts (geometric mean rate ratio $\mathrm{RR} = 0.32$; 95\% CI $[0.20, 0.54]$;  = -1.06$).
Norwegian danger levels and Utah occurrence counts decrease concordantly.
Five-country European data reveal a predictable geographic gradient ( = 0.69$), demonstrating systematic, non-stochastic surface coupling.
A count--rating dissociation (danger \emph{ratings} increase while occurrence \emph{counts} decrease) and a differential trigger response (human-triggered avalanches increase relative to natural) reveal a `loaded gun'' mechanism: modest SSW-associated cold regimes suppress surface-energy triggers ($-57\%$ warming, $-29\%$ rain) while preserving snowpack instability.
This `threshold amplification'' transforms modest stratospheric-origin weather shifts into disproportionately large hazard responses.
Hemispheric verification confirms that every link in the mechanism chain is independently significant, with planetary wave forcing identified as the common cause.'''

text = text.replace(old_abstract, '\n' + new_abstract + '\n')

pls_match = re.search(r'\\noindent\\textbf\{Plain Language Summary\.\}(.*?)(?=\n\n%|$)', text, re.DOTALL)
if pls_match:
    pls_text = pls_match.group(0)
    text = text.replace(pls_text, '')
else:
    print("PLS not found")
    pls_text = ""

start_marker = r'\\subsection\*\{Out-of-sample geographic gradient: ALBINA and expanded Norwegian validation\}'
end_marker = r'\\section\*\{Discussion\}'

sections_match = re.search(f'({start_marker}.*?)(?={end_marker})', text, re.DOTALL)
sections_text = ""
if sections_match:
    sections_text = sections_match.group(1)
    text = text.replace(sections_text, '')
else:
    print("Sections not found")

with open('paper/main.tex', 'w', encoding='utf-8') as f:
    f.write(text)

with open('paper/supplementary_information.tex', 'r', encoding='utf-8') as f:
    si_text = f.read()

insertion_idx = si_text.find('\\clearpage\n') + len('\\clearpage\n')

si_insert = '''
% ============================================================
%  EXTENDED ABSTRACT AND PLAIN LANGUAGE SUMMARY
% ============================================================
\\section{Extended Abstract and Plain Language Summary}
\\subsection{Extended Abstract}
''' + old_abstract + '''

\\subsection{Plain Language Summary}
''' + pls_text.replace('\\noindent\\textbf{Plain Language Summary.}', '').strip() + '''

% ============================================================
%  ADDITIONAL VALIDATIONS AND SPECIFICITY ANALYSES
% ============================================================
\\section{Additional Validations and Specificity Analyses}
''' + sections_text.replace('\\subsection*{', '\\subsection{') + '\n'

si_text = si_text[:insertion_idx] + si_insert + si_text[insertion_idx:]

with open('paper/supplementary_information.tex', 'w', encoding='utf-8') as f:
    f.write(si_text)

print("Done. Text moved to SI.")
