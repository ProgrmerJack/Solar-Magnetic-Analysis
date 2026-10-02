#!/bin/sh
# Build main.pdf and supplementary.pdf (pdflatex; cross-references via xr-hyper need both, twice).
cd "$(dirname "$0")"
for d in supplementary main; do pdflatex -interaction=nonstopmode $d.tex > /dev/null; done
bibtex main > /dev/null; bibtex supplementary > /dev/null
for i in 1 2; do for d in supplementary main; do pdflatex -interaction=nonstopmode $d.tex > /dev/null; done; done
for d in main supplementary; do
  echo "$d: $(grep -c '^!' $d.log) errors, $(grep -c 'undefined' $d.log) undefined, $(grep -c 'Overfull' $d.log) overfull"
done
