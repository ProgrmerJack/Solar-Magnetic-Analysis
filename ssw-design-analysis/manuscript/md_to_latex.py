#!/usr/bin/env python3
r"""
md_to_latex.py -- one-time conversion of main.md and supplementary.md to LaTeX.

After this conversion the .tex files in manuscript/latex/ are the manuscript source.
Citations ^{key1,key2}^ become \cite{key1,key2}; refs.bib is written from
number_citations.REFS (each entry's formatted text as its note), and natbib/unsrt
number them as superscripts in order of first citation. Figures are included from 09_figures/out/ after the legends.
The .tex files were later restyled onto latex/preamble.tex (pdflatex); build with latex/build.sh.
"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from number_citations import REFS  # noqa: E402

OUT = HERE / "latex"
FIG = Path("../../09_figures/out")
CIT = re.compile(r"\^\{([A-Za-z0-9,]+)\}\^")
MAIN_FIGS = ["fig1_summary", "fig2_selection", "fig3_regional", "fig4_predictability", "fig5_operational"]
ED_FIGS = [f"figED{i}_{n}" for i, n in enumerate(["shift", "archetypes", "loeffel", "sweep", "heldout", "continuity",
                                                  "shiftrule", "translation", "forecasts"], start=1)]

PREAMBLE = r"""\documentclass[11pt]{article}
\usepackage[a4paper,margin=2.5cm]{geometry}
\usepackage{fontspec}
\setmainfont{Noto Serif}
\usepackage{setspace}\onehalfspacing
\usepackage{graphicx}
\usepackage[hidelinks]{hyperref}
\usepackage{enumitem}
\usepackage[super,sort&compress]{natbib}
\usepackage{xcolor}
\setlength{\parskip}{4pt}
\newcommand{\authors}[1]{\textcolor{red}{[#1]}}
"""


def esc(t):
    t = t.replace("\\", r"\textbackslash{}")
    for a, b in (("&", r"\&"), ("%", r"\%"), ("#", r"\#"), ("_", r"\_"), ("$", r"\$"), ("~", r"\textasciitilde{}")):
        t = t.replace(a, b)
    return t


class Numberer:
    def __init__(self):
        self.order = []

    def cite(self, m):
        nums = []
        for k in m.group(1).split(","):
            if k not in REFS:
                raise KeyError(f"citation key {k!r} not in REFS")
            if k not in self.order:
                self.order.append(k)
            nums.append(self.order.index(k) + 1)
        nums = sorted(set(nums))
        parts, i = [], 0
        while i < len(nums):                               # compress 3,4,5 -> 3-5
            j = i
            while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
                j += 1
            parts.append(f"{nums[i]}--{nums[j]}" if j - i >= 2 else ",".join(map(str, nums[i:j + 1])))
            i = j + 1
        return "@@CITE{" + m.group(1) + "}@@"


def inline(t, num):
    t = CIT.sub(num.cite, t)
    code = []
    t = re.sub(r"`([^`]+)`", lambda m: code.append(m.group(1)) or f"@@CODE{len(code) - 1}@@", t)
    t = re.sub(r"\[\[(.*?)\]\]", lambda m: f"@@AUTH{{{m.group(1)}}}@@", t)
    t = esc(t)
    t = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", t)
    t = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"\\emph{\1}", t)
    t = re.sub(r"@@CITE\{(.*?)\}@@", r"\\cite{\1}", t)
    t = re.sub(r"@@AUTH\{(.*?)\}@@", r"\\authors{\1}", t)
    t = re.sub(r"@@CODE(\d+)@@", lambda m: r"\texttt{" + esc(code[int(m.group(1))]).replace("/", r"/\allowbreak{}") + "}", t)
    t = re.sub(r'"([^"]*)"', "\u201c\\1\u201d", t)
    t = t.replace("\u2264", r"$\le$").replace("\u2265", r"$\ge$").replace("\u221a", r"$\surd$")
    return t


def blocks(md, num, top_level_map):
    out, para, items = [], [], []

    def flush():
        if para:
            out.append(inline(" ".join(para), num) + "\n")
            para.clear()
        if items:
            out.append("\\begin{itemize}[leftmargin=*,itemsep=1pt]\n" + "".join(f"  \\item {inline(i, num)}\n" for i in items)
                       + "\\end{itemize}\n")
            items.clear()
    for ln in md.splitlines():
        s = ln.rstrip()
        if not s.strip():
            flush(); continue
        if s.strip() == "---":
            flush(); continue
        m = re.match(r"^(#{1,3})\s+(.*)$", s)
        if m:
            flush()
            lvl, title = len(m.group(1)), m.group(2)
            cmd = top_level_map.get(lvl, "subsection*")
            out.append(f"\\{cmd}{{{inline(title, num)}}}\n")
            continue
        if re.match(r"^\s*[-*]\s+", s):
            if para:
                flush()
            items.append(re.sub(r"^\s*[-*]\s+", "", s)); continue
        if items and s.startswith("  "):
            items[-1] += " " + s.strip(); continue
        para.append(s.strip())
    flush()
    return "\n".join(out)


def main():
    if (OUT / "main.tex").exists() and "--force" not in sys.argv:
        raise SystemExit("latex/main.tex exists and is now the manuscript source; refusing to overwrite (use --force)")
    OUT.mkdir(exist_ok=True)
    md = (HERE / "main.md").read_text()
    title = re.match(r"#\s+(.*)", md).group(1)
    md = md.split("\n", 1)[1]
    md = re.sub(r"^\*Draft for Nature Geoscience.*?\*\n", "", md, flags=re.S | re.M)
    num = Numberer()
    # split off references placeholder and legends so figures can follow the legends
    md = re.sub(r"## References\s*\n\s*\[\[REFERENCES\]\]\s*\n", "@@REFLIST@@\n", md)
    body = blocks(md, num, {2: "section*", 3: "subsection*"})
    bib = "".join(f"@misc{{{k},\n  note = {{{inline(REFS[k], Numberer())}}}\n}}\n\n" for k in REFS)
    (OUT / "refs.bib").write_text(bib, encoding="utf8", newline="\n")
    body = body.replace("@@REFLIST@@", "{\\small\\bibliographystyle{unsrt}\\bibliography{refs}}\n")
    figs = "\\clearpage\n\\section*{Figures}\n" + "".join(
        f"\\begin{{figure}}[p]\\centering\\includegraphics[width=\\textwidth]{{{FIG / f}.pdf}}\n"
        f"\\caption*{{{lab}}}\\end{{figure}}\n" for f, lab in
        [(f, f"Fig. {i}") for i, f in enumerate(MAIN_FIGS, 1)] + [(f, f"Extended Data Fig. {i}") for i, f in enumerate(ED_FIGS, 1)])
    tex = (PREAMBLE + "\\title{" + inline(title, Numberer()) + "}\n\\author{\\authors{AUTHORS}}\\date{}\n"
           "\\begin{document}\n\\maketitle\n" + body + figs + "\\end{document}\n")
    (OUT / "main.tex").write_text(tex, encoding="utf8", newline="\n")
    si = (HERE / "supplementary.md").read_text()
    si_title = re.match(r"#\s+(.*)", si).group(1); si = si.split("\n", 1)[1]
    sib = blocks(si, num, {2: "section*", 3: "subsection*"})        # continues the numbering
    (OUT / "supplementary.tex").write_text(PREAMBLE + "\\title{" + inline(si_title, Numberer()) + "}\\date{}\n\\begin{document}\n\\maketitle\n"
                                           + sib + "\\end{document}\n", encoding="utf8", newline="\n")
    print(f"{len(num.order)} references; -> {OUT.relative_to(HERE)}/main.tex, supplementary.tex")


if __name__ == "__main__":
    main()
