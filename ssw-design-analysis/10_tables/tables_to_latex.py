#!/usr/bin/env python3
"""
tables_to_latex.py -- render the display tables for the LaTeX manuscript.

Reads the Markdown tables written by tableED1_tests.py and tableS1_shiftrule.py
(09_figures/out/*.md) and writes booktabs LaTeX next to them:
  tableED1_tests.tex      Extended Data Table 1 (landscape longtable, \\input by main.tex)
  tableS_shiftrule.tex    Supplementary Tables 1-4 (\\input by supplementary.tex)
Recomputes nothing: cells are copied, only typeset.

texify() is also the rule set used to convert the manuscript prose to LaTeX
(Unicode minus, sigma, degrees, units, p/r/n statements -> math mode).
Run after the two table scripts:  python 10_tables/tables_to_latex.py
"""
import re
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "09_figures" / "out"
REGION = {"NEURASIA": "N. Eurasia", "HI_EUROPE": "High-lat. Europe", "MID_EASIA": "Mid-lat. E. Asia",
          "MID_NAMER": "Mid-lat. N. America"}
HEAD = {"q": "$q$", "rule p": "rule $p$", "p (rule)": "$p$ (rule)", "p (clim.)": "$p$ (clim.)", "f [95%]": "$f$ [95\\%]",
        "RR obs": "RR obs.", "region": "Region", "window": "Days", "observed": "Observed", "RR rule": "RR rule",
        "pooled with V1": "Pooled with 1959--2022"}
SYM = {"σ²": r"$\sigma^2$", "σ": r"$\sigma$", "ρ": r"$\rho$", "Δ": r"$\Delta$", "ε": r"$\varepsilon$",
       "φ": r"$\varphi$", "Φ": r"$\Phi$", "π": r"$\pi$", "×": r"$\times$", "≤": r"$\le$", "≥": r"$\ge$",
       "≈": r"$\approx$", "±": r"$\pm$", "R²": r"$R^2$", "²": r"$^2$", "⁻¹": r"$^{-1}$", "−": r"$-$"}


def _protected(s):
    """Split s into (is_protected, text) runs: math, and commands whose argument must not change."""
    out, i, buf = [], 0, []
    cmds = ("\\cite", "\\texttt", "\\label", "\\ref", "\\url", "\\includegraphics", "\\href", "\\bibliography",
            "\\bibliographystyle", "\\input", "\\externaldocument")
    while i < len(s):
        if s[i] == "$":
            j = s.index("$", i + 1)
            out.append((False, "".join(buf))); buf = []
            out.append((True, s[i:j + 1])); i = j + 1
            continue
        hit = next((c for c in cmds if s.startswith(c, i) and not s[i + len(c):i + len(c) + 1].isalpha()), None)
        if hit:
            j = i + len(hit)
            while j < len(s) and s[j] in "[{":                  # [opt]{arg}{arg}
                close = "]" if s[j] == "[" else "}"
                depth, k = 0, j
                while True:
                    if s[k] in "[{" and (s[k] == s[j]):
                        depth += 1
                    elif s[k] == close:
                        depth -= 1
                        if depth == 0:
                            break
                    k += 1
                j = k + 1
            out.append((False, "".join(buf))); buf = []
            out.append((True, s[i:j])); i = j
            continue
        buf.append(s[i]); i += 1
    out.append((False, "".join(buf)))
    return out


def _plain(t):
    num = r"\d[\d,]*(?:\.\d+)?"
    t = re.sub(r"([+−]?\d+)\.\.([+−]?\d+)", r"\1 to \2", t)                       # days +8..+25
    t = re.sub(r"(?<![\w)\]$])([−+±])(" + num + r")(σ²|σ)?",
               lambda m: "$" + {"−": "-", "+": "+", "±": r"\pm"}[m.group(1)] + m.group(2).replace(",", "{,}")
               + {"σ²": r"\sigma^2", "σ": r"\sigma", None: ""}[m.group(3)] + "$", t)
    t = re.sub(r"(?<![\w.$])(\d+(?:\.\d+)?)(σ²|σ)", lambda m: "$" + m.group(1) + (r"\sigma^2" if m.group(2) == "σ²" else r"\sigma") + "$", t)
    t = re.sub(r"m s⁻¹", r"m\\,s$^{-1}$", t)
    t = re.sub(r"°\s?([NSEW])\b", r"$^\\circ$\\,\1", t)
    t = t.replace("°", r"$^\circ$")
    for k in sorted(SYM, key=len, reverse=True):
        t = t.replace(k, SYM[k])
    t = re.sub(r"(?<![\w\\$])([prnqzG]) (=|<|>|\$\\le\$|\$\\ge\$) (\$?[-+]?\d[\d.,]*\$?)",
               lambda m: "$" + m.group(1) + " " + m.group(2).strip("$") + " " + m.group(3).strip("$") + "$", t)
    t = re.sub(r"(\d|\$) (hPa|K|Pa|km|gpm)\b", r"\1\\,\2", t)
    t = t.replace("$$", "")                                                       # adjacent inline math
    return t


def texify(s):
    return "".join(x if prot else _plain(x) for prot, x in _protected(s))


def cell(x):
    x = x.strip()
    for k, v in REGION.items():
        x = x.replace(k, v)
    x = (x.replace("\\", r"\textbackslash{}").replace("&", r"\&").replace("%", r"\%").replace("_", r"\_")
         .replace("#", r"\#"))
    x = re.sub(r"(?<![\w.])-(\d)", "−\\1", x)                                     # ASCII minus in numbers
    x = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", x)
    return texify(x)


def parse(md):
    """Blocks of a Markdown file: ('title', s), ('text', s) or ('table', [rows])."""
    blocks, rows = [], []
    for ln in md.splitlines() + [""]:
        if ln.startswith("|"):
            if not re.match(r"^\|[\s\-:|]+\|$", ln):
                rows.append([c for c in ln.strip().strip("|").split("|")])
            continue
        if rows:
            blocks.append(("table", rows)); rows = []
        if ln.strip().startswith("**") and ln.strip().endswith("**"):
            blocks.append(("title", ln.strip().strip("*")))
        elif ln.strip():
            blocks.append(("text", ln.strip()))
    return blocks


def tabular(rows, spec=None, size="\\footnotesize"):
    n = len(rows[0])
    spec = spec or "l" + "r" * (n - 1)
    head = " & ".join(r"\textbf{" + (HEAD[c.strip()] if c.strip() in HEAD else cell(c)) + "}" for c in rows[0]) + r" \\"
    body = "\n".join(" & ".join(cell(c) for c in r + [""] * (n - len(r))) + r" \\" for r in rows[1:])
    return (f"{{{size}\\setlength{{\\tabcolsep}}{{4pt}}\\begin{{adjustbox}}{{max width=\\textwidth}}\\begin{{tabular}}{{{spec}}}\n\\toprule\n{head}\n\\midrule\n"
            f"{body}\n\\bottomrule\n\\end{{tabular}}\\end{{adjustbox}}}}")


def supplementary_tables():
    blocks = parse((OUT / "tableS1_shiftrule.md").read_text(encoding="utf8"))
    out, cur = [], None
    for kind, x in blocks:
        if kind == "title":
            if cur:
                out.append(cur)
            m = re.match(r"Supplementary Table (\d+) \| (.*)", x)
            cur = {"n": m.group(1), "title": m.group(2), "parts": []}
        else:
            cur["parts"].append((kind, x))
    out.append(cur)
    tex = ["% generated by 10_tables/tables_to_latex.py from tableS1_shiftrule.md -- do not edit"]
    for t in out:
        tex.append("\\begin{table}[p]\\centering")
        notes = [cell(x) for k, x in t["parts"] if k == "text"]
        tex.append(f"\\caption{{\\textbf{{{cell(t['title'])}}}{(' ' + ' '.join(notes)) if notes else ''}}}\\label{{tab:{t['n']}}}")
        for k, x in t["parts"]:
            if k == "table":
                tex.append(tabular(x)); tex.append("\\par\\medskip")
        tex.append("\\end{table}")
    (OUT / "tableS_shiftrule.tex").write_text("\n".join(tex) + "\n", encoding="utf8", newline="\n")
    print("->", (OUT / "tableS_shiftrule.tex").name, f"({len(out)} tables)")


def ed_table1():
    rows = next(x for k, x in parse((OUT / "tableED1_tests.md").read_text(encoding="utf8")) if k == "table")
    spec = (r"@{}>{\raggedright\arraybackslash}p{1.5cm}>{\raggedright\arraybackslash}p{6.8cm}"
            r">{\raggedright\arraybackslash}p{2.3cm}>{\raggedright\arraybackslash}p{3.0cm}p{1.3cm}"
            r">{\raggedright\arraybackslash}p{5.6cm}p{1.1cm}p{1.1cm}@{}")
    head = " & ".join(r"\textbf{" + c + "}" for c in ("ID", "Test", "Role", "Registration", "Commit", "Statistic", "$p$",
                                                       r"$p_{\mathrm{Holm}}$")) + r" \\"
    body = "\n".join(" & ".join([r"\texttt{" + cell(r[0]) + "}"] + [cell(c) for c in r[1:4]] + [r"\texttt{" + r[4].strip() + "}"]
                               + [cell(c) for c in r[5:]]) + r" \\" for r in rows[1:])
    tex = ("% generated by 10_tables/tables_to_latex.py from tableED1_tests.md -- do not edit\n"
           "\\begin{landscape}\n{\\scriptsize\\setlength{\\tabcolsep}{3pt}\\renewcommand{\\arraystretch}{1.15}\n"
           f"\\begin{{longtable}}{{{spec}}}\n"
           "\\caption{\\textbf{Test register.} \\edtableonelegend}\\label{edtab:1}\\\\\n"
           f"\\toprule\n{head}\n\\midrule\n\\endfirsthead\n"
           f"\\multicolumn{{8}}{{l}}{{\\textit{{Extended Data Table 1, continued}}}}\\\\\n\\toprule\n{head}\n\\midrule\n\\endhead\n"
           "\\bottomrule\n\\endlastfoot\n"
           f"{body}\n\\end{{longtable}}}}\n\\end{{landscape}}\n")
    (OUT / "tableED1_tests.tex").write_text(tex, encoding="utf8", newline="\n")
    print("->", (OUT / "tableED1_tests.tex").name, f"({len(rows) - 1} rows)")


if __name__ == "__main__":
    ed_table1()
    supplementary_tables()
