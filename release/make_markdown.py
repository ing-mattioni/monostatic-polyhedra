#!/usr/bin/env python3
"""Generate the readable Markdown companion from main.tex with Pandoc 3.14."""
import argparse
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
LINKS = {
 "DKLRV20": "https://arxiv.org/abs/1810.05382v3",
 "Reg25": "https://doi.org/10.1556/096.2025.00142",
 "Public16": "https://github.com/shuhualin1023/mono-unstable-16/tree/cf21399a241e48a1ca872b20846be46f6a7b8025",
 "BDPR24": "https://arxiv.org/abs/2401.17906",
 "PRDB23": "https://doi.org/10.1016/j.ejor.2023.04.028",
 "Res14": "https://doi.org/10.1142/S0218195914500022",
 "Hum26": "https://github.com/ohresearchresults/gomboc-research/releases/tag/research-2026-10-12.3",
}
REFS = {
 "sec:intro":"1", "sec:mechanism":"2", "sec:proof":"3", "sec:census":"3.4",
 "sec:family":"4", "sec:verification":"5", "app:ledger":"A",
 "thm:main":"1", "cor:family":"1", "lem:vertex":"1", "fig:plate":"1",
 "tab:comparison":"1", "tab:points":"2", "tab:witnesses":"3", "tab:faceledger":"4", "tab:edgeledger":"5",
 "eq:thickness":"1", "eq:flatten":"2", "eq:facets":"3", "eq:support":"4",
 "eq:mass":"5", "eq:centroid":"6", "eq:facepredicate":"7", "eq:edgepredicate":"8",
}


def convert(pandoc):
    tex = (ROOT/"paper/latex/main.tex").read_text(encoding="utf-8")
    tex = tex.replace(r"\input{equilibrium_ledger.tex}", (ROOT/"paper/latex/equilibrium_ledger.tex").read_text())
    tex = re.sub(r"\\path\{([^}]+)\}", lambda m: r"\texttt{"+m[1].replace("_", r"\_")+"}", tex)
    tex = re.sub(r"\\cite(?:\[([^]]*)\])?\{([^}]+)\}",
                 lambda m: r"\href{"+LINKS[m[2]]+"}{["+str(list(LINKS).index(m[2])+1)+
                 (", "+m[1] if m[1] else "")+"]}", tex)
    tex = re.sub(r"\\(eqref|ref)\{([^}]+)\}",
                 lambda m: "("+REFS[m[2]]+")" if m[1] == "eqref" else REFS[m[2]], tex)
    tex = tex.replace(r"\begin{thebibliography}{10}", r"\section*{References}")
    tex = tex.replace(r"\end{thebibliography}", "")
    tex = re.sub(r"\\bibitem\{([^}]+)\}",
                 lambda m: r"\noindent\textbf{["+str(list(LINKS).index(m[1])+1)+"]} ", tex)
    abstract = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", tex, re.S)[1].strip()
    tex = re.sub(r"\\begin\{abstract\}.*?\\end\{abstract\}", "", tex, flags=re.S)
    result = subprocess.run([pandoc, "--from=latex", "--to=gfm", "--wrap=none"],
                            input=tex, text=True, encoding="utf-8", capture_output=True,
                            cwd=ROOT/"paper/latex", check=True)
    md = result.stdout
    md = re.sub(r"<figure.*?</figure>", 
                "![Projected outline, signed heights and escape witnesses](../examples/V14_plate.svg)\n\n"
                "*Figure 1. Projection onto the (x,y)-plane; labels give vertex indices and signed heights. "
                "Arrows are escape witnesses, and the cross marks the projected volume centroid.*", md, flags=re.S)
    md = re.sub(r" \{#[^}]+\}", "", md)
    md = md.replace('<div class="proof">', "").replace("</div>", "")
    md = re.sub(r'<div[^>]*>', '', md)
    # All internal references are literal numbers; citations remain direct links.
    md = re.sub(r"\n{3,}", "\n\n", md).strip()
    header = ("# A homogeneous mono-unstable polyhedron with fourteen vertices\n\n"
              "**Claudio Mattioni** — claudio.mattioni@gmail.com\n\n"
              "17 November 2026 · Version 1.0 · Preprint, not peer reviewed\n\n"
              "DOI reserved for this deposit: [10.5281/zenodo.24286668](https://doi.org/10.5281/zenodo.24286668).\n\n"
              "This companion is generated from [the LaTeX source](latex/main.tex). "
              "The [PDF](latex/main.pdf) fixes the equation, table and theorem numbering.\n\n"
              "## Abstract\n\n"+abstract+"\n\n")
    return header+md+"\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pandoc", default="pandoc")
    args = parser.parse_args()
    (ROOT/"paper/paper.md").write_text(convert(args.pandoc), encoding="utf-8", newline="\n")
    print("Generated paper/paper.md from the C-only LaTeX manuscript.")
