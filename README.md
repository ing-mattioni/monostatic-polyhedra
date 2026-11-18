# A homogeneous mono-unstable polyhedron with fourteen vertices

Claudio Mattioni · Version 1.0 · Prepared 17 November 2026 · Preprint, not peer reviewed

An explicit convex polyhedron with uniform **volume** density has **14 vertices,
36 edges and 24 triangular facets**, with **S = 2, H = 1, U = 1** about its actual
volume centroid. Every positive vertical contraction by a factor at most one
also has exactly one unstable equilibrium. The full S/H census is asserted
for the supplied integer example, not for every contraction.

Read the [paper PDF](paper/latex/main.pdf) or open the
[interactive demo (English)](https://claude.ai/artifact/FLTsDzvhxB97iQ1TL2zFbP)
to see why V0 is the only vertex on which the solid can balance.
The [Markdown companion](paper/paper.md) and [review guide](REVIEW_GUIDE.md) are also available.

## Reproduce the proof

Python 3.10 or later; Python 3.11.9 has been tested. No pip packages, virtual
environment or network access are required.

```sh
python reproduce.py
```

Eight checks cover file integrity, three exact checking paths, regenerated
reports, numerical manuscript tables and the complete equilibrium sign ledger.
The script writes reproduce.log and exits nonzero if any check fails.
For explicit isolation from installed packages:

```sh
python -I -S reproduce.py
```

The supplied transcript is [release/reproduce_extracted.log](release/reproduce_extracted.log).
Do not use that filename for a new log. The SHA-256 manifest normalizes CRLF
to LF for text files; binary files are hashed unchanged. It detects accidental
changes, but is not a digital signature or a substitute for reading the proof.

## Contents

- [examples/V14_integer.json](examples/V14_integer.json): exact coordinates.
- [examples/V14_short_check.py](examples/V14_short_check.py): short integer existence proof.
- [verify/exact_verify.py](verify/exact_verify.py): general rational hull and equilibrium checker.
- [verify/check_c.py](verify/check_c.py): separate integer hull, complete census and family certificate.
- [examples/V14_exact_report.json](examples/V14_exact_report.json) and
  [examples/V14_family_report.json](examples/V14_family_report.json): full exact reports.
- [paper/latex/](paper/latex/): manuscript, sign ledger and figure.
- [release/BUILD.md](release/BUILD.md): manuscript and package build instructions.
- [release/PROVENANCE.md](release/PROVENANCE.md): scope, prior work and provenance.

All three proof checkers use arbitrary-precision integers or exact fractions.
Decimal values in the general checker's report are display aids only.
The figure generator optionally uses matplotlib, and PDF compilation requires
a TeX engine; neither is needed to verify the mathematics.

## Scientific scope

This is a constructive upper bound V_U <= 14, not a minimality proof.
A public sixteen-vertex construction is explicitly cited. The construction
has more edges and facets than that example; it is not an improvement in
every size measure. The body has two stable equilibria and is not a polyhedral
Gomboc. The literature comparison is limited to the sources found and checked.

The draft DOI is [10.5281/zenodo.24286668](https://doi.org/10.5281/zenodo.24286668).
It was reserved during preparation and will resolve after the record is published.
The private Git repository is not required for reproduction.

## Earlier repository material

The complete development repository also retains the historical unpublished
A/B/C manuscript in paper/history/abc-94a1f1e and its workflow reproduce_abc.py.
Those files and the A/B certificates are deliberately outside this C-only
deposit. Old theorem labels in the unchanged short checker refer to that
historical manuscript: former Theorem C is now Theorem 1, and its cited
vertex lemma is now Lemma 1.

## Author, AI assistance and licensing

Claudio Mattioni · claudio.mattioni@gmail.com.
Searches, code and writing received substantial assistance from Claude
(Anthropic) and Codex (OpenAI). No external human expert review is claimed.

Code (*.py and *.sh): MIT. All other original material: CC BY 4.0.
See [LICENSE.md](LICENSE.md) and the full texts in LICENSES/.
