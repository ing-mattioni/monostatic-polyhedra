# Preparation checks

Prepared 17 November 2026.

- Python 3.11.9, Windows; verification commands also run with -I -S.
- Three exact paths confirm the original fourteen-vertex example.
- The separate integer checker also verifies the interval 0 < t <= 1 for U = 1.
- The development repository's existing verify/test_exact_verify.py suite passed
  all 23 tests, including degenerate features, hull robustness, rational input,
  CLI failure modes and independent NumPy/SciPy comparisons. That broader suite
  and its unrelated seeds are outside this minimal archive.
- Tectonic 0.17.0 produced the 8-page PDF. Its TeX log has no overfull/underfull
  boxes, unresolved references or LaTeX warnings. The process emitted a
  Fontconfig configuration notice; fonts were embedded and all rendered pages
  were visually inspected.
- Numerical coordinate, facet and witness tables are checked against exact
  data; the 24-facet and 36-edge sign ledger is regenerated and compared.
- Pandoc 3.14 produces the Markdown companion from the LaTeX source, with
  citations and local file paths preserved.
- The supplied clean-extraction transcript is release/reproduce_extracted.log.
  The final archive is extracted and verified again by the package builder.

These checks are computational and editorial validation, not external human
peer review or a guarantee of priority.
