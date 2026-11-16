# Licences

Copyright © 2026 Claudio Mattioni. This repository is distributed under two licences, chosen by file type.

## Code: MIT

The MIT licence ([`LICENSES/MIT.txt`](LICENSES/MIT.txt)) covers **all program code**: every file ending in
`.py` or `.sh`. This includes:

* the verifiers: `lb_joint/verify_joint.py`, `lb_joint/enumerate_spheres.py`, `lb_joint/flip_closure.py`,
  `penta/verify.py`, `penta/core.py`, `penta/shapes.py`, `verify/exact_verify.py` and
  `examples/V14_short_check.py`;
* the audit checkers in `audit/`;
* the certificate generators and searches;
* `reproduce.py` and `paper/latex/make_figure.py`.

## Everything else: CC BY 4.0

The Creative Commons Attribution 4.0 International licence
([`LICENSES/CC-BY-4.0.txt`](LICENSES/CC-BY-4.0.txt), <https://creativecommons.org/licenses/by/4.0/>) covers
**every other file**. This includes:

* **The manuscript and supplement:** `paper/paper.md`, `paper/supplement.md`, `paper/latex/*.tex` and
  the compiled PDFs in `paper/latex/`.
* **The figures:** `examples/V14_plate.svg` and `paper/latex/V14_plate.pdf`.
* **The certificates:** `lb_joint/certs_V*.jsonl.gz`, `lb_joint/exp1_V8_certG.jsonl.gz`,
  `penta/certs_Z.json.gz` and `penta/certs_P.json.gz`.
* **The other data files:** `*.json`, `*.jsonl.gz`, `*.json.gz` and `*.csv`, among them the example
  polyhedra and their exact reports.
* **The logs** (`*.log`).
* **The documentation:** `README.md`, `REVIEW_GUIDE.md`, `LICENSE.md`, the audit reports and the other
  `*.md` files, `CITATION.cff`, `CHECKSUMS.sha256` and the requirements files.

Attribution: Claudio Mattioni, *Monostatic polyhedra: no mono-stable pentahedron, and 9 ≤ V_U ≤ 14*
(2026), together with the repository or archive in which the material was obtained.

## Third-party material

The files `search_faces/certified/reshetov1*.json` contain the integer plane data of the four
mono-stable polyhedra published by A. Reshetov, *A unistable polyhedron with 14 faces*, Int. J. Comput.
Geom. Appl. 24 (2014) 39–59. They are included only so that those bodies can be re-verified, and should
be attributed to that source. The exact verification reports next to them are our own.

The licence texts in `LICENSES/` are reproduced unchanged.
