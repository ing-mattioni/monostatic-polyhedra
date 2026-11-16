"""Extract (as DATA) the polynomial vectors nu[f], eta[(i,j)] that the certificates refer to.
This is the only place where penta/shapes.py is imported; recheck.py verifies the dumped polynomials
independently (exactly) against covectors derived from the vertices."""
import sys, json
sys.path.insert(0, r"C:\Projects\Claude4\penta")
import sympy as sp
import shapes
out = {}
for name, fn in (("Z", shapes.chart_Z), ("P", shapes.chart_P)):
    ch = fn()
    xs = ch["xs"]
    def pd(e):
        p = sp.Poly(sp.expand(e), *xs)
        return [[list(m), int(c)] for m, c in zip(p.monoms(), p.coeffs())]
    out[name] = {"xs": [str(x) for x in xs],
                 "nu": {str(f): [pd(ch["nu"][f][c]) for c in range(3)] for f in range(5)},
                 "eta": {f"{i},{j}": [pd(ch["eta"][(i, j)][c]) for c in range(3)] for (i, j) in shapes.ADJ_PAIRS}}
json.dump(out, open("shapes_data.json", "w"))
print("dumped")
