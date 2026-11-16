"""
derive_check.py -- INDEPENDENT derivation (audit): face covectors l_k of the pentahedron of each chart,
computed from the VERTICES only (own centroid formula via the divergence theorem on a consistently
oriented triangulated boundary, own plane computation), then an exact symbolic comparison with the
polynomial vectors nu/eta that the certificates refer to (shapes_data.json, dumped from shapes.py).

Prints, for every vector, the ratio  nu_k / l_k  resp.  eta_ij / (l_j - l_i)  in factored form.
"""
import json
import sympy as sp

u, v, w = sp.symbols("u v w")
rho, p, q = sp.symbols("rho p q")


def vertices(chart):
    if chart == "Z":
        t = [u * v * w, v * w, w]
        B = [sp.Matrix([1, 0, 0]), sp.Matrix([0, 1, 0]), sp.Matrix([0, 0, 1])]
        A = [t[i] * B[i] for i in range(3)]
        xs = (u, v, w)
    else:
        sig = [rho, rho * p, rho * p * q]
        s = [1, p, p * q]
        Bxy = [(0, 0), (1, 0), (0, 1)]
        G = (sp.Rational(1, 3), sp.Rational(1, 3))
        B = [sp.Matrix([a, b, 0]) for a, b in Bxy]
        A = [sp.Matrix([(1 - sig[i]) * Bxy[i][0] + sig[i] * G[0], (1 - sig[i]) * Bxy[i][1] + sig[i] * G[1], s[i]])
             for i in range(3)]
        xs = (rho, p, q)
    V = {"A1": A[0], "A2": A[1], "A3": A[2], "B1": B[0], "B2": B[1], "B3": B[2]}
    return V, xs


# consistently oriented boundary (every edge once in each direction); face index = 0..4 = Q1,Q2,Q3,A,B
FACES = [("B3", "B2", "A2", "A3"),      # Q1 (lateral edges 2,3)
         ("B1", "B3", "A3", "A1"),      # Q2 (lateral edges 1,3)
         ("B2", "B1", "A1", "A2"),      # Q3 (lateral edges 1,2)
         ("A1", "A3", "A2"),            # A
         ("B1", "B2", "B3")]            # B


def check_orientation():
    edges = {}
    for f in FACES:
        for k in range(len(f)):
            e = (f[k], f[(k + 1) % len(f)])
            assert e not in edges
            edges[e] = 1
    for (a, b) in edges:
        assert (b, a) in edges
    return True


def centroid(V):
    """divergence theorem: signed cone volumes from the origin over the fan triangulation of each
    oriented face; an overall orientation sign cancels in mom/vol."""
    vol = 0
    mom = sp.zeros(3, 1)
    for f in FACES:
        for k in range(1, len(f) - 1):
            a, b, c = V[f[0]], V[f[k]], V[f[k + 1]]
            d = sp.Matrix.hstack(a, b, c).det()
            vol += d
            mom += d * (a + b + c)
    vol = sp.expand(vol)          # = 6 * signed volume
    mom = mom.applyfunc(sp.expand)  # = 24 * signed first moment
    return vol, mom, mom / (4 * vol)


def covectors(chart):
    check_orientation()
    V, xs = vertices(chart)
    vol, mom, g = centroid(V)
    L = []
    for f in FACES:
        P = [V[nm] for nm in f]
        if len(P) == 3:
            n = (P[1] - P[0]).cross(P[2] - P[0])
        else:
            n = (P[2] - P[0]).cross(P[3] - P[1])
        n = n.applyfunc(sp.expand)
        d = sp.expand((n.T * P[0])[0])
        for X in P:   # planarity
            assert sp.expand((n.T * X)[0] - d) == 0
        e = sp.cancel(d - (n.T * g)[0])
        L.append(n.applyfunc(lambda c: sp.cancel(c / e)))
    return L, xs, vol, g


def load_data(chart, xs):
    D = json.load(open(r"C:\Projects\Claude4\audit\shapes_data.json"))[chart]
    def poly(terms):
        return sp.Add(*[c * sp.Mul(*[x ** k for x, k in zip(xs, m)]) for m, c in terms])
    nu = {int(f): sp.Matrix([poly(t) for t in D["nu"][f]]) for f in D["nu"]}
    eta = {tuple(int(a) for a in k.split(",")): sp.Matrix([poly(t) for t in D["eta"][k]]) for k in D["eta"]}
    return nu, eta


def ratio(vec, l):
    """exact check vec x l == 0, return the scalar r with vec = r l (rational function)"""
    cr = vec.cross(l).applyfunc(sp.cancel)
    assert all(c == 0 for c in cr), "not parallel"
    for c in range(3):
        if sp.cancel(l[c]) != 0:
            r = sp.cancel(vec[c] / l[c])
            assert all(sp.cancel(vec[k] - r * l[k]) == 0 for k in range(3))
            return r
    raise AssertionError("zero covector")


if __name__ == "__main__":
    names = ["Q1", "Q2", "Q3", "A", "B"]
    for chart in ("Z", "P"):
        L, xs, vol, g = covectors(chart)
        print("=" * 20, "chart", chart, " 6*vol =", sp.factor(vol))
        nu, eta = load_data(chart, xs)
        for f in range(5):
            r = ratio(nu[f], L[f])
            print(f"  nu[{names[f]}] / l = {sp.factor(r)}")
        for (i, j) in sorted(eta):
            r = ratio(eta[(i, j)], (L[j] - L[i]).applyfunc(sp.cancel))
            print(f"  eta[{names[i]},{names[j]}] / (l_j - l_i) = {sp.factor(r)}")
