"""
shapes.py -- affine normal forms ("charts") of convex pentahedra and the polynomial data
used in the proof that no homogeneous convex pentahedron is mono-stable.

Faces are numbered 0,1,2 = Q1,Q2,Q3 (lateral faces; Q_i does NOT contain lateral edge i),
3 = A, 4 = B (the two triangles).  All pairs of faces are adjacent except {A, B}.

For a chart with parameters x = (x0, x1, x2) this module provides, as polynomials in x,
    nu[f]        (f = 0..4)   a POSITIVE multiple of the covector l_f of face f
                              (plane l_f . y = 1 in coordinates centred at the centroid),
    eta[(i, j)]  (i < j adjacent)  a POSITIVE multiple of  l_j - l_i,
together with all first partial derivatives.  "Positive multiple" holds at every genuine
(non-degenerate) shape of the chart; the multipliers are products of
    * Delta_k = d_k * vol - n_k . mom  > 0   (the centroid is an interior point),
    * vol > 0,
    * monomials in the chart variables (> 0 at genuine shapes).
Derivation (sympy, exact):  plane k is  n_k . y <= d_k  in the chart frame, the centroid is
g = mom / vol, hence l_k = n_k / (d_k - n_k . g) = vol * n_k / Delta_k  and
l_j - l_i = vol * (Delta_i n_j - Delta_j n_i) / (Delta_i Delta_j).

Chart Z ("cone frame", variables u, v, w in [0,1]):  apex O = 0, B_i = e_i, A_i = t_i e_i,
    t_3 = w, t_2 = w v, t_1 = w v u   (0 < t_1 <= t_2 <= t_3 <= 1;  t_3 = 1 is a pyramid).
Chart P ("prism frame", variables rho, p, q):  B_1 = (0,0,0), B_2 = (1,0,0), B_3 = (0,1,0),
    G = (1/3, 1/3),  sigma = rho (1, p, pq),  s = (1, p, pq),
    A_i = ((1 - sigma_i) B_i + sigma_i G,  s_i)     (rho = 0: prisms with parallel edges;
    rho > 0: affine image of the cone frame with t_i = 1 - sigma_i, apex (G, 1/rho)).
"""
import itertools
import sympy as sp


def require(cond, msg="shapes.py: internal identity check failed"):
    """explicit check (unlike `assert`, not removed by `python -O`)"""
    if not cond:
        raise RuntimeError(msg)

FACES = ["Q1", "Q2", "Q3", "A", "B"]
ADJ_PAIRS = [(i, j) for i in range(5) for j in range(i + 1, 5) if (i, j) != (3, 4)]   # 9 pairs


def _tet_data(verts, tets, sample):
    """vol6 = sum of |det| (6 x volume), mom = sum |det| * (sum of the 4 vertices) (24 x first moment).
    The sign of each determinant is fixed by evaluating at the sample point (must be != 0)."""
    vol6 = 0
    mom = sp.Matrix([0, 0, 0])
    for (a, b, c, d) in tets:
        M = sp.Matrix.hstack(verts[b] - verts[a], verts[c] - verts[a], verts[d] - verts[a])
        det = sp.expand(M.det())
        val = det.subs(sample)
        require(val != 0, "degenerate tetrahedron at the sample point")
        sgn = 1 if val > 0 else -1
        vol6 += sgn * det
        mom += sgn * det * (verts[a] + verts[b] + verts[c] + verts[d])
    return sp.expand(vol6), mom.applyfunc(sp.expand)


def _monomial_content(vec, xs):
    """largest monomial in the chart variables dividing all components of vec"""
    polys = [sp.Poly(c, *xs) for c in vec if sp.expand(c) != 0]
    exps = [min(min(m[k] for m in p.monoms()) for p in polys) for k in range(len(xs))]
    mono = sp.Integer(1)
    for x, e in zip(xs, exps):
        mono *= x ** e
    return mono


def _integer_scale(vec, xs):
    """multiply the polynomial vector by the POSITIVE integer lcm of its coefficient denominators,
    so that all coefficients are integers (then evaluation needs no divisions: exact in ball arithmetic)"""
    den = 1
    for c in vec:
        c = sp.expand(c)
        if c != 0:
            for coeff in sp.Poly(c, *xs).coeffs():
                den = sp.ilcm(den, sp.Rational(coeff).q)
    out = vec.applyfunc(lambda c: sp.expand(c * den))
    for c in out:
        if c != 0:
            require(all(sp.Rational(k).q == 1 for k in sp.Poly(c, *xs).coeffs()), "non-integer coefficient")
    return out


def _finish(name, xs, planes, vol, mom):
    """planes: list of (n_k (3-vector), d_k) for faces 0..4;  centroid = mom / vol_den with
    vol_den = vol (already including the factor 4: g = mom / vol)."""
    Delta = [sp.expand(d * vol - (n.T * mom)[0]) for n, d in planes]
    nu = [_integer_scale(n, xs) for n, d in planes]
    eta, mono = {}, {}
    for (i, j) in ADJ_PAIRS:
        raw = (Delta[i] * planes[j][0] - Delta[j] * planes[i][0]).applyfunc(sp.expand)
        m = _monomial_content(raw, xs)
        red = raw.applyfunc(lambda c: sp.expand(sp.cancel(c / m)))
        # exact identity check:  red * m == raw
        require(all(sp.expand(red[k] * m - raw[k]) == 0 for k in range(3)), "monomial-content identity failed")
        red = _integer_scale(red, xs)          # positive rational multiple with integer coefficients
        eta[(i, j)] = red
        mono[(i, j)] = m
    return dict(name=name, xs=xs, planes=planes, vol=vol, mom=mom, Delta=Delta, nu=nu, eta=eta, mono=mono)


def chart_Z():
    u, v, w = sp.symbols("u v w")
    t = [u * v * w, v * w, w]
    e = [sp.Matrix([1, 0, 0]), sp.Matrix([0, 1, 0]), sp.Matrix([0, 0, 1])]
    O = sp.Matrix([0, 0, 0])
    A = [t[i] * e[i] for i in range(3)]
    B = e
    verts = {"B1": B[0], "B2": B[1], "B3": B[2], "A1": A[0], "A2": A[1], "A3": A[2]}
    # frustum = cone from A1 over the faces not containing A1:  B, and the quad (B2,B3,A3,A2)
    tets = [("B1", "B2", "B3", "A1"), ("A1", "B2", "B3", "A3"), ("A1", "B2", "A3", "A2")]
    sample = {u: sp.Rational(1, 2), v: sp.Rational(2, 3), w: sp.Rational(3, 4)}
    vol6, mom = _tet_data(verts, tets, sample)
    vol = sp.expand(4 * vol6)                       # centroid g = mom / (4 vol6)
    # closed form check:  g_i = (1 - tau t_i) / (4 (1 - tau))
    tau = t[0] * t[1] * t[2]
    for i in range(3):
        require(sp.simplify(mom[i] / vol - (1 - tau * t[i]) / (4 * (1 - tau))) == 0, "centroid closed form mismatch")
    planes = [(-e[0], sp.Integer(0)), (-e[1], sp.Integer(0)), (-e[2], sp.Integer(0)),
              (-sp.Matrix([1, u, u * v]), -u * v * w),            # A:  x1 + u x2 + uv x3 >= uvw
              (sp.Matrix([1, 1, 1]), sp.Integer(1))]              # B:  x1 + x2 + x3 <= 1
    # the planes really contain the vertices they should
    for k, names in enumerate([("A2", "A3", "B2", "B3"), ("A1", "A3", "B1", "B3"), ("A1", "A2", "B1", "B2"),
                               ("A1", "A2", "A3"), ("B1", "B2", "B3")]):
        n, d = planes[k]
        require(all(sp.expand((n.T * verts[nm])[0] - d) == 0 for nm in names), f"plane {k} does not contain {names}")
    return _finish("Z", (u, v, w), planes, vol, mom)


def chart_P():
    rho, p, q = sp.symbols("rho p q")
    sig = [rho, rho * p, rho * p * q]
    s = [sp.Integer(1), p, p * q]
    Bxy = [(sp.Integer(0), sp.Integer(0)), (sp.Integer(1), sp.Integer(0)), (sp.Integer(0), sp.Integer(1))]
    G = (sp.Rational(1, 3), sp.Rational(1, 3))
    B = [sp.Matrix([Bxy[i][0], Bxy[i][1], 0]) for i in range(3)]
    A = [sp.Matrix([(1 - sig[i]) * Bxy[i][0] + sig[i] * G[0], (1 - sig[i]) * Bxy[i][1] + sig[i] * G[1], s[i]])
         for i in range(3)]
    verts = {"B1": B[0], "B2": B[1], "B3": B[2], "A1": A[0], "A2": A[1], "A3": A[2]}
    tets = [("B1", "B2", "B3", "A1"), ("A1", "B2", "B3", "A3"), ("A1", "B2", "A3", "A2")]
    sample = {rho: sp.Rational(1, 5), p: sp.Rational(2, 3), q: sp.Rational(3, 4)}
    vol6, mom = _tet_data(verts, tets, sample)
    vol = sp.expand(4 * vol6)
    nA = (A[1] - A[0]).cross(A[2] - A[0]).applyfunc(sp.expand)
    require(nA[2].subs(sample) > 0, "normal of A not oriented upwards")
    planes = [(sp.Matrix([3, 3, rho]), sp.Integer(3)),            # Q1 through B2, B3
              (sp.Matrix([-3, 0, rho]), sp.Integer(0)),           # Q2 through B1, B3
              (sp.Matrix([0, -3, rho]), sp.Integer(0)),           # Q3 through B1, B2
              (nA, sp.expand((nA.T * A[0])[0])),                  # A (top)
              (sp.Matrix([0, 0, -1]), sp.Integer(0))]             # B (bottom, z >= 0)
    for k, names in enumerate([("A2", "A3", "B2", "B3"), ("A1", "A3", "B1", "B3"), ("A1", "A2", "B1", "B2"),
                               ("A1", "A2", "A3"), ("B1", "B2", "B3")]):
        n, d = planes[k]
        require(all(sp.expand((n.T * verts[nm])[0] - d) == 0 for nm in names), f"plane {k} does not contain {names}")
    ch = _finish("P", (rho, p, q), planes, vol, mom)
    ch["verts"] = verts
    return ch


def build_evaluator(ch):
    """Returns f(x0, x1, x2) -> flat list of values, generic in the number type (float, Fraction, arb):
    for every face f: nu[f] (3 values) and its 3 partial derivatives (9 values);
    for every adjacent pair: eta (3 values) and its partial derivatives (9 values).
    Layout: see `unpack`."""
    xs = ch["xs"]
    exprs = []
    for f in range(5):
        vec = ch["nu"][f]
        exprs += [vec[c] for c in range(3)]
        exprs += [sp.diff(vec[c], xs[k]) for k in range(3) for c in range(3)]
    for pr in ADJ_PAIRS:
        vec = ch["eta"][pr]
        exprs += [vec[c] for c in range(3)]
        exprs += [sp.diff(vec[c], xs[k]) for k in range(3) for c in range(3)]
    exprs = [sp.horner(sp.expand(e), *xs) if sp.expand(e).free_symbols else sp.expand(e) for e in exprs]
    fn = sp.lambdify(xs, exprs, modules=[], cse=True)
    return fn


def unpack(vals):
    """-> nu[f] = (vec, [d/dx0 vec, d/dx1 vec, d/dx2 vec]),  eta[(i,j)] likewise"""
    nu, eta = {}, {}
    pos = 0
    for f in range(5):
        vec = vals[pos:pos + 3]
        der = [vals[pos + 3 + 3 * k: pos + 6 + 3 * k] for k in range(3)]
        nu[f] = (vec, der)
        pos += 12
    for pr in ADJ_PAIRS:
        vec = vals[pos:pos + 3]
        der = [vals[pos + 3 + 3 * k: pos + 6 + 3 * k] for k in range(3)]
        eta[pr] = (vec, der)
        pos += 12
    return nu, eta


_CACHE = {}


def get_chart(name):
    if name not in _CACHE:
        ch = chart_Z() if name == "Z" else chart_P()
        ch["fn"] = build_evaluator(ch)
        _CACHE[name] = ch
    return _CACHE[name]


if __name__ == "__main__":
    for name in ("Z", "P"):
        ch = get_chart(name)
        print("=" * 30, "chart", name, "variables", ch["xs"])
        print("vol =", sp.factor(ch["vol"]))
        for k in range(5):
            print(f"  Delta[{FACES[k]}] =", sp.factor(ch["Delta"][k]), "   nu =", list(ch["nu"][k]))
        for pr in ADJ_PAIRS:
            print(f"  eta[{FACES[pr[0]]}->{FACES[pr[1]]}] (divided by {ch['mono'][pr]}) =",
                  [sp.factor(c) for c in ch["eta"][pr]])
