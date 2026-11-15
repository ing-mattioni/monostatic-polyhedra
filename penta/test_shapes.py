"""Independent numerical cross-check of the chart data of shapes.py against direct geometry:
for random parameters, build the vertices of the pentahedron in the chart frame, compute the centroid by
a tetrahedral decomposition and the face covectors l_f from the vertices, and check that
nu[f] is a positive multiple of l_f and eta[(i,j)] a positive multiple of l_j - l_i."""
import numpy as np, sys
from shapes import get_chart, unpack, ADJ_PAIRS, FACES

def verts(chart, x):
    if chart == "Z":
        u, v, w = x; t = [u*v*w, v*w, w]
        B = np.eye(3); A = np.array([t[i]*B[i] for i in range(3)])
    else:
        rho, p, q = x; sig = [rho, rho*p, rho*p*q]; s = [1.0, p, p*q]
        Bxy = np.array([[0,0],[1,0],[0,1]], float); G = np.array([1/3, 1/3])
        B = np.c_[Bxy, np.zeros(3)]
        A = np.array([[*( (1-sig[i])*Bxy[i] + sig[i]*G ), s[i]] for i in range(3)])
    return A, B

def direct(A, B):
    A1, A2, A3 = A; B1, B2, B3 = B
    vol = 0; mom = np.zeros(3)
    for (a, b, c, d) in [(B1,B2,B3,A1), (A1,B2,B3,A3), (A1,B2,A3,A2)]:
        det = abs(np.dot(np.cross(b-a, c-a), d-a)); vol += det; mom += det*(a+b+c+d)/4
    g = mom/vol
    faces = [(A2,A3,B3,B2), (A3,A1,B1,B3), (A1,A2,B2,B1), (A1,A2,A3), (B1,B2,B3)]
    L = []
    for poly in faces:
        n = np.zeros(3)
        for k in range(len(poly)): n += np.cross(poly[k], poly[(k+1) % len(poly)])
        if np.dot(n, poly[0] - g) < 0: n = -n
        L.append(n/np.dot(n, poly[0] - g))
    return np.array(L)

rng = np.random.default_rng(1)
worst = 0
for chart in ("Z", "P"):
    ch = get_chart(chart)
    for it in range(3000):
        x = rng.uniform(0.02, 0.98, size=3)
        if chart == "P": x[0] *= 0.6
        A, B = verts(chart, x)
        L = direct(A, B)
        nu, eta = unpack([float(t) for t in ch["fn"](*x)])
        for f in range(5):
            a = np.array(nu[f][0]); b = L[f]
            cosang = a@b/np.linalg.norm(a)/np.linalg.norm(b)
            assert cosang > 0, (chart, x, FACES[f], "nu not a positive multiple")
            worst = max(worst, 1 - cosang)
        for (i, j) in ADJ_PAIRS:
            a = np.array(eta[(i, j)][0]); b = L[j] - L[i]
            cosang = a@b/np.linalg.norm(a)/np.linalg.norm(b)
            assert cosang > 0, (chart, x, FACES[i], FACES[j], "eta not a positive multiple")
            worst = max(worst, 1 - cosang)
        # derivative check by central differences (first component only, cheap)
        if it < 200:
            eps = 1e-6
            for k in range(3):
                xp = x.copy(); xp[k] += eps; xm = x.copy(); xm[k] -= eps
                fp = [float(t) for t in ch["fn"](*xp)]; fm = [float(t) for t in ch["fn"](*xm)]
                nup, etap = unpack(fp); num, etam = unpack(fm)
                for f in range(5):
                    fd = (np.array(nup[f][0]) - np.array(num[f][0]))/(2*eps)
                    assert np.allclose(fd, nu[f][1][k], rtol=1e-5, atol=1e-6), ("dnu", chart, f, k)
                for pr in ADJ_PAIRS:
                    fd = (np.array(etap[pr][0]) - np.array(etam[pr][0]))/(2*eps)
                    assert np.allclose(fd, eta[pr][1][k], rtol=1e-5, atol=1e-5), ("deta", chart, pr, k)
    print(f"chart {chart}: 3000 random shapes OK (covector directions and derivatives agree with direct geometry)")
print("max 1 - cos(angle) over all comparisons:", worst)
