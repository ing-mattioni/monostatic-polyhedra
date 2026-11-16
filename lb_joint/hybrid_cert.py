"""
hybrid_cert.py -- "hybrid" two-cone certificates (record kind "certH" in verify_joint.py).

Setting as in certify_joint.py: (G, v, tree p), shadow matrices S_i = S_{i->p(i)}, Gram matrix R of the
vertices (PSD, nonzero).  Let C_v and C_q be the cone decompositions from v and from a second vertex q
(tetrahedra conv(v, F), F a triangle of G not containing v; resp. for q).  Split C_v = Good u Bad.

CERTIFICATE: vectors mu, mu' in Z^V and multipliers y >= 0 (y_v = 0) such that
    for every T in Good:              2 M_T   = sum_i y_i (2 S_i) + mu 1_T^T + 1_T mu^T                       > 0,
    for every T in Bad, T' in C_q:    2 M_TT' = sum_i y_i (2 S_i) + mu 1_T^T + 1_T mu^T
                                                                    + mu' 1_T'^T + 1_T' mu'^T             > 0
(each block with its own y).
SOUNDNESS: by the centroid identity for C_v there is T* in C_v of positive volume with
<sym(mu 1_T*^T), R> <= 0, and by the identity for C_q there is T'* in C_q of positive volume with
<sym(mu' 1_T'*^T), R> <= 0.  If T* is in Good, <M_T*, R> <= 0; if T* is in Bad, <M_T*T'*, R> <= 0.  Either
way a positive definite matrix has a non-positive inner product with R >= 0, R != 0: contradiction.
(Good = C_v, Bad = {} is the single-cone certificate; Good = {} is the full pair-product certificate.)

Compared with the two-cone certificate with a coupling matrix W ("certG"), this has no shared matrix
variable, so the SDP is far cheaper; the Bad set is chosen as the triangles deepest in the tree.
"""
import math

import numpy as np

import certify_joint as cj
import strengthen_joint as sj


def depths(parent, v):
    d = []
    for i in range(len(parent)):
        x, s = i, 0
        while x != v:
            x = parent[x]
            s += 1
        d.append(s)
    return d


def hybrid_blocks(V, v, q, faces, bad):
    """blocks for strengthen_joint.solve_multi (D = 2) and their (face_v, face_q) labels"""
    Cv = [tuple(f) for f in faces if v not in f]
    Cq = [tuple(f) for f in faces if q not in f]
    blocks, labels = [], []
    for f in Cv:
        if f in bad:
            for g in Cq:
                blocks.append([cj.tet_indicator(V, v, f), cj.tet_indicator(V, q, g)])
                labels.append((f, g))
        else:
            blocks.append([cj.tet_indicator(V, v, f), None])
            labels.append((f, None))
    return blocks, labels


def round_hybrid(V, v, parent, blocks, res, max_bits=46):
    """integer rounding + exact positive-definiteness check (Bareiss) of every block"""
    z = res["z"]
    if z <= 1e-9:
        return None
    start = max(4, int(math.ceil(math.log2(16.0 / z))))
    for bits in range(start, max_bits, 2):
        D = float(2 ** bits)
        mus = [[int(t) for t in np.rint(mu * D)] for mu in res["mus"]]
        Yi = np.maximum(np.rint(res["Y"] * D), 0.0)
        Yi[:, v] = 0
        ok = True
        for blk, y in zip(blocks, Yi):
            N = [[0] * V for _ in range(V)]
            for i in range(V):
                if i == v or y[i] == 0:
                    continue
                j = parent[i]
                yi = int(y[i])
                N[i][i] += 2 * yi
                N[i][j] -= yi
                N[j][i] -= yi
            for d, one in enumerate(blk):
                if one is None:
                    continue
                idx = [k for k in range(V) if one[k]]
                for a in range(V):
                    for b in idx:
                        N[a][b] += mus[d][a]
                        N[b][a] += mus[d][a]
            if not cj.bareiss_pd(N):
                ok = False
                break
        if ok:
            return mus, [[int(t) for t in y] for y in Yi], bits
    return None


def try_hybrid(V, v, parent, faces, q, bad):
    blocks, labels = hybrid_blocks(V, v, q, faces, bad)
    res = sj.solve_multi(V, v, parent, blocks, 2)
    rc = round_hybrid(V, v, parent, blocks, res) if res["z"] > 1e-7 else None
    return res["z"], rc, labels


def certify_tree(V, v, parent, faces, max_bad=3):
    """Search: second apex q (shallow vertices first) x Bad = the k deepest triangles (k = 1..max_bad).
    Returns (record dict, z) or (None, best z)."""
    d = depths(parent, v)
    Cv = [tuple(f) for f in faces if v not in f]
    ranked = sorted(Cv, key=lambda f: (-min(d[x] for x in f), -sum(d[x] for x in f), f))
    qs = sorted((u for u in range(V) if u != v), key=lambda u: (d[u], u))
    best = -1.0
    for k in range(1, max_bad + 1):
        bad = set(ranked[:k])
        for q in qs[:3]:
            z, rc, labels = try_hybrid(V, v, parent, faces, q, bad)
            best = max(best, z)
            if rc is not None:
                mus, Y, bits = rc
                single = [[list(f), y] for (f, g), y in zip(labels, Y) if g is None]
                pairs = [[list(f), list(g), y] for (f, g), y in zip(labels, Y) if g is not None]
                return dict(kind="certH", root=v, parent=list(parent), q=q, mu=mus[0], mu2=mus[1],
                            z=float(f"{z:.6g}"), single=single, pairs=pairs), z
    return None, best
