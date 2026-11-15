"""Diagnose a near-boundary (G -> 0-) solution: for each stable face other than i*, report its margin,
the distance of the foot q_i to the nearest vertex / edge, and the smallest dihedral bend
(angle between normals) to an adjacent face -> detects 'keel pairs' (nearly coplanar split faces)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from core import evaluate_p

for f in sys.argv[1:]:
    P = np.load(f)
    r = evaluate_p(P)
    V, faces, s = r["vertices"], r["faces"], r["s"]
    sc = r["vol"] ** (1 / 3)
    U = P / np.linalg.norm(P, axis=1)[:, None]
    print(f"{os.path.basename(f)}: F={r['nfaces']} G={r['G']:.3e} nstable={r['nstable']} istar={r['istar']}")
    for i in np.argsort(-s)[:r["nstable"] + 1]:
        if not r["present"][i]:
            continue
        h = (P[i] ** 2).sum() / np.linalg.norm(P[i])
        q = r["c"] - (U[i] @ r["c"] - h) * U[i]
        dv = min(np.linalg.norm(V[v] - q) for v in faces[i]) / sc
        nb = [j for j in range(len(P)) if j != i and r["present"][j] and len(set(faces[i]) & set(faces[j])) >= 2]
        bend = min(np.degrees(np.arccos(np.clip(U[i] @ U[j], -1, 1))) for j in nb)
        jb = min(nb, key=lambda j: -U[i] @ U[j])
        print(f"   face {i:2d} s={s[i]:+.3e} foot-to-nearest-vertex={dv:.3e} min bend to neighbour={bend:6.2f} deg "
              f"(face {jb}, s={s[jb]:+.2e}) area={r['area'][i]:.3f}")
