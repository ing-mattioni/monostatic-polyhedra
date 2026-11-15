import numpy as np, sys
sys.path.insert(0, 'explore')
from mu3d import f, info
from flat_pyramid_2d import polygon, area_centroid

def flat_pyramid(dphi_deg, shrink=0.0, eps=0.02, lower=0.5):
    """Build 3D flat pyramid from symmetric spiral angles. Spiral radii r_i = r_{i-1} cos(dphi)(1-shrink).
    Apex placed at A' = on axis between -3g and bottom (fraction 'lower' of the slack)."""
    d = np.radians(dphi_deg)
    phis = np.cumsum(d); r = np.cumprod(np.cos(d)*(1-shrink))
    xs = r*np.sin(phis); zs = r*np.cos(phis)
    P = np.array([(0.0, 1.0)] + list(zip(xs, zs)) + list(zip(-xs[::-1], zs[::-1])))
    A, cx, cz = area_centroid(P)
    m = len(d); zbot = P[m,1]
    # want reference c' = (3g + A')/4 below origin slightly: choose A' between bottom and -3g
    zA_needed = -3*cz           # gives c' exactly at origin
    zA = zbot + (1-lower)*(zA_needed - zbot)   # lower => c' below origin
    X = np.array([(x, 0.0, z) for x, z in P] + [(0.0, eps, zA)])
    return X

ang8 = [31.49,19.14,17.19,16.56,16.59,17.13,18.26,20.18]
ang7 = [34.85,20.71,18.86,18.47,18.88,20.05,22.21]
for name, ang in (("m=8 V=18", ang8), ("m=7 V=16", ang7)):
    for eps in (0.05, 0.01, 0.002):
        for lower in (0.3, 0.6, 0.9):
            X = flat_pyramid(ang, eps=eps, lower=lower)
            print(name, "eps", eps, "lower", lower, "f=%.3e" % f(X), info(X)["nvert"], info(X)["unstable"])
np.save('explore/seed_V16.npy', flat_pyramid(ang7, eps=0.01, lower=0.6))
np.save('explore/seed_V18.npy', flat_pyramid(ang8, eps=0.01, lower=0.6))
