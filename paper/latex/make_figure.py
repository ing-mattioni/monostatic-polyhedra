"""make_figure.py -- top view of the 14-vertex mono-unstable polyhedron of Theorem C (figure for main.tex).

    python paper/latex/make_figure.py      ->  paper/latex/V14_plate.pdf

The outline, the shadowing tree i -> p(i) and the centre of mass are computed exactly from
examples/V14_integer.json and the witness table of examples/V14_short_check.py; the drawing is matplotlib.
"""
import json
import os
import sys
from fractions import Fraction

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "examples"))
import V14_short_check as chk  # noqa: E402  (data only: POINTS, TRIANGLES, WITNESS, UNSTABLE)

INK, MUTED, FILL, ACCENT = "#1a1a1a", "#6b6a66", "#efeeea", "#b3261e"

pts = [tuple(int(t) for t in p) for p in json.load(open(os.path.join(ROOT, "examples", "V14_integer.json")))["vertices"]]
assert pts == chk.POINTS
D, C = 0, [0, 0, 0]
for a, b, c in chk.TRIANGLES:
    A, B, Cc = pts[a], pts[b], pts[c]
    d = (A[0] * (B[1] * Cc[2] - B[2] * Cc[1]) - A[1] * (B[0] * Cc[2] - B[2] * Cc[0])
         + A[2] * (B[0] * Cc[1] - B[1] * Cc[0]))
    D += d
    C = [C[m] + d * (A[m] + B[m] + Cc[m]) for m in range(3)]
cx, cy = (float(Fraction(C[m], 4 * D)) for m in range(2))

fig, ax = plt.subplots(figsize=(3.9, 3.9))
xy = [(p[0], p[1]) for p in pts]
ax.fill([q[0] for q in xy], [q[1] for q in xy], color=FILL, zorder=0)
ax.plot([q[0] for q in xy] + [xy[0][0]], [q[1] for q in xy] + [xy[0][1]], color=MUTED, lw=0.9, zorder=1)
for i, j in chk.WITNESS.items():
    (x0, y0), (x1, y1) = xy[i], xy[j]
    s0, s1 = 0.12, 0.12
    ax.annotate("", xy=(x1 - s1 * (x1 - x0), y1 - s1 * (y1 - y0)), xytext=(x0 + s0 * (x1 - x0), y0 + s0 * (y1 - y0)),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.8, mutation_scale=7), zorder=2)
for k, (x, y) in enumerate(xy):
    if k == chk.UNSTABLE:
        ax.plot(x, y, marker="o", ms=6, color=ACCENT, zorder=3)
    else:
        ax.plot(x, y, marker="o", ms=3.2, color=INK, zorder=3)
    dx, dy = x - cx, y - cy
    n = (dx * dx + dy * dy) ** 0.5
    ax.text(x + 75 * dx / n, y + 75 * dy / n, f"{k}\n{pts[k][2]:+d}".replace("-", "−"),
            ha="center", va="center", fontsize=7.5, color=INK if k else ACCENT, linespacing=0.9)
ax.plot(cx, cy, marker="+", ms=8, mew=1.1, color=INK, zorder=3)
ax.text(cx + 25, cy - 45, r"$c$", fontsize=8, color=INK)
ax.set_aspect("equal")
ax.set_xlim(-760, 780)
ax.set_ylim(-600, 1110)
ax.axis("off")
fig.tight_layout(pad=0.1)
out = os.path.join(HERE, "V14_plate.pdf")
fig.savefig(out)
print("wrote", out, " centre of mass (x, y) =", (round(cx, 3), round(cy, 3)))
