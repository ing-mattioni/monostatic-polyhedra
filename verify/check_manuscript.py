#!/usr/bin/env python3
"""Compare the paper's numerical tables and generated sign ledger with exact data."""
import argparse
import importlib.util
import json
from pathlib import Path
import re
import runpy
import sys

ROOT = Path(__file__).resolve().parents[1]


def need(ok, msg):
    if not ok:
        raise ValueError(msg)


def table_rows(tex, label):
    section = tex.split(r"\label{"+label+"}", 1)[1].split(r"\end{tabular}", 1)[0]
    rows = []
    for line in section.splitlines():
        fields = [v.strip() for v in line.replace(r"\\", "").split("&")]
        if fields and all(re.fullmatch(r"-?\d+", v) for v in fields):
            rows.append([int(v) for v in fields])
    return rows


def check_tables():
    short = runpy.run_path(str(ROOT/"examples/V14_short_check.py"))
    points = [tuple(map(int, p)) for p in json.loads((ROOT/"examples/V14_integer.json").read_text())["vertices"]]
    need(points == short["POINTS"], "JSON and short checker coordinates differ")
    report = json.loads((ROOT/"examples/V14_family_report.json").read_text())
    tex = (ROOT/"paper/latex/main.tex").read_text(encoding="utf-8")
    rows = table_rows(tex, "tab:points")
    coordinates = {}
    for row in rows:
        need(len(row) == 8, "Coordinate row length")
        for group in (row[:4], row[4:]):
            need(group[0] not in coordinates, "Repeated coordinate row")
            coordinates[group[0]] = tuple(group[1:])
    need(coordinates == dict(enumerate(points)), "Printed coordinate table differs")
    expected = [[int(r[k]) for k in ("i", "j", "A", "B", "A_plus_B")] for r in report["family_witnesses"]]
    need(table_rows(tex, "tab:witnesses") == expected, "Printed family/witness table differs")
    need({r[0]: r[1] for r in expected} == short["WITNESS"], "Witness maps differ")
    facet_section = tex.split(r"\label{eq:facets}", 1)[1].split(r"\end{equation}", 1)[0]
    facets = [tuple(map(int, f)) for f in re.findall(r"\((\d+),(\d+),(\d+)\)", facet_section)]
    need(facets == short["TRIANGLES"], "Printed facet list differs")
    canonical = lambda face: min(tuple(face[i:]+face[:i]) for i in range(3))
    need(sorted(canonical(list(f)) for f in facets) == sorted(canonical(f) for f in report["faces"]),
         "Printed faces differ from independently rebuilt hull")
    # Check the displayed centroid and volume numerators in addition to the tables.
    for value in ("135509395", "51679332", "87837359", "90438831", "1717399", "25839666"):
        need(value in tex, "Missing exact centroid value: "+value)
    print("Coordinates, 24 oriented facets, 13 coefficient rows and witnesses agree.")


def check_ledger():
    spec = importlib.util.spec_from_file_location("c_ledger", ROOT/"paper/latex/make_ledger.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    stored = (ROOT/"paper/latex/equilibrium_ledger.tex").read_text(encoding="utf-8")
    need(stored == module.ledger(), "Regenerated sign ledger differs")
    print("All 24 facet and 36 edge sign rows agree with regenerated integer predicates.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["tables", "ledger"])
    args = parser.parse_args()
    try:
        (check_tables if args.mode == "tables" else check_ledger)()
    except (ValueError, OSError) as exc:
        print("FAILED:", exc, file=sys.stderr)
        raise SystemExit(1)
