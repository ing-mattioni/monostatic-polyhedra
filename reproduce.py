#!/usr/bin/env python3
"""Reproduce the fourteen-vertex paper using Python 3.10+ and its standard library.

    python reproduce.py
    python -I -S reproduce.py --log reproduce.log

The historical A/B/C workflow is preserved separately as reproduce_abc.py.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent
REQUIRED = {
    "reproduce.py", "examples/V14_integer.json", "examples/V14_short_check.py",
    "verify/exact_verify.py", "verify/check_c.py", "verify/check_manuscript.py",
    "examples/V14_exact_report.json", "examples/V14_family_report.json",
    "paper/latex/main.tex", "paper/latex/equilibrium_ledger.tex", "paper/latex/make_ledger.py",
}


def need(condition, message):
    if not condition:
        raise ValueError(message)


def integrity():
    seen = set()
    for line in (ROOT/"CHECKSUMS_C.sha256").read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        digest, rel = line.split("  ", 1)
        path = (ROOT/rel).resolve()
        need(path.is_relative_to(ROOT) and rel not in seen, "Invalid manifest path: "+rel)
        seen.add(rel)
        data = path.read_bytes()
        # Match Git LF text on every platform. PDF and other binary files are raw.
        if path.suffix.lower() in {".py", ".md", ".tex", ".json", ".cff", ".txt", ".svg", ".log"}:
            data = data.replace(b"\r\n", b"\n")
        need(hashlib.sha256(data).hexdigest() == digest, "Hash mismatch: "+rel)
    need(REQUIRED <= seen, "Manifest omits required verification inputs")
    return len(seen)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", type=Path, default=ROOT/"reproduce.log")
    args = parser.parse_args()
    log = args.log.resolve()
    # Never overwrite an input or the supplied frozen transcript.
    need(log.suffix == ".log" and log.name not in {"reproduce_extracted.log", "validation.log"},
         "Choose a new .log filename, not an archived transcript")
    lines = []
    def emit(text):
        print(text)
        lines.append(text)
    started = time.monotonic()
    passed = 0
    try:
        emit("Fourteen-vertex construction: version 1.0")
        emit("Python "+platform.python_version()+"; "+platform.system()+"; standard-library verification")
        emit("PASS 1/8: integrity of "+str(integrity())+" archived files (LF-normalized text)")
        passed += 1
        with tempfile.TemporaryDirectory(prefix="mono14-") as temp:
            temp = Path(temp)
            commands = [
                ("short integer existence proof", ["examples/V14_short_check.py"]),
                ("rational hull and complete equilibrium classification",
                 ["verify/exact_verify.py", "examples/V14_integer.json", "--report", str(temp/"exact.json")]),
                ("independent integer hull, equilibrium census and contraction family",
                 ["verify/check_c.py", "--report", str(temp/"family.json")]),
            ]
            for label, command in commands:
                proc = subprocess.run([sys.executable, "-I", "-S", "-B", *command], cwd=ROOT,
                                      capture_output=True, text=True, encoding="utf-8")
                # Keep transcripts portable and free of temporary absolute paths.
                output = (proc.stdout+proc.stderr).replace(str(temp), "<temporary>")
                emit(output.rstrip())
                need(proc.returncode == 0, label+" failed")
                passed += 1
                emit(f"PASS {passed}/8: {label}")
            for generated, stored in [("exact.json", "V14_exact_report.json"), ("family.json", "V14_family_report.json")]:
                actual = json.loads((temp/generated).read_text(encoding="utf-8"))
                expected = json.loads((ROOT/"examples"/stored).read_text(encoding="utf-8"))
                need(actual == expected, "Regenerated report differs: "+stored)
                passed += 1
                emit(f"PASS {passed}/8: regenerated {stored} agrees")
            for mode in ("tables", "ledger"):
                proc = subprocess.run([sys.executable, "-I", "-S", "-B", "verify/check_manuscript.py", mode],
                                      cwd=ROOT, capture_output=True, text=True, encoding="utf-8")
                emit((proc.stdout+proc.stderr).rstrip())
                need(proc.returncode == 0, "Manuscript "+mode+" mismatch")
                passed += 1
                emit(f"PASS {passed}/8: manuscript {mode}")
        emit(f"SUCCESS: {passed}/8 checks in {time.monotonic()-started:.3f} seconds.")
        return 0
    except (OSError, ValueError) as exc:
        emit(f"FAILED after {passed}/8 checks: {exc}")
        return 1
    finally:
        log.write_text("\n".join(lines)+"\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
