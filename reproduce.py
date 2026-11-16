#!/usr/bin/env python3
"""
reproduce.py -- one command that re-checks all three theorems from the files in this repository.

    python reproduce.py                 # everything (about 10 minutes with 8+ cores)
    python reproduce.py --quick         # Theorem C, Theorem B for V <= 7, sphere-list completeness (seconds)
    python reproduce.py --workers 4     # number of processes for penta/verify.py (default: all cores)
    python reproduce.py --with-tests    # also run the tamper tests of the verifiers
    python reproduce.py --with-audit    # also run the independent Theorem B checker (audit/theoremB)

Setup from a clean checkout (Python 3.11 tested):
    python -m venv venv-verify
    venv-verify/bin/python -m pip install -r requirements-verify.txt      # Windows: venv-verify\Scripts\python
    venv-verify/bin/python reproduce.py

Steps:
  0. environment report (Python, platform, versions of sympy / python-flint);
  1. SHA-256 of the certificate files and of the trusted checking code against CHECKSUMS.sha256
     (text files are hashed with LF line endings, as stored in git, so the hashes do not depend on the
     platform's checkout settings);
  2. Theorem C: examples/V14_short_check.py and verify/exact_verify.py on the 14-vertex example;
  3. Theorem B: lb_joint/verify_joint.py for V = 4..8, and completeness of the sphere lists twice
     (exhaustive enumeration, lb_joint/enumerate_spheres.py, and Wagner's flip theorem,
     lb_joint/flip_closure.py);
  4. Theorem A: penta/verify.py, run with `python -O` (the verifier contains no `assert`).
Every step is a separate process; the script exits with status 0 iff every step exits with status 0.
A transcript is written to reproduce.log (or --log FILE).
"""
import argparse
import hashlib
import os
import platform
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
CHECKSUMS = os.path.join(ROOT, "CHECKSUMS.sha256")

TRUSTED_CODE = ["lb_joint/verify_joint.py", "lb_joint/enumerate_spheres.py", "lb_joint/flip_closure.py",
                "examples/V14_short_check.py", "verify/exact_verify.py",
                "penta/verify.py", "penta/core.py", "penta/shapes.py"]
CERTIFICATES = ["lb_joint/certs_V4.jsonl.gz", "lb_joint/certs_V5.jsonl.gz", "lb_joint/certs_V6.jsonl.gz",
                "lb_joint/certs_V7.jsonl.gz", "lb_joint/certs_V8.jsonl.gz", "lb_joint/exp1_V8_certG.jsonl.gz",
                "examples/V14_integer.json", "penta/certs_Z.json.gz", "penta/certs_P.json.gz"]


def sha256(rel):
    with open(os.path.join(ROOT, rel), "rb") as fh:
        data = fh.read()
    if not rel.endswith(".gz"):
        data = data.replace(b"\r\n", b"\n")
    return hashlib.sha256(data).hexdigest()


def steps(quick, workers):
    py = [sys.executable]
    out = [("Theorem C: short exact check of the 14-vertex example", py + ["examples/V14_short_check.py"]),
           ("Theorem C: full exact equilibrium classification", py + ["verify/exact_verify.py",
                                                                     "examples/V14_integer.json"])]
    for V in range(4, 9):
        if quick and V == 8:
            continue
        files = [f"lb_joint/certs_V{V}.jsonl.gz"] + (["lb_joint/exp1_V8_certG.jsonl.gz"] if V == 8 else [])
        out.append((f"Theorem B: certificates, V = {V}", py + ["lb_joint/verify_joint.py", str(V)] + files))
    for V in range(4, 9):
        if quick and V == 8:
            continue
        out.append((f"Theorem B: sphere list complete (enumeration), V = {V}",
                    py + ["lb_joint/enumerate_spheres.py", str(V), f"lb_joint/certs_V{V}.jsonl.gz"]))
        out.append((f"Theorem B: sphere list complete (flip closure), V = {V}",
                    py + ["lb_joint/flip_closure.py", str(V), f"lb_joint/certs_V{V}.jsonl.gz"]))
    if not quick:
        out.append(("Theorem A: rigorous cover of all pentahedra (python -O)",
                    py + ["-O", "penta/verify.py", "--workers", str(workers)]))
    return out


def tests():
    py = [sys.executable]
    return [("tamper tests of verify_joint.py", py + ["lb_joint/test_verifier.py"])]


def audit_steps(quick, workers):
    """the independent re-implementation of the Theorem B check (audit/theoremB, separate code)"""
    py = [sys.executable]
    out = []
    for V in range(4, 8 if quick else 9):
        files = [f"lb_joint/certs_V{V}.jsonl.gz"] + (["lb_joint/exp1_V8_certG.jsonl.gz"] if V == 8 else [])
        out.append((f"Theorem B, independent audit checker, V = {V}",
                    py + ["audit/theoremB/check_theoremB.py", str(V)] + files + [f"--procs={workers}", "--quiet"]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--with-tests", action="store_true")
    ap.add_argument("--with-audit", action="store_true", help="also run the independent Theorem B checker")
    ap.add_argument("--log", default=os.path.join(ROOT, "reproduce.log"))
    ap.add_argument("--write-checksums", action="store_true", help="(maintainers) regenerate CHECKSUMS.sha256")
    a = ap.parse_args()
    if a.write_checksums:
        with open(CHECKSUMS, "w", encoding="ascii", newline="\n") as fh:
            for rel in CERTIFICATES + TRUSTED_CODE:
                fh.write(f"{sha256(rel)}  {rel}\n")
        print(f"wrote {CHECKSUMS}")
        return 0

    log = open(a.log, "w", encoding="utf-8", newline="\n")

    def say(text=""):
        print(text, flush=True)
        log.write(text + "\n")
        log.flush()

    t_start = time.time()
    say(f"reproduce.py  {time.strftime('%Y-%m-%d %H:%M:%S')}  {'(quick)' if a.quick else '(full)'}")
    say(f"python {platform.python_version()} ({platform.python_implementation()}), {platform.platform()}, "
        f"{os.cpu_count()} logical CPUs")
    for mod in ("sympy", "flint"):
        try:
            m = __import__(mod)
            say(f"  {mod} {getattr(m, '__version__', '?')}")
        except ImportError:
            say(f"  {mod}: NOT INSTALLED" + ("" if a.quick else "  (needed for Theorem A)"))
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=ROOT,
                               capture_output=True, text=True).stdout.strip()
        say(f"git commit {commit or '?'}{' (with local modifications)' if dirty else ''}")
    except OSError:
        say("git: not available")

    results = []
    # 1. checksums
    say("\n== SHA-256 of certificate files and trusted code (CHECKSUMS.sha256)")
    ok_hash = True
    with open(CHECKSUMS, encoding="ascii") as fh:
        listed = [line.split() for line in fh if line.strip()]
    for digest, rel in listed:
        if not os.path.exists(os.path.join(ROOT, rel)):
            say(f"  MISSING  {rel}")
            ok_hash = False
            continue
        got = sha256(rel)
        say(f"  {'ok      ' if got == digest else 'MISMATCH'} {got[:16]}...  {rel}")
        ok_hash &= got == digest
    if {rel for _, rel in listed} != set(CERTIFICATES + TRUSTED_CODE):
        say("  the list of files in CHECKSUMS.sha256 is not the expected one")
        ok_hash = False
    results.append(("SHA-256 checksums", ok_hash, 0.0))

    # 2.-4. verifiers
    todo = (steps(a.quick, a.workers) + (tests() if a.with_tests else [])
            + (audit_steps(a.quick, a.workers) if a.with_audit else []))
    env = dict(os.environ, OMP_NUM_THREADS="1", PYTHONIOENCODING="utf-8")
    for title, cmd in todo:
        say(f"\n== {title}\n$ {' '.join(['python'] + cmd[1:])}")
        t0 = time.time()
        p = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True, encoding="utf-8",
                           errors="replace")
        dt = time.time() - t0
        lines = (p.stdout + p.stderr).rstrip().splitlines()
        for line in lines[-12:] if len(lines) > 12 else lines:
            say("  | " + line)
        say(f"  -> exit status {p.returncode}  [{dt:.1f}s]")
        results.append((title, p.returncode == 0, dt))

    say("\n== SUMMARY")
    for title, ok, dt in results:
        say(f"  {'PASS' if ok else 'FAIL'}  {title}  [{dt:.0f}s]")
    all_ok = all(ok for _, ok, _ in results)
    say(f"\n{'ALL STEPS PASSED' if all_ok else 'SOME STEPS FAILED'}"
        f"{' (quick mode: V = 8 and Theorem A skipped)' if a.quick else ''}; total {time.time() - t_start:.0f}s")
    log.close()
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
