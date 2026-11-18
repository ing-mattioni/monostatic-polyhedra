#!/usr/bin/env python3
"""Build the curated C-only deposit and verify both candidate and final extractions."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
SLUG = "mono-unstable-14-v1.0"
TEXT_SUFFIXES = {".py", ".md", ".tex", ".json", ".cff", ".txt", ".svg", ".log"}


def need(ok, message):
    if not ok:
        raise ValueError(message)


def inventory():
    names = [s.strip() for s in (ROOT/"release/files.txt").read_text().splitlines()
             if s.strip() and not s.startswith("#")]
    need(len(names) == len(set(names)), "Repeated allowlist entry")
    for rel in names:
        path = ROOT/rel
        need(path.resolve().is_relative_to(ROOT) and not path.is_symlink(), "Unsafe path: "+rel)
        need(path.is_file(), "Missing file: "+rel)
        need(not any(p in {"outreach", ".git", "history", "tmp", "output"} for p in Path(rel).parts),
             "Forbidden package path: "+rel)
    return sorted(names)


def normalized(path):
    content = path.read_bytes()
    return content.replace(b"\r\n", b"\n") if path.suffix in TEXT_SUFFIXES else content


def prepare():
    transcript = ROOT/"release/reproduce_extracted.log"
    if not transcript.exists():
        transcript.write_text("Pending candidate extraction; replaced by the build command.\n", encoding="utf-8")
    rows = [hashlib.sha256(normalized(ROOT/rel)).hexdigest()+"  "+rel for rel in inventory()]
    (ROOT/"CHECKSUMS_C.sha256").write_text(
        "# C-only version 1.0; SHA-256, LF-normalized text, raw binary files.\n"+"\n".join(rows)+"\n",
        encoding="utf-8", newline="\n")
    print("Prepared manifest for", len(rows), "files.")


def pack(path):
    names = inventory()+["CHECKSUMS_C.sha256"]
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for rel in sorted(names):
            info = zipfile.ZipInfo(SLUG+"/"+rel, date_time=(2026,11,17,0,0,0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            archive.writestr(info, normalized(ROOT/rel))
    return len(names)


def check_extraction(archive_path, work, log):
    expected = {SLUG+"/"+rel for rel in inventory()+["CHECKSUMS_C.sha256"]}
    with zipfile.ZipFile(archive_path) as archive:
        need(set(archive.namelist()) == expected and len(archive.namelist()) == len(expected),
             "Archive inventory differs")
        for entry in archive.namelist():
            need((work/entry).resolve().is_relative_to(work.resolve()), "Unsafe ZIP entry")
        need(archive.testzip() is None, "Corrupted ZIP member")
        archive.extractall(work)
    extracted = work/SLUG
    proc = subprocess.run([sys.executable, "-I", "-S", "-B", "reproduce.py", "--log", str(log.resolve())],
                          cwd=extracted, capture_output=True, text=True, encoding="utf-8")
    need(proc.returncode == 0, "Extracted verification failed:\n"+proc.stdout+proc.stderr)
    need("SUCCESS: 8/8" in proc.stdout, "Incomplete verification transcript")
    need((extracted/"paper/latex/main.pdf").read_bytes() == (ROOT/"paper/latex/main.pdf").read_bytes(),
         "Embedded PDF differs")
    return proc.stdout


def build():
    output = ROOT/"output/zenodo"
    output.mkdir(parents=True, exist_ok=True)
    scratch = ROOT/"tmp"
    scratch.mkdir(exist_ok=True)
    prepare()
    with tempfile.TemporaryDirectory(prefix="c-release-", dir=scratch) as dirname:
        work = Path(dirname).resolve()
        need(work.is_relative_to(scratch.resolve()), "Temporary directory escaped project")
        candidate = work/"candidate.zip"
        pack(candidate)
        candidate_log = work/"candidate.log"
        check_extraction(candidate, work/"candidate", candidate_log)
        shutil.copyfile(candidate_log, ROOT/"release/reproduce_extracted.log")
        prepare()
        final_zip = output/(SLUG+".zip")
        count = pack(final_zip)
        final_log = output/"verification_final.log"
        stdout = check_extraction(final_zip, work/"final", final_log)
        pdf = output/"mono-unstable-14.pdf"
        shutil.copyfile(ROOT/"paper/latex/main.pdf", pdf)
        hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (final_zip, pdf)}
        (output/"SHA256SUMS.txt").write_text(
            "".join(digest+"  "+name+"\n" for name,digest in hashes.items()), encoding="utf-8", newline="\n")
        report = {"version":"1.0", "archive":final_zip.name, "archive_bytes":final_zip.stat().st_size,
                  "archive_files":count, "pdf_bytes":pdf.stat().st_size, "sha256":hashes,
                  "candidate_extraction":"8/8 checks passed", "final_extraction":"8/8 checks passed",
                  "verification_summary":stdout.strip().splitlines()[-1]}
        (output/"BUILD_REPORT.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
        print(json.dumps(report,indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--prepare", action="store_true")
    group.add_argument("--build", action="store_true")
    args = parser.parse_args()
    (prepare if args.prepare else build)()
