# Building and verifying this deposit

## Verify the mathematics

Run python reproduce.py with Python 3.10+ (3.11.9 tested).
No third-party packages or network access are used. The command also works
with python -I -S reproduce.py.

## Regenerate presentation artifacts (optional)

From the archive root:

```sh
python paper/latex/make_ledger.py
python paper/latex/make_figure.py
tectonic -X compile paper/latex/main.tex --keep-logs
python release/make_markdown.py --pandoc pandoc
```

Only make_figure.py needs matplotlib (3.11.2 used during preparation).
The supplied figure is already compiled, so neither verification nor compiling
the manuscript requires regenerating it. Tectonic 0.17.0 and Pandoc 3.14 were
used. Tectonic may fetch fonts/packages on first use. Standard pdfLaTeX-compatible
packages are used, but a separate pdfLaTeX build has not been tested.

The PDF was visually checked on every page. Presentation-source or PDF changes
invalidate the corresponding manifest entries: regenerate the release
intentionally, rather than suppressing an integrity failure.

## Build the curated archive

release/files.txt is the explicit allowlist. The builder never scans the full
repository or includes local outreach files.

```sh
python release/build_archive.py --prepare
python reproduce.py
python release/build_archive.py --build
```

The first command regenerates CHECKSUMS_C.sha256 for the allowlisted files.
The build command creates a candidate ZIP, extracts it to a fresh temporary
directory, runs python -I -S reproduce.py there and stores the transcript
as release/reproduce_extracted.log. It then regenerates the manifest and
creates the final ZIP. Finally it extracts and verifies the final ZIP again,
checks the inventory, and supplies a separate copy of the identical PDF and
raw SHA-256 hashes for upload.

The transcript is included in the manifest. Verification logs from subsequent
runs are written outside the frozen archive. ZIP members use fixed timestamps
and ordering; byte-for-byte repeatability assumes identical inputs, including
the recorded execution transcript and compiled PDF. PDF rebuilds can differ
because TeX includes creation metadata.

The two upload files are output/zenodo/mono-unstable-14-v1.0.zip and
output/zenodo/mono-unstable-14.pdf. Metadata-entry instructions and checksum
files sit beside them. Upload the two research files; the form instructions
are for the depositor.
