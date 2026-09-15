#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python prepare_report.py
command -v pdflatex >/dev/null || { echo 'Missing pdflatex' >&2; exit 1; }
pdflatex -interaction=nonstopmode -halt-on-error report.tex > build.log
if command -v bibtex >/dev/null 2>&1; then
  bibtex report >> build.log
elif command -v bibtex8 >/dev/null 2>&1; then
  bibtex8 --8bit report >> build.log
else
  echo 'Missing BibTeX' >&2; exit 1
fi
for pass in 1 2 3; do pdflatex -interaction=nonstopmode -halt-on-error report.tex >> build.log; done
if grep -Eq 'undefined citations|undefined references|Label\(s\) may have changed|Rerun to get cross-references right' report.log; then
  echo 'Unresolved report references' >&2; exit 1
fi
cp report.pdf AE4_parameter_provenance.pdf
pdfinfo AE4_parameter_provenance.pdf | grep -E '^(Pages|Page size|File size):'
