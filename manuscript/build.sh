#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python manuscript/scripts/build_figures.py
cd manuscript
pdflatex -interaction=nonstopmode -halt-on-error main.tex
if command -v bibtex >/dev/null 2>&1 && bibtex --version >/dev/null 2>&1; then
  bibtex main
elif command -v bibtex.original >/dev/null 2>&1; then
  bibtex.original main
else
  echo "BibTeX is required but was not found" >&2
  exit 1
fi
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
