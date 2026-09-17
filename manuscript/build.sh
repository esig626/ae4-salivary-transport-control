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
pdflatex -interaction=nonstopmode -halt-on-error main.tex

if grep -Eq "undefined citations|undefined references|Label\(s\) may have changed|Rerun to get cross-references right" main.log; then
  echo "Manuscript build did not converge cleanly" >&2
  grep -E "undefined citations|undefined references|Label\(s\) may have changed|Rerun to get cross-references right" main.log >&2 || true
  exit 1
fi

cp main.pdf AE4_manuscript_current.pdf
python - <<'PYFLAT'
from pathlib import Path
s = Path('main.tex').read_text()
for name in ('model', 'results', 'discussion'):
    s = s.replace('\\input{' + name + '}', Path(name + '.tex').read_text())
Path('AE4_manuscript_standalone.tex').write_text(s)
PYFLAT
