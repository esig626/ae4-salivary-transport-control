#!/usr/bin/env bash
set -euo pipefail
report_dir="$(cd "$(dirname "$0")" && pwd)"
cd "$report_dir"
export SOURCE_DATE_EPOCH="$(git show -s --format=%ct f783785a863440df459e9c5530beba4c7eb59110)"
export FORCE_SOURCE_DATE=1
report_build_dir="$(mktemp -d)"
trap 'rm -rf "$report_build_dir"' EXIT
pdflatex -interaction=nonstopmode -halt-on-error -output-directory="$report_build_dir" report.tex > "$report_build_dir/pass1.txt"
(cd "$report_build_dir" && BIBINPUTS="$report_dir:" bibtex report > bibtex.txt)
pdflatex -interaction=nonstopmode -halt-on-error -output-directory="$report_build_dir" report.tex > "$report_build_dir/pass2.txt"
pdflatex -interaction=nonstopmode -halt-on-error -output-directory="$report_build_dir" report.tex > "$report_build_dir/pass3.txt"
python - "$report_build_dir/report.log" <<'PY'
from pathlib import Path
import sys
source = Path(sys.argv[1])
log = source.read_text().replace(str(source.parent), 'BUILD_DIRECTORY')
Path('build.log').write_text('\n'.join(line.rstrip() for line in log.splitlines()).rstrip()+'\n')
PY
if rg -n 'undefined|multiply defined|Overfull' build.log > "$report_build_dir/problems.txt"; then
  cat "$report_build_dir/problems.txt"
  exit 1
fi
cp "$report_build_dir/report.pdf" AE4_full_system_mathematics.pdf
pdfinfo AE4_full_system_mathematics.pdf | sed -n '/Pages:/p;/Page size:/p'
