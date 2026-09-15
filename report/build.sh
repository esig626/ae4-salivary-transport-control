#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 assemble.py "$@"
mkdir -p build
export SOURCE_DATE_EPOCH=1789430400
export FORCE_SOURCE_DATE=1
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build report.tex > build/console.log 2>&1 || {
  tail -n 100 build/console.log >&2
  exit 1
}
cp build/report.pdf AE4_mathematical_analysis_report.pdf
python3 verify_report.py
