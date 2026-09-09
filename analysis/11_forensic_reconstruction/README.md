# Phase 11: full forensic reconstruction

## Outcome

The archived MATLAB snapshot is graded **`UNLIKELY`** to be the exact
implementation that generated the 2018 article. It is a strongly lineaged
calibration/reduction branch, but its equations, dependencies, calcium/time
protocol, and AE4 phenotype do not match the publication strongly enough to
open the clean-reconstruction gate.

The project classification remains **`STOP`**. No clean 2018 implementation,
physiological identifiability analysis, or manuscript was created.

## Reading order

1. `archive_audit.md` — file-level ledger for all 67 archive items and
   integrity evidence;
2. `code_lineage.md` — MATLAB dependency graph, state reductions, calibration
   lineage, and complete `Par.mat` interpretation;
3. `paper_code_concordance.md` — equation-by-equation and
   parameter-by-parameter comparison;
4. `figure_and_manuscript_provenance.md` — visual, TeX, bibliography, and
   chronology evidence;
5. `reproduction.md` — literal-runtime limits, exact calibration checks, and
   translated historical trajectories;
6. `blocker_resolution.md` — disposition and recomputation of every Task 10
   blocker;
7. `literature_crosscheck.md` — published-source evidence kept separate from
   historical implementation evidence;
8. `identity_evidence.md` — supplementary graded identity matrix and
   conditional-gate detail; and
9. `final_assessment.md` — controlling identity decision, conditional gates, and direct
   answers to the ten required questions.

## Reproduce the machine checks

From the repository root:

```bash
python analysis/11_forensic_reconstruction/tools/dump_par_mat.py --check
python analysis/11_forensic_reconstruction/tools/reproduce_historical_matlab.py --check
PYTHONPATH=src python -m unittest discover -s tests -p 'test*.py' -v
```

MATLAB and Octave were unavailable. Exact calibration and parameter-lineage
checks are direct evaluations of archived expressions; trajectories are
explicitly labelled SciPy-BDF translations of those expressions, not native
MATLAB execution.

Machine-readable evidence is in `results/11_forensic_reconstruction/`.
Temporary renders, public-paper downloads, and writable execution copies are
outside the committed analysis. Nothing under `archive/` was edited.
