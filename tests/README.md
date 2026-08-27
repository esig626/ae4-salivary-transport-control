# Tests

Testing should be built alongside the model, not added at the end.

Minimum planned test classes

- published baseline regression tests
- steady state residual tests
- mass balance and stoichiometry checks
- unit and sign convention checks
- finite difference versus implicit sensitivity checks
- continuation reproducibility tests
- mechanism variant limiting case tests

Scientific claims depending on failed tests must not enter the results ledger.

The current Phase 10 tests certify the exact two-flux algebra, the regular
implicit-sensitivity identity on an independent nonlinear system, and the
material inconsistencies found by direct evaluation of the printed baseline.
They do not certify a full AE4-model equilibrium.

The Phase 11 forensic tests additionally certify source/copy hashes, the
constructed historical calibration, the calibration-to-`Par.mat` lineage,
historical WT residuals, deterministic translated BDF outputs, and the reported
nonmatching AE2/AE4 flow ratios. They do not claim native MATLAB execution or a
validated published-2018 physiological model.

Run them with

```bash
PYTHONPATH=src python -m unittest discover -s tests -p 'test*.py' -v
```
