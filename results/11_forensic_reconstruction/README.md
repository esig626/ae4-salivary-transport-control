# Phase 11 machine-readable evidence

| File | Role |
| --- | --- |
| `archive_hashes.csv` | 67-file integrity ledger; two connector transport exceptions use authoritative remote Git blob evidence |
| `figure_mapping.csv` | Historical figure source, build, manuscript, chronology, and published-figure comparison map |
| `par_mat_inventory.json` | Complete top-level and 45-field `Par.mat` value/shape/type dump with exact `Par.m` comparison |
| `reproduction_summary.json` | Runtime limits, calibration values, parameter lineage, solver metadata, endpoints, residuals, and phenotype ratios |
| `reproduction_residuals.csv` | Initial RHS values before and after each historical global time multiplier |
| `reproduction_trajectories.csv` | Six translated archived-equation trajectories with state, flow, voltage, and flux diagnostics |

Regenerate or verify these outputs using the commands in
`analysis/11_forensic_reconstruction/README.md`. The trajectories are forensic
SciPy-BDF translations because MATLAB/Octave was unavailable; they are not a
new scientific model or a native-execution claim.
