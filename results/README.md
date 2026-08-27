# Results

Machine readable outputs from completed analyses belong here. Results should be regenerable from version controlled code and documented inputs.

Suggested naming is `phase##_short_description.ext`.

Temporary outputs belong in `results/tmp/` and are ignored by Git.

`10_identifiability_discrimination/` contains the deterministic Phase 10
printed-equation audit and certified stoichiometric geometry. Its local README
defines the scope and prevents these outputs from being mistaken for a
reproduced physiological simulation.

The directory also retains outputs from the affine software fixture introduced
by the preceding branch commit. Files using `Q_star`, diagnostic restoring
matrices, or proxy panel metrics are implementation checks only. The local
Phase 10 README distinguishes them from the certified aggregate map and
printed-equation audit.
