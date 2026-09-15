# Full system mathematics: Tasks 44 and 45

Continue the recovered Task 44 analysis without changing the production model.
All scripts load the frozen Task 40 production objects and verify their hashes.
All new results remain in this directory. The parameter audit lives on the
separate Task 43 branch.

## Recompute

From the repository root:

```sh
python analysis/44_full_system_mathematics/run_analysis.py --skip-inverse
python analysis/44_full_system_mathematics/verify_exact.py
python analysis/44_full_system_mathematics/equilibrium_review.py
python analysis/44_full_system_mathematics/sensitivity_analysis.py
python analysis/44_full_system_mathematics/inverse_analysis.py --phase all --recompute
python analysis/44_full_system_mathematics/inverse_verify.py
python analysis/44_full_system_mathematics/prepare_inverse_report.py
python analysis/44_full_system_mathematics/verify.py
python analysis/44_full_system_mathematics/verify_completion.py
python analysis/44_full_system_mathematics/prepare_report.py
bash analysis/44_full_system_mathematics/build.sh
```

The inverse sensitivity phase can instead be distributed over the eight named
parameters after the baseline and long phases. Run `--phase combine` after all
parameter jobs finish. For an independent replay always use `--recompute`.
`inverse_detail/trials/` is an optional local speed cache, excluded from version
control. The published top level inverse files retain the ordered scan, refined
crossings, gates, selected reproduction, long trajectories, equilibrium checks,
all local derivative validations and the independent quadrature verification.

## Interpretation

Exact balance identities do not imply zero total AE4 sensitivity. Constant
stimulus equilibria differ from the declared 600 s production experiment.
Mathematical roots outside the physiological domain and thermodynamically
opposed fluxes are retained explicitly. No model parameter was fitted.
The Task 41 comparison threshold is one reported standard error below the
experimental mean, not a confidence bound. Its construction remains target
selected. No global uncertainty interval is inferred from local difference
steps. See the final report and status for the complete scope.
