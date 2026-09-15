# Task 46G local solver verification

The required local solver verification passes. This checkpoint establishes the
computational prerequisite for the Task 46 CBM screen. It does not establish an
AE4 transport network, a phenotype result or an accepted replacement model.

## Authority and isolation

The user invoked `prompts/46g_use_local_cbm_core_and_isolate_carbonscope.md`,
read at AE4 main commit `1148397ca7add328b1e3da71864c56e5382ae324`, with
prompt blob `7ba7922cd63fbe4fd6f461a8e3a104866679e9ab`. The existing Task 46
branch started at `3598bd866c33d41d4a6914bda56d17c73900e82a`.

All solver imports come from `vendor/carbonscope_flux_core/` in this Task 46
directory. Its README, provenance and third party notices were read first.
The adapted local port cites CarbonScope source commit
`a11e25f176a2cc70de22c3f36a1c0e623248393a`; this is source provenance, not
a runtime dependency or a claim that the port is an identical copy.

No CarbonScope repository was cloned, read, modified, branched, committed,
pushed or otherwise accessed during this work. No CarbonScope pull request
was opened. The script checks every loaded core module's resolved path and
rejects an imported `fluxemu` package. It records all core file SHA256 hashes
and verifies that they remain unchanged. No generic solver extension was needed.

## Runtime and smoke test

The declared requirements were installed with:

```bash
python -m pip install -r analysis/46_physiology_constrained_model_reconstruction/vendor/carbonscope_flux_core/requirements.txt
```

The environment used Python 3.12.14, NumPy 2.3.5, pandas 2.2.3 and highspy
1.12.0. NumPy and pandas were already present; only highspy was installed.
The supplied smoke test passed both as a module and in the recorded verification:

```bash
PYTHONDONTWRITEBYTECODE=1 python -m analysis.46_physiology_constrained_model_reconstruction.vendor.carbonscope_flux_core.smoke_test
```

## Independent analytical control

This control was derived for Task 46G. It does not reuse the smoke test's single
pool chain or any AE4 model. Two abstract pools A and B have separate inputs,
separate drains and a reversible transfer T from A to B. Units are arbitrary.
The reaction order is `(U_A, U_B, T, D_A, D_B)` and balance order is `(A, B)`.

$$
S=\begin{pmatrix}1&0&-1&-1&0\\0&1&1&0&-1\end{pmatrix},\qquad Sv=0.
$$

| Flux | Lower bound | Upper bound | Interpretation |
| :--- | ---: | ---: | :--- |
| U_A | 0 | 6 | Input to A |
| U_B | 0 | 4 | Input to B |
| T | −2 | 3 | Transfer from A to B; negative values reverse it |
| D_A | 0 | 5 | Drain from A |
| D_B | 0 | 7 | Drain from B |

Maximise `z = D_A + D_B`. Adding the balances gives `z = U_A + U_B <= 10`.
Equality is attained at `(6, 4, 1, 5, 5)`. All maximisers have the form
`(6, 4, t, 6-t, 4+t)` with `1 <= t <= 3`; the optimum deliberately has
several feasible solutions. A second objective, minimising T on the same
feasible region, has optimum −2, attained at `(0, 2, −2, 2, 0)`.

The exact FVA ranges are:

| Flux | No objective constraint | Retain z >= 8 | Retain z = 10 |
| :--- | :--- | :--- | :--- |
| U_A | [0, 6] | [4, 6] | [6, 6] |
| U_B | [0, 4] | [2, 4] | [4, 4] |
| T | [−2, 3] | [−1, 3] | [1, 3] |
| D_A | [0, 5] | [1, 5] | [3, 5] |
| D_B | [0, 7] | [3, 7] | [5, 7] |

For retained level q of 8 or 10, the balances and bounds give `U_A >= q-4`,
`U_B >= q-6`, `T = U_A-D_A >= q-9`, `D_A = U_A-T >= q-7`, and
`D_B = z-D_A >= q-5`. Combined with the original bounds these give the table.
For the unconstrained region all original bounds are attainable. The script
supplies an independently specified complete feasible witness for each of the
30 endpoints. All witness entries are integers, so exact substitution proves
attainability without relying on an LP solver.

Each region is solved with reference FVA, VFFVA using two workers and chunk size
one, and prepared VFFVA using one worker and chunk size two. Both fast paths
are compared with the table and with reference FVA. Both FBA states are checked
against the balances, bounds and objectives. The compiled ordering, matrix,
bounds and objective are checked against independent arrays.

FVA endpoints remain marginal diagnostics. Combining all unconstrained maxima
gives `(6, 4, 3, 5, 7)`, with balance residual `(−2, 0)`. The script verifies
this incompatibility. Complete states support all attainability claims.

## Recorded result and reproduction

`output/checkpoint_46g_local_core_verification.json` records one supplied smoke
test, two correct analytical FBA objectives, three FVA regions, 30 independently
specified witnesses, 90 comparisons with analytical endpoints, and 60 comparisons
between reference and fast paths. There were zero failures and zero absolute
endpoint error in this run.

Analytical checks use absolute tolerance `1e-8` and zero relative tolerance.
The unchanged core uses HiGHS primal and dual tolerance `1e-9`, primal validation
tolerance `1e-7` and objective validation tolerance `1e-7`. The record includes
versions, file hashes, imports, full FBA states, diagnostics, witnesses and ranges.

Reproduce with a new output filename; the script refuses to replace an existing
record:

```bash
PYTHONDONTWRITEBYTECODE=1 python analysis/46_physiology_constrained_model_reconstruction/verify_local_cbm_core.py --output /tmp/ae4_task46g_verification.json
```

This control verifies the features exercised here. It is not a proof for every
possible LP or a validation of future AE4 constraints.

## Next scientific stage

The evidence split and Task 46E network construction remain pending. Document
the transport ledger, stoichiometry, storage assumptions, charge audit, directions,
bounds and shared WT/KO capacities before scientific optimisation. Keep the
isolated NKCC assay's protocol limitations visible; do not invent a numerical
tolerance. Publish the network checkpoint before phenotype optimisation.
The preserved Task 44 replay does not need to be repeated for this purpose.
