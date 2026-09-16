# Task 50B implementation verification

**PASS: ten tests, zero failed scientific comparisons.** Classification:
**TARGET-CALIBRATED CONSTRUCTION**.

The public model is `Task50EffectiveCouplingModel` in
`src/modern_full_model/task50_effective_coupling.py`. Its only new constant is
`FIXED_LAMBDA = 0.89488127156712`. It inherits the existing Task 41 conductance
reconstruction, RHS and integration dispatch, and overrides the multiplier and
diagnostic labels. It does not create the old recruitment parameter object.
There is no new model state, signalling law, timing, fit or conductance input.

`effective_coupling_factor(a,beta,e)` evaluates
`1-FIXED_LAMBDA*a*beta*(1-e)`. Exact off case branches return `1.0`; the inherited
`effective_model(1.0)` then returns the original parent object. Active cases
change only the existing apical chloride conductance before full current
closure. The inherited calcium Hill gate remains intact.

| Check | Result |
| --- | --- |
| WT, REST, CCh only null and AE2 knockout parent nesting | 21 fixed state comparisons, bit identical full RHS and physical diagnostics |
| Beta one factor versus the full frozen Task 41 coefficient | Maximum absolute difference `4.5102810375396984e-15`; fixed allowance `8e-15` |
| Full beta one RHS at the shared rest and each saved endpoint | Six comparisons; largest native component difference `8.090750291955828e-15` |
| Membrane voltage comparison | Maximum difference `2.8449465006019636e-16 V` |
| Membrane currents | Maximum difference `7.237830359838992e-25 A` |
| Membrane amount sources and fluxes | Maximum difference `8.090750291955828e-15 fmol/s` |
| Charge/current/conservation checks | Same source identities; all inherited conservation gates pass |
| Unrelated fixed state evaluations | Water, AE4, homeostasis, speciation, CO2, outflow and state charge unchanged |
| Parameters, states and protocols | Only the declared conductance differs; parent state names, AE4/NKCC/NBC/regulation and protocol settings unchanged |
| Public integration helper | Mock confirms the solver receives the amended RHS; no integration performed |
| Frozen Task 41 execution inputs | All seven original SHA-256 values match |

All numerical allowances were published in D50-07 before tests. The current
closure tolerances were not changed. The two coefficient literals are still
distinct; these are numerical comparisons under the user's authorisation.
They are not proof of exact mutant trajectory identity or a global solution
error bound.

Verification attempt 01 performed no tests and no model evaluations: the old
Task 40 fixture checked the historical AGENTS hash. The full audit found only
that governance file differed among 150 inputs. Attempt 02 pins the current
binding AGENTS blob, verifies the other 149 original inputs and the Task 40
execution/rest/parameter hashes, then uses the unchanged parent factory. No
historical source, output or instruction file was restored or altered.

Attempt 02 performed 80 core evaluations including constructor checks. Across
both attempts there were zero production trajectories, stationary solves,
parameter fits, optimisation calls or mechanism searches. Both logs are kept
under `output/`. `run_verification.py` prevents an optional repeat after success.
The tested source hashes are recorded in `output/verification_attempt_02.json`.
