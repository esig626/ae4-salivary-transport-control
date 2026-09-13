# Active phase: Task 29

Work only on `codex/task-29-exact-null-feasibility-audit`.
Read and execute `prompts/29_exact_null_feasibility_audit.md`.

This task is an algebraic and constitutive feasibility audit. It is NOT an
optimisation task and it must not resume Task 28 numerical closure attempts or
Task 27 phenotype optimisation.

## Hard restrictions

- No combinatorial, factorial, exhaustive pair, grid, random, global,
  evolutionary, Bayesian, high-dimensional or phenotype-directed search.
- No multidimensional nonlinear optimisation or calibration.
- No stimulated ODE trajectories.
- No subagents, parallel scientific workers or automatic task decomposition.
- Do not add a new transporter, signalling cascade, hidden acid source,
  bicarbonate sink, pH clamp, genotype-specific switch or output-level repair.
- Do not alter the production model equations during the audit. A small
  standalone analysis helper is allowed if needed to evaluate existing laws
  and the Task 28 NHE1 candidate.
- Maximum numerical budget: 50 fixed-state transporter/RHS evaluations,
  12 one-dimensional bracketed root solves, one numerical worker, one BLAS
  thread, and 10 minutes of numerical execution. Stop before any limit is
  exceeded. No retries disguised as a new method.
- The scientific result may be feasible, infeasible under specified fixed
  parameters/target conditions, or unresolved. Do not force a contradiction.
- Distinguish sign/topology impossibility from insufficient capacity and from
  mere solver failure.
- Git authentication is not a scientific prerequisite. Test ordinary Git
  access once at task start. If unavailable, record that fact and continue.
  Do not attempt credential recovery, repeated login, manual Git objects,
  Base64, archives or alternative publication machinery.
- Do not merge or modify `main`, archive branches, Task 28, or other remote
  task branches.

The exact Task 28 partial findings supplied by the user are summarized in
`reference/task28_partial_handoff.md`. They are handoff evidence, not a new
validated repository result. Re-derive every equation used for the Task 29
conclusion from the seed model code.
