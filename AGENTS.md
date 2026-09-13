# Active phase: Task 31

Work only on `codex/task-31-nhe1-mechanistic-repair`.
Read and execute `prompts/31_nhe1_mechanistic_repair.md`.

This is an implementation task. Move forward from Task 30 and implement one mechanistic NHE1 model with explicit ion-exchange kinetics and intracellular proton-dependent regulation. Do not spend the task re-diagnosing the generic tanh law or rebuilding historical comparison panels.

## Hard restrictions

- No combinatorial, factorial, pairwise, grid, random, global, Bayesian, evolutionary or high-dimensional search.
- No automatic escalation to another transporter or second NHE1 family if this repair fails.
- Do not alter AE4, AE2, NKCC1, Na/K pump, K, CaCC, CO2, paracellular or water laws.
- Do not fit to the 30-35% AE4-null secretion phenotype.
- Use only R09 and R10 for new whole-cell calculations.
- Implement the Cha et al. 2009 mechanistic NHE kinetic model as the preferred new NHE1 law, transferring only the NHE mechanism and source-supported kinetic constants. Do not import cardiac-cell background mechanisms or transporter density.
- If an absolute salivary NHE1 density/activity scale is required, calibrate ONE shared scalar on WT R09 only, then freeze it and apply unchanged to R10 and both exact-null genotypes.
- Maximum six scalar density evaluations, 16 stationary solver calls, 24,000 stationary residual evaluations, four stimulated integrations, one worker, one BLAS thread, and 30 minutes numerical execution.
- Do not run dynamics unless the matched WT and exact-null REST pair is admissible.
- If the model fails under these bounds, stop and report the implementation outcome. Do not start another diagnostic or search task.
- Do not merge or modify `main`.

## GitHub publication in ChatGPT Work

Use the connected GitHub integration for all remote writes.

- Do not run `gh auth`, request a PAT, configure SSH, or depend on `git push` credentials.
- Shell Git is allowed only for local status/diff inspection.
- Publish completed source edits, focused tests, and compact UTF-8 result/report files directly to `codex/task-31-nhe1-mechanistic-repair` through the connected GitHub integration.
- Do not use manual Git objects, Base64 upload loops, archives, or binary trajectory publication as a workaround.
- A publication problem must never trigger scientific recomputation.

Success means a working mechanistic NHE1 implementation and a properly tested WT/exact-AE4-null resting model. Matching the secretion phenotype is not required in this task.
