# Active phase: Task 30

Work only on `codex/task-30-nhe1-model-repair`.
Read and execute `prompts/30_nhe1_model_repair.md`.

This task is an implementation task. Stop diagnosing the existing NHE1 defect and replace it with a transporter-specific NHE1 model based on the published 2018 salivary secretion model.

## Hard restrictions

- No combinatorial, factorial, pairwise, grid, random, global, Bayesian, evolutionary or high-dimensional search.
- No automatic escalation to other transporters if NHE1 repair fails.
- Do not alter AE4, AE2, NKCC1, pump, K, CaCC, CO2, paracellular or water laws in this task.
- Do not fit to the 30-35% AE4-null secretion phenotype.
- Use only R09 and R10 for new calculations.
- Implement one new NHE1 law: the published 2018 Eq. 26 law, with its stated parameters and units. Preserve the current law as a selectable historical comparator, not as the active repair.
- One optional WT-only scalar calibration of NHE1 activity per background is allowed only if the published activity value prevents a physiologically reasonable WT rest. No knockout information may enter that calibration.
- Maximum numerical budget: 12 stationary solver calls total, 10,000 stationary residual evaluations, 8 scalar NHE1-activity evaluations total, 6 stimulated integrations total, one worker, one BLAS thread, 30 minutes numerical execution.
- If an admissible WT and exact-null REST are obtained, run only the central Ca=0.25 uM stimulation for WT/null in R09/R10. Do not tune after seeing secretion.
- If the new NHE1 law fails under this bounded implementation, stop and report failure. Do not invent a second NHE1 family or start another diagnostic/search task.
- Do not merge or modify `main`.

## GitHub publication in ChatGPT Work

Use the connected GitHub integration for publication. Do not depend on ordinary shell Git authentication.

- Do not run `gh auth`, request a PAT, configure SSH, or stop because `git push` asks for credentials.
- Shell Git is allowed for local `status` and `diff` inspection only.
- At the end of the task, publish all completed source edits, focused tests, and compact UTF-8 result/report files directly to `codex/task-30-nhe1-model-repair` through the connected GitHub integration.
- Open or update a PR only through the connected GitHub integration if the task requests one.
- Do not use manual Git objects, Base64 upload loops, archives, or binary trajectory publication as a workaround.
- A publication problem must never trigger scientific recomputation.

Success means a working NHE1 implementation plus properly re-equilibrated WT/null states and, if those rests are admissible, a small matched dynamic check. It does not require matching the secretion phenotype.
