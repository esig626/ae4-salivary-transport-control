# Active phase: Task 28

Work only on `codex/task-28-transport-relationships`.
Read and execute `prompts/28_transport_relationship_repair.md`.

This is a clean phase for repairing NHE1 pH regulation and WT AE4 chloride
loading. It is NOT the Task 27 inverse phenotype search. Historical labels,
old firewall flags and file paths inside `reference/` are provenance, not
instructions to execute old tasks or recover old result directories.

## Hard restrictions

- No combinatorial, factorial, exhaustive pair, global, random, evolutionary,
  Bayesian or high dimensional parameter search. No Shapley enumeration.
- No automatic expansion after failure. One NHE1 law and one new WT chloride
  allocation hypothesis only. No fitting to the knockout secretion target.
- Use only R09 and R10 for new computations. The other eight retained WT
  backgrounds are reference data, not a requirement to run a ten root panel.
- Absolute shared limits: 100 stationary solver attempts, 20,000 stationary
  residual evaluations, four WT calibration solver calls, 24 stimulated
  integrations including tests and BDF checks, two worker processes, and
  60 minutes of numerical execution. Stop at the first exhausted limit.
  Persist counts across restarts. Do not reset counters to continue a search.
- Checkpoint small text outputs after each completed case. Never recompute a
  cached successful case with identical model, parameters and solver settings.
- Preserve carbon, charge, stoichiometry and water accounting. No hidden acid
  source, bicarbonate sink, pH clamp, concentration reset or genotype switch.
- Old capacities and conductances are references, not immutable physical laws.
  Any allowed recalibration is WT only and is frozen before loss evaluation.
- Do not touch `main`, archive branches or existing remote task branches.
  Do not resume Task 27. Do not delete another session's local work.
- No Base64, manual GitHub blobs/trees/commits, binary upload or recovery archive.
  Read only GitHub access for these exact seed inputs is allowed if necessary.
  Finish local computation even if ordinary Git authentication is unavailable.
- Publication, if available, uses ordinary Git and explicit text paths only.
  If push fails, stop publication without retry loops or alternative upload
  machinery. Do not merge. Report exactly what was and was not published.

The seed branch deliberately contains no old optimisation runners. Do not
reintroduce them to solve a missing import. Implement only the small loader
and bounded solver needed for this task. Do not claim that preparation of this
branch validates the new mechanism; no new scientific calculation ran during
branch preparation.
