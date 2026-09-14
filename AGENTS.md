# Active phase: Task 40 equal Na/K AE4 cation routing

Work only on `codex/task-40-ae4-equal-cation-routing`.
Read `analysis/40_ae4_equal_cation_routing/source_contract.md` first, then execute `prompts/40_ae4_equal_cation_routing.md`.

Task 40 tests exactly one scientific change to the completed Task 39 model: replace donor-concentration-weighted AE4 cation routing with a fixed 50:50 Na/K split while preserving the inherited net AE4 cycle law and every other mechanism.

## Frozen parent

Parent scientific lineage is merged Task 39 at:

`afd101448763439f369f2682c467f069ea18a442`

Keep fixed:

- Palk/Benjamin NKCC1 law and `alpha_eff`;
- inherited NKCC1 stimulation multiplier;
- Cha NHE1 kinetics/carrier amount/stimulation;
- electrogenic 1 Na : 2 HCO3 NBC and its capacity/current/recruitment;
- total AE4 net cycle law, capacity, expression scaling, 1 Cl : 1 cation : 2 HCO3 stoichiometry, and beta/cAMP/PKA regulation;
- AE2, pump, K channels, CaCC, CO2/acid-base, paracellular transport, water, bath, geometry and standard stimulus.

## Fixed Task 40 AE4 routing

Let `J4` be the inherited AE4 net chloride source, positive inward.

Use exactly:

- `S_Na = -0.5 J4`
- `S_K = -0.5 J4`
- `S_Cl = +J4`
- `S_TIC = -2 J4`
- `S_TA = -2 J4`.

Do not change `J4` itself. Do not route by donor concentrations. Do not test another Na/K fraction.

The charge-equivalent source must remain zero and total cation source must remain `-J4`.

## Algebraic preflight

At the accepted Task 31 R09 state:

- donor-weighted Na fraction = `0.09062450161183984`
- `J4 = 0.004937067441546469 fmol/s`
- expected change in Na source = `-0.0020211144444590447 fmol/s`
- expected change in K source = `+0.0020211144444590447 fmol/s`
- direct Cl/TIC/TA AE4 sources unchanged
- charge perturbation zero.

Exact Task 31 REST nesting is therefore not expected.

## Execution order

1. Implement and source-test the 50:50 routed AE4 evaluator.
2. Reproduce the algebraic preflight independently.
3. Perform exactly one intended WT resting-state solve from the accepted Task 31 R09 state as the sole initial guess. One purely numerical retry is allowed; no calibration, multistart or root search family.
4. If the WT rest is admissible, freeze it.
5. From that same frozen WT rest run exactly three standard 600 s trajectories: WT, AE4=0.05, AE4=0.00.
6. Compare Na/K/Cl, NKCC1 compensation, AE2, NBC/NHE1/pump/channel responses, and secretion against Task 40 WT and Task 39 compensation.
7. Report held-out phenotype context without fitting.

## No search

No optimisation, routing-fraction sweep, alternate split, R10, genotype-specific rest, alternate NKCC1/NBC/NHE1/AE4 mechanism, parameter retuning, stimulus changes, new transporter, or phenotype fitting.

The 50:50 split is the only routing hypothesis in this task.

## Compute limits

- one intended WT stationary solve;
- at most one purely numerical stationary retry;
- three intended production integrations if WT rest passes;
- at most one purely numerical integration retry across the task;
- zero optimisation calls;
- zero parameter sweeps;
- one worker;
- one BLAS thread;
- maximum 15 minutes numerical execution.

## Publication

Use the connected GitHub integration for remote writes.
Do not use shell Git authentication, PATs, SSH configuration or `gh auth`.
Do not merge or modify `main`.
Publish compact UTF-8 source, tests and result files only.

Stop at the first declared scientific failure or after the Task 40 final report.
