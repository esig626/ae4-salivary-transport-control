# Task 31: implement a mechanistic NHE1 repair and test exact AE4 null

Repository: `esig626/ae4-salivary-transport-control`
Branch: `codex/task-31-nhe1-mechanistic-repair`

Read `AGENTS.md` first. This is an implementation task, not another diagnostic audit.

## Goal

Replace the still-unsatisfactory NHE1 representation with one biologically credible, transporter-specific kinetic model and test whether that repaired NHE1 permits a physiologically reasonable WT and exact-AE4-null resting state.

Do not spend this task comparing against the historical generic `tanh` law. Task 30 already established that the generic law is too large at ordinary pH/Na and implemented Vera-Siguenza et al. (2018) Eq. 26. Task 31 moves forward from that result.

The preferred mechanistic model is the published NHE kinetic formulation of Cha et al., Biophysical Journal 97 (2009) 2674-2683, DOI 10.1016/j.bpj.2009.08.053, which separates the ion-exchange cycle from intracellular proton-dependent modifier regulation and was fitted to forward/reverse mode, ion dependence and turnover data. Use the published paper as the primary source.

This is a cross-cell-type kinetic transfer. Preserve that caveat. Do not import cardiac cell geometry, cardiac background fluxes, cardiac Na/K pump, NCX, Ca handling or cardiac transporter density. Import only the NHE kinetic mechanism and its source-supported kinetic constants. Map the resulting exchange rate into the existing whole-cell salivary model with one explicit salivary-cell NHE1 density/activity scale if an absolute scale is required.

## 1. Implement one NHE1 model

Implement ONE new NHE1 model family only.

Requirements:

- obligatory electroneutral 1 Na inward : 1 H outward exchange;
- forward and reverse transport supported by the same cycle;
- correct dependence on intracellular/extracellular Na and H;
- explicit intracellular proton-dependent modifier/allosteric regulation as in the source model;
- finite rate at all positive concentrations;
- zero carbon source and zero net charge source;
- identical law and parameters in WT and exact AE4 null;
- no genotype switch, pH clamp, hidden acid source or bicarbonate sink.

Implement the Cha et al. equations and kinetic constants directly where the paper identifies them. If the paper leaves an absolute cell-specific expression/density factor, represent that as ONE scalar NHE1 density/activity parameter. Do not reinterpret a cardiac expression density as a salivary acinar measurement.

Retain the Task 30 Vera-Siguenza Eq. 26 evaluator as a selectable reference implementation, but the Task 31 production model must explicitly select the new mechanistic NHE1 law. Do not run a historical-law comparison panel.

Before whole-cell solves, add focused tests for:

- 1:1 Na/H source stoichiometry;
- zero carbon and charge source;
- forward and reverse signs under imposed gradients;
- finite saturation;
- stronger activation with intracellular acidification over the source-supported range;
- suppression at alkaline intracellular pH;
- source-model fixed-state checks or figure/value reconstruction where the paper provides enough information.

## 2. Calibrate only the salivary NHE1 amount/density

Do NOT calibrate multiple NHE1 shape parameters.

Use the source-model kinetic shape and constants unchanged. If an overall salivary-cell NHE1 density/activity factor is required, calibrate that ONE scalar using WT information only.

Primary calibration background: R09.

Use the measured WT resting pH centre 6.91 as the scalar calibration target. The diagnostic interval is 6.84 to 6.98. Do not use AE4-null data, secretion, or chloride phenotype in calibration.

Use a bounded one-dimensional root solve or secant/bracketing procedure. No grid. No multistart. Maximum six accepted scalar density evaluations. Freeze the single calibrated NHE1 density before any AE4-null solve.

After calibration, apply exactly the same NHE1 kinetic parameters and density to R10. Do not recalibrate separately for R10. R10 is a validation background.

Report the calibrated density relative to the source model's normalization and explain its units. A large density factor is not automatically a failure, but do not hide it.

## 3. Re-equilibrate WT and exact AE4 null

Use only R09 and R10.

For each background:

1. solve WT REST with the new NHE1 model;
2. freeze the NHE1 model and density;
3. set AE4 expression exactly to zero;
4. solve the exact-null REST with every other production parameter unchanged.

Do not alter AE4, AE2, NKCC1, Na/K pump, K channels, CaCC, CO2 permeability, paracellular pathways, water laws, cAMP/PKA regulation, chloride allocation or geometry.

The old fitted WT coordinates are warm starts, not immutable targets.

A resting state must pass the inherited numerical, conservation, current, positivity and volume gates.

Physiological checks to report, not all to fit:

WT context:
- pH: 6.91 +/- 0.07;
- Cl_i: 50.10 +/- 1.50 mM.

Exact AE4-null context:
- pH: 6.89 +/- 0.02;
- Cl_i: 36.50 +/- 1.60 mM.

Additional sanity checks:
- Na_i should remain in an ordinary cellular range; use <= 30 mM as the inherited diagnostic ceiling, not as an optimisation target;
- K_i positive and broadly cell-like;
- cell volume <= 3 pL;
- no huge TIC/HCO3 accumulation;
- no transporter capacity violation.

Success for this task is not defined by matching null chloride or secretion. The key success question is whether the repaired NHE1 permits an admissible exact-null resting cell with near-normal pH and ordinary sodium without changing another transporter.

## 4. Small dynamic check only after resting success

Only if a matched WT and exact-null REST pair is admissible for a background, run the existing CCh + IPR protocol at Ca = 0.25 uM for that pair.

Do not run Ca = 0.10 or 0.50 in this task.

Do not tune anything after seeing the dynamic result.

Report:

- integrated secretion;
- null/WT secretion ratio as an observation only;
- NHE1, AE2, AE4, NKCC1 and CaCC fluxes;
- Na, K, Cl, TIC, HCO3, pH and volume through the trajectory;
- luminal osmolarity and water fluxes.

The 30-35% experimental secretion deficit is context only. Do not fit to it.

## 5. Hard prohibition on search expansion

No combinatorial search.
No factorial or pairwise search.
No parameter grid.
No random/global/Bayesian/evolutionary optimisation.
No multidimensional phenotype fitting.
No chloride-allocation search.
No second NHE1 family if the Cha model fails.
No automatic move to another transporter.

If this one model cannot produce an admissible exact-null REST under the permitted scalar density calibration, stop and report that implementation outcome. Do not launch another diagnostic task inside this task.

## 6. Compute budget

This is a bounded implementation test, not a large computation.

- one scientific worker;
- one BLAS thread;
- at most 16 stationary solver calls total;
- at most 24,000 stationary residual evaluations total;
- at most 6 scalar NHE1-density evaluations;
- at most 4 stimulated integrations total;
- at most 30 minutes numerical execution.

Each individual stationary solve must have its own residual/evaluation ceiling. A single hard case must not consume the entire global budget.

Checkpoint completed cases. Never rerun an accepted identical case.

## 7. Required outputs

Write compact UTF-8 outputs under:

- `results/31_nhe1_mechanistic_repair/`
- `analysis/31_nhe1_mechanistic_repair/`

Include:

- `contract.json`
- `budget.json`
- `nhe1_source_and_equations.md`
- `nhe1_parameters.json`
- `resting_states.csv`
- `flux_ledger.csv`
- `trajectory_summary.csv` if trajectories run
- `verification.json`
- focused test log
- `final_answer.md`

The final answer must state plainly:

1. Was the mechanistic NHE1 model implemented correctly?
2. What single salivary NHE1 density/activity scale was obtained from WT R09?
3. Did R09 WT close physiologically?
4. Did R09 exact AE4 null close with ordinary Na, near-normal pH and cell volume <= 3 pL?
5. Did the same frozen NHE1 model work in R10 without recalibration?
6. If dynamics ran, what happened to secretion and chloride transport?
7. What remains unresolved?

## 8. GitHub publication

Use the connected GitHub integration for all remote writes.

Do NOT depend on shell `git push`, `gh auth`, SSH credentials or a PAT. Shell Git may be used only for local status/diff inspection.

Publish completed source edits, tests and compact text results to `codex/task-31-nhe1-mechanistic-repair` through the connected GitHub integration. Do not create binary result archives. Do not merge to main.

Publication failure must not trigger recomputation.
