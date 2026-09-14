# Active phase: Task 36 manual mechanism verification

Work only on `codex/task-36-minimal-nahco3-alkalinity`.

This branch is being verified interactively before any Codex execution. Do not launch a search programme from this branch.

## Fixed scientific decisions

- Start from the accepted Task 31 R09 resting model and mechanistic Cha et al. NHE1 carrier amount `2.339370005697548e-05 fmol`.
- Preserve the existing AE4 routed topology and 1 Cl : 1 monovalent cation : 2 HCO3 working stoichiometry.
- Preserve AE2, NKCC1, pump, K channels, CaCC, CO2, paracellular, water, bath, geometry and existing AE4 beta/PKA regulation unless a literal wiring defect is found.
- Do not force the stimulated 70/30 chloride partition at REST.
- The accepted Task 31 REST must remain an exact nesting limit.

## Task 36 structural intervention

The discarded 1 Na : 1 HCO3 electroneutral NBC draft is not valid.

The active Task 36 surrogate is a stimulus-recruited electrogenic NBCe-like cycle:

`Na_o + 2 HCO3_o <-> Na_i + 2 HCO3_i`.

Positive inward cycle source signature in conserved cell coordinates is

`(+1 Na, 0 K, 0 Cl, +2 TIC, +2 TA)`.

It carries net `-1` charge-equivalent into the cell per inward cycle and therefore must participate explicitly in basolateral current closure.

The inward affinity is

`A_NBC = log((Na_o HCO3_o^2)/(Na_i HCO3_i^2)) + V_b/(RT/F)`.

The flux law is

`J_NBC = u_sec * G_NBC * tanh(A_NBC/2)`.

Use

`u_sec = clip((Ca - 0.058)/(0.25 - 0.058), 0, 1)`.

Thus NBC is exactly zero at Task 31 REST and fully recruited at the standard central `Ca=0.25 uM` condition.

The same secretory coordinate applies the source-fixed NHE1 stimulation

`G_NHE = 1 + 1.3 u_sec`,

so NHE1 is exactly Task 31 at REST and 2.3-fold active at the central stimulated endpoint. Do not change any Cha kinetic constant or add an NHE state/time constant.

The current fixed structural capacity is

`G_NBC = 0.11570913197464398 fmol/s`.

This is derived from the WT stimulated 70/30 chloride-loading feasibility algebra plus self-consistent current closure at the Task 31 reference chemical state. It is not a measured NBC abundance and not a knockout-phenotype fit.

## No search

- No Cartesian grids.
- No parameter sweeps.
- No pairwise or factorial searches.
- No random, global, evolutionary or Bayesian optimisation.
- No alternative NBC stoichiometries in the same task.
- No automatic addition of another transporter.
- No widening or retuning the NBC capacity after seeing genotype results.
- No fitting to the AE4-null secretion phenotype.

Before any genotype dynamics, verify the stoichiometry, charge/current closure, carbon accounting and exact Task 31 REST nesting. Then run at most the direct WT stimulated trajectory. Only after WT is physically admissible should AE4 5% and 0% be evaluated.

Do not merge or modify `main`.
