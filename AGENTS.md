# AGENTS.md

These instructions apply to all work in this repository.

## Current scientific objective

The active task is Task 18. Execute:

`prompts/18_fixed_wt_inverse_transporter_rebalance.md`

Work only on branch:

`codex/task-18-fixed-wt-inverse-rebalance`

Tracking issue: #13.

Task 18 supersedes Task 17's inherited-capacity feasibility experiment.

## Scientific inheritance

Repository inheritance is completed Task 17 commit:

`fb2afcf3e23af1d19b297280ad67953370f7c9a9`

The unchanged production scientific model baseline remains completed Task 14B commit:

`4d4403f279fc36aec940c9e4759486d02f07e3f6`.

Task 17 correctly kept the WT conserved states fixed, but it still treated inherited capacity assumptions such as total K conductance 14 nS as hard feasibility constraints. Task 18 removes that mistake.

## Fixed WT means fixed conserved state

For every inherited accepted WT root, keep all conserved intracellular and luminal amounts and both volumes exactly unchanged, together with the inherited resting regulatory coordinate.

Do not perform a new WT state search, state continuation, state fitting, rounding or reconstruction.

Algebraic membrane voltages may reclose self-consistently under the rebalanced capacities. Report their changes. Do not treat historical voltage targets as immutable laws.

## Chloride allocation

Use exactly three conditions:

- inherited baseline;
- AE4 positive basolateral chloride-loading share 0.10;
- AE4 positive basolateral chloride-loading share 0.30.

No additional share values or threshold searches are allowed.

For each root, identify the actual positive basolateral intracellular chloride-loading pool from the production flux ledger.

At each target, increase AE4 by scaling the complete coupled transporter so it carries the declared share of that inherited positive-loading pool.

Reduce every other positive basolateral chloride-loading pathway proportionally so the total positive-loading pool remains equal to the inherited WT value.

Do not scale AE4 chloride independently of its Na/K/HCO3/TIC/alkalinity coupling.

Negative chloride counterfluxes are not part of the positive-loading denominator.

## Rebalance the rest of the network

The changed AE4 contribution perturbs Na, K, carbon, alkalinity and current balance even though the WT conserved state is fixed.

Restore all complete production steady-state equations by adjusting existing uncertain capacity-like parameters in the current topology.

At minimum allow the solver to consider:

- total K conductance and its membrane allocation;
- CaCC conductance/capacity;
- apical and basolateral Na/K pump capacities;
- NHE1 capacity;
- neutral CO2 exchange capacities;
- non-Cl paracellular Na, K and HCO3 capacities;
- paracellular Cl capacity only if required to preserve fixed chloride balance under voltage reclosure;
- any other already-existing active capacity multiplier required for full source rank.

Do not add a transporter, background current or new mechanism.

## Hard constraints versus reference values

Hard constraints are structural physics and the fixed WT state:

- exact fixed WT conserved state;
- full resting RHS equal to zero within inherited numerical tolerances;
- target AE4 positive-loading share;
- preserved total positive basolateral chloride-loading pool;
- existing stoichiometry, topology, source signs and thermodynamic structure;
- nonnegative physical capacities/permeabilities;
- fractions in [0,1];
- bath, chemistry and geometry;
- charge, current, carbon, alkalinity, water and lumen accounting;
- finite physical model quantities.

The following are **reference values, not hard feasibility walls**:

- 14 nS total K conductance;
- 31.4 nS CaCC maximum conductance;
- old NHE1 parameter bounds;
- old AE4 carrier bounds;
- old common-conductance bounds;
- previous pump capacities;
- previous paracellular capacities;
- previous CO2 exchange capacities;
- capacity values imported from parotid or another preparation as modeling assumptions.

Report departures from those values; do not reject a mathematically and physically valid fixed-WT solution merely for exceeding them.

## Canonical solution rule

Before genotype evaluation, solve the fixed-WT inverse problem for every root/share using no genotype information.

If multiple parameterizations satisfy the hard constraints, choose one canonical solution by:

1. minimizing the largest absolute log-fold change among adjustable positive capacities relative to baseline;
2. then minimizing the sum of squared log-fold changes.

Use independent positive membrane capacities where this avoids arbitrary penalties on partition fractions.

Do not use the known AE4-loss phenotype or Task 14B secretion result in feasibility, optimization, ranking or root selection.

A local optimizer failure is not an infeasibility proof. Preserve full rank/nullspace information and distinguish structural infeasibility from numerical nonconvergence.

## Freeze before genotype evaluation

First complete, verify, commit and push the entire WT inverse solution set.

Only after that frozen WT checkpoint may genotype results be accessed.

For every feasible modified parameterization:

- run matched WT dynamics at calcium 0.10, 0.25 and 0.50 uM;
- reduce AE4 expression to 0.05 with all Task 18 parameters otherwise frozen;
- evaluate AE2 loss with the same frozen parameters;
- do not refit either genotype;
- do not attempt exact AE4 zero.

Reuse hash-valid baseline trajectories wherever scientifically identical.

Absolute whole-gland mapping remains nonblocking. Detailed AE4-loss time-course shape remains diagnostic only.

## General discipline

Do not edit `archive/`.

Do not merge to `main`.

Do not draft manuscript text.

Do not perform unrelated model reduction, GSPT or identifiability work.

Commit and push completed Task 18 work only to `codex/task-18-fixed-wt-inverse-rebalance`.
