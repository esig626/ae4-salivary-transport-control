# AGENTS.md

These instructions apply to all work in this repository.

## Current scientific objective

The active task is Task 13C. The executable prompt is `prompts/13C_calcium_input_sweep.md`.

Task 13B produced a modern conservation-explicit WT model family with ten retained native-source WT resting roots. Those roots survived the required root confirmation and geometry checks. Their WT 600-second dynamics passed the numerical, conservation, solver, sustainment, co-stimulation, regulatory, and nearby-state gates. They failed only the absolute one-SMG WT secretion-scale gate. The strict AE4-null secretion magnitude and time course were not revealed.

Task 13C asks one narrow question only.

> Is the remaining absolute WT secretion deficit explained by the assumed stimulated calcium amplitude?

This is a one-variable sensitivity experiment on the frozen Task 13B WT models. It is not a reconstruction task.

## Hard freeze

Do not refit or recalibrate anything in Task 13C.

Keep frozen exactly as inherited from Task 13B

- all ten final WT resting roots;
- all whole-cell equations;
- all transporter capacities and source-class multipliers;
- NKCC1 scale 4 in the retained `N_ABS_NKCC` family;
- AE4 transport parameters and Na/K routing;
- NHE1 and AE2 parameters;
- pump and K-channel topology/fractions;
- CaCC and K-channel maximum conductances;
- hydraulic coefficients and lumen/outflow parameters;
- bath composition and geometry;
- cAMP/PKA regulatory equations, gains, and kinetic sensitivity members;
- solver tolerances and numerical gates;
- the one-SMG observation geometry and WT 9--10 uL/min gate.

The only variable that may change is the stimulated calcium input supplied to the existing calcium gate.

Do not compensate for calcium by changing any other parameter.

## Calcium sweep

Use the predeclared stimulated calcium grid

`0.10, 0.15, 0.18, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50 uM`.

Keep resting calcium at the existing `0.058 uM`.

The primary arm is the matched Task 13B `CCH_IPR` WT stimulation protocol with beta input unchanged. A `CCH_ONLY` arm may be run as a contextual diagnostic because the later parotid measurements used 300 nM carbachol alone, but it is not a substitute for the SMG CCh+IPR WT gate.

Do not introduce a calcium rise time, oscillation law, spatial gradient, or new calcium ODE in this task. Task 13C tests amplitude only.

## Source context for the sweep

Use the following published results only to justify the sensitivity domain and interpret it. Preserve gland and protocol context.

- Vera-Sigüenza et al. 2019, Bull Math Biol, DOI `10.1007/s11538-018-0534-z`: mouse parotid lobules stimulated with 300 nM CCh showed a cell-wide calcium response rising from about 72 nM to an oscillatory mean around 200 nM. The corresponding secretion model produced about 78 um^3/s mean flow from about 22 um^3/s at rest. The paper also showed that the spatio-temporal mean calcium was sufficient to reproduce mean secretion.
- Vera-Sigüenza et al. 2020, Bull Math Biol, DOI `10.1007/s11538-020-00712-3`: the seven-cell parotid acinus at mean calcium around 0.18--0.20 uM produced about 481 um^3/s total, with representative individual cells around 68.7--69.7 um^3/s. The paper states that earlier dual-agonist CCh+IPR studies corresponded to calcium around 0.5 uM and that raising calcium to 0.5 uM reproduces the larger secretion increase in that model.

These parotid values are not direct matched-SMG calibration targets. They justify testing the calcium range. Do not claim that 0.20 or 0.50 uM is the measured SMG calcium for the exact 2015 protocol unless a primary source directly establishes that.

## Holdout firewall

The AE4-null stimulated secretion magnitude and time course remain sealed.

Do not open, parse, digitize, evaluate, hash-decode, or use the held-out target values. Do not run AE4-null or AE2-null genotype trajectories in Task 13C.

Task 13C ends at the WT calcium result.

## Numerical standard

Reuse the Task 13B production implementation and gates. Do not rewrite the model unless a minimal runner change is required to expose calcium as the sweep variable.

At minimum retain

- positive physical states;
- charge/current/carbon/water accounting;
- the same 600-second physical-time protocol;
- the same one-SMG absolute flow-scale gate;
- the same minute-wise flow-shape logic;
- the same sustainment gate;
- Radau/BDF crosscheck for confirmatory cases;
- no first-passer pruning across the ten roots.

Task 13C should be computationally small. Do not launch a broad multi-agent reconstruction. One implementation/reproduction stream plus one independent audit is sufficient.

## Required outputs

Create

- `analysis/13C_calcium_input_sweep/evidence_and_freeze.md`
- `analysis/13C_calcium_input_sweep/final_answer.md`
- `results/13C_calcium_input_sweep/calcium_screen.csv`
- `results/13C_calcium_input_sweep/calcium_confirmatory.csv`
- `results/13C_calcium_input_sweep/calcium_threshold.json`
- `results/13C_calcium_input_sweep/frozen_manifest.json`
- `results/13C_calcium_input_sweep/artifact_hashes.csv`

The final report must state for every tested calcium amplitude

- per-cell flow range across the ten roots;
- mean and minute-wise flow behavior;
- existing one-SMG scale interval and whether it passes;
- implied effective cell count where meaningful;
- endpoint Cl, Na, K, pH, volume, and membrane potentials;
- calcium-gate activation fraction;
- numerical/conservation status.

## Completion classifications

End with exactly one of

1. `WT ABSOLUTE FLOW GATE PASSES WITH CALCIUM-ONLY CHANGE`
2. `WT ABSOLUTE FLOW GATE PASSES ONLY AT HIGH CALCIUM EDGE — MATCHED-SMG CALCIUM REMAINS DECISION-CRITICAL`
3. `CALCIUM-ONLY CHANGE INSUFFICIENT TO CLOSE WT ABSOLUTE FLOW DEFICIT`
4. `CALCIUM SWEEP EXPOSES PHYSIOLOGICAL OR NUMERICAL FAILURE BEFORE FLOW CLOSURE`

Do not convert a WT pass into an AE4 phenotype claim. A genotype test is a later task.

## General discipline

Do not edit `archive/`.

Do not merge to `main`.

Do not draft a manuscript.

Do not perform model reduction, GSPT, bifurcation analysis, or identifiability analysis.

Do not state that the 2018 paper was wrong or that AE4 itself is insufficient.
