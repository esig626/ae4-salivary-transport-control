# AGENTS.md

These instructions apply to all work in this repository.

## Current scientific objective

The active task is Task 13B. The executable prompt is `prompts/13B_modern_full_model_reconstruction_multiagent.md`.

Task 13 established that modern state-resolved AE4 transport structure can add genuine transporter-level behavior, but no tested candidate produced a valid whole-cell WT calibration on the inherited seven-state chassis. It also showed that the current evidence does not identify a unique whole-cell correction from AE4-null Cl and pH alone.

Task 13B therefore has one priority only.

> Construct and validate a physiologically credible modern full salivary acinar-cell model before attempting any model reduction, GSPT, bifurcation analysis, identifiability analysis, or manuscript drafting.

Do not work on a three-ODE reduction in this task. Do not run continuation merely because it is available. The full model must work first.

## Required interpretation of earlier tasks

Do not state or imply that AE4 itself is insufficient.

The supported statement is narrower.

> The tested coarse-grained and state-resolved AE4 implementations did not yield a valid WT-calibrated whole-cell model on the inherited historical chassis over the tested source-supported domains.

The biological evidence directly supports an important role for AE4 in intracellular chloride accumulation and stimulated saliva secretion. The failure is a model-reconstruction problem, not evidence against AE4 biology.

## Evidence hierarchy

When sources disagree, use the following order.

1. Conservation laws, charge balance, thermodynamic feasibility, dimensional consistency, and membrane orientation.
2. Direct primary experimental measurements in salivary acinar cells or glands.
3. Direct transporter-level measurements in heterologous expression systems, with assay context preserved.
4. Structural and mutagenesis evidence that constrains molecular mechanism.
5. Reproducible behavior of historical implementations.
6. Published historical equations and parameter tables.
7. New modeling decisions, which must be labelled and stress-tested.

The 2018 article and archived code are historical model records, not ground truth. The old article must not be framed as wrong. This project is a source-informed revisit using later biology.

## Primary evidence that must be incorporated

At minimum use and audit the following primary sources.

- Peña-Münzenmayer et al. 2015, JBC, DOI `10.1074/jbc.M114.612895`, for AE4-null saliva, resting chloride and pH phenotypes, AE2 comparisons, and the reported early-versus-sustained secretion pattern.
- Peña-Münzenmayer et al. 2016, JGP, DOI `10.1085/jgp.201611571`, for direct Na and K transport, electroneutrality, reversal, cation dose response, Hill summaries, and source-discussed stoichiometries.
- Peña-Münzenmayer et al. 2021, AJP Gastrointestinal and Liver Physiology, DOI `10.1152/ajpgi.00145.2021`, for beta-adrenergic activation of AE4, H89 sensitivity, PKAc activation, and the S173/S273 mutant evidence.
- Catalán et al. 2025, AJP Cell Physiology, DOI `10.1152/ajpcell.00346.2024`, for cation coordination, Na/K-specific mutant behavior, and later transport-cycle hypotheses.
- Almássy et al. 2018, Pflügers Archiv, DOI `10.1007/s00424-018-2109-0`, for apical as well as basolateral Na/K-ATPase localization, apical Ca-activated K current, and the associated whole-cell transport topology.
- The 2018 AE4 mathematical model, DOI `10.1007/s11538-017-0370-6`, as historical mathematical lineage only.

Do not promote proposed molecular cartoons, MD trajectories, or speculative stoichiometries to measured fact.

## Dynamic cAMP and PKA regulation is mandatory

Task 13 treated PKA mainly through imposed activity changes and did not establish a physical activation time scale. Task 13B must explicitly investigate dynamic regulation.

The biological chain to represent is

`beta-adrenergic stimulation -> adenylate cyclase/cAMP -> PKA -> AE4 regulatory state -> AE4 transport`.

The 2021 experiments establish regulation and S173 dependence, but they do not automatically establish a unique kinetic scheme, a direct phosphorylation rate, a physiological time constant, or a universal multiplicative gain.

Therefore build the smallest dynamic regulatory models justified by the data and label all unmeasured kinetic choices. Compare nested alternatives such as

- an effective cAMP state driving a PKA activation state;
- a single effective PKA/AE4 regulatory state if the upstream cAMP kinetics are not separately identifiable;
- phosphorylation/dephosphorylation of an AE4 regulatory fraction;
- state-specific modulation of the AE4 transport cycle versus a common capacity modulation.

Do not choose a delay or time constant because it makes the AE4-null secretion curve look right.

## Full-model reconstruction scope

The model may need to exceed the inherited seven-state chassis. That is allowed and expected if required by conservation or new biology.

The full modern model should explicitly decide how to represent

- intracellular Na, K, Cl, HCO3, H/ pH, and CO2 or total inorganic carbon;
- cell volume and osmotic water balance;
- luminal Na, K, Cl, and any required bicarbonate/carbon species;
- apical and basolateral membrane potentials or equivalent current constraints;
- NKCC1;
- NHE1;
- AE2;
- state-resolved or reduced-but-source-consistent AE4;
- Na/K-ATPase with experimentally supported apical and basolateral topology;
- Ca-activated K channels with experimentally supported apical and basolateral topology;
- apical Cl current;
- tight-junction/paracellular transport;
- CO2/HCO3 chemistry and buffering;
- dynamic beta/cAMP/PKA regulation of AE4;
- the calcium stimulus/input used to drive secretion.

Do not add complexity merely to improve fit. Every added state or flux must solve a documented closure, conservation, timing, or evidence problem.

## Calibration and holdout firewall

The AE4-null stimulated secretion magnitude and time course remain the decisive held-out validation.

Do not fit to

- the approximately 35 percent reduction in 10-minute saliva;
- the 0.65 KO/WT total saliva ratio;
- the observed early similarity followed by later divergence;
- any digitized AE4-null secretion curve.

The reconstruction should use WT data, transporter assays, resting ionic phenotypes, source-supported geometry/topology, and where necessary released knockout ionic data only according to the staged rules in the Task 13B prompt.

A model that is tuned until it produces the desired KO secretion ratio has not explained that ratio.

## Multi-agent standard

Task 13B must be run as a real multi-agent scientific investigation. Independent roles must include at least

1. primary experimental evidence and protocol extraction;
2. cAMP/PKA regulatory biology and kinetic-model audit;
3. AE4 transport-cycle and thermodynamic modeling;
4. whole-cell epithelial transport and acid-base reconstruction;
5. numerical calibration and root finding;
6. dynamic simulation and timing calibration;
7. adversarial scientific audit;
8. independent numerical reproduction;
9. lead synthesis and decision making.

Agents must preserve disagreements and uncertainty rather than averaging them away.

## Iteration requirement

Do not stop after the first failed reconstruction.

Proceed through successive model generations. Each generation must state

- which defect from the previous generation is being addressed;
- which equations or topology changed;
- which parameters were newly introduced;
- what independent evidence constrains them;
- which gates improved or worsened;
- whether the change is retained or reverted.

Continue until a source-supported full model passes the declared physiology gates or until one experimentally unmeasured quantity is shown to be decision-critical and cannot be bounded from existing sources.

A generic `STOP` is not acceptable.

## Numerical standard

Every accepted model must satisfy

- species and charge conservation where applicable;
- thermodynamic sign and reversal checks;
- positive physical states;
- stable deterministic root finding with multistart or continuation as needed;
- WT resting physiology before knockout predictions;
- independent solver reproduction for dynamics;
- sensitivity to tolerances and initial conditions;
- explicit dimensional and time-scale ledgers.

Do not hide numerical failure behind normalized ratios.

## Claim discipline

Separate throughout

1. direct experimental fact;
2. source-supported mechanistic constraint;
3. exact model deduction;
4. numerical model result;
5. historical implementation evidence;
6. new modeling assumption.

Do not claim the old article was wrong. Do not claim AE4 itself is insufficient. Do not call a parameterized extension identified unless the data actually identify it.

## Completion standard

Task 13B is complete only with one of these substantive outcomes.

1. `FULL MODEL VALIDATED — HELD-OUT AE4 PHENOTYPE RECONSTRUCTED`
2. `FULL MODEL VALIDATED — PHENOTYPE PARTIALLY RECONSTRUCTED; ONE SPECIFIC DYNAMIC OR FLUX DEFICIT REMAINS`
3. `MULTIPLE FULL MODELS VALIDATED — PHENOTYPE MECHANISM NONUNIQUE`
4. `FULL MODEL NOT IDENTIFIABLE FROM EXISTING DATA — ONE DECISION-CRITICAL MEASUREMENT SPECIFIED`

Outcomes 1 to 3 require a valid WT model and held-out knockout evaluation. Outcome 4 requires proof that the missing measurement changes the model-class decision, not merely that more data would be useful.

## Manuscript and reduction discipline

Do not draft a paper, reduce the system to three ODEs, apply GSPT, or perform a bifurcation study during Task 13B.

First freeze a validated full model. Reduction and dynamics are later tasks.