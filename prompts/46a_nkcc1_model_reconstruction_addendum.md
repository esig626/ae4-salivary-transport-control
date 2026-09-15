# Task 46A addendum: reconstruct NKCC1 rather than merely cap compensation

This addendum is mandatory when executing Task 46.

## Scientific premise

Treat the current Palk/Benjamin-derived NKCC1 implementation as a candidate model defect, not as a trusted fixed law.

The completed Task 44 model permits very strong state-driven NKCC1 compensation after AE4 loss. This is in tension with the salivary experiments of Peña-Münzenmayer et al. (2015), in which isolated NKCC1-dependent chloride uptake under HCO3-free conditions was unchanged in AE4-null and AE2-null acinar cells relative to controls, while AE4-null cells showed reduced resting intracellular chloride and slower stimulated chloride uptake. The 2018 AE4 model also predicted only minimal NKCC1 upregulation after AE4 loss.

Do not merely impose an arbitrary ceiling on the existing NKCC1 law. Investigate whether the law itself, its state dependence, capacity calibration, stimulation multiplier, reversal structure, saturation, carrier representation, or coupling to sodium/potassium/chloride and pump activity is inappropriate for salivary acinar physiology.

## Required NKCC1 reconstruction programme

Before accepting any Task 46 model:

1. Audit the present NKCC1 equation term by term and identify which parts are:
   - directly supported by salivary acinar data;
   - inherited from Palk/Benjamin or another preparation;
   - phenomenological;
   - calibrated from the present model;
   - responsible for large AE4-null compensation.

2. Perform a literature search for alternative NKCC1 kinetic formulations, including mechanistic carrier-state, thermodynamic, saturating and regulatory models. Prefer formulations with direct NKCC1 evidence and assess whether they are appropriate for salivary acinar cells rather than importing a law merely because it exists.

3. Construct a finite candidate panel of NKCC1 models. At minimum include:
   - the current Task 44 law as the reference failure case;
   - a thermodynamically reversible saturating NKCC1 law with explicit 1 Na : 1 K : 2 Cl stoichiometry;
   - any better-supported mechanistic carrier model found in the literature;
   - a parsimonious salivary-specific effective law if the literature does not support a full carrier-state model.

4. Every candidate must preserve exact 1 Na : 1 K : 2 Cl stoichiometry, electroneutrality, reversal direction, finite capacity and physically sensible saturation. Do not use a law that can generate unlimited compensatory flux merely because intracellular chloride falls.

5. Calibrate WT NKCC1 only from WT evidence and explicitly declared construction constraints. Do not use AE4-null NKCC1 behaviour to tune WT parameters unless this is clearly labelled as consuming that observation as calibration evidence.

6. Treat the Peña isolated NKCC1 assay as a principal holdout wherever possible. Under a faithful or honestly approximated HCO3-free/AE-blocked/CaCC-blocked protocol, the accepted model should not predict a material genotype-dependent increase in NKCC1 activity after AE4 loss.

7. The accepted whole-cell model must not recover most missing AE4 chloride by increasing NKCC1 activity. NKCC1 may remain the dominant WT chloride loader and may change secondarily through normal state dependence, but strong compensatory replacement of the missing AE4 chloride is a rejection criterion unless directly supported by protocol-matched experimental evidence.

8. Quantify for each candidate:
   - WT NKCC1 flux and chloride-loading fraction;
   - AE4-null NKCC1 flux under physiological stimulation;
   - AE4-null NKCC1 flux under the isolated assay mimic;
   - local chloride compensation gain;
   - finite-time replacement fraction;
   - integrated NKCC1 increase;
   - effects on Na_i, K_i, Cl_i, pH, volume, membrane potential and secretion.

9. Reject candidate NKCC1 laws that reproduce the AE4 secretion phenotype only by violating WT physiology, thermodynamics, isolated NKCC evidence, intracellular chloride behaviour or other declared holdouts.

10. Once an NKCC1 law is accepted, freeze it before evaluating downstream architectural additions. Do not use later CaCC, NBC, NHE1 or pump modifications to conceal a still-unphysical NKCC1 response.

## Checkpoint discipline

Publish NKCC1 work incrementally to the Task 46 branch. Required additional checkpoints are:

1. current NKCC1 law forensic audit;
2. literature/model candidate ledger;
3. candidate equations and provenance before simulation;
4. results for each NKCC1 candidate family;
5. accepted/rejected NKCC1 decision and frozen accepted law;
6. only then continue the remaining Task 46 whole-cell reconstruction.

At every checkpoint fetch remote state, pull with `--ff-only`, verify, commit and push immediately. Never force-push.

## Final reporting

The Task 46 final report must state explicitly whether the current Palk/Benjamin-derived NKCC1 implementation was retained, modified or replaced, why, which experimental facts discriminate the models, and how the accepted NKCC1 reconstruction changes the excessive compensation identified by Task 44.
