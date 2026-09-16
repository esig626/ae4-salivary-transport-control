# Task 49 exclusion ledger: do not repeat completed work

This file exists because the repository, not conversational recollection, is the scientific memory of the project. Every item below has already been executed sufficiently to exclude it from Task 49 unless Task 49 discovers a literal implementation defect in the cited work.

## AE4-only reconstruction

### Task 12: `analysis/12_ae4_mechanism_reconstruction/`

Controlling result: `AE4 ALONE INSUFFICIENT; MISSING PATHWAY LOCALIZED`.

The C1-C8 AE4 candidate family, reversible mixed-cation variants, saturation/cooperativity variants, source-compatible moving roots, positive PKA folds and genotype resting-state tests did not produce a source-admissible AE4-only explanation. The AE4-null vector field is independent of AE4 candidate kinetics and already showed the chronic low-Cl/pH failure.

**Do not:** reopen an AE4-only kinetic family, Hill/cooperativity search, arbitrary cation selectivity search or another source-vector grid.

## Calcium amplitude

### Task 13C: `analysis/13C_calcium_fast_screen/`

Increasing the calcium stimulus through the frozen 0.10/0.25/0.50 uM panel did not fix the WT secretion problem and is not an AE4-loss mechanism.

**Do not:** alter calcium amplitude as the Task 49 repair.

## Pump coupling

### Task 19: `analysis/19_ae4_loss_nak_pump_coupling/`

At 50% pump capacity, pump inhibition itself reduced secretion strongly but AE4-low cells were relatively less affected; the interaction opposed the proposed harmful AE4-pump mechanism. The work does not justify a new AE4-pump coupling law.

**Do not:** use pump capacity reduction or AE4-dependent pump scaling.

## NHE1

### Tasks 30-31, especially `analysis/31_nhe1_mechanistic_repair/`

The Cha mechanistic NHE1 model was reconstructed and tested. It did not resolve the genotype resting-state problem; exact-null resting closure remained problematic and WT chloride remained wrong.

**Do not:** replace NHE1 again or tune NHE capacity to the AE4 phenotype.

## Stimulated bicarbonate pathway

### Task 36: `analysis/36_minimal_nahco3_alkalinity/`

A stimulus-recruited electrogenic 1 Na:2 HCO3 pathway was introduced as an NBCe-like structural surrogate to supply alkalinity for productive stimulated AE4 flux.

Critical provenance: the Task 36 design explicitly states that the implemented NBC is **the stimulus-recruited increment, not a claim that native basal NBC transport is absent**. Its activation is exactly zero at rest solely to nest the accepted Task 31 rest.

### Task 37: `analysis/37_wt_nbc_validation/`

The stimulated NBC architecture produced a valid WT trajectory. Positive Cl loading was about 79.9% NKCC1 / 20.1% AE4 / 0% AE2; the NKCC fraction remained above the 65-75% context band.

**Task 49 rule:** do not alter NBC. The zero-basal-NBC construction is retained as an explicitly unresolved resting-homeostasis issue for a later task if necessary. Do not pretend it was experimentally established.

## Acute AE4 loss and NKCC compensation

### Task 38: `analysis/38_ae4_perturbation_validation/`

Acute/shared-rest AE4 deletion produced only a small secretion deficit because NKCC compensation was large.

### Task 39: `analysis/39_palk_nkcc1_full_validation/`

The source-fixed Palk/Benjamin NKCC1 concentration law was implemented and validated. It did not suppress compensation; compensation became at least as problematic. The active reduced law is bath-conditioned and lacks extracellular substrate arguments.

**Do not:** replace the NKCC core again, add a hard NKCC cap, or claim the isolated low-Cl NKCC assay is represented by the current core.

Task 49 may add the independently documented **beta/cAMP regulatory input** to the existing NKCC activity multiplier. That is a missing stimulus input, not another NKCC kinetic-law replacement.

## AE4 cation routing

### Task 40: `analysis/40_ae4_equal_cation_routing/`

Equal Na/K routing produced the frozen modern reference but only about a 3.86% acute null deficit.

**Do not:** reopen the Na/K split as the phenotype mechanism.

## Target-selected apical conductance construction

### Task 41: `analysis/41_ae4_loss_algebraic_design/`

An explicit AE4-expression-dependent reduction in stimulated TMEM16A/CaCC recruitment was deliberately selected against the desired secretion range and produced about a 30.26% null deficit. It was not independent validation. Its temporal pattern was wrong, mutant chloride was high, pH/volume were abnormal, and the specific AE4-CaCC coupling is unsupported.

**Do not:** resurrect AE4-dependent CaCC recruitment.

Task 49 is different: independent adult mouse SMG experiments demonstrate that beta/cAMP secretion uses a distinct TMEM16A-independent apical anion pathway. The new branch must have no direct AE4-expression multiplier.

## AE4 stoichiometric source classes

### Task 42: `analysis/42_catalan_2025_ae4_mechanism_classes/`

Seven predeclared Catalan-2025-compatible source classes failed to recover the phenotype; valid classes gave at most a small loss or even increased secretion, and carbonate classes damaged WT pH.

**Do not:** test another AE4 source-vector class in Task 49.

## CBM structural screening

### Task 46: `analysis/46_physiology_constrained_model_reconstruction/`

Conservation-explicit CBM showed full reference chloride export remains structurally feasible after AE4 deletion. Resource/coupling analysis did not supply a defensible numeric cap that blocks rescue.

**Do not:** return to generic CBM resource caps.

## Potassium-recycling ceiling

### Task 47: `analysis/47_dynamic_potassium_recycling_reconstruction/`

The frozen dynamic model carried the additional K burden through increased K-channel efflux, storage and voltage changes. K efflux increased about 18.1% while pump turnover increased only about 1.86%. No binding K/pump recycling ceiling emerged.

**Do not:** reduce K conductance or pump capacity to force the phenotype.

## Joint experimental reconstruction

### Task 48: `analysis/48_joint_experimental_constraint_reconstruction/`

The full data/protocol audit exposed an exact structural defect: in AE4 knockout, beta input disappears from the core equations, so CCh and CCh+IPR trajectories are identical for every admissible shared parameter vector. The source reports AE4-KO uptake `2.30 +/- 0.10` under CCh and `0.90 +/- 0.09` under CCh+IPR. The contrast has zero sensitivity to all 33 provenance-eligible parameter directions.

Genotype-specific resting solves also fail badly: AE4 KO remains too chloride-rich and becomes strongly alkaline; secretion is slightly increased rather than reduced. Several isolated assays remain unrepresentable because of bath-conditioned NKCC and instantaneous carbonate chemistry.

**Do not:** try to parameter-fit around the exact beta-invariance obstruction.

## Task 49 is allowed to test exactly one previously unrepresented architecture

Repository search at the Task 48 final state found no implementation or prior test under `VRAC`, `LRRC8`, `DCPIB` or `volume-regulated anion channel`.

Independent adult mouse SMG physiology shows that beta/cAMP secretion is TMEM16A-independent, CFTR-independent and ClC-2-independent, is strongly sensitive to VRAC blockers, and has a slow swelling-associated time course.

Therefore Task 49 is restricted to reconstructing that missing cAMP branch, with the existing Ca/TMEM16A branch retained, plus the independently documented beta/cAMP input to the existing NKCC regulatory multiplier.

If this source-backed branch fails, record the failure. Do not substitute another mechanism in the same task.
