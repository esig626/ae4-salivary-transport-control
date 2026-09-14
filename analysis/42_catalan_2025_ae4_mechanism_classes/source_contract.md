# Task 42 source contract: Catalán et al. 2025 AE4 mechanism classes

## Purpose

Task 42 asks whether transport stoichiometries explicitly permitted by Catalán et al. 2025 can change the AE4 loss phenotype on the completed Task 40 model, without the target selected AE4 dependent CaCC recruitment introduced in Task 41.

This task is a source constrained mechanism discrimination exercise. It is not a search over arbitrary cation fractions and it must not fit the approximately 35% secretion phenotype.

## Immutable repository parents

The Task 42 branch was created directly from current `main` at:

`fe51473d7cc201e2e04007aa62fbc1b45098732e`

The completed Task 40 scientific result is already contained in this lineage. Task 41 is preserved separately in draft PR #33 at publication commit:

`86f135e075aa91a70a752b650d7bb894353fee04`

Task 41 must not be merged, cherry picked, force pushed, rewritten or modified in Task 42. Its target selected CaCC recruitment law is forbidden in Task 42. Task 41 may be used only as a post freeze comparison result.

## Primary source

Catalán MA et al. Molecular determinants of HCO3− and cation transport in the human cation dependent Cl−/HCO3− exchanger AE4. American Journal of Physiology Cell Physiology. 2025;328:C2070–C2084. DOI: 10.1152/ajpcell.00346.2024.

The source establishes or proposes the following points relevant here.

1. Human AE4 supports bicarbonate transport in the presence of either Na+ or K+ as the principal extracellular monovalent cation, while replacement by impermeant NMDG+ abolishes transport.
2. Na+ and K+ do not appear to use identical coordination. The T756A T448I double mutant abolishes transport in Na+ but retains substantial transport in K+.
3. The authors propose that extracellular Cl− may be stabilised by extracellular Na+ and that both may move inward, while intracellular HCO3−, or CO3^2−, may be stabilised with intracellular K+ and move outward.
4. For that proposed NaCl inward and K base outward cycle, the explicitly proposed electroneutral stoichiometries are:
   * 1:1:1:1 for Cl− : Na+ : HCO3− : K+;
   * 1:1:1:2 for Cl− : Na+ : CO3^2− : K+.
5. The authors also propose an alternative in which Cl− is stabilised by the positive charge of K879 rather than requiring cotransported Na+ at the external site. For that class they give:
   * 1:1:2 for Cl− : cation : HCO3−;
   * 1:1:1 for Cl− : cation : CO3^2−.
6. The paper explicitly states that further experiments are required. These are proposed mechanism classes, not experimentally resolved AE4 stoichiometries.

## Task 40 control

Task 40 retains the inherited total AE4 cycle `J4`, positive for inward Cl−, but replaces donor concentration weighted cation routing with an equal split. Its cellular source vector per positive cycle is:

`Na = -0.5 J4`

`K = -0.5 J4`

`Cl = +1 J4`

`TIC = -2 J4`

`TA = -2 J4`

This is the legacy 1 Cl inward : 1 monovalent cation outward : 2 HCO3 outward class with an imposed 50:50 Na/K source split. It is the control, not a Catalán 2025 conclusion.

## Source vector classes for the first discrimination stage

The first stage must keep the inherited Task 40 scalar AE4 cycle law `J4(state)` unchanged and alter only the transported species source vector. This isolates stoichiometry and direction from unknown class specific kinetics.

### C0: Task 40 legacy control

`Na = -0.5 J4`

`K = -0.5 J4`

`Cl = +J4`

`TIC = -2 J4`

`TA = -2 J4`

### C1: NaCl inward, K plus HCO3 outward

Catalán proposed 1:1:1:1 Cl : Na : HCO3 : K.

`Na = +J4`

`K = -J4`

`Cl = +J4`

`TIC = -J4`

`TA = -J4`

### C2: NaCl inward, 2K plus CO3 outward

Catalán proposed 1:1:1:2 Cl : Na : CO3 : K.

`Na = +J4`

`K = -2 J4`

`Cl = +J4`

`TIC = -J4`

`TA = -2 J4`

The existing model does not track carbonate as an independent species. One outward CO3^2− removes one unit of total inorganic carbon and two equivalents of total alkalinity. Do not silently represent it as bicarbonate.

### C3a and C3b: K879 stabilised Cl, cation plus 2 HCO3 outward

Catalán proposed 1:1:2 Cl : cation : HCO3. Because the paper supports both Na+ and K+ dependent transport but does not resolve which cation participates in this specific whole cycle in salivary acinar cells, use the two pure endpoints as predeclared bracketing subvariants. Do not fit or sweep a mixture.

C3a, K coupled:

`Na = 0`

`K = -J4`

`Cl = +J4`

`TIC = -2 J4`

`TA = -2 J4`

C3b, Na coupled:

`Na = -J4`

`K = 0`

`Cl = +J4`

`TIC = -2 J4`

`TA = -2 J4`

### C4a and C4b: K879 stabilised Cl, cation plus CO3 outward

Catalán proposed 1:1:1 Cl : cation : CO3. Again use pure K and pure Na endpoints only.

C4a, K coupled:

`Na = 0`

`K = -J4`

`Cl = +J4`

`TIC = -J4`

`TA = -2 J4`

C4b, Na coupled:

`Na = -J4`

`K = 0`

`Cl = +J4`

`TIC = -J4`

`TA = -2 J4`

## Required algebraic checks

Before any resting solve or production trajectory:

1. Verify each source vector is electroneutral in charge equivalents.
2. Verify reversal of `J4` reverses every transported source consistently.
3. Derive the Na, K, Cl, TIC and TA whole cell balance identities for each class.
4. Compute the chemical gradient contribution to the forward transport affinity at the accepted Task 40 WT state for each proposed stoichiometry. Include membrane electrical work explicitly and show its cancellation only where electroneutrality justifies it.
5. For carbonate classes, derive CO3^2− consistently from the model acid base state. Do not substitute HCO3− for CO3^2−.
6. State clearly that keeping `J4(state)` fixed is an intentional first stage assumption because Catalán 2025 proposes stoichiometries and coordination mechanisms but does not provide a whole cell kinetic law suitable for direct insertion into this salivary model.

## Scientific firewall

The approximately 35% experimental AE4 loss secretion deficit is held out during mechanism implementation, WT assessment and genotype execution.

Do not select classes, cation identity, capacities, equilibrium constants, stimulus parameters, roots or solver settings based on closeness to the phenotype.

All predeclared classes that pass their WT gates must be run through the same genotype protocol. Freeze and commit all class predictions before comparing them with the experimental secretion phenotype or with Task 41.

## Interpretation limits

If none of these classes recovers the phenotype, the allowed conclusion is:

`CATALAN 2025 STOICHIOMETRIC SOURCE CLASSES ALONE DO NOT RECOVER THE AE4 LOSS PHENOTYPE UNDER THE INHERITED TASK 40 CYCLE LAW.`

That does not rule out the Catalán mechanisms themselves, because class specific kinetic laws and regulation remain unresolved experimentally.

If one or more classes does recover the phenotype without target fitting, preserve every other class and report the result as a source constrained prediction requiring independent validation.
