# Task 48 protocol audit

Advisory worker only. Frozen input: `5c945c4a7857657b5f7c76d5bf9e3d3abc7d5721`. No model modification, inference, commit or push was performed.

The inherited model already contains AE2 and an explicit AE4 activation state driven by beta input. Several assay protocols cannot be represented exactly by the inherited equations and bath coordinates. A programme which silently changes their bath values and calls the resulting output the source assay would be scientifically misleading.

## Primary source and experimental settings

[Peña Münzenmayer et al. (2015), PMC4409235](https://pmc.ncbi.nlm.nih.gov/articles/PMC4409235/), Methods, Table 1 and Figures 1 to 5. Local primary HTML is available at `task48_advisory/primary_data/jbc.html`. Source recipes below retain their chemical species rather than replacing omitted ions with chloride.

All listed concentrations are mM. All experiments were at 37 °C. Bicarbonate solutions were gassed with 95% O2 and 5% CO2 for at least 30 minutes; bicarbonate absent solutions used 100% O2. pH was 7.4.

| Bath | NaCl | Sodium gluconate | NaHCO3 | KCl | Potassium gluconate | CaCl2 | MgCl2 | Glucose | HEPES |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| B+ high Cl | 120 | 0 | 25 | 4.3 | 0 | 1 | 1 | 5 | 10 |
| B+ low Cl | 0 | 120 | 25 | 0 | 4.3 | 1 | 1 | 5 | 10 |
| B− high Cl | 145 | 0 | 0 | 4.3 | 0 | 1 | 1 | 5 | 10 |
| B− low Cl | 0 | 145 | 0 | 0 | 4.3 | 1 | 1 | 5 | 10 |

Thus B+ high Cl has Na 145, K 4.3 and Cl 128.3; B− high Cl has Cl 153.3. Both low Cl baths retain Cl 4 from the divalent salts. The Na absent B+ high bath replaces NaCl with NMDG chloride and NaHCO3 with choline bicarbonate. The source describes the corresponding low Cl bath using NMDG glutamate; its sentence refers to replacing NaCl although the low Cl parent contains sodium gluconate. Preserve that source wording discrepancy rather than silently asserting a corrected recipe.

## Protocol map

| Protocol | Genotypes | Documented intervention sequence | Direct implementation and present limitation |
| --- | --- | --- | --- |
| Rest | WT, established AE4 KO, established AE2 KO | B+ high Cl, no agonist | Set `StimulusArm.REST`; solve each genotype with only the relevant expression deleted. Source B+ bath cannot yet be represented exactly. |
| CCh | WT, AE4 KO; AE2 control and KO | B+ high Cl, 0.3 µM CCh | `CCH_ONLY`; calcium drive changes while beta remains zero. Frozen calcium step is an effective input, not a measured calcium trace. |
| IPR | WT, AE4 KO; AE2 control and KO | B+ high Cl, 5 µM IPR | `IPR_ONLY`; beta input changes while calcium remains basal. Do not manufacture an initial Cl exit: source describes none in this arm. |
| CCh plus IPR | WT, AE4 KO; AE2 control and KO | B+ high Cl, 0.3 µM CCh plus 5 µM IPR | `CCH_IPR`; gland stimulation lasts 600 s and collection is at 60 s intervals. Cell fluorescence and gland collection use different observation operators. |
| NKCC isolation | Pooled matched controls, AE4 KO, AE2 KO | B−, EZA 30 µM, T16Ainh A01 50 µM; low Cl depletion then high Cl readdition | Tmem16A block maps to `g_cl_apical_S=0`. EZA has no target in equilibrium carbonate chemistry. Exact bicarbonate absence is outside positive concentration domain. Active NKCC law lacks extracellular substrate dependence. |
| NKCC plus bumetanide | Control | Same as above; bumetanide 50 µM added after initial depletion | Set NKCC activity to zero only in the inhibitor stages. A 100% block is an idealisation; do not set it to 95.6% solely to reproduce observed inhibition. |
| NHE response | Control, AE4 KO, AE2 KO | CCh 0.3 µM plus IPR 5 µM, BCECF measurement | Compare pH trajectory or initial alkalinisation via a documented operator, not raw NHE flux. EIPA 10 µM maps to NHE activity zero. This is an acute inhibitor after genotype rest, not an NHE KO resting state. |
| Isolated exchanger | Control, AE4 KO, AE2 KO, double KO | B+ low Cl depletion; then bumetanide 50 µM plus T16Ainh A01 50 µM during low Cl; then B+ high Cl with both inhibitors retained | Delete only intended pathways; both AE2 and AE4 remain active in WT. Exact low Cl and spectator ions require expanded bath bookkeeping. |
| Bicarbonate dependence of exchanger | AE4 KO, AE2 KO, control | Corresponding depletion and readdition in B− plus EZA 30 µM and the two blockers | EZA and absence of inorganic carbon cannot be represented by changing unrelated capacities. This is blocked in inherited architecture. |
| Isolated exchanger with IPR | Control and AE4 KO | Exchanger isolation as above, with or without IPR | Use existing AE4 beta driven activation. The numerical time of IPR preincubation and individual slope windows are not reported in text; do not infer an activation time constant from the uptake slope. |
| Isolated exchanger Na dependence | AE2 KO | B+ isolation with normal Na versus substituted Na absent solutions | Existing extracellular Na domain excludes zero. NMDG/choline charges are omitted. Substitution is not AE4 deletion and does not alter K. |

The low Cl and inhibitor stage lengths are illustrated by figure traces but not numerically specified in the text. Without a checked figure based timing record, retain them as unknown. The word initial refers to the regression segment drawn in each figure, not universally to a derivative at stimulus time zero. A programme may expose stage lengths as required protocol inputs, but must not give fabricated defaults and then use its results for inference.

## Exact inherited code lineage

`analysis/47_dynamic_potassium_recycling_reconstruction/common.py:load_model` reconstructs the frozen object from `input/active_parameters.json`. The actual scientific model is:

`MinimalNbcModel(StimulatedNkcc1Model(ModernFullModel(...)))`.

Its AE4 evaluator is `equal_routing_adapter`, retaining the Task 40 shared carrier flux and assigning half the transported cation to Na and half to K. It uses `Nkcc1Kinetics(law='palk_benjamin_2010_eq17')`, Cha 2009 NHE and the existing AE2 law. The older `model/equations.md` file is a historical equation skeleton, not the implementation being reconstructed.

`CompositeAe4NkccRegulation` contains `R1EffectiveActivation` with frozen `tau_activation_s=30`, basal multiplier 1, fully activated increment 0.25 and coupling scale 1. Its state obeys `da/dt=(beta-a)/tau`; AE4 capacity is multiplied by `1+0.25 a`. The current model therefore already represents the dependency which Task 48 permits adding only if absent. The gain is not yet identified from the current joint assay data.

AE2 is implemented in `membranes.evaluate_homeostasis` using a reversible tanh of the Cl/HCO3 log affinity, capacity `0.005 fmol/s`. `Genotype.ae2_expression=0` removes its sources. AE4 KO analogously removes the AE4 sources. A double deletion should be represented by both expression flags, never by reassigning capacities of surviving proteins.

## Structural representability blockers

1. **Exact bath composition.** `ModernFullModel._validate_reference_electroneutrality` requires bath Na + K − Cl − total alkalinity = 0. The source recipes contain Ca, Mg, gluconate, HEPES and replacement cations which that relation does not track. Merely replacing source Cl by a value derived from the incomplete relation would alter the experiment. The frozen bath is Na 145, K 5, Cl 126.16159, TIC 25 with no HEPES, whereas the source B+ recipe uses K 4.3 and bicarbonate 25. TIC 25 is not bicarbonate 25.
2. **Zero substrate domains.** `FullModelParameters` requires positive bath Na, Cl and TIC. `evaluate_homeostasis` requires positive HCO3 and Na before computing all log affinities, even for deleted or inhibited transporters. Membrane Nernst expressions, the NBC affinity and AE4 environment also assume positive substrate concentrations. Small positive floors are numerical approximations and require documented limiting calculations, not silent substitution for exact zero.
3. **Carbonic anhydrase inhibition.** `acid_base.py` imposes instantaneous carbonate equilibrium. It has no hydration kinetic state or carbonic anhydrase activity parameter. Setting AE4, AE2 or NBC capacity to zero is not EZA inhibition. Modelling EZA accurately would require an additional chemical kinetic assumption which is beyond a pure protocol wrapper.
4. **NKCC bath changes.** The active `palk_cycle_flux_fmol_s(na_i,k_i,cl_i,...)` has no extracellular arguments. Changing bath Cl from 4 to 153.3 cannot directly change NKCC activity at fixed intracellular state. A generic bath sensitive tanh law exists elsewhere but switching to it changes the inherited equations and must not be called R0.
5. **Acinar assay geometry.** The integrated model has a finite secretory lumen with outflow. Isolated acinar clumps are not a whole gland and their apical solution access is not specified quantitatively. The same basolateral bath command need not reproduce every membrane exposure in an isolated assay. Do not invent a bath clamp of the lumen without stating it as an assay boundary assumption.
6. **Dye mapping and timing.** SPQ slopes and BCECF ratios need observation operators; no single transporter flux equals these measurements. Source regression windows and the intracellular dye calibration cannot be supplied by optimisation.
7. **Residual uptake after double deletion.** With ideal complete blockade of NKCC and Tmem16A and deletion of AE2 plus AE4, the represented cellular chloride sources vanish. Chloride concentration can still change through volume; this must not be labelled residual exchanger activity. The observed remaining uptake does not uniquely identify a new transporter, since incomplete block and unrepresented pathways are alternatives.

These are representational limitations, not biochemical claims that Na absent or bicarbonate absent experiments are impossible in real cells.

### Permissible positive B+ bath bookkeeping

An external spectator bookkeeping adapter can reproduce the listed transported ions of positive B+ without changing intracellular balances. Set Na 145, K 4.3, Cl 128.3, HCO3 25 and pH 7.4; convert bicarbonate to TIC using the inherited carbonate fractions. Represent Ca and Mg as spectator charge +4 mEq/L and assign Ca, Mg, glucose and HEPES a combined spectator osmolarity of 17 mOsm/L. If HEPES is treated as an inert external osmolyte under fixed bath pH, explicitly record that effective reservoir assumption. The residual charge beyond Ca/Mg arises from carbonate speciation and unreported titration counterions. It is a bookkeeping closure, not an experimentally measured spectator concentration. Source titrant quantities, HEPES protonation and their counterions are not reported, so this is an effective representation of source B+ with exact listed transported ions, not a complete chemical reconstruction. It is defensible for standard positive bath diagnostics and does not fix the blocked isolation assays.

### Source figure timing check

Figures 1 to 3 were visually inspected from the primary worker's saved JPGs. Figure 1 CCh plus IPR cell traces show a short baseline, rapid exit and subsequent reuptake, with separate dashed regression segments. Figure 2 NKCC low Cl bar is approximately 60 to 250 s; Figure 3 exchanger low Cl bar is approximately 60 to 360 s, with blockers introduced approximately 280 s before readdition. These are visual orientation estimates only: they are not calibrated digitisation, have no established timing error estimate and must not become authoritative protocol durations. For quantitative use save axis calibration, measured bar positions and timing uncertainty first.

## Minimal executable layer design

Keep the scientific parameter vector immutable and shared. Build a `ProtocolSpec` with an ordered tuple of `StageSpec` objects. Each stage should carry source bath recipe, stimulus arm, exact inhibitor targets, stage duration and timing provenance. Every protocol declares `representable`, `partial` or `blocked`, with machine readable reasons. A blocked protocol must produce an explicit unavailable prediction rather than a numerical surrogate entering the objective.

The model builder must rebuild the stimulus in both `ModernFullModel` and `CompositeAe4NkccRegulation`; changing only the outer stimulus risks stale NKCC regulation. The AE4 beta input and calcium drive must refer to the same stage clock. State vectors remain continuous between stages, including AE4 activation, while extracellular inputs and acute inhibitor activities change only at recorded boundaries.

For each parameter vector, solve rest independently for WT, AE4 KO and AE2 KO under the resting protocol. The previous WT rest is a numerical starting guess only. Preserve total charge, fixed intracellular material and all shared kinetic constants; do not use fitted genotype specific ion offsets. Store scaled RHS residuals, positivity checks, stability/relaxation evidence, parameter hash and genotype. Root finding does not by itself establish a stable established animal state. Each assay begins at its own genotype rest and subsequently undergoes the documented conditioning stages.

Intervention regression checks should show exact zero of the targeted source, unchanged unrelated parameters, charge consistent source removal, and unchanged state across a bath switch. Check that zero dose or absent CCh/IPR arms keep their corresponding inputs basal. Test the rest cache key includes parameter hash, genotype, bath and stimulus. None of these checks is evidence that the blocked assays were successfully reconstructed.

The honest minimum for Task 48B is an executable inventory and supported intervention layer, plus explicit failures for unrepresentable source baths and assays. Resolving those failures needs a scientific decision by the orchestrator. They must not be hidden by enlarging an optimisation vector.


Orchestrator integration: reviewed and accepted for the 48A evidence freeze. This document supplies evidence and restrictions, not a completed fit.
