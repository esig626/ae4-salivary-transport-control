# Mandatory research ledger: AE4 salivary transport project

This file is the authoritative cumulative scientific memory for this repository.

**Every scientific task, Codex run, agent, mechanism proposal, model edit, numerical solve, fit, or interpretation decision must read this file before acting.** The repository, not conversational memory, is the project memory.

A scientific checkpoint is not complete unless this ledger is updated in the same dependency boundary when the checkpoint establishes a new result, failure, exclusion, correction, or interpretation rule. Historical entries are never silently rewritten: later work may supersede or qualify an entry, but the earlier result must remain traceable.

## 0. Current project state and remaining budget

- Canonical repository: `esig626/ae4-salivary-transport-control`.
- Canonical `main` before Task 50 staging includes merged Tasks 47-48 plus the repository-level research guardrails.
- Task 49 is complete on `analysis/task-49-camp-vrac-secretory-branch-reconstruction`, final closeout `ab7f1cb054b8c8150343653727773f1ca59ebbcf`. Its scientific conclusions are part of the cumulative ledger even though its branch is not the production parent for Task 50.
- **Three funded/scientific Codex execution shots remain.** Do not spend one rediscovering anything in this ledger.
- The primary modelling target is the magnitude and direction of the AE4-loss secretion phenotype, not exact experimental timing. `docs/PHENOTYPE_TARGET_CONVENTION.md` is binding.

Planned budget, unless new evidence forces a revision:

1. **Shot 1 / Task 50:** one-parameter beta-conditioned effective secretory coupling, no search.
2. **Shot 2 / Task 51:** robustness/uncertainty of the accepted reduced model only; no new mechanism family unless Task 50 falsifies the algebraic construction.
3. **Shot 3 / Task 52:** final reproducibility, uncertainty, article-facing lock; no mechanism fishing.

## 1. Source-backed experimental facts that constrain every model

The primary phenotype dataset is Peña-Münzenmayer et al., JBC 2015, DOI `10.1074/jbc.M114.612895`, PMCID `PMC4409235`. Values below were source-checked in Task 48.

### 1.1 AE4 line: resting state and chloride uptake

| Quantity | WT | AE4 KO | Notes |
| --- | ---: | ---: | --- |
| Resting intracellular Cl | `50.10 +/- 1.50 mM`, n=7 | `36.50 +/- 1.60 mM`, n=6 | Established-genotype measurement |
| Resting pH | `6.91 +/- 0.07`, n=4 | `6.89 +/- 0.02`, n=4 | pH approximately preserved despite large Cl loss |
| Initial Cl exit, CCh+IPR | `-0.21 +/- 0.01`, n=10 | `-0.18 +/- 0.02`, n=9 | SPQ-derived signal |
| Initial Cl exit, CCh | `-0.18 +/- 0.03`, n=7 | `-0.15 +/- 0.01`, n=6 | SPQ-derived signal |
| Initial Cl uptake, CCh+IPR | `2.02 +/- 0.10 x10^-3 s^-1`, n=10 | `0.90 +/- 0.09 x10^-3 s^-1`, n=9 | Strong KO reduction |
| Initial Cl uptake, CCh | `2.18 +/- 0.20 x10^-3 s^-1`, n=7 | `2.30 +/- 0.10 x10^-3 s^-1`, n=6 | Near-preserved KO response |
| Initial Cl uptake, IPR | `0.40 +/- 0.07 x10^-3 s^-1`, n=4 | `0.20 +/- 0.03 x10^-3 s^-1`, n=7 | Positive response remains in KO |

The exact SPQ observation map is not fully identified. Do not silently equate transporter flux with fluorescence slope. Concentration derivatives must include dilution: `(dn_Cl/dt - [Cl] dV/dt)/V`.

### 1.2 Whole-gland secretion phenotype

- Under CCh+IPR, AE4 KO total saliva over 10 min is reported approximately **35 +/- 4.7% lower than WT**, n=6 per genotype.
- The source describes relatively preserved early secretion and a larger sustained deficit later. **This timing pattern is qualitative context, not a hard fitting target.** Do not reject a model because it does not reproduce minute-by-minute timing.
- The model target is a substantial, robust reduction of the same order under the stimulated condition, not an exact `35.000%` or exact 600 s trace.

### 1.3 AE2 control line

| Quantity | Control | AE2 KO |
| --- | ---: | ---: |
| Resting Cl | `53.40 +/- 1.80 mM`, n=6 | `54.50 +/- 1.80 mM`, n=6 |
| Resting pH | `6.87 +/- 0.01`, n=4 | `6.95 +/- 0.05`, n=6 |
| CCh+IPR uptake | `2.11 +/- 0.50 x10^-3 s^-1`, n=13 | `2.39 +/- 0.20 x10^-3 s^-1`, n=18 |
| CCh uptake | `2.20 +/- 0.20 x10^-3 s^-1`, n=8 | `2.30 +/- 0.20 x10^-3 s^-1`, n=9 |
| IPR uptake | `0.70 +/- 0.10 x10^-3 s^-1`, n=7 | `0.58 +/- 0.05 x10^-3 s^-1`, n=6 |

AE2 deletion is essentially secretion-neutral in the primary whole-gland experiment. Any accepted AE4 mechanism should preserve this specificity unless new evidence says otherwise.

### 1.4 Isolated NKCC and NHE assays

- Isolated NKCC assay: HCO3-free, `30 uM` ethoxyzolamide, `50 uM` T16Ainh-A01, chloride depletion/readdition. No statistically detectable genotype difference among control, AE4 KO and AE2 KO; n=16/8/10. `50 uM` bumetanide reduces uptake by `95.6 +/- 0.9%`, n=9.
- **Do not convert a nonsignificant difference into exact equality or an invented tolerance.**
- The current Palk-reduced NKCC law lacks extracellular substrate dependence and cannot faithfully reproduce the depletion/readdition assay. This assay therefore remains a qualitative/protocol-specific constraint unless a future task supplies a faithful observation model.
- NHE assay: `0.3 uM` CCh + `5 uM` IPR, BCECF. No detectable genotype difference; n=13/8/6. `10 uM` EIPA reverses alkalinisation to acidification, n=7.
- These assays constrain gross genotype-specific regulation but do not justify arbitrary flux equality under normal whole-cell stimulation.

### 1.5 Exchanger and beta/PKA evidence

- Isolated exchanger, NKCC/TMEM16A inhibited: AE4 KO reduces uptake `38.2 +/- 6.6%`; AE2 KO `55.1 +/- 3.4%`; double KO `87.7 +/- 3.1%`.
- IPR increases isolated exchanger uptake in control from `1.56 +/- 0.19` to `4.33 +/- 0.91 x10^-3 s^-1`, p=.016; AE4 KO shows `0.95 +/- 0.10` to `1.08 +/- 0.09`, p=.366.
- Peña-Münzenmayer et al. AJP GI 2021, DOI `10.1152/ajpgi.00145.2021`, establishes beta/cAMP/PKA activation of AE4 and S173 dependence. It does **not** identify acute signalling time constants. Treat signalling as an effective multiplier unless dynamics are independently necessary.
- The 2021 source contains two reasonable gain readings already preserved by the model: approximately `+25%` from prose and approximately `+70%` from visual figure interpretation. They are sensitivities, not exact native gland constants.

### 1.6 AE4 biophysics and stoichiometry limits

Peña-Münzenmayer et al. JGP 2016, DOI `10.1085/jgp.201611571`:

- AE4 is electroneutral and associated with Cl/HCO3/Na/K transport.
- External Na EC50 approximately 49 mM, Hill about 2.0; K EC50 approximately 62 mM, Hill about 1.8.
- Na and K both participate.
- The model's equal Na/K source routing is a frozen effective allocation, **not a directly measured microscopic stoichiometry**.
- A working `1 Cl : 1 monovalent cation : 2 HCO3` cycle is an inferred/effective representation, not a direct stoichiometric measurement.

## 2. Current production-model architecture and important parameter status

The modern model uses twelve conserved core states: intracellular Na, K, Cl, total inorganic carbon (TIC), total alkalinity (TA), cell volume; luminal Na, K, Cl, TIC, TA, lumen volume. Membrane voltages are algebraic current-closure variables.

Current accepted/frozen ingredients inherited by the Task 40/41 lineage:

- NKCC1: reversible source-fixed Palk/Benjamin concentration-response core; exact `1 Na : 1 K : 2 Cl`; stimulation represented by an effective scalar activity multiplier. Current `N1_CCH_ALGEBRAIC` multiplier is `1` at rest and `1.75` at the declared stimulated endpoint. The gain is an effective WT construction, not a measured universal constant.
- NHE1: mechanistic Cha-family law from Task 31; stimulated multiplier inherited from the accepted construction.
- NBC-like pathway: electrogenic `1 Na : 2 HCO3`, introduced in Task 36/37 to provide productive alkalinity supply during stimulation. Its zero-rest recruitment was a nesting choice, not evidence that basal NBC is absent.
- AE4: current Task 40 source routing is `(-A/2, -A/2, +A, -2A, -2A)` in Na/K/Cl/TIC/TA order, with the inherited scalar AE4 cycle law and effective beta/PKA capacity regulation.
- CaCC/TMEM16A: apical chloride conductance driven by the existing calcium activation.
- Pump/K channels, paracellular paths, CO2 exchange, water laws and geometry are inherited from the modern reconstruction unless a task explicitly states otherwise.

Task 40 WT resting checkpoint values include approximately Na `11.505 mM`, K `116.974 mM`, Cl `60.572 mM`, pH `6.91098`, volume `1.4533 pL`. These are valid model equilibria but do not match every experimental WT observable exactly.

## 3. Exact algebraic facts that must guide future reasoning

### 3.1 Task 41 transient and stationary identity

With equal AE4 routing, NKCC cycles `N`, AE4 cycles `A`, NBC cycles `B`, NHE1 flux `H`, AE2 chloride-in flux `E`, pump cycles `P`, and apical chloride export `J_CaCC`:

`dNa_i/dt = N + H + B - A/2 - 3P`

`dTA_i/dt = H + 2B - 2A - E`

`dCl_i/dt = 2N + A + E - J_CaCC`

Eliminating `N`, `B`, and `E` gives the exact transient identity:

`J_CaCC = 6P - H + 2 dNa_i/dt - dTA_i/dt - dCl_i/dt`.

At stationary state:

`J_CaCC = 6P - H`.

**AE4 disappears explicitly from the sustained export identity.** This is why the network can compensate strongly after AE4 removal. Changing AE4 stoichiometry/capacity alone is not expected to create a large sustained secretion phenotype unless it changes the rest of the network enough to alter pump/NHE/export support.

### 3.2 CBM structural identities from Task 46

For the zero-storage coarse balance model:

- `Db + Da = H`.
- `E = H + 2B - 2A`.
- `J = 2N + H + 2B - A`.
- `P = (H + J)/6`.
- `Kout = N - A/2 + 2P`.
- Under AE4 KO, `Kout = J/3 + N + H/3 >= J/3`.

Task 46 showed full reference chloride export is structurally feasible without AE4 under loose shared capacities. Therefore a generic structural inability to secrete without AE4 is false; the phenotype must come from regulation/effective coupling or quantitatively constrained capacities, not conservation alone.

## 4. Chronological result ledger and exclusions

### R00-R11: inventory, reproduction, identifiability and forensic lineage

- Tasks 00-08 established baseline reproduction, knockout behaviour, dimensionless/Jacobian analyses, mechanism variants, identifiability and robustness on earlier chassis versions.
- Task 10 formalised local identifiability and exact stoichiometric reductions. The printed historical Table 1 row is not simultaneously steady under all printed equations; Palk unit corrections matter.
- Task 11 reconstructed the historical MATLAB lineage. The archived snapshot shares ancestry with the publication but is `UNLIKELY` to be the exact publication-generating implementation. Do not treat old archive code as an authoritative hidden solution.

### R12: AE4-only mechanism reconstruction

`analysis/12_ae4_mechanism_reconstruction/`

- C1-C8 tested; no source-admissible AE4-only mechanism reproduced the phenotype.
- The AE4-null vector field on that chassis predicted resting `Cl_i=48.7822 mM` and `pH_i=7.3785`, already demonstrating the chronic genotype-rest failure.
- Rest-compatible `2:1:3` variants still predicted AE4 KO/WT secretion >1 in the screened domain.
- Fixed-chassis algebra showed a large target secretion loss requires changes beyond AE4 source coefficients; however its pump/K localisation was a frozen-state diagnostic, not proof of a changed pump or K channel.
- Controlling conclusion: **AE4 ALONE INSUFFICIENT** on the tested chassis.

### R13B: modern conservation-explicit full model

`analysis/13B_modern_full_model/`

- Established the modern 12-state conservation architecture, explicit carbon/alkalinity, finite compartments, current closure, mixed-cation AE4 and effective beta regulation.
- The full model remained weakly identified in absolute scale and several transporter/regulatory parameters.
- The exact positive-swelling `Vbeta`/VRAC-style gate was already explored diagnostically here. Do not claim Task 49 invented it.
- No validated minimal AE4 phenotype mechanism emerged.

### R13C: calcium amplitude

`analysis/13C_calcium_fast_screen/`

- Higher calcium amplitudes did not repair the central model/phenotype problem. Do not reopen calcium-amplitude sweeps.

### R16-R18: WT chloride allocation/rebalancing

- Task 16 tried increasing AE4's WT positive chloride-loading share to 10% or 30% by capacity scaling while conserving the positive loading pool. Modified endpoints became physiologically inadmissible before genotype testing; inherited baseline AE4 near-loss increased secretion.
- Task 17 proved the fixed-WT repartition was incompatible with existing source/capacity bounds, especially K handling, at the unchanged state.
- Task 18 relaxed uncertain capacities enough to realise 10% and 30% AE4 WT shares, but required extreme parameter shifts: up to hundreds-fold and implausibly large K conductance in some cases. Across modified models AE4 near-loss still almost always increased secretion; the approximately 35% loss was not recovered.
- Therefore simply giving AE4 a larger WT chloride share is not sufficient.

### R19: pump coupling

`analysis/19_ae4_loss_nak_pump_coupling/`

- Available 50% pump perturbation results opposed the proposed harmful AE4-pump interaction. No evidence supports an explicit AE4-pump coupling as the missing mechanism.

### R20: pooled/mixed cation work

- Alternative pooled cation/no-slip handling did not provide the missing phenotype. Later equal routing became the cleaner fixed control. Do not restart routing searches without new direct evidence.

### R28-R31: NHE1 / acid-base repair programme

- Task 28 explicitly targeted excessive NHE proton extrusion, alkalinisation and compensatory chloride loading; it also explored a single WT chloride-allocation hypothesis.
- Task 29 derived the exact AE4-null stationary identities `J_AE2=J_NHE1` and `J_CO2=J_NHE1` and audited sign/capacity feasibility near the measured knockout Cl/pH.
- Task 30 implemented the published 2018 NHE1 law. WT required scalar calibration; the exact-null solve exhausted the numerical budget and produced no accepted null state.
- Task 31 implemented the mechanistic Cha NHE1 law. WT calibration succeeded; exact-null resting closure still failed and candidate null states were too alkaline.
- `analysis/ae4_carbon_nhe1_compensation/` independently varied CO2 permeability up to x10 and NHE1 down to x0.25 across retained roots. Alkalinisation and compensation persisted; no tested condition produced the required secretion decrease.
- **Closed:** generic NHE scaling, another NHE-law replacement, CO2-permeability scaling, or rediscovery that the acute acid-base laws cannot generate the chronic knockout rest.

### R36-R37: minimal NBC-like alkalinity supply

- Task 36 introduced electrogenic `1 Na : 2 HCO3` NBC-like transport because the prior alkalinity balance could not support productive AE4 transport.
- Task 37 validated the WT stimulation architecture. Positive stimulated chloride loading was approximately `79.9% NKCC1 / 20.1% AE4 / 0% AE2` in the accepted reference.
- Resting NBC was set to zero to preserve nesting, not because basal NBC was experimentally proven absent.

### R38: acute AE4 perturbation on frozen WT rest

`analysis/38_ae4_perturbation_validation/`

- AE4=5% cumulative secretion ratio `0.984315` (1.57% deficit).
- AE4=null cumulative ratio `0.981605` (1.84% deficit).
- 60-600 s WT positive chloride loading: NKCC1 `155.723 fmol`, AE4 `39.114 fmol`.
- Null NKCC1 loading rises to `197.604 fmol`, replacing most missing AE4 chloride.
- Conclusion: excessive network compensation almost erases AE4 loss.

### R39: source-fixed Palk NKCC1 core

`analysis/39_palk_nkcc1_full_validation/`

- Replacing the generic NKCC law with the source-fixed Palk/Benjamin core **increased** relative compensation.
- Null NKCC compensation `+33.1013%`.
- Null cumulative secretion deficit only `0.6273%`; endpoint flow exceeds WT.
- The current effective NKCC stimulation multiplier remains `1.75`; it is not a measured universal gain.
- **Closed:** Palk-law replacement as the phenotype solution, arbitrary NKCC hard caps, or assuming state-driven compensation will disappear automatically.

### R40: equal AE4 cation routing control

`analysis/40_ae4_equal_cation_routing/`

- Fixed 50:50 Na/K source routing reduced but did not remove NKCC compensation.
- Null 60-600 s NKCC compensation `+23.16345%`.
- Null cumulative secretion ratio `0.9614251370`, deficit **3.857486%**.
- 5% deficit `3.417578%`.
- At 600 s null Cl remains lower than WT by about `2.196 mM`, but compensation remains large.
- This is the primary parent/control for the Task 41 constructive result.

### R41: target-selected algebraic construction — IMPORTANT CONSTRUCTIVE RESULT

`analysis/41_ae4_loss_algebraic_design/`

Task 41 deliberately used the desired secretion range as a design target. It is **not independent validation**.

The selected construction multiplied stimulated CaCC conductance by an AE4-expression-dependent factor:

`g_CaCC = g_parent * {(1-a) + a[b + (1-b)e_AE4]}`

with selected `b = 0.10511872843288`.

Results:

- WT exactly nests Task 40.
- AE4=5% cumulative secretion deficit: **23.1634%**.
- AE4=null cumulative secretion deficit: **30.2612%**.
- Null instantaneous flow deficit at 600 s: **27.04%**.
- Null NKCC compensation falls to **+2.79%**, compared with +23.16% in Task 40.
- Null CaCC chloride export falls about **26.88%** relative to WT.
- Mutant Cl is retained/high rather than catastrophically depleted; pH/volume remain imperfect.

Interpretation boundary:

- This proves that **one extra network coupling affecting effective apical chloride secretory capacity is sufficient to generate the right order of phenotype**.
- It does **not** prove that AE4 physically recruits TMEM16A/CaCC or that 89.5% of TMEM16A recruitment is molecularly AE4-dependent.
- The parameter is phenotype-calibrated and must be labelled as such.
- Earlier complaints about its exact temporal pattern are no longer hard rejection criteria under `PHENOTYPE_TARGET_CONVENTION.md`.

### R42: Catalán 2025 source/stoichiometric classes

`analysis/42_catalan_2025_ae4_mechanism_classes/`

- Seven predeclared source classes tested without phenotype-informed tuning.
- Largest valid null loss was `7.67%` (C3b); C1 predicted `10.76%` more secretion; C3a approximately neutral/slight increase; carbonate classes failed WT pH.
- Changing AE4 source stoichiometry alone does not recover the phenotype under inherited kinetics.
- Task 42 explicitly concludes that an additional network coupling can create the magnitude (Task 41) but the biological identity remains unresolved.

### R43-R45: parameter provenance / mathematical audit

- These tasks refined provenance, uncertainty and model bookkeeping. Treat unmeasured conductances, routing fractions, regulation gains and scale factors as model assumptions/sensitivities rather than source facts.
- Do not promote a historical fitted number to a biological constraint without primary evidence.

### R46: physiology-constrained structural reconstruction / CBM

`analysis/46_physiology_constrained_model_reconstruction/`

- Zero-storage CBM showed WT feasible and AE4 KO can structurally support full reference export under loose shared capacities.
- Pure NKCC full rescue required increased NKCC and K efflux; structural feasibility is not evidence that this compensation occurs physiologically.
- The isolated NKCC assay cannot be converted into an exact sustained WT=KO flux constraint because the protocol is transient, inhibitor-specific and lacks a numerical equivalence tolerance.
- Task 46H stopped after H04; later coupling stages H05-H07 were not completed.

### R47: dynamic potassium recycling

`analysis/47_dynamic_potassium_recycling_reconstruction/`

- No binding K-recycling ceiling emerged.
- AE4-null K-channel efflux increased about `18.07%`, pump turnover only about `1.86%`; voltage/K storage absorbed the burden.
- Acute/shared-rest null secretion deficit remained about `3.8575%`.
- **Closed:** generic K-channel/pump recycling ceiling as the missing phenotype mechanism.

### R48: joint experimental constraint reconstruction

`analysis/48_joint_experimental_constraint_reconstruction/`

Exact structural contradiction:

- In the inherited model beta input affects the core through AE4 transport only.
- With AE4 deleted, the knockout core vector field is independent of beta.
- Therefore AE4-KO CCh-only and CCh+IPR core trajectories are mathematically identical under a common observation operator for every shared parameter vector.
- Experiment reports KO CCh uptake `2.30 +/- 0.10` versus CCh+IPR `0.90 +/- 0.09 x10^-3 s^-1`; the relaxed common-mean chi contribution is `108.2873`.
- The contrast Jacobian is zero across all 33 provenance-eligible shared parameter directions.
- IPR-only AE4 KO is also predicted to have zero uptake from exact rest, while experiment reports `0.20 +/- 0.03`.

Resting-state numerical failure at unchanged shared parameters:

- WT Cl predicted `57.855` vs `50.10 +/- 1.50`.
- AE4 KO Cl predicted `56.954` vs `36.50 +/- 1.60`.
- AE4 KO pH predicted `7.2246` vs `6.89 +/- 0.02`.
- Established-genotype model predicted essentially no secretion loss (slight increase), not 35%.

Conclusion: the inherited architecture fails both chronic rest and beta-containing knockout response. No unique correction was identified.

### R49: swelling-gated VRAC and resting source geometry

Branch `analysis/task-49-camp-vrac-secretory-branch-reconstruction`, final `ab7f1cb054b8c8150343653727773f1ca59ebbcf`.

- The positive-swelling VRAC gate duplicated an earlier Task 13B diagnostic that had been overlooked.
- Task 49 integrated the separate current into full electrical/chloride closure and passed nine implementation tests.
- Independent data did not identify a usable conductance for the chosen one-scale law.
- More decisively, AE4-KO IPR from exact rest has no chemical perturbation to start positive swelling, so the gate remains off and predicts zero uptake; experiment is `0.20 +/- 0.03`.
- Resting source geometry gave a cleaner SVD/Farkas certificate for the already-known acid-base/resting incompatibility. It did **not** discover a new scientific problem or uniquely identify a correction.
- Zero production trajectories and zero parameter fits were run in Task 49.
- **Closed:** another strictly positive-swelling gate on the unchanged beta-blind core; rediscovery of NHE/AE2/CO2 resting incompatibility.

## 5. Hard no-repeat list

Do not spend a remaining shot on any of the following without genuinely new independent evidence that changes the equations:

- AE4-only stoichiometry, Na/K routing, Hill/cooperativity, scalar capacity or PKA-family searches;
- calcium-amplitude escalation;
- generic NHE1 scaling/replacement, CO2 permeability scaling, or another proof of the chronic acid-base mismatch;
- simply switching on basal NBC as an assumed cure;
- Palk NKCC replacement, arbitrary NKCC hard cap, or broad NKCC candidate-family search;
- Na/K pump capacity/coupling as the phenotype mechanism;
- K-recycling ceiling;
- AE4-expression-dependent TMEM16A recruitment interpreted as established molecular biology;
- positive-swelling-only VRAC gate on the unchanged knockout core;
- transporter subsets, all-pairs/all-triples, Cartesian grids, best-subset/L0, stochastic/global/evolutionary searches;
- exact experimental timing, onset, half-time or minute-wise fitting.

## 6. Phenotype-target convention for all remaining work

Binding rules from `docs/PHENOTYPE_TARGET_CONVENTION.md`:

- Primary target: a substantial, robust AE4-null secretion reduction of the experimentally observed order.
- Do not optimise exact time course, first 2-3 min, onset, half-time or minute landmarks.
- A reduction that exists only at one cherry-picked time is unacceptable; a robust cumulative/sustained reduction over a reasonable window is sufficient.
- Collapse CCh/beta/cAMP/PKA/Ca signalling to the smallest effective scalar/algebraic input needed. Do not add signalling ODEs merely to manufacture timing.
- Experimental timing is approximate biology, not an exact mathematical constraint.

## 7. Current reduced hypothesis for Task 50 — NOT YET AN ESTABLISHED RESULT

The remaining constructive opportunity is to **recast Task 41 as a reduced effective network coupling**, rather than claim a molecular AE4->TMEM16A mechanism.

Proposed one-parameter law:

`g_Cl,eff = g_parent * [1 - lambda * a_Ca * beta * (1 - e_AE4)]`

where:

- `a_Ca` is the existing normalized calcium/channel activation already in the parent model;
- `beta` is the existing effective protocol input (`0` without IPR/beta stimulation, `1` in the standard beta-containing stimulated condition); it is not a cAMP/PKA dynamic state;
- `e_AE4` is AE4 expression/activity fraction used only as the existing genotype coordinate;
- `lambda = 1 - b = 0.89488127156712`, inherited exactly from Task 41's phenotype-calibrated selected value `b=0.10511872843288`.

Algebraic nesting expected before any numerical work:

1. WT: `e_AE4=1` -> multiplier exactly 1 under every protocol.
2. Rest: `a_Ca=0` -> multiplier exactly 1.
3. AE4 KO + CCh only: `beta=0` -> multiplier exactly 1; exact parent CCh-only model.
4. AE2 KO: AE4 remains present (`e_AE4=1`) -> multiplier exactly 1; exact parent model.
5. AE4 KO + CCh+IPR at the standard fully beta-stimulated condition: `beta=1`; when the same calcium activation used by Task 41 applies, the law reduces to the selected Task 41 conductance factor and should reproduce its already frozen combined-stimulus phenotype.

Interpretation if implemented:

- `lambda` is **an effective phenotype-calibrated network-coupling parameter**, not a biochemical constant and not evidence that AE4 physically recruits TMEM16A.
- The scientific statement is only that a beta-conditioned AE4-dependent contribution to effective apical chloride secretory capacity is sufficient to generate the required phenotype while preserving WT/CCh-only/AE2 nesting.
- Task 50 must not search over `lambda`; it is fixed from the already published Task 41 construction. If algebra/nesting fails, stop and report the contradiction.

## 8. Mandatory decision protocol for Codex and all agents

Before **every new scientific decision**, not merely at task start:

1. Read this ledger and `docs/PHENOTYPE_TARGET_CONVENTION.md`.
2. Identify the ledger entries bearing on the decision.
3. Check `docs/RESULTS_LEDGER.md`, `analysis/README.md`, relevant final reports, `prompts/`, `analysis/`, `results/`, and remote branches for overlap/synonyms.
4. Record the decision and its ledger basis in the task's `NOVELTY_CHECK.md` or `DECISION_LOG.md` **before** running dependent compute.
5. If the decision reopens a closed axis without genuinely new independent evidence, stop. Do not compute.
6. If a checkpoint establishes a new result/failure/exclusion/qualification, update this master ledger in the same checkpoint commit/push. A checkpoint without its ledger update is incomplete.
7. Never silently overwrite a prior result. Append a superseding entry with the exact old/new relationship.

## 9. Required classification language

Every result must be labelled as one of:

- **SOURCE FACT:** directly supported by the cited experiment/source.
- **FORMAL DEDUCTION:** exact consequence of declared equations/assumptions.
- **NUMERICAL PREDICTION:** output of a frozen model/parameter set.
- **TARGET-CALIBRATED CONSTRUCTION:** parameter/model selected using the phenotype target; demonstrates sufficiency, not independent validation.
- **ASSUMPTION/SENSITIVITY:** unmeasured modelling choice.
- **NEGATIVE/EXCLUSION RESULT:** tested hypothesis failed within its declared domain.
- **UNRESOLVED:** data/equations do not identify a unique conclusion.

Task 41 and the proposed Task 50 effective coupling are explicitly **TARGET-CALIBRATED CONSTRUCTIONS** unless future independent evidence changes that status.

## 10. What is still genuinely unresolved

- We do not have an independently identified microscopic mechanism explaining why AE4 loss yields the large secretion phenotype.
- The acute model does not explain the established knockout's chronic low-Cl / near-WT-pH resting adaptation.
- The precise beta-associated pathway producing the positive IPR-only KO uptake is not identified by the current model.
- The isolated NKCC assay cannot be faithfully reconstructed by the current reduced NKCC law.
- The current model's state compensation strongly masks AE4 loss unless an additional network coupling is introduced.
- We **do** have a constructive one-dimensional direction (Task 41) that produces the right secretion magnitude. Task 50's purpose is to reduce that construction to the smallest protocol-conditioned effective law and remove the unjustified molecular interpretation, not to discover another mechanism family.

This ledger is mandatory project memory. Future agents are not permitted to replace it with recollection, a short conversation summary, or an ad hoc mechanism list.

## 11. R50A: exact coefficient provenance qualification and audit stop

Task 50 started at `8f5fefd9fc0b563a6bd963895211b4acf76972b0` on
`analysis/task-50-minimal-beta-conditioned-effective-coupling`.

**Classification:** FORMAL DEDUCTION / provenance qualification concerning a
**TARGET-CALIBRATED CONSTRUCTION**. **Disposition: STOP at 50A, before scientific
implementation.** This is not a biological rejection of the effective coupling.

The generic identity
`(1-a)+a[b+(1-b)e] = 1-(1-b)a(1-e)` is correct. The proposed Task 50 law also
has exact formal WT, REST, AE4 null CCh only, and AE2 knockout parent nesting.

However, R41 and section 7 above quote a shortened report value for `b`.
The actual Task 41 candidate 08 design, frozen input manifest and public driver
all store `b=0.10511872843288446`, preserved from commit
`86f135e075aa91a70a752b650d7bb894353fee04`. All seven execution file hashes in
that frozen manifest match at the Task 50 start. The mandated Task 50 literal
`lambda=0.89488127156712` instead has exact decimal complement
`0.10511872843288`.

Consequently the Task 50 factor minus the actual frozen Task 41 factor at
beta one is exactly `-4.46e-15*a*(1-e)` for those decimal literals. The
fully stimulated null case is a counterexample to the required exact
multiplier identity. The difference is tiny; no biological significance,
numerical secretion difference or full RHS error bound has been established.

**Relationship to prior entries:** this append qualifies the claim in section 7
that the mandated literal was inherited *exactly* from the frozen computation.
It does not invalidate Task 41's output or the algebraic construction and does
not silently rewrite its historical coefficient or files. Exact identity to
the rounded report equation must not be confused with identity to the frozen
production model. Passing a numerical tolerance would be a different claim.

The repository and all 53 fetched remote branch heads were searched for the
beta conditioned Task 41 form and synonyms; no prior implementation was found
in that scoped search. Existing Task 49 and earlier exclusions remain closed.

The Task 50 prompt's stop rule applies before 50B. No scientific source edit,
RHS evaluation, model test, trajectory, stationary solve, parameter search,
fit or new phenotype prediction was run. No frozen mutant trajectory was
relabelled as an exact Task 50 result. 50B, 50C and 50D are unexecuted; Task 51
has not started. This execution is the designated Task 50 shot; the planned
Task 51 and Task 52 shots remain unstarted.

Evidence and the precise code mapping are in
`analysis/50_minimal_beta_conditioned_effective_coupling/EQUIVALENCE_PROOF.md`,
`DECISION_LOG.md`, `NOVELTY_CHECK.md` and the task's `output/` audit records.
Continuing requires explicit resolution of the fixed literal versus exact
frozen inheritance conflict; this audit selects neither alteration.

R50A publication: audit commit `e4c72143fe7d71ae0034d15887312062626b4636` was published
on the designated Task 50 branch, independently fetched and remotely verified.
All 12 checkpoint artefact hashes match at that commit. The verification record
is `analysis/50_minimal_beta_conditioned_effective_coupling/output/publication_50A.json`.
No dependent scientific work followed the failed gate.

## 12. R50A continuation: user authorised numerical treatment of rounding

After the verified audit closeout `9bb9c9513a5816b872d75fd7ab86ee621effdb73`,
the user explicitly authorised retaining exactly `lambda=0.89488127156712`
and treating the full frozen Task 41 coefficient difference through numerical
equivalence. This supersedes R50A's operational stop only for that precision
issue. The exact difference and historical audit remain valid and preserved.

Task 50 may now proceed through 50B, 50C and 50D with the original scientific
scope, fixed coefficient, exact parent nesting and publication requirements.
The numerical software allowances and saved state fixtures are fixed in
D50-07 before computation. Reuse of mutant Task 41 summaries must explicitly
state the numerical precision qualification; no bit identity or global
trajectory error bound is implied. The model remains a
TARGET-CALIBRATED CONSTRUCTION. No fit, mechanism search, signalling ODE or
new production trajectory is authorised merely by this precision clarification.
This is continuation of Task 50; Tasks 51 and 52 remain unstarted.

### R50B setup qualification

The first software verification attempt stopped before tests or model evaluation
because the old Task 40 fixture also hashes its historical `AGENTS.md`. A full
150 file audit found exactly one difference: today's binding instructions.
All 149 other inputs match. D50-08 authorises a task local loader that pins the
current instructions while verifying unchanged scientific inputs, rather than
changing old files or restoring obsolete instructions. The failed setup is
preserved; no scientific tolerance or model parameter is changed.

### R50B: fixed reduced coupling implemented and verified

The new public `Task50EffectiveCouplingModel` implements exactly
`1-0.89488127156712*a_Ca*beta*(1-e_AE4)` through the unchanged Task 41 conductance
machinery before electrical closure. It adds no state or adjustable parameter
and changes no other scientific source. Classification:
**TARGET-CALIBRATED CONSTRUCTION**.

Ten tests pass. Twenty one parent comparisons give bit identical complete RHS
and physical diagnostics for WT, REST, CCh only null and AE2 knockout cases.
Six beta one comparisons at frozen initial/final states have maximum RHS
component difference `8.090750291955828e-15` in native units relative to the
full stored Task 41 coefficient; the maximum factor difference is
`4.5102810375396984e-15`. The predeclared numerical allowances and unchanged
conservation gates pass. Currents, charge and chloride source bookkeeping are
covered, and unrelated fixed state evaluations remain unchanged. Integration
dispatch is tested with a mock, without an ODE integration.

Attempt 02 used 80 core evaluations including constructor checks. Attempt 01
was the governance hash setup failure recorded above and performed none.
There are zero new production trajectories, stationary solves or fits. All
seven frozen Task 41 execution hashes remain unchanged. Numerical equivalence
under the authorised rounding qualification is established on these fixtures;
no exact mutant solution identity or global trajectory error bound is claimed.
See `IMPLEMENTATION_VERIFICATION.md` and `output/verification_attempt_02.json`
under the Task 50 directory. D50-09 requires remote 50B publication before reuse.

### R50C: inherited phenotype retained through qualified frozen reuse

50B was published and remotely verified at
`537ce3204d11a1a6cc51b1dd39032d12cb42df5e`. The file only 50C audit verifies
28 artefacts byte for byte against original Task 41 publication
`86f135e075aa91a70a752b650d7bb894353fee04`, including saved trajectories,
summaries, states, protocols and source. Shared parameter, initial state, onset
and solver records match. The common rest hash remains
`a0b96ef16e5967de4dfe2bdabb882f2a28d478e211f920bc19cd735a5a143acb`.

**TARGET-CALIBRATED CONSTRUCTION:** the reused cumulative deficits are
`23.163365263893244%` at 5% AE4 and `30.26115064458511%` at null AE4. Saved
endpoint flow deficits are `20.108927720908575%` and `27.038018389118477%`.
The effect is substantial across the saved 60 s cumulative observations through
600 s, not only at one selected instant. This is finite window persistence,
not a new uncertainty study, timing fit or long duration validation.

The mutant outputs remain inherited Task 41 results under the user's explicit
numerical precision authorisation. They are not newly simulated Task 50 outputs,
bit identical solutions for unequal coefficients, independent validation or a
rigorous global trajectory error bound. WT, REST, null CCh only and AE2
knockout remain exact parent initial value problems. That identity alone does
not establish agreement with the corresponding experimental assays.

Chronic knockout resting Cl/pH adaptation, positive IPR only null uptake from
genotype rest and the microscopic identity remain unresolved. The new factor
is off when the calcium coordinate is zero, so it cannot repair the IPR only
discrepancy. No source, bath, genotype rest or observation rule was changed.
New 50C model evaluations, trajectories, stationary solves and fits: zero.
One reporter filename error was corrected without scientific changes and is
preserved in `output/reuse_attempt_01.json`. See Task 50 `FROZEN_RESULT_REUSE.md`,
`output/inherited_phenotype.json` and `output/frozen_reuse_manifest.json`.

## 13. R50D: final minimal beta conditioned construction

**Controlling Task 50 result: COMPLETE WITH AUTHORISED PRECISION QUALIFICATION.**
**Classification: TARGET-CALIBRATED CONSTRUCTION.** The fixed coefficient remains
exactly `lambda=0.89488127156712`. One new wrapper changes only apical chloride
conductance through the declared factor before the inherited current closure.
No earlier scientific source, state, protocol, transporter parameter or frozen
result was changed. No signalling ODE or adjustable parameter was added.

Exact WT, REST, AE4 null CCh only and AE2 knockout parent identities hold.
At beta one the algebra is identical to the displayed Task 41 equation, while
comparison to its full frozen coefficient retains the R50A precision
qualification. The user's explicit authorisation, recorded in D50-06, supersedes
the original audit's operational stop only. It does not supersede the exact
coefficient inequality or turn fixed state comparisons into a trajectory bound.

Ten software tests pass, including 21 exact parent comparisons and six active
comparisons. All 28 reused artefacts match the original Task 41 publication.
The inherited cumulative deficits are 23.1634% at 5% AE4 and 30.2612% at null;
the saved observations establish finite window persistence. These magnitudes
were already phenotype calibrated and are not independent validation or newly
simulated Task 50 outputs. No molecular AE4 to TMEM16A interaction is identified.

The final report is
`analysis/50_minimal_beta_conditioned_effective_coupling/MINIMAL_EFFECTIVE_COUPLING_REPORT.md`.
The unresolved chronic resting adaptation, IPR only null uptake, microscopic
identity, pH/volume limitations and long duration validity remain explicit.
All earlier exclusions remain closed.

The original 50A audit was verified at
`e4c72143fe7d71ae0034d15887312062626b4636`; the authorised 50A continuation
was verified at `cc69edfb88fbcf17b6d60e0cbc70c0130c0d978a` before source edits.
50B was verified at `537ce3204d11a1a6cc51b1dd39032d12cb42df5e` before reuse.
50C was verified at `88997ef6888498e6ea6bdd9bd4acd7e5ec108dbc`, including
all ten receipt hashes, before final reporting. The final 50D publication
receipt records its exact remotely verified commit after publication.

Across Task 50 there were 80 core evaluations, all in software verification
and constructor checks; zero production trajectories, stationary solves,
parameter fits or mechanism searches. The historical setup and reporter errors
remain recorded. The two original verification logs are included explicitly in
50D after the repository's generic ignore rule excluded them from 50B; the JSON
test records were already published. Final integrity verification uses file
and Git checks only. Stop after remotely verified 50D. Tasks 51 and 52 and the
two remaining scientific shots remain unstarted.

R50D final integrity audit: all 39 artefact hashes across the original 50A,
authorised 50A continuation, 50B and 50C receipts match at their own published
commits. The three tested source hashes and all 28 frozen artefacts also match.
Earlier scientific source, analyses, results, archive and manuscript files are
unchanged. Binding instructions, phenotype convention and execution prompt are
unchanged; the earlier ledger contents are preserved as an exact prefix. The
audit passed without any model evaluation. Its record is
`analysis/50_minimal_beta_conditioned_effective_coupling/output/final_integrity_audit.json`.
