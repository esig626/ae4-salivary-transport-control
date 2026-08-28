# Task 13B final answer

## Classification

`FULL MODEL NOT IDENTIFIABLE FROM EXISTING DATA — ONE DECISION-CRITICAL MEASUREMENT SPECIFIED`

The modern model is mathematically closed and its native-source WT resting
roots and dynamics are numerically robust, but no tested member passes the
absolute one-SMG WT secretion-scale gate.  The frozen AE4-null secretion
holdout was therefore not revealed.  This is a model-reconstruction result,
not evidence that AE4 is biologically insufficient and not a claim that the
2018 article was wrong.

## Direct answers

1. **State vector.** Twelve conserved core states are intracellular amounts of
   Na, K, Cl, total inorganic carbon and total alkalinity; cell volume; luminal
   Na, K, Cl, total inorganic carbon and total alkalinity; and lumen volume.
   R1 adds one effective regulatory state, R2 two effective cAMP/downstream
   states, and R3 PKA plus an AE4 regulatory fraction. Apical and basolateral
   voltages are algebraic current-closure variables; AE4 carrier occupancy is
   QSS.
2. **Differences from 2018.** The reconstruction uses amount coordinates,
   explicit finite compartments, conserved carbon/alkalinity chemistry,
   finite buffer speciation, two membrane-current closures, reversible
   NKCC1/NHE1/AE2, mixed Na/K QSS AE4, apical and basolateral pump/K topology,
   paracellular transport, corrected dimensional water permeabilities,
   physical seconds, and dynamic beta/cAMP/PKA regulation.
3. **Task 13 inheritance.** Shared-pool mixed-cation AE4 structure,
   electroneutrality, reversal, and mutation asymmetry were retained. Explicit
   fast carrier states were simplified to QSS. The inherited seven-state
   chassis, Na-only AE4 closure, and common membrane-conductance multiplier
   were rejected for this task.
4. **Smallest defensible regulation.** R1, one effective beta/PKA/AE4
   activation state, is the smallest dynamic subsystem. R2 and symmetric R3
   are observationally equivalent under the available data; R4 is not
   justified.
5. **What 2021 identifies.** It identifies beta/adenylate-cyclase/PKA-dependent
   AE4 activation and S173 dependence, with S273A not showing the same loss.
   It does not identify cAMP, PKA, phosphorylation, dephosphorylation, or acute
   whole-cell kinetic constants, nor prove direct phosphorylation of S173.
6. **Slow scale.** No unique source-supported slow time scale is identified.
   The frozen 10/30/90 s panel is an assumption sensitivity learned without KO
   secretion.
7. **Regulatory coupling.** Common capacity modulation is sufficient for all
   available transporter/regulatory constraints; state-specific modulation is
   not distinguishable and was rejected as unnecessary complexity.
8. **Acid-base closure.** Total inorganic carbon plus total alkalinity with
   finite intracellular/luminal buffers and explicit bath/outflow boundary
   fluxes is required. It conserves internal carbon and avoids an
   infinite-buffer shortcut.
9. **Cation topology.** Apical and basolateral pumps and Ca-activated K
   channels are allowed and required for current/source closure. Total K
   capacity is conserved; apical fractions remain assumption panels rather
   than measured localization fractions.
10. **WT rest.** Ten native `N_ABS_NKCC` scale-4 roots pass the declared Cl,
    pH, Na, K, positivity, current, charge, carbon, rank, and boundary gates.
    All five topology panels survive in two routing/pump contexts.
11. **WT stimulation.** All 800 production dynamic cases integrate and pass
    numerical, solver, nearby-state, sustainment, co-stimulation and regulatory
    checks on a 600 s physical axis. None passes absolute one-SMG scaling:
    required scales are `6362.31--8210.72`, versus ceiling `3088.15386`.
12. **AE2 deletion.** The earlier superseded ensemble predicted only a small
    effect, but no replacement native model passed G5; a final native AE2
    genotype prediction is therefore not licensed.
13. **AE4 deletion.** Not evaluated. Exact deletion code exists and is tested,
    but staged release forbids its scientific use without a valid WT G5 model.
14. **Ten-minute deficit.** Not evaluated; the strict target remains sealed.
15. **Early/sustained pattern.** Not evaluated; the strict timing target
    remains sealed.
16. **Claim types.** Conservation, charge identities, nesting and reversal are
    exact deductions; source protocols and phenotype statements are measured;
    roots and trajectories are numerical; routing, kinetics and assay-to-cell
    capacity transfer are assumptions.
17. **Weakly identified parameters.** Absolute NKCC1 capacity transfer,
    effective calcium, pump/K partition, AE4 Na/K routing, regulatory kinetics,
    effective secretory-unit count, and the whole-cell amplitude of the
    beta/volume-sensitive apical anion exit remain weak or unmeasured.
18. **Minimal successful mechanism.** None is yet validated. Resting closure
    specifically requires the tested NKCC1 absolute-capacity direction, but
    that direction does not close absolute stimulated secretion.
19. **Decision-critical measurement.** Measure 600 s absolute lumen-volume
    secretion per imaged WT acinar cell under `0.3 uM` CCh + `5 uM` IPR at
    10 s resolution. The frozen discriminator is
    `0.00291436--0.00323818 pL/s` per `1.30 pL` modeled cell. This chooses
    whether the remaining failure is the cellular flux model or the
    cell-to-gland observation map.
20. **Ready for reduction/GSPT?** **NO.** The absolute WT dynamic gate has not
    passed, so reduction would preserve an unresolved scale defect.

