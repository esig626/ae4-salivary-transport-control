# Task 51 decision log

This log is binding. Every new scientific decision must cite the relevant `docs/MANDATORY_RESEARCH_LEDGER.md` entry before dependent compute.

## D51-01 — canonical start

**Decision:** Task 51 starts only from merged `main` commit `4bb2c89870fe3f4d4bb887cf9d40ba7de377a725` on branch `analysis/task-51-beta-nkcc-swelling-vrac-mechanistic-reconstruction`.

**Ledger basis:** Tasks 49 and 50 are complete and merged; two scientific shots remain.

## D51-02 — preserve Task 50 independently

**Decision:** the Task 50 effective-coupling model is a frozen benchmark, not the Task 51 parent mechanism. Do not edit or repoint `archive/task-50-working-effective-coupling-benchmark` and do not import its multiplier into Task 51 equations.

**Benchmark:** `docs/TASK50_EFFECTIVE_COUPLING_BENCHMARK.md`.

## D51-03 — one permitted mechanistic chain

**Decision:** test exactly

`beta/IPR -> beta/cAMP NKCC1 activation -> solute loading -> swelling -> Task49 VRAC-like current -> secretion`

with existing

`beta/PKA -> AE4 activation`.

No alternate mechanism family may be introduced during Task 51.

**Ledger basis:** R48 exact beta-blind KO defect; R49 swelling-only self-start failure; R50 target-calibrated proof that a substantial beta-conditioned network effect is sufficient; section 1.7/R50E source evidence.

## D51-04 — NKCC core frozen

**Decision:** retain the source-fixed Palk/Benjamin NKCC concentration-response core, exact 1Na:1K:2Cl stoichiometry and all thermodynamic/reversal properties. Add only the minimum beta contribution to the existing activity multiplier. At beta=0 the Task 51 NKCC law must exactly nest the parent regulator.

No NKCC cap, transporter-law replacement or mechanism search is authorised.

## D51-05 — minimum beta-regulation law

**Decision:** use the smallest no-interaction algebraic combination of the existing Ca/CCh activity arm and one beta/cAMP arm. Do not introduce Ca×beta synergy, a signalling ODE, phosphorylation state, delay or Hill coefficient unless a primary source makes that extra degree unavoidable.

The preferred structure is additive above baseline so beta=0 exactly preserves the inherited Ca arm and IPR-only beta stimulation can increase NKCC from rest.

The exact algebraic form must be written and source-audited before implementation.

## D51-06 — beta-NKCC quantitative treatment

**Decision:** direct rat salivary NKCC data establish sign and plausible scale but are not exact mouse-SMG constants.

Primary source anchors:

- rat submandibular IPR response approximately 2.5-fold (PMID `9880083`);
- rat parotid approximately 3-fold functional response and K1/2 21.5 nM (PMID `1313447`);
- approximately 20 nM phosphorylation half-maximal response and stronger recruitment context from PMID `7559664`;
- PKA/cAMP regulatory support from PMIDs `10600770`, `11880270`.

Use mouse-SMG IPR swelling/current observations to identify/check the effective mouse strength wherever possible. Do not select the beta gain from AE4 secretion.

## D51-07 — VRAC law and conductance identification hierarchy

**Decision:** reuse the Task 49 conservation-tested positive-swelling VRAC-like current law. Do not search alternate gate families.

Parameter-identification hierarchy, fixed before AE4 reveal:

1. First retrieve/use the Catalán supplementary information and derive a conductance map from the IPR-induced/DCPIB-sensitive current only if the intracellular/extracellular chloride conditions and current decomposition support that derivation.
2. If a direct current-to-conductance map remains non-identifiable, the predeclared fallback is one scalar VRAC conductance calibration to the independent mouse-SMG IPR swelling magnitude (approximately 12.5%) while the beta-NKCC gain is independently source-fixed/bounded. The current and blocker data then serve as validation/consistency checks.
3. Do not use whole-gland AE4 phenotype or Task 50 output to identify VRAC.

This hierarchy is an observation-map fallback, not a mechanism search.

## D51-08 — independent calibration data only

**Decision:** before the pre-AE4 freeze, calibration/constraint data may include direct beta-NKCC salivary studies, Catalán mouse-SMG IPR swelling, current/voltage and blocker information, and non-AE4 beta-pathway controls. The JBC 2015 AE4 phenotype and Task 50 benchmark are held out.

## D51-09 — immutable reveal boundary

**Decision:** after the mechanistic parameters and all source/protocol mappings are fixed, publish a remote immutable prediction checkpoint containing WT, AE4 5%, AE4 null and AE2/control predictions. Only then reveal/compare:

- AE4-null 10-minute secretion deficit;
- AE4-null CCh vs CCh+IPR uptake contrast;
- AE4-null IPR-only uptake;
- Task 50 benchmark outputs.

No post-reveal retuning.

## D51-10 — success interpretation

**Decision:** success does not require exact 35.000% or exact Task 50 numbers. A successful source-backed mechanism must prospectively generate a substantial persistent AE4-loss deficit of the same order, improve/remove the beta-blind IPR-only structural failure, preserve CCh-only nesting, and not destroy AE2 specificity or physical/conservation constraints.

Task 50 is the quantitative reference result, not a fitting target.

## D51-11 — crash-safe publication

**Decision:** publish and remotely verify 51A through 51F before dependent work. Every checkpoint that establishes a new fact/failure/qualification must update `docs/MANDATORY_RESEARCH_LEDGER.md` in the same dependency boundary.

Do not start Task 52 automatically.

## D51-12 — accepted source/architecture and operational start (51A)

Ledger reread: §§0,1.7,2,3.1,R39–R50E,5,6,8,9,14. The user-specified operational head ad0e7344 supersedes the older D51-01 wording only for checkout; scientific parent remains4bb2c898. Two-shot budget applies. Read-only primary-source advisers verified the beta-NKCC anchors; Paulais is the correct first author of JCI1992. SI access failed at official routes, recorded in SOURCE_ARCHITECTURE_FREEZE.md, so no direct conductance identification is claimed. The prescribed swelling fallback is accepted; this does not reopen a closed axis. All old source/results/archive files remain frozen.

## D51-13 — exact implementation and early-time implication (51A)

Ledger §§1.7,2,3.1,R48,R49,14.2 were reread. Accept the exact additive regulator adapter and unchanged Task49 current. Both NKCC-recomputing layers must receive the changed multiplier. The conditional V-second-derivative proof in SOURCE_ARCHITECTURE_FREEZE.md removes the old exact deadlock; sustained behaviour awaits computation. Task50/41 conductance wrappers are absent. No new signalling state, synergy, core, stoichiometry, cap or other pathway change is authorised.

## D51-14 — fixed fallback/protocol/observation plan before compute (51A)

Ledger §§1.1,1.7,2,R40,R48,R49,6,9,14.3 and phenotype convention reread. Central M_beta2.5, sensitivity3 with identical g. One scalar g solve against WT Track2 mean300–600s relative swelling0.125, bracket0 to1e-6S and Brent tolerance1e-13S/max40 iterations, as detailed in the freeze. This is an explicit amplitude-map assumption, not a source clock fit, and narrative/figure discrepancy remains. Two tracks, cached rests, solver settings, diagnostics, numerical budget and conditional SPQ limitations are fixed in that document. No AE4/Task50 value is a calibration residual. Non-identification is reportable through51F, without substituting the Task49 theorem.
