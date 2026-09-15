# Task 46H04: paired WT/KO evidence protocol

Parent checkpoint: `cb18477adc393d5a18ea6f2ab8f41f4dc42376fc` (46H03).

## Question

Can the preserved isolated NKCC1 assay be translated into a quantitative cross-genotype constraint for the sustained physiological CBM before the six declared coupling diagnostics are run?

## Evidence result

No.

The preserved assay reports no detectable genotype difference in bumetanide-sensitive initial chloride uptake under bicarbonate-free solution with carbonic anhydrase and CaCC inhibited. The preserved audit supplies no numerical equivalence tolerance. It therefore does not justify an exact equality, percentage band, or upper bound relating WT and AE4-null NKCC1 flux during normal physiological stimulation.

There is also a model-class mismatch. The assay observable is an initial uptake rate and therefore permits transient intracellular chloride storage. The frozen CBM imposes zero sustained amount storage. In addition, the assay removes bicarbonate-dependent conditions and inhibits carbonic anhydrase and CaCC, whereas the CBM represents a conditional sustained secretion sector rather than that inhibited transient protocol. A protocol-matched reconstruction would require new transient state equations and inhibitor-specific assumptions, which Task 46H forbids.

Accordingly the assay remains an independent qualitative, protocol-specific diagnostic. It is not converted into a numerical physiological constraint.

## Smallest paired LP

`cbm/paired_lp.py` builds two copies of the frozen compiled LP above the unchanged local vendor core, fixes WT `J=1`, optionally accepts an explicit cross-genotype NKCC equality row, and maximises KO `J`. The evidence-approved H04 problem does **not** activate that optional equality row because no transferable numerical equality/tolerance is supported.

The paired maximum is exactly `KO J=1`. No new optimisation is needed to establish this value: the reporting region gives the upper bound `J<=1`, while H02 supplies a validated WT `J=1` witness and H03 supplies a validated KO `J=1` witness. With no evidence-supported cross-genotype numerical row added at H04, those two states are jointly admissible.

At fixed KO `J=1`, NKCC flux can range from `N=0` to `N=0.5`. Thus strong NKCC recruitment is optional, not structurally required, even at maximal reported KO output.

## Assay-matched structural diagnostic

The frozen cellular chloride identity is

`2*N + A + E - J = 0`.

For the isolated uptake protocol, the bicarbonate-dependent AE4/AE2 routes are unavailable in the intended structural check and CaCC is inhibited, so setting `A=0`, `E=0`, `J=0` reduces the sustained identity to

`2*N = 0`.

Hence the zero-storage CBM cannot represent positive initial NKCC uptake under that assay protocol. Positive initial uptake requires transient intracellular chloride storage, which is precisely what the sustained CBM omits. This is why the assay cannot be converted into a normal-stimulation NKCC flux equality without changing model class.

## Frozen paired construction

H04 freezes the cross-genotype rules before H05/H06 diagnostics:

1. WT and KO use the same audited stoichiometric matrix and the same non-genotype-specific computational bounds.
2. KO differs only by `A = 0` for AE4.
3. The normalisation `0 <= J <= 1` is a reporting restriction, not a measured capacity.
4. Computational `M=100` remains the primary numerical bound with the single already-declared `M=1000` sensitivity interpretation.
5. No realised WT and KO flux is forced equal merely because a capacity assumption is shared.
6. No NKCC1 WT/KO equality, tolerance band, ratio, or stimulation cap is imposed from the isolated uptake assay.
7. No reserved secretion magnitude, trajectory, resting chloride phenotype, pH, CaCC behaviour, Task 42 outcome, or old fitted capacity is used to choose a paired restriction.
8. Later single-resource tests may share the same algebraically derived cap between genotypes only as predeclared limiting diagnostics. Such a cap is not called a measured physiological capacity.
9. H05 is restricted to the six coupling questions C1-C6 frozen in `output/46h/constraint_budget.json`.
10. H06 is restricted to the six single resource tests G1-G6 frozen in that same budget. No pairs, subsets, grids, or target-selected tuning are permitted.

## Interpretation of old-model style strong NKCC rescue

An old-model style roughly 90% replacement of lost AE4 chloride loading is **merely allowed by loose constraints**. It is not structurally required because KO `J=1` remains feasible with `N=0`, and it cannot be declared incompatible with the preserved isolated NKCC evidence because that assay supplies no numerical normal-stimulation tolerance that can be transferred into this sustained CBM.

This is not a claim that physiological NKCC compensation is unlimited. It is the narrower statement that the preserved evidence does not provide the numerical sustained-flux restriction required to exclude the H03 rescue witness inside this CBM.

## Decision

**Accepted evidence translation:** qualitative protocol-specific NKCC diagnostic only.

**Rejected translations:** exact WT=KO NKCC flux; arbitrary epsilon-equivalence band; ratio bound; transfer of the isolated uptake result to normal stimulation; import of old model NKCC capacities.

**Next milestone:** H05, six predeclared targeted coupling questions, one question per independent worker where tooling permits. No new family may be introduced.
