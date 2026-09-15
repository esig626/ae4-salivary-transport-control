# Task 46H04: paired WT/KO evidence protocol

Parent checkpoint: `cb18477adc393d5a18ea6f2ab8f41f4dc42376fc` (46H03).

## Question

Can the preserved isolated NKCC1 assay be translated into a quantitative cross-genotype constraint for the sustained physiological CBM before the six declared coupling diagnostics are run?

## Evidence result

No.

The preserved assay reports no detectable genotype difference in bumetanide-sensitive initial chloride uptake under bicarbonate-free solution with carbonic anhydrase and CaCC inhibited. The preserved audit supplies no numerical equivalence tolerance. It therefore does not justify an exact equality, percentage band, or upper bound relating WT and AE4-null NKCC1 flux during normal physiological stimulation.

There is also a model-class mismatch. The assay observable is an initial uptake rate and therefore permits transient intracellular chloride storage. The frozen CBM imposes zero sustained amount storage. In addition, the assay removes bicarbonate-dependent conditions and inhibits carbonic anhydrase and CaCC, whereas the CBM represents a conditional sustained secretion sector rather than that inhibited transient protocol. A protocol-matched reconstruction would require new transient state equations and inhibitor-specific assumptions, which Task 46H forbids.

Accordingly the assay remains an independent qualitative, protocol-specific diagnostic. It is not converted into a numerical physiological constraint.

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

## Mathematical consequence at H04

With no defensible quantitative cross-genotype measurement to add, the admissible paired physiological region contains independent WT and KO feasible copies under the common construction assumptions. H03 already established that the KO copy contains a valid `J=1` witness with AE4 deleted. H04 therefore does not, by evidence alone, remove that witness.

This is not a claim that physiological NKCC compensation is unlimited. It is the narrower statement that the preserved assay does not provide the numerical sustained-flux restriction required to exclude the H03 witness inside this CBM.

## Decision

**Accepted evidence translation:** qualitative protocol-specific NKCC diagnostic only.

**Rejected translations:** exact WT=KO NKCC flux; arbitrary epsilon-equivalence band; ratio bound; transfer of the isolated uptake result to normal stimulation; import of old model NKCC capacities.

**Next milestone:** H05, six predeclared targeted coupling questions, one question per independent worker where tooling permits. No new family may be introduced.
