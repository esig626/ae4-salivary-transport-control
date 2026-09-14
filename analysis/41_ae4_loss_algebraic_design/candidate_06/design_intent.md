Candidate 5 did not converge. No new trajectory has run after the failed
capacity-only set. The static shared-capacity designs have not yielded an
accepted full model meeting the target with chloride/pH persistence.

Test a single alternative feedback hypothesis on the unchanged Task 40
model: intracellular alkalinization reduces CaCC conductance. Task 40
already predicts a large WT/null pH separation. A shared state-dependent
outlet gate can make that separation reduce chloride export instead of
allowing unrestricted export to deplete chloride while NKCC1 compensates.
This is an explicit phenomenological hypothesis, not a measured TMEM16A
pH law. None of the preceding NKCC1, pump or NBC/AE4 changes is carried over.

Use g_CaCC(pH)=g_CaCC_parent*expit((pH_half-pH_i)/0.02).
The width 0.02 pH units is an engineering shape assumption, fixed before
inverse selection. Solve stationary WT/null ion balances and a 30% null
sustained flow deficit to determine pH_half. The gate enters the existing
electrical closure, with no secretion-output or genotype multiplier.
Both forward and reverse channel currents use the same positive conductance.

The sole genotype parameter remains AE4 expression; AE4 stays exactly 50:50
with the Task 40 total cycle law. If an admissible algebraic root is found,
implement the feedback, solve one WT rest, and test three 600 s trajectories
against the unchanged physiological/conservation gates. The primary target
remains a 20-35% reduction in cumulative 0-600 s secretion, to be measured
in trajectories, not inferred from stationary flow.
