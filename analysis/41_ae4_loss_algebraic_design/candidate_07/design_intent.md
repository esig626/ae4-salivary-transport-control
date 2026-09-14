Candidate 6's pH feedback stationary solution gives a 30% loss but null
pH=7.3191, outside the inherited long-run box. Its gate would also activate
late relative to Task 40's measured finite-time pH. It is not simulated.

Construct a direct AE4-dependent stimulated CaCC recruitment hypothesis:

  g_CaCC = g_parent * [(1-a) + a * (b + (1-b)*e_AE4)]

where a is the inherited normalized calcium activation, e_AE4 is AE4
expression, and b is the residual recruitment fraction in a stimulated
null. The same law/parameter b applies to all genotypes. WT is exactly
Task 40 for every state and stimulus, and at rest the conductance is the
parent value for every genotype. AE4 expression explicitly affects channel
recruitment in this NEW hypothesis; it no longer affects transport only.
This dependency is not deduced or validated by AE4 stoichiometry. It is
a transparent constructive model selected using the requested phenotype.

Derive b from a sustained null/WT flow ratio of 0.70 by solving the full
stationary ion/water/electrical balances. Use the inherited WT stationary
state and one null algebraic state, then freeze b before production. A
stationary pH limitation must be retained in the report; finite 600 s
trajectories retain all original pH/Cl/conservation gates. They determine
the actual primary 20-35% cumulative secretion target.

Only this recruitment dependency is added to Task 40. No finite NKCC1 cap,
pump affinity change, NBC/AE4 capacity increase, or pH conductance law from
earlier candidates is carried into this model. The channel factor enters
the original electrical closure, including NBC; secretion is computed by
the unchanged water and lumen-outflow equations.
