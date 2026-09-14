# Task 28 partial handoff for Task 29

This note records the Task 28 partial result supplied by the user after the
Codex run stopped at its predeclared per-solve residual limit. It is included
only so Task 29 can reason from the same candidate law and observed failure
mode. Re-derive all equations from the seed code before using them.

## Candidate NHE1 correction used in Task 28

Task 28 retained the existing reversible NHE1 chemical-affinity factor and
multiplied it by one bounded intracellular-pH activation factor, normalized at
the saved WT resting pH:

    J_H = C_H * tanh(log((Na_e * H_i)/(Na_i * H_e)) / w)
          * g(pH_i) / g(pH_WT_ref)

    g(pH) = 1 / (1 + 10**(1.40 * (pH - 6.57)))

The pH-response shape was borrowed from Hisamitsu et al. (2007), Table 1,
human NHE1 expressed in PS120 cells. It is not an identified salivary-acinar
kinetic law. Task 28 treated the extra pH dependence as a modelling hypothesis,
not as transported-proton stoichiometry. The same law was to apply to WT and
AE4 null.

Constitutive checks reported by Task 28 passed: positivity/bounds, fixed WT
normalization, 1 Na inward per 1 H outward, zero transported charge, zero
carbon source, genotype independence, chemical reversal, and declining proton
extrusion as the cell alkalinises.

## Completed Task 28 WT observations

Only R09 and R10 were used.

Imported/reference WT REST:

- R09 pH 6.88759, Cl 60.39895 mM.
- R10 pH 7.14663, Cl 66.60712 mM.
- Positive AE4 loading shares were 1.7374% and 1.3766% respectively.

With the corrected NHE1 law at stationary WT stimulation, Ca = 0.25 uM:

- R09 pH 6.86685, Cl 54.79236 mM, AE4 positive loading share 1.2457%,
  ductal outflow 0.00156594 pL/s.
- R10 pH 7.13083, Cl 63.33491 mM, AE4 positive loading share 0.9207%,
  ductal outflow 0.00193493 pL/s.

Signed stationary fluxes, fmol/s:

- R09: AE4 Cl +0.00460616; NKCC1 Cl +0.36514990; AE2 Cl -0.00194788;
  NHE1 H extrusion +0.00726445; apical CaCC Cl export +0.36780819.
- R10: AE4 Cl +0.00406234; NKCC1 Cl +0.43715788; AE2 Cl -0.00091214;
  NHE1 H extrusion +0.00721254; apical CaCC Cl export +0.44030808.

Negative AE2 means chloride export/bicarbonate import under the Task 28 sign
convention.

## Failed 20% allocation construction

Task 28 attempted one prescribed WT allocation in which AE4 supplied 20% of
the positive basal chloride-loading pool at stationary Ca = 0.25 uM. The
construction failed full stationary balance closure in both backgrounds.

The algebraically required AE2 capacities became approximately:

- R09: 6,816,633 fmol/s.
- R10: 21,596,316 fmol/s.

Reference AE2 capacity was 0.005 fmol/s. The candidate states also failed the
full independent residual gate by many orders of magnitude. Task 28 therefore
did NOT establish a 20% productive AE4 allocation. Do not resume or enlarge
that allocation experiment in Task 29.

## Exact-null partial iterate

The first R09 corrected-NHE1 exact-null REST solve stopped after the predeclared
1,800 residual evaluations for that solver call. Its best retained iterate had:

- pH 7.01334;
- cell volume 1.54031 pL;
- squared scaled residual 0.00410439.

This was not an equilibrium and must not be used as a solved null state. R10
null was not attempted. No stimulated trajectories were run.

The notable qualitative observation is that the corrected NHE1 candidate moved
the partial R09 null iterate far away from the previous catastrophic pH near 8
and multi-pL swelling, but Task 28 did not establish closure.

## Balance identities reported in Task 28

Task 28 reported, for positive directions defined in the seed model:

    d nCl_i/dt = 2 J_N + J_2 + J_4 - J_a
    d A_i/dt   = J_H - J_2 - 2 J_4
    d T_i/dt   = J_CO2,b + J_CO2,a - J_2 - 2 J_4

Thus at exact AE4 null and stationary intracellular chloride/base/carbon:

    J_2 = J_H
    J_CO2,b + J_CO2,a = J_H
    J_a = 2 J_N + J_H

Task 29 must verify these identities directly from the model code, then ask
whether they can be satisfied near the experimental AE4-null pH and chloride
without a numerical search.

## Experimental context to use, not fit

Peña-Münzenmayer et al. (2015) reports approximately:

- WT resting Cl_i 50.10 +/- 1.50 mM.
- AE4-null resting Cl_i 36.50 +/- 1.60 mM.
- WT resting pH 6.91 +/- 0.07.
- AE4-null resting pH 6.89 +/- 0.02.

The near-unchanged resting pH is the relevant target. Task 29 must not impose
acidification as a requirement and must not use the approximately 30-35%
secretion reduction as an optimisation objective.
