# Task 42 source, charge and thermodynamic audit

This execution uses exactly C0, C1, C2, C3a, C3b, C4a and C4b. All seven focused tests pass before any stationary solve or production integration. The frozen reference is the accepted Task 40 WT rest. Every inherited scientific input is checked against the remote preparation tree.

The source layer replaces only the AE4 contribution object. The complete model already consumes separate TIC and TA properties, so carbonate contributes one carbon and two alkalinity equivalents without duplicating or editing the whole cell equations. C0 gives a bitwise identical complete RHS to Task 40 in the tested resting and stimulated states and all three expressions. Tests also check every other class against the literal declared source difference, invariant scalar J4, unchanged non AE4 diagnostics and independent cellular balances.

The table gives cellular source coefficients per positive inward chloride cycle and forward affinity divided by RT. Concentrations approximate activities.

| Class | Na | K | Cl | TIC | TA | A/RT | Forward direction | J4 agrees |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| C0 | -0.5 | -0.5 | 1 | -2 | -2 | -2.106668764 | opposed | False |
| C1 | 1 | -1 | 1 | -1 | -1 | 4.845371436 | supported | True |
| C2 | 1 | -2 | 1 | -1 | -2 | 6.871880179 | supported | True |
| C3a | 0 | -1 | 1 | -2 | -2 | 0.7365707451 | supported | True |
| C3b | -1 | 0 | 1 | -2 | -2 | -4.949908274 | opposed | False |
| C4a | 0 | -1 | 1 | -1 | -2 | 1.18540116 | supported | True |
| C4b | -1 | 0 | 1 | -1 | -2 | -4.501077859 | opposed | False |

For each row, charge is exactly `nu_Na + nu_K - nu_Cl - nu_TA = 0`, using rational arithmetic. For every transported quantity, `S(-J4) = -S(J4)`. Source coefficients are fixed; there is no adjustable cation fraction.

## Whole cell amount balances

Let N be inward NKCC1 cycles, H inward NHE1 cycles, B inward NBC cycles, E inward AE2 chloride cycles and P the sum of both Na/K pump cycle rates. Let K_a and K_b be outward apical and basolateral potassium fluxes, C_a outward CaCC chloride flux and D the sum of inward CO2 fluxes. All rates below use fmol/s; TA counts equivalents. For each class, insert its five coefficients from the table:

```
dNa_i/dt  = N + H + B - 3P + nu_Na J4
dK_i/dt   = N + 2P - K_a - K_b + nu_K J4
dCl_i/dt  = 2N + E - C_a + nu_Cl J4
dTIC_i/dt = D + 2B - E + nu_TIC J4
dTA_i/dt  = H + 2B - E + nu_TA J4
```

Volume evolves under the unchanged bath to cell and cell to lumen water fluxes. Concentration changes follow `dc/dt = (dn/dt - c dV/dt)/V`. Paracellular transport enters the inherited lumen balances and current closure, not these cellular transport sources. Multiplying the amount balances by charge coefficients gives `d(Na_i+K_i-Cl_i-TA_i)/dt = -B-P-K_a-K_b+C_a`; AE4, NKCC1, NHE1 and AE2 cancel exactly. The unchanged membrane current solve closes the remaining electrogenic terms.

## Thermodynamic calculation

For signed cellular species coefficients nu_s, the forward free energy divided by RT is `sum_s nu_s log(c_i,s/c_o,s) + (V_b/V_T) sum_s nu_s z_s`. Forward affinity is its negative. Each individual chemical and electrical term is saved in preflight.json. The electrical contribution cancels because every class is electroneutral, including carbonate valence minus two. The nonzero membrane potential is included before this cancellation.

Carbonate is taken from the inherited acid base speciation in the cell and bath, with `CO3/HCO3 = 10^(pH-pKa2)`. The inherited pKa2, temperature, reference concentrations and potential are recorded in preflight.json. Carbonate is never replaced by bicarbonate in the audit or cellular sources.

Positive A supports the proposed inward chloride direction, negative A opposes it, and values within 1e-12 of zero are labelled indeterminate. This direction diagnostic does not remove any algebraically valid class from the prescribed panel. Keeping J4(state) fixed is an intentional source comparison assumption. A mismatch between its direction and a class affinity is recorded, not repaired by changing rates. These models do not claim a thermodynamically derived kinetic law for every new class.

## Execution and publication boundaries

The seven source classes, parameter values, solver tolerances, physiology gates and stimulus are fixed before execution. Each class receives one WT root solve from the accepted Task 40 state. Every class passing rest receives WT, 5% and null trajectories from that identical class rest, with the inherited 600 s protocol. No root search, parameter fitting, sweep, added class or genotype rest is permitted. Numerical retry allowances remain those in the supplied task; no retry follows a physiology failure.

The source audit is published first, then all seven rests, then each class immediately after its three trajectories. The complete prediction checkpoint is published last and execution stops before phenotype reveal. Source and execution hashes are frozen in preflight.json and each class manifest. Publication history records confirmed remote commits.

Previous Task 42 source scripts were reused after review; previous numerical outputs were not copied into this execution. Changes before the audit comprise publication stage controls, a verifier that does not assume previous outcomes, compact reporting and accurate class completion status when a trajectory fails. The scientific source layer and all inherited mechanisms remain unchanged.

The required instructions and Task 40 report contain a phenotype reference. No experimental phenotype source or Task 41 content is inspected for this execution. No target value is used by source code, gates, execution or reports. This checkpoint makes no claim about agreement with the held out phenotype.
