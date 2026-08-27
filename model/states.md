# Baseline model states

## Conventions

The published model has three homogeneous compartments: fixed interstitium (`e`), dynamic cytoplasm (`i`), and fixed-volume lumen (`l`). Concentration units are mM unless a cited table states otherwise; volumes are in pL. Intracellular balance equations are written for amounts, `[X]_i ω_i`, rather than concentrations alone.

This file records the published state inventory. It does not choose between inconsistent published sign conventions or add missing dynamics.

## Differential states in the QSS formulation

| State | Compartment | Stored mathematical quantity | Published location | Notes |
|---|---|---|---|---|
| Intracellular Cl− | Cytoplasm | `[Cl−]_i ω_i` | §2.3 Eq. (1) | Amount balance |
| Intracellular Na+ | Cytoplasm | `[Na+]_i ω_i` | §2.3 Eq. (2) | Amount balance |
| Intracellular K+ | Cytoplasm | `[K+]_i ω_i` | §2.3 Eq. (3) | Amount balance |
| Intracellular HCO3− | Cytoplasm | `[HCO3−]_i ω_i` | §2.3 Eq. (4) | Amount balance |
| Intracellular H+ | Cytoplasm | `[H+]_i ω_i` | §2.3 Eq. (5) | Amount balance; pH conversion is not specified explicitly |
| Intracellular CO2 | Cytoplasm | `[CO2]_i ω_i` | §2.3 Eq. (6) | Amount balance |
| Luminal Na+ | Lumen | `[Na+]_l` with fixed `ω_l` | §2.6 Eq. (14) | Concentration balance with convective removal |
| Luminal K+ | Lumen | `[K+]_l` with fixed `ω_l` | §2.6 Eq. (15) | Concentration balance with convective removal |
| Luminal Cl− | Lumen | `[Cl−]_l` with fixed `ω_l` | §2.6 Eq. (16) | Concentration balance with convective removal |
| Cell volume | Cytoplasm | `ω_i` | §2.5 Eq. (9), §2.7 | **Blocked:** opposite signs are printed in these two locations |

The paper identifies this reduction as ten differential equations.

## Algebraic variables and constraints

| Variable | Definition/source | Role |
|---|---|---|
| `V_b` | §2.4 Eq. (7), QSS constraint in §2.9 | Basolateral membrane potential |
| `V_a` | §2.4 Eq. (8), QSS constraint in §2.9 | Apical membrane potential |
| `V_t` | Appendix 10, `V_t = V_a - V_b` | Transepithelial/tight-junction potential |
| Channel and tight-junction reversal potentials | Appendices 3, 4, and 10 | Nernst relations evaluated from compartment concentrations |
| `x_i` | Appendix 9 Eqs. (32)–(33) | Fixed moles of impermeant intracellular negative charge, derived from electroneutrality |
| Luminal/interstitial electroneutrality | Appendix 9 Eqs. (34)–(35) | Algebraic compartment constraints |
| CO2/HCO3−/H+ equilibrium relations | Appendix 8 | Algebraic relations for lumen and cytoplasm; their exact use with dynamic Eq. (6) needs clarification |

## Prescribed inputs and fixed quantities

| Quantity | Published status |
|---|---|
| `[Ca2+]_i(t)` | Fixed input, §2.8 Fig. 2; resting value 58 nM and stimulation from minute 6 to 12, but no numerical function or data table |
| Interstitial concentrations | Constant infinite-bath assumption, §2.1; nominal values in Table 1 |
| `ω_l` | Constant lumen volume; ratio `ω_l/ω_i0 = 0.02` in Table 2 |
| `Ψ_l` | Fixed untracked luminal osmolyte term calibrated at steady state, §2.5 and Table 1 |
| Temperature and physical constants | Table 2 |

## Variables not established as independent differential states

The paper reports or uses luminal/interstitial H+, HCO3−, and CO2, but does not include them among the ten QSS differential states. They must therefore be fixed or algebraically derived in the published formulation. The exact closure is not sufficiently explicit for implementation and remains open.
