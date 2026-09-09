# Baseline observables

## Primary observable

The model’s primary output is total primary-saliva flow,

\[
q_{tot}=q_a+q_t,
\]

defined in §2.5 Eq. (13). With constant lumen volume, the paper treats this as both total water entering the lumen and water leaving toward the duct. Its dimensional unit follows the water-permeability convention and must be confirmed before implementation.

Published comparison locations:

- Fig. 3: control fluid flow during stimulation.
- Fig. 8: control, Ae2-knockout, and Ae4-knockout flow; caption reports about 24% lower flow for Ae4 knockout and no significant Ae2 change.
- §3.2: the paper notes that the simulated Ae4 response reaches its new state faster than the experiment in Peña-Münzenmayer et al. (2015).

## State-derived observables

| Observable | Definition/location | Baseline target location |
|---|---|---|
| Intracellular ion concentrations | Divide intracellular amounts by `ω_i`; Eqs. (1)–(6) | Table 1; stimulated time courses in Figs. 4 and 6 |
| Luminal Na+, K+, Cl− | Direct concentration states; Eqs. (14)–(16) | Table 1; stimulated time courses in Figs. 5 and 6 |
| Cell volume | `ω_i`, governed by Eq. (9) subject to unresolved sign | Table 1; stimulated time course in Fig. 3 |
| Apical potential | Algebraic `V_a` in QSS model | Table 1; stimulated/knockout time courses in Fig. 7 |
| Basolateral potential | Algebraic `V_b` in QSS model | Table 1; stimulated/knockout time courses in Fig. 7 |
| Transepithelial potential | `V_t = V_a - V_b`, Appendix 10 | Derived; no separate numeric table row |
| Intracellular pH | Presumably derived from `[H+]_i` | Table 1; exact conversion/units not stated |
| Compartment osmolarity | Sums defined after Eq. (11), plus `x_i/ω_i` or `Ψ_l` where applicable | §2.5 reports 292.6 mM interstitium, 296.6 mM cell, 297.4 mM lumen |

## Mechanistic observables

The individual transporter, pump, channel, buffer, and tight-junction fluxes are model outputs once the states and algebraic potentials are known. Fig. 9 reports selected Ae2, Ae4, Nkcc1, and Nhe1 flux comparisons under control and knockout conditions. These are useful regression outputs but are not independent experimental validation targets in the paper.

## Required reproducibility metadata for Phase 01

Any machine-readable result should include:

- explicit units and sign convention for each flow/current/flux;
- parameter-set version and source status (direct, inherited, or paper-derived);
- complete initial state and algebraic-consistency residuals;
- numerical definition of `[Ca2+]_i(t)`;
- solver, tolerances, time grid, and DAE initialization;
- mass-balance and electroneutrality residuals;
- definition of knockout (which density/activity is set to zero);
- normalization used for flow-rate comparisons.

No baseline value is accepted into `docs/RESULTS_LEDGER.md` in Phase 00; this phase records sources and blockers only.
