# Phase 00 open questions

These issues block a faithful, independently testable baseline implementation. They must be resolved from published sources, author clarification, or an explicit new modelling decision. Historical MATLAB or unpublished derivations are not acceptable shortcuts.

## Critical equation and sign ambiguities

1. **Cell-volume sign.** Section 2.5 Eq. (9) prints `dω_i/dt = q_a - q_b`, while the collected system in §2.7 prints `dω_i/dt = q_b - q_a`. Which is intended?
2. **Water-flux labels.** Eqs. (10)–(11) make `q_a` depend on lumen/cytoplasm and `q_b` on cytoplasm/interstitium, but the prose immediately after Eq. (9) describes the membrane associations inconsistently. Define directions and membrane labels before coding.
3. **CO2 net-flux sign.** Appendix 8 prints `J_CO2 = P_CO2(2[CO2]_i - [CO2]_l - [CO2]_e)`, yet describes diffusion into the cytoplasm and Eq. (6) adds `J_CO2`. Confirm whether the printed sign or the verbal direction is intended.
4. **Luminal chloride loss term.** Eq. (16) prints `q_tot[Cl−]` without the lumen subscript; §2.7 supplies `[Cl−]_l`. Record any correction explicitly.
5. **NaK reaction stoichiometry typo.** Appendix 2’s net-reaction display appears to produce one intracellular K+, whereas the text and balances use two. Confirm against Palk et al. (2010) and Smith and Crampin (2004).
6. **Ae4 intermediate-state equations.** Check Eqs. (27)–(30) against Fig. 10 and an independent Markov derivation before relying on them; the simplified Eq. (31) may conceal typographical state-transition errors.

## Units, scaling, and geometry

7. **Current versus molar flux.** The main balances use `I/(Fz)`, appendices sometimes name expressions `J` after already dividing by `F`, and the potential balances use `I`. Establish one dimensional ledger for every term.
8. **Ae4 units.** Table 7 gives `G_Ae4` in amol/μm^3, prose after Eq. (31) says fmol, and describes `J_Ae4` as concentration/time, while Eqs. (1)–(4) require amount/time.
9. **Transporter density scaling.** `α_Nkcc1` and `α_NaK` are reported in amol/μm^3 and called membrane densities. Determine whether they multiply cell volume, membrane area, or are already whole-cell amounts.
10. **Tight-junction notation.** Appendix 10 uses `J` and lowercase `g`; Eqs. (7)–(8), (14)–(15), and Table 2 use `I` and uppercase `G`. Confirm equivalence and units.
11. **CO2 reverse-rate units.** Table 8 reports `k_-1` as s^-1 although `J_Buffer` multiplies it by two concentrations. Check Sharp et al. (2015) for the intended units and any concentration normalization.
12. **Membrane capacitance.** `C_m` appears in Eqs. (7)–(8) but is absent from Table 2. It is unnecessary only after adopting QSS; the full dynamic model cannot be reconstructed.
13. **Cell geometry.** The paper reports cell and lumen volumes but no explicit membrane areas, although several parameters are called densities. Determine the geometry/scaling used to obtain whole-cell fluxes.

## State closure and initial conditions

14. **Lumen acid-base closure.** The ten-state QSS system has no lumen H+, HCO3−, or CO2 differential state. Specify which are fixed and which are algebraically derived from electroneutrality/equilibrium.
15. **Interstitial acid-base values.** Table 1 marks `[HCO3−]_e` and `[CO2]_e` as model-derived despite treating the interstitium as a fixed bath. Identify an independent bath-composition source or label the values as new assumptions.
16. **pH definition and units.** The paper does not explicitly state the conversion from `[H+]` to pH in its unit system. Table 1 also appends “mM” to the model `pH_l` entry, apparently a typographical error.
17. **Impermeant charge `x_i`.** It is solved from resting electroneutrality and then enters osmotic flux. Confirm whether it remains fixed in moles under changing cell volume and how its valence is represented when `z_x ≤ -1`.
18. **Luminal osmolyte term `Ψ_l`.** It is selected by a steady-state volume/flow condition, but the target, fitting equation, and uniqueness are not fully reported.
19. **Initial algebraic consistency.** A full initial state and the method used to make `V_a`, `V_b`, electroneutrality, and acid-base relations consistent are not published as a machine-readable set.

## Calibration and numerical reproducibility

20. **Model-derived parameters.** Channel/tight-junction conductances, transporter densities/activities, all Ae4 parameters, `x_i`, `Ψ_l`, and several concentrations are calibrated in the paper. The fitting objective, weights, bounds, parameter correlations, and uniqueness are not reported.
21. **Calcium input.** Fig. 2 gives a qualitative prescribed `[Ca2+]_i(t)` with rest at 58 nM and stimulation during minutes 6–12, but no equation or data table. Obtain a citable numerical input.
22. **DAE details.** The paper names MATLAB `ode15s` but does not report tolerances, maximum step, consistent-initialization method, Jacobian treatment, or output grid.
23. **Baseline validation tolerances.** Table 1 mixes exact reused values, approximate experimental values, and paper-derived values. Define acceptance tolerances without treating fitted outputs as independent validation.
24. **Knockout definition.** Confirm exactly which parameter is set to zero for Ae2 and Ae4 knockout and whether any compensatory expression change is imposed or emerges only from state changes.
25. **Flow units and normalization.** The precise units used for `q_tot`, and the normalization behind plotted flow comparisons, must be recovered before regression testing.

## Literature gaps

26. Obtain and audit the upstream published model sources used for Nkcc1/NaK (Palk et al. 2010), Ae2/Nhe1 (Falkenberg and Jakobsson 2010), buffer chemistry (Sharp et al. 2015), and channel gating (Takahata et al. 2003; Arreola et al. 2002).
27. Obtain the primary experimental sources behind Table 1, especially those for intracellular ions, membrane potentials, lumen composition, and volume, and determine whether gland/cell type and stimulation conditions match the baseline.
28. Audit the 2016 Ae4 transport paper against the 1 Cl− : 2 HCO3− : 1 monovalent-cation cycle used in the model, including whether Na+ and K+ competition should be represented by simple concentration fractions.
29. Identify whether a published correction, supplement, repository, or later paper resolves the sign, unit, or calibration ambiguities above.

## Phase boundary

Phase 00 does not resolve these by inspecting historical code, inventing values, or choosing preferred equations. Phase 01 may begin only after the critical signs, dimensions, state closure, and numerical input are documented in `model/` with published provenance or explicitly approved new assumptions.
