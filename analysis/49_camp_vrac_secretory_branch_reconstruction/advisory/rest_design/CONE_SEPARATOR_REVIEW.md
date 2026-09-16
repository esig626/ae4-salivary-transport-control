# Advisory algebra review of the proposed cone certificate

The proposed certificate is correct **conditional on the lifted AE4-KO state having `H=J_NHE1>0` and `J_2=J_AE2<0`**, with AE4 deleted and resting NBC/VRAC unavailable. No model evaluation, source matrix construction, rank calculation or solver was run for this review.

Let `ell` select the intracellular TA rate from the full conserved vector, let `C` be the single predeclared source matrix oriented in the allowed addition directions, and let `b=-F(x)` be the required correction. Then

`ell^T C >= 0`, while `ell^T b=-(H-J_2)<0`.

This is an exact Farkas separating functional: no `z>=0` can satisfy `C z=b`. It remains valid after imposing current closure or tying coefficients across genotypes, since those restrictions only shrink the admissible cone. On the four-cohort stacked matrix, use the same functional on the AE4-KO TA row and zero on every other row. No second matrix construction or optimisation is needed. Positive rescaling of rows does not change the sign certificate.

Be precise about two distinct claims:

1. **Addition cone:** no nonnegative addition of the allowed, correctly oriented native fluxes cancels the positive TA residual.
2. **Complete capacity reweighting at this fixed chemical state:** with nonnegative capacities, TA balance reads `a_H H+a_2 |J_2|=0`. It is impossible whenever NHE remains nonzero (or AE2 remains nonzero). Setting both activities to zero is the trivial degenerate escape, not a physiological reconstruction. This stronger claim does not follow merely by calling negative corrections forbidden; it follows from the signs of the *total* reweighted fluxes. It assumes the kinetic affinity and state are fixed.

An inward NBC addition contributes `+2J_B` to the same TA balance and fails the same certificate. If the computed NBC affinity is inward at this state, increasing any nonnegative basal NBC capacity therefore cannot fix the discrepancy. A desired outward NBC term would require a changed state or driving force; it is not supplied by the same law with positive capacity at an inward affinity. TMEM16A/VRAC, CO2, pump, K, NKCC and water have zero **direct** cell-TA signatures, though changing their dynamics can move the state and hence the affinities.

## Useful exact necessary relations

For AE4 KO the native resting equations are

`H+2B=J_2`,

`J_CO2,b+J_CO2,a=H`,

where `B=0` in the inherited gate. With the frozen AE2 law `J_2=G_2 tanh(A_2/w)`, finite positive concentrations imply

`0<H=J_2<G_2` when `B=0`,

and `H+2B<G_2` for inward `B>=0`.

Thus `H<=G_2` is a necessary weak bound. Its use must state whether `G_2` is the frozen capacity; if AE2 capacity is itself reweighted, the corresponding bound uses that new capacity. The negative-affinity sign obstruction is stronger and cannot be repaired by any positive capacity increase at the fixed state.

The AE2 sign condition can also be written entirely explicitly:

`A_2>0` iff `HCO3_i > Cl_i HCO3_o/Cl_o`.

Since the implemented carbonate equilibrium gives `HCO3_i=CO2_i 10^(pH_i-pKa1)`, define

`Ccrit=(Cl_i HCO3_o/Cl_o) 10^(pKa1-pH_i)`,

`D=P_CO2,b+P_CO2,a`,

`Cbar=(P_CO2,b CO2_b+P_CO2,a CO2_l)/D`.

Stationary TIC minus TA balance fixes `CO2_i=Cbar-H/D`, so native `B=0` requires the scalar compatibility equation

`H = G_2 tanh(log((Cbar-H/D)/Ccrit)/w)`.

In particular, `Cbar>Ccrit` and

`0<H<min(G_2, D(Cbar-Ccrit))`

are necessary. This is an equation-level test, not another pathway search. It also makes the latent-state dependence explicit: the reference lumen CO2 enters `Cbar`, and NHE turnover depends on unmeasured intracellular Na. No parameter values or signs were evaluated here.

## Scope of the conclusion

The source *span* can contain `b` while the admissible *cone* does not. Report that the active native flux directions at the declared state cannot provide the required negative intracellular-TA rate. Do **not** infer a new chemical species, a unique transporter, or a globally impossible shared model solely from this cone result.

Cl/pH do not determine intracellular Na/TIC/volume or the lumen state. Changing these unmeasured coordinates can change AE2/NBC affinity while keeping the measured Cl/pH unchanged. Therefore the separator rules out fixed-state capacity-only rescue and any joint shared solution constrained to those same lifts; it does not rule out every observation-compatible resting state. Its rigorous deliverable is the explicit TA sign certificate and the compatibility equation above, with measured and nuisance uncertainty attached. Only a sign certificate that survives all admissible nuisance states could support a global impossibility conclusion.

