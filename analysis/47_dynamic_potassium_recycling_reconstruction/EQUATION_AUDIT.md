# Potassium recycling equation audit

The active system is Task 40 equal AE4 routing, the Palk NKCC law, Cha NHE1 and the minimal electrogenic NBC wrapper. The generic tanh NKCC law and the older membrane closure without NBC are not the production laws. The Task 44 parameter payload is copied unchanged into `input/active_parameters.json`; `common.py` checks all six parameter hashes and the WT initial state against Task 40 before evaluation.

## Units, signs and exact balances

Amounts are fmol, volumes pL, concentrations mM, time s, voltages V, conductances S and currents A. Thus `n/V` is mM and `j = 1e15 I/(z F)` is fmol/s. Positive K channel flux is cell to lumen or bath; positive pump flux exports 3 Na and imports 2 K; positive NKCC imports 1 Na, 1 K and 2 Cl. Positive AE4 exports half a Na and half a K and imports one Cl per cycle, with two TIC and two TA equivalents exported. Charge is `Na+K-Cl-TA`; TIC is not itself a charge coordinate.

Write N for NKCC cycles, H for NHE1, B for NBC, E for AE2, A for AE4, P for both pumps, J for apical chloride export, and K for both K channel effluxes. The active cellular rows are

```
dn_Na/dt = N + H + B - A/2 - 3P
dn_K/dt  = N - A/2 + 2P - K
dn_Cl/dt = 2N + E + A - J
dn_TA/dt = H + 2B - E - 2A.
```

Concentration dynamics include dilution: `d[K_i]/dt = (dn_K/dt - [K_i] dV_i/dt)/V_i`, likewise for Na. A concentration change is not a change in the same amount when volume changes. There is no fixed intracellular K reservoir or hidden K source. The lumen has `dn_K,l/dt = -2P_a + K_a - K_para - Q_out[K_l]`, and analogous Na export and washout. Luminal K therefore feeds back on the apical reversal potential and pump.

Eliminating the other rows gives the exact dynamic extension of the CBM clue:

```
K = J/3 + N + H/3 - A/2 + S
S = (-2 dn_Na/dt + dn_TA/dt + dn_Cl/dt)/3 - dn_K/dt
  = -(dn_Na/dt + 2 dn_K/dt)/3
```

The last equality uses charge conservation. S need not vanish during stimulation. Only at zero storage does this reduce to the Task 46 identity. The runner will check both forms against independently assembled channel fluxes. Neither form changes the ODE.

## NKCC1

`N = alpha*m*(a1-a2 X)/(a3+a4 X)`, `X=[Na_i][K_i][Cl_i]^2`.
The coefficients are `a1=157.5`, `a2=2.0096e-5 mM^-4`, `a3=1.0306`, `a4=1.3852e-6 mM^-4`, with the source convention treating the ratio as dimensionless. `alpha=0.017334746894096052 fmol/s` matches an inherited WT resting flux; it is not a measured capacity. The calcium recruitment m is 1 at rest and 1.75 during the declared stimulus. The generic capacity 0.32 fmol/s is inactive.

The exact derivative `dN/dX = -alpha*m*(a2*a3+a1*a4)/(a3+a4*X)^2` is negative. K feedback is already present. Na has the same log concentration derivative, and Cl has twice that derivative. NKCC is electroneutral, so no explicit voltage term belongs in its ideal transport affinity.

The diagnostic called affinity in this selected law is `log(a1/(a2 X))`. The ideal bath affinity is separately `log([Na_e][K_e][Cl_e]^2/X)`. These are not identical at the active bath: the source coefficients encode a different fixed bath reversal product. This is a real limitation of transferring the inherited fixed bath law; agreement of the current law with its own diagnostic is not proof of full bath thermodynamics. Record both signs along the trajectories. A wholesale kinetic rederivation cannot be replaced by changing one numerator coefficient without independent rate information. In particular, this issue does not identify a K conductance or pump ceiling.

## K channels and membrane voltage

For membrane r, `E_K,r=(RT/F) log([K_r]/[K_i])`, `I_K,r=g_K,r*(V_r-E_K,r)` and `K_r=1e15 I_K,r/F`. This has the correct reversal, current sign and nonnegative passive dissipation `I_K,r*(V_r-E_K,r)`.

`g_K,a = f_K*g_K_total*h(Ca) + g_background,a` and `g_K,b = (1-f_K)*g_K_total*h(Ca) + g_background,b`, where `h=Ca^p/(Ca^p+K_Ca^p)`. Active values are total 14 nS, f_K=0.3, K_Ca=0.26 micromolar, p=1.46 and zero backgrounds. These are effective inherited settings, not measured reserve. Fractions partition conductance, not realised flux. The law lumps channel populations and has no explicit BK voltage gating or dynamic inactivation; that alone does not prescribe a unique replacement for the mixed BK/IK current.

There is no membrane voltage ODE or stored capacitive charge. Both voltages are algebraic functions of the evolving amounts, volume and calcium. With conventional outward membrane current, the closure is `I_a-I_para=0`, `I_b+I_para=0`, where `I_a=I_K,a+I_Cl,a+F*1e-15*P_a`, and `I_b=I_K,b+F*1e-15*(P_b+B)`. NBC current is essential here. Summing yields the additional instantaneous identity `K+P+B=J`. Omitting NBC would assign its electrical effect incorrectly to the pump or K channels.

The voltage Jacobian is `[[g_K,a+g_Cl+g_para,-g_para],[-g_para,g_K,b+g_NBC+g_para]]`, with `g_NBC=F*1e-15*capacity*u*sech(affinity/width)^2/(width*RT/F) >= 0`. Positive channel conductances make the closure locally nonsingular. The executable bracket of plus or minus 0.5 V is numerical and is not a physiological range.

## Na/K pump

`P_r = C_r * [Na_i]^3/([Na_i]^3+K_Na^3) * [K_r]^2/([K_r]^2+K_K^2)`.
The total nominal C is 0.08 fmol/s; `C_a=f_P*C`, `C_b=(1-f_P)*C`, with f_P=0.075075, K_Na=10 mM and K_K=1.5 mM. Both pumps share intracellular Na and use their own extracellular K. For a fixed instantaneous external K state the upper limit as Na increases is `C_eff=sum C_r*g_Kout,r`; utilisation is `P/C_eff=g_Na`, not just P/C. Na clearance is 3P and pump K influx is 2P. The local elasticity is `d log P/d log Na = 3*(1-g_Na)`.

The pump is bounded but voltage independent and has no ATP, ADP, phosphate or reversal law. Thus the proposed voltage to pump arrow has no direct implementation; indirect coupling remains through Na, luminal K and the other transporters. This is a model limitation, not an independently measured correction. The retained topology supports apical pumps, but the functional fraction is uncertain.

## Provenance and admissible conclusions

`parameter_provenance.csv` retains Task 43 classifications and limitations for the recycling block and its state, electrical and water dependencies. Published kinetic coefficients remain literature derived; physical constants remain constants; inherited model values, WT constrained scales and unconstrained effective parameters remain distinct. No parameter in this audit is promoted to a direct measurement. Frozen initial coordinates are inherited solved states, not observations. There is no justified numerical uncertainty interval for the effective pump capacity, K conductance, localisation or gating parameters.

Primary experimental context supports both BK and IK channels in the apical domain, but does not identify the present 30% whole cell functional fraction or a matched WT/null sustained capacity. See [Almassy et al. (2012)](https://doi.org/10.1085/jgp.201110718). Regional channel density is not the same quantity as integrated whole membrane conductance.

The audit finds no missing K amount, sign, pump stoichiometry or channel driving force term in the active equations. It does identify incomplete kinetic and thermodynamic information, especially the fixed bath NKCC transfer and the effective channel and pump laws. Whether the current reserve is biologically too large remains an evidence question. Baseline diagnostics must be published before deciding whether a specific correction is justified.
