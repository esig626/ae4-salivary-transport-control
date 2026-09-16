# Mandatory research ledger — pre-Task52 algebraic appendix

**Binding status:** this is an append-only continuation of `docs/MANDATORY_RESEARCH_LEDGER.md` created after verified Task51 closeout `75b50cdc056876ce7c8e83e0a1e2ba502a80c904`. It must be read together with the master ledger before Task52. Historical ledger entries remain unchanged and traceable.

**Current budget:** Task51 is complete. **Exactly one funded/scientific Codex shot remains: Task52.** It is reserved for numerics. No new mechanism family is authorised.

Full derivation: `analysis/pre_task52_algebraic_gate/PRE_TASK52_ALGEBRAIC_GATE.md`.

## R51G — pre-Task52 algebraic gate; no scientific compute

**Classification:** FORMAL DEDUCTION + source-scale calculations from existing frozen results. No model integration, fit, optimisation, stationary solve or parameter search was performed.

### Exact conserved-coordinate identities

With `N` NKCC, `A` AE4, `B` NBC, `H` NHE1, `E` AE2, `P` pump, `K` K-efflux, `C` CO2 influx and `J=J_TMEM16A+J_VRAC` total apical Cl export:

`Na' = N + H + B - A/2 - 3P`

`K' = N - A/2 + 2P - K`

`Cl' = 2N + A + E - J`

`TA' = H + 2B - E - 2A`

`TIC' = C + 2B - E - 2A`.

Therefore exactly

`J = 2N + H + 2B - A - TA' - Cl'`

and

`J = 6P - H + 2Na' - TA' - Cl'`.

At stationarity:

`J = 2N + H + 2B - A = 6P - H`.

For stationary AE4 KO versus WT:

`Delta J = 2 Delta N + Delta H + 2 Delta B + A_WT = 6 Delta P - Delta H`.

Thus loss of WT AE4 enters the KO-minus-WT export balance with a **positive** term. Increased NKCC compensation also pushes export upward. A sustained KO deficit requires sufficiently large opposing changes in NKCC/NBC/NHE/pump, not merely an extra apical Cl channel.

### NKCC suppression scale required for the phenotype

Using frozen Task40 WT source scales and assuming forward KO NHE/NBC (`H_KO,B_KO>=0`), the stationary lower bound `J_KO>=2N_KO` implies approximately:

- for a 35% secretion deficit: `N_KO/N_WT <= 0.859` on the broad 60–600 s scale, i.e. at least **14.1% below WT**;
- for the Task50 null ratio `0.6974`: `N_KO/N_WT <= 0.922`, at least **7.8% below WT**.

At the 600 s WT endpoint scale the requirements strengthen to approximately **17.4%** and **11.4%** below WT respectively.

Task40 instead has `N_KO/N_WT=1.232`: NKCC rises ~23.2%.

### Exact finite-window decomposition

For any window:

`Integral J = 2 Integral N + Integral H + 2 Integral B - Integral A - Delta_state(TA+Cl)`.

Task40 60–600 s KO-WT source-flux contribution is `+7.6654 fmol`, favouring KO export, while actual export difference is `-8.6279 fmol`. The sign reversal is supplied by `16.2934 fmol` of differential `(TA+Cl)` storage. Thus the inherited small deficit is primarily transient storage, not a sustained source-flux disadvantage.

A 35% deficit at the same WT broad-window scale, with the same source term, would require ~`80.23 fmol` of differential storage. This is an illustrative scale, not a hard bound under new parameters.

### Palk NKCC monotonicity

The active NKCC core is `N=alpha*M*f(x)`, `x=Na_i*K_i*Cl_i^2`, with

`f(x)=(A1-A2*x)/(A3+A4*x)`.

Exactly

`f'(x)=-(A2*A3+A4*A1)/(A3+A4*x)^2 < 0`.

So lower intracellular `Na*K*Cl^2` increases inward NKCC. At the Task40 WT endpoint, the favourable 35% threshold would require `x_KO` ~4.17% above WT; the Task50-ratio threshold ~2.28% above WT. Actual Task40 KO `x` is below WT.

The measured resting Cl phenotype (`36.5 mM` KO versus `50.1 mM` WT) reduces the Cl-squared part of `x` to ~0.531 if Na/K were unchanged. Around the Task40 operating point that would increase Palk inward flux roughly 3.8-fold. To overcome that Cl effect and raise `x` 4.17% above WT would require `Na_i*K_i` nearly 1.96-fold WT. This is a scale calculation, not a predicted KO state.

### Pump/NHE coupled sign constraint

The pump Na gate is strictly increasing in intracellular Na. The mechanistic Cha NHE1 law is strictly decreasing in intracellular Na at fixed pH/external concentrations. Therefore raising Na to suppress Palk NKCC simultaneously raises `P` and lowers `H`. Since stationary export is `J=6P-H`, **both effects increase export**, fighting the desired KO deficit.

### Osmotic sign and VRAC genotype direction

The water model counts intracellular `Na+K+Cl+TIC+fixed osmoles`. Summing source stoichiometries gives

`M' = 4N + H + 3B - 2A - P - K - J + C`.

Therefore NKCC adds `+4` tracked osmoles/cycle, but the current equal-routing AE4 positive cycle contributes **-2 osmoles/cycle**. Within this model AE4 is a net osmotic sink. At a common state, deleting AE4 therefore favours swelling by `+2A_WT`; upward NKCC compensation adds `4 Delta N` more. A genotype-independent positive-swelling VRAC gate consequently has no natural reason to be weaker in KO and can have the wrong genotype direction.

VRAC itself removes intracellular Cl and is negative feedback on swelling. Task51 confirms the sign numerically: increasing `g_scale` from 0 to `1e-6 S` reduced mean late swelling from 5.75% to 2.44%.

### Electrical redundancy of TMEM16A and VRAC

Task51 gives both apical Cl currents the same `E_Cl`. At any fixed state after common voltage closure:

`I_VRAC/(I_TMEM16A+I_VRAC)=g_V/(g_CaCC+g_V)`.

They are parallel conductances, not independent driving forces. In the Task51 high-g trial, ~87.8% of added VRAC export replaced lost TMEM16A export; total Cl export rose only ~12.1%. This is a topology diagnostic, not a physiological calibration.

### Catalan Supporting Information recovered after Task51

Task51 correctly recorded that official SI routes were inaccessible during execution. A web mirror of the Supporting Information was subsequently retrieved. It reports nystatin-perforated IPR recordings with internal `58 mM NMDG-Cl + 60 mM NMDG-glutamate` and external `90 mM NaCl +25 mM NaHCO3 +4.3 mM KCl +1 mM CaCl2 +1 mM MgCl2`, isotonic osmolality `315 mOsm/kg`.

Using `Cl_i~58 mM`, `Cl_o~98.3 mM` gives an approximate room-temperature `E_Cl~-13.5 mV`; at -60 mV the driving-force magnitude is ~46.5 mV. Combined with the previously preserved Task49 current extraction:

- total IPR current `208.6+/-35.2 pA` -> local whole-cell ohmic-equivalent ~4.49 nS;
- DCPIB-blocked component `107.8+/-17.2 pA` -> ~2.32 nS.

Mapping these onto 12.5% swelling gives Task49/51 conductance-scale sensitivities of roughly `3.6e-8 S` and `1.9e-8 S` respectively.

These are not identified apical conductances: the patch current is whole-cell, rectifying, has unknown apical fraction, and the current/swelling measurements are not guaranteed to be matched-cell simultaneous measurements. They nevertheless show that Task51's numerical upper bracket `1e-6 S` was ~30–50 times larger than the source-scale whole-cell mapping.

### Scale relative to Task50

The preserved Task50 WT combined-stimulus effective CaCC conductance is ~31.4 nS. Source-mapped VRAC at 12.5% swelling is only ~2.3–4.5 nS (~7–14% of CaCC). Even an artificially favourable WT-only VRAC comparison gives raw conductance ratios ~0.93–0.875, far weaker than the Task50 null secretion ratio 0.697 or the experimental ~0.65 ratio. Matching those ratios by conductance alone would require ~13.6–16.9 nS of WT-only extra conductance. This is a scale warning, not a secretion-ratio theorem.

### Binding Task52 numerical gates

Task52 must use Codex only for the numerical pieces the algebra cannot settle:

1. constrain `g_scale` from electrophysiological source mappings/sensitivities, never from AE4 secretion;
2. with `g_scale` fixed, infer at most one beta-NKCC effective gain from independent mouse WT IPR swelling, and report its relation to the 2.5–3-fold salivary evidence / sixfold context;
3. freeze before AE4 phenotype comparison;
4. report `N,H,B,P,A,J`, `TA'+Cl'`, osmotic-source terms and VRAC/TMEM16A conductances for WT and KO;
5. decompose any deficit with the exact transient/integrated identities above;
6. if KO NKCC is not suppressed by the required order (~8–17% below WT on relevant scales), or if the deficit is primarily finite storage, the literature mechanism does not explain the sustained phenotype in this chassis;
7. if KO swelling/VRAC activation is at least WT, record the wrong genotype direction explicitly;
8. no new mechanism, no Task50 multiplier, no phenotype retuning.

If these gates fail, Task50 remains the quantitative proof that an additional unidentified beta-conditioned network coupling is required.

## References

- Catalan MA et al. *PNAS* 2015. DOI `10.1073/pnas.1415739112`; PMID `25646474`; PMCID `PMC4343136`.
- Catalan Supporting Information mirror retrieved 2026-09-16: `https://esdocs.com/doc/477676/download-supporting-information--pdf-`.
- Paulais M, Turner RJ. *J Clin Invest* 1992. DOI `10.1172/JCI115695`; PMID `1313447`; PMCID `PMC442971`.
- Chaib N et al. 1998. DOI `10.1016/S0196-9781(98)00134-X`; PMID `9880083`.
- Tanimura A et al. 1995. DOI `10.1074/jbc.270.42.25252`; PMID `7559664`.
- Kurihara K et al. 1999. DOI `10.1152/ajpcell.1999.277.6.C1184`; PMID `10600770`.
- Kurihara K et al. 2002. DOI `10.1152/ajpcell.00352.2001`; PMID `11880270`.
- Pena-Munzenmayer G et al. *JBC* 2015. DOI `10.1074/jbc.M114.612895`; PMID `25745107`; PMCID `PMC4409235`.
- Pena-Munzenmayer G et al. *J Gen Physiol* 2016. DOI `10.1085/jgp.201611571`.
- Pena-Munzenmayer G et al. *AJP GI* 2021. DOI `10.1152/ajpgi.00145.2021`; PMID `34585968`; PMCID `PMC8887885`.

Repository provenance: Task40 frozen flux/state outputs; Task49 conserved VRAC implementation; Task51 verified result at `75b50cdc056876ce7c8e83e0a1e2ba502a80c904`; Task50 preserved benchmark branch `archive/task-50-working-effective-coupling-benchmark` at `b16c30094b95f61a73d8f1977cd58e79c7bb50f6`.
