# Published-paper / historical-code concordance

## Scope, evidence levels, and status vocabulary

This is a forensic comparison, not a clean implementation. The published
source is Vera-Sigüenza et al. (2018), DOI
[`10.1007/s11538-017-0370-6`](https://doi.org/10.1007/s11538-017-0370-6),
using the complete public author manuscript at
[PMCID PMC5792321](https://pmc.ncbi.nlm.nih.gov/articles/PMC5792321/).
Historical implementation evidence comes from the seven immutable files under
`archive/legacy-2017/Ae4_Dynamics_Project/`. Numerical statements below marked
as recomputations were independently evaluated from those archived expressions
and the values stored in `Par.mat`; they are not published observations. The
corresponding machine-readable evidence is
`results/11_forensic_reconstruction/reproduction_summary.json`,
`reproduction_residuals.csv`, and `par_mat_inventory.json`.

The tables use the required categories literally:

- **exact match** — the mathematical expression and numerical convention agree;
- **equivalent after algebra** — the code eliminates a variable or changes a
  sign representation without changing the relevant equation;
- **publication typo likely** — code plus conservation or an upstream
  correction supplies specific evidence for a printed error;
- **implementation convention** — a coherent historical choice absent from the
  published executable specification;
- **calibration choice** — a fitted historical constant or transformation makes
  the archived baseline close, but is not transparently reported in the paper;
- **substantive model difference** — the archived code implements different
  dynamics or transport stoichiometry; and
- **unresolved** — neither record establishes a unique commensurate convention.

## Executive finding

The archived programs are directly related to the published project, but no
archived MATLAB variant is equation-identical to the 2018 model. The strongest
links are the exact resting voltages and major ion concentrations, the geometry
that gives exactly 1.3 pL, the same transporter functional families, and several
near-exact transformations from `Par.mat` back to published fitted parameters.
The decisive non-identity is Ae4: every archived balance and flux law is
**Na-specific**, while the published model transports a nonselective Na/K
cation and partitions the Ae4 cycle between both intracellular cation
balances. At the archived resting point, merely restoring the published
partition creates equal-and-opposite Na and K residuals of
`0.0321793147` in the archived internal flux units.

The outer program also has seven differential states, eliminates intracellular
H⁺ and CO₂, and uses a one-step calcium protocol. The paper reports ten
differential equations and a multi-stage calcium input. H⁺ elimination can
be derived from electroneutrality and QSS charge balance, but replacing the
CO₂ differential balance by an algebraic steady relation changes the
transient model. The historical code nevertheless explains several otherwise
hidden fitting and unit conventions. It is therefore strong ancestry and
implementation-provenance evidence, not a drop-in executable specification of
the published model.

## Historical variants compared

| Historical file | Mathematical role relevant to this concordance |
| --- | --- |
| `Original/Parameters.m` | Defines geometry, the archived resting state, raw kinetic constants, and hand scaling factors. |
| `Original/Saliva_Ae4.m` | Calibrates conductances/activities so its chosen resting state has essentially zero residual; records a 12-row diagnostic but calls a reduced seven-state integrator. |
| `Par.m` / `Par.mat` | Stores the calibrated parameter structure used by all outer functions. The MAT values agree with `Par.m`. |
| `Salivary.m` | Seven-state outer model: luminal Na/K, cell height, intracellular Na/K/Cl/HCO₃; potentials, H⁺, CO₂, and luminal Cl are algebraic. |
| `Salivary2.m` | Five-state reduction with fixed luminal ions and a different global right-hand-side multiplier. |
| `Salivary_ex.m` | Runs an Ae2 knockout (`g2=0`) and a reduced WT simulation; no archived script runs all three published WT/Ae2-KO/Ae4-KO cases. |

Line references below are to the immutable archived sources.

## Equation-by-equation concordance

### Balance equations, constraints, and water transport

| Published item | Historical implementation | Assessment | Status |
| --- | --- | --- | --- |
| Eq. (1), intracellular Cl amount | `Salivary.m:80` and `Original/Saliva_Ae4.m:117` contain `2 Nkcc1 + Ae2 + Ae4 + signed CaCC`, followed by the product-rule volume dilution term. | With historical `jcl=I_CaCC/F` and `z_Cl=-1`, adding the negative electrical-equivalent Cl flux is algebraically the same as `-I_CaCC/(F z_Cl)`. | **equivalent after algebra** |
| Eq. (2), intracellular Na amount | `Salivary.m:78` and `Original/Saliva_Ae4.m:115` subtract the **entire** Ae4 cycle from Na. | The paper subtracts `Na_i/(Na_i+K_i)` of Ae4. At the archived baseline the missing K fraction is `120/145`, giving a `+0.0321793147` Na residual if the published partition is imposed. | **substantive model difference** |
| Eq. (3), intracellular K amount | `Salivary.m:79`, `Salivary2.m:72`, and `Original/Saliva_Ae4.m:116` have no Ae4 term. | The paper subtracts `K_i/(Na_i+K_i) J_Ae4`; the omitted archived term is `0.0321793147` at rest. This is not an algebraic reduction. | **substantive model difference** |
| Eq. (4), intracellular HCO₃ amount | `Salivary.m:82` and `Salivary2.m:75` use `Buffer - 2 Ae4 - Ae2` and the product-rule volume term. | Stoichiometry agrees, subject to the different Ae4 law and historical units. | **equivalent after algebra** |
| Eq. (5), intracellular H⁺ amount | No H⁺ ODE in either outer function. H⁺ is reconstructed from electroneutrality at `Salivary.m:15`. | Differentiating cell electroneutrality and using the other ion balances plus QSS current balance recovers `d(H omega)/dt=Buffer-Nhe1`; it is a redundant charge row under those constraints. | **equivalent after algebra** |
| Eq. (6), intracellular CO₂ amount | No CO₂ ODE. `Salivary.m:16` solves `J_CO2=J_Buffer` algebraically. | This is the steady zero set of the historical CO₂ row, not the published dynamic equation. The coefficient `1.99997` is also undocumented. | **implementation convention** |
| Eqs. (7)–(8), two membrane potentials | `Salivary.m:28-40` solves the two linear QSS current equations analytically at every RHS evaluation. | The paper explicitly adopts QSS potentials for simulations. Given the historical current convention, the analytic elimination is exact. | **equivalent after algebra** |
| Eq. (9), cell volume | Every archived dynamics function uses `q_b-q_a` (`Salivary.m:77`; `Salivary2.m:70`; diagnostic `Original/Saliva_Ae4.m:114`). | This agrees with the collected published system and conservation, not the opposite sign printed in Eq. (9). | **publication typo likely** |
| Eq. (10), apical water flux | After cell/lumen electroneutrality, `Salivary.m:53` has the same reduced ionic structure. | The code omits the tiny luminal H⁺/HCO₃ pair and uses a fitted `b1`; structurally equivalent, numerically not the published `P_a`. | **calibration choice** |
| Eq. (11), basolateral water flux | `Salivary.m:54` uses `2(Na_i+K_i+H_i)+CO2_i-Ie`. | Published algebra gives the same cell term, but its bath sum includes H⁺ and CO₂. Historical `Ie=288.1` omits them and uses HCO₃ₑ=40 mM rather than 42.9 mM. | **substantive model difference** |
| Eq. (12), tight-junction water flux | `Salivary.m:55` uses the corresponding reduced lumen-minus-bath expression. | It shares the historical bath omission and its own fitted `b3`. | **calibration choice** |
| Eq. (13), total flow | `Salivary.m:56` and `Salivary_ex.m:118` use `qT=qa+qt`. | Definition agrees exactly; absolute units are not recorded in the code. | **exact match** for the sum; **unresolved** for units |
| Eqs. (14)–(15), luminal Na/K | `Salivary.m:74,76` use tight-junction particle flux minus `qT` times lumen concentration. | The code has no explicit `omega_l` divisor. Its fluxes are pre-scaled and the entire RHS is multiplied by `10^4`; the algebraic form matches but the volume/time scale is hidden. | **implementation convention** |
| Eq. (16), luminal Cl | No independent Cl-lumen state; `cll=nal+kl` at `Salivary.m:17`. | With lumen electroneutrality, equal H⁺/HCO₃, the two lumen cation balances, and apical current QSS, the Cl equation is redundant. The paper's missing `_l` on the loss term is restored by its collected system. | **equivalent after algebra**; printed subscript is **publication typo likely** |
| Ten-state DAE claim in §2.9 | `Salivary.m` integrates seven ODEs; `Salivary2.m` integrates five. Potentials are analytically eliminated, H⁺ and CO₂ are algebraic, and luminal Cl follows electroneutrality. | H⁺ and luminal Cl can be eliminated as redundant charge rows under QSS. The extra CO₂ QSS assumption is not an algebraic identity and changes transients. | **implementation convention** with a CO₂ **substantive model difference** |
| Cell electroneutrality, Eqs. (32)–(33) | `h=cl+hco+Xi/hi-na-k` at `Salivary.m:15`. | Exact formula when `Xi/hi=x_i/omega_i`. Historical `Xi/hi=85.0001230` mM differs from the published 82.8 mM because the code uses HCO₃ᵢ=10 rather than 12.1 mM. | **equivalent after algebra** plus **calibration choice** |
| Lumen electroneutrality, Eq. (34) | `cll=nal+kl`. | Exact only after neglecting or cancelling H⁺ₗ and HCO₃ₗ. At the printed resting values their mismatch is only `4.88e-6` mM. | **equivalent after algebra** |
| Bath electroneutrality, Eq. (35) | Not enforced. `Parameters.m:47` uses HCO₃ₑ=40 mM; exchanger laws hard-code 21 mM. | Published electroneutrality requires approximately 42.90004 mM. | **substantive model difference** |
| Product rule for variable volume | Intracellular code rows subtract `(q_b-q_a)[X]_i` and divide by height. | This is the concentration form of `d(omega_i[X]_i)/dt`, provided cell area is constant and historical fluxes are areal. | **equivalent after algebra** |

### Constitutive fluxes and numerical method

| Published item | Historical implementation | Assessment | Status |
| --- | --- | --- | --- |
| Nkcc1, Eq. (17) | Same rational function at `Salivary.m:60`. `Par.mat` stores `a2=2.0096e-5`, `a4=1.3852e-6` for mM states. | These are exactly the Palk-corrigendum M-based coefficients converted for mM numerical states. Literal use of the paper's `2.0096e7` and `1.3852e6` with mM states reverses the flux. | **publication typo likely** for the printed units/numerical pairing; density scale is **calibration choice** |
| NaK, Eq. (18) | Same rational topology, with `s=0.001` inserted explicitly (`Original/Saliva_Ae4.m:43`; `Salivary.m:19-20`). | This is an M-based concentration convention. The paper labels `r` and `alpha1` as mM-based. There is no published NaK corrigendum and the fitted density conversion remains hidden. | **implementation convention** and **unresolved** scale |
| NaK reaction display | All archived K balances add `2 J_NaK`, and the voltage solve treats one pump cycle as one net outward positive charge. | The paper's Appendix reaction display shows only one intracellular K product, while its prose, Eqs. (2)–(3), and historical code require two. | **publication typo likely** |
| CaKC gate, Eq. (19) | Code uses `1/(1+(K/Ca)^eta)` rather than `(Ca/(Ca+K))^eta`; `K=100.3`, `eta=1.7`. | For `eta != 1` the gates are not algebraically equal. Paper uses `K=0.182 microM`, `eta=2.54`. | **substantive model difference** |
| CaKC current/Nernst, Eqs. (20)–(21) | `vk=RTF log(K_e/K_i)` and `jk=gk Pk(V_b-vk)/F`. | Functional sign and Nernst relation agree. Conductance and gate do not. | **equivalent after algebra** for sign; **calibration choice** for magnitude |
| CaCC gate, Eq. (22) | Same alternative Hill form; `K=100.26`, `eta=1.49`. | Paper uses `(Ca/(Ca+K))^eta`, `K=0.26 microM`, `eta=1.46`. | **substantive model difference** |
| CaCC current/Nernst, Eqs. (23)–(24) | Code stores the positive magnitude `vcl=RTF log(Cl_l/Cl_i)` and evaluates `V_a+vcl`. | Since the published Cl Nernst potential has `z_Cl=-1`, `V_a+vcl=V_a-V_Cl`; this sign transformation is exact. | **equivalent after algebra** |
| Ae2, Eq. (25) | Same bidirectional Michaelis-Menten difference at `Salivary.m:65-66`. | Historical `KB=1e-4` mM versus printed `10^4` mM, and historical external HCO₃ is hard-coded at 21 mM versus published 42.9 mM. Both code layers use `1e-4`, supporting a missing minus sign in the paper, but the external concentration remains different. | `K_B`: **publication typo likely**; bath value: **substantive model difference** |
| Nhe1, Eq. (26) | `Salivary.m:68-69` uses first powers of both proton saturation terms. | The paper squares both proton terms. The values `KH` and `KNa` match, but the kinetic law does not. | **substantive model difference** |
| Ae4 Markov model, Eqs. (27)–(30) | No archived code contains the four Markov states or their denominator. | Neither archived full nor reduced program implements the published Markov system. The paper's own Eq. (29)/Eq. (31) reduction remains questionable, as recorded in Task 10. | **substantive model difference**; printed derivation remains **unresolved** |
| Simplified Ae4, Eq. (31) | Code uses only Na on each side: forward proportional to `Cl_e HCO3_i^2 Na_i`, reverse to `Cl_i 21^2 Na_e`. | The paper uses `(Na+K)` on each side. Code also uses effective `kminus=0.0015930588`, while the paper gives `1.3e-5`; the original code constructs the former by multiplying `1.3418407e-5` by `118.7219`. | **substantive model difference** and hidden **calibration choice** |
| Ae4 cation balance partition | Full Ae4 leaves through Na only; K has no Ae4 term in all archived variants. | This contradicts the central published nonselective-cation mechanism at baseline `K_i=120` mM. | **substantive model difference** |
| CO₂ membrane transport, Appendix 8 | Code follows the printed aggregate sign `P(about 2 CO2_i - CO2_l - CO2_e)` and replaces 2 by 1.99997. | It does not follow the sum of the two individually printed inward Fickian fluxes. The historical implementation confirms which printed branch it used but does not make that branch conservation-consistent. | **implementation convention**; conservation issue **unresolved** |
| Buffer, Appendix 8 | Same mass-action difference, but both rates are multiplied by `0.012`, then by fitted `gb=0.19965256`. | The equilibrium relation is not imposed, and the historical baseline violates it by a factor 2.2697 in CO₂. | **calibration choice** and **substantive model difference** |
| Tight junction, Eqs. (36)–(37) | `jt=g(Vt-E)/F` with the published Nernst ratios. | Sign and conversion agree for monovalent cations. Historical conductances are pre-scaled. | **equivalent after algebra** plus **implementation convention** |
| Current versus particle flux | Channels and tight junctions are stored as `I/F`; Cl is negative for outward motion, while the Cl mass balance adds that signed value. NaK current is `F` times its net cycle flux in the voltage solve. | This is internally coherent and closes both QSS current laws. It supplies the missing sign convention, but not a published dimensional scale for all `g` values. | **implementation convention**; absolute scale **unresolved** |
| Calcium input, §2.8/Fig. 2 | Both outer functions use 0.05 microM until `t>100`, then 0.55 microM. | The paper states 58 nM at rest and a qualitative multi-stage input during minutes 6–12. The archived step is not that curve. | **substantive model difference** |
| Solver, §2.9 | `ode15s`, `RelTol=AbsTol=1e-6`. | Solver family agrees and the historical tolerance fills one omission. The code solves an ODE after analytic/algebraic eliminations, not the reported ten-state DAE. | Solver: **exact match**; formulation: **substantive model difference** |
| Time scaling | `Salivary.m` multiplies its complete RHS by `10^4`; `Salivary2.m` uses `10^3`. | Neither multiplier is stated in the paper and the variants disagree with each other. Transient seconds/minutes cannot be mapped uniquely. | **unresolved** |
| Knockout definition | `Salivary_ex.m:7` sets `g2=0`; function arguments permit `g4=0`. | This agrees with the paper's definition of setting the respective activity/density to zero. Only Ae2-KO is explicitly run in the archived script. | Definition: **exact match**; archived coverage: **unresolved** |

## Parameter-by-parameter concordance

### Resting state and fixed compositions

| Quantity | Published value | Historical value/use | Status and evidence |
| --- | ---: | ---: | --- |
| Caᵢ | 0.058 microM | 0.05 microM | **substantive model difference**; stimulation is 0.55 microM step. |
| Clᵢ, Kᵢ, Naᵢ | 50, 120, 25 mM | 50, 120, 25 mM | **exact match**. |
| HCO₃ᵢ | 12.1 mM | 10 mM | **substantive model difference**. |
| pHᵢ | 6.91 | 6.91 (`H=1.23026877e-4` mM) | **exact match**. |
| CO₂ᵢ | 6.6 mM | 6.60012370 mM | Numerical match is a **calibration choice**, not chemical equilibrium. |
| `x_i/omega_i` | 82.8 mM | `Xi/hi=85.00012303` mM | **calibration choice**; the archived value makes its HCO₃=10 state exactly electroneutral. |
| Vₐ, Vᵇ | -50.24, -62.8 mV | -50.24, -62.8 mV to machine precision | **exact match** obtained by calibration. |
| Naₗ, Kₗ, Clₗ | 118.7, 5.6, 124.3 mM | 118.7, 5.6, `Na_l+K_l=124.3` mM | **exact match/equivalent after algebra**. |
| HCO₃ₗ, pHₗ | `1.5e-4` mM, 6.81 | HCO₃ₗ is initialized as Hₗ=`1.5488e-4` mM; pH 6.81 | **exact match** to reported precision. |
| Psiₗ | 48.8 mM | `xl=48.80013572` | **calibration choice**, numerically exact to rounding. |
| Clₑ, Kₑ, Naₑ | 102.6, 5.3, 140.2 mM | same | **exact match**. |
| HCO₃ₑ | 42.9 mM | 40 mM in the water sum; 21 mM in Ae2/Ae4 | **substantive model difference**. |
| pHₑ | 7.4 | 7.41 | Small **calibration choice**. |
| CO₂ₑ | 1.9 mM | 1.6 mM; CO₂ₗ=11.6 mM, sum 13.2 | **substantive model difference**. |
| Cell volume | 1.3 pL | `delta*Hi0=45.29616725*28.7=1300 um^3=1.3 pL` | **exact match** with a recovered geometry convention. |
| Lumen/cell volume ratio | 0.02 | No explicit lumen volume in outer equations | **unresolved**. |

The printed Table 1 row is itself rounded: its cell electroneutrality
residual is `+0.100123027` mM. The historical row eliminates that residual
by changing HCO₃ᵢ and impermeant charge; this is evidence of a different
calibration row, not an exact reproduction of Table 1.

### Geometry, conductances, water, and transporters

The historical geometry is `Aa=36.09256354 um^2`,
`Ab=54.49977095 um^2`, mean area `delta=45.29616725 um^2`, and height
`Hi0=28.7 um`. Several fitted values transform back to the published table
with less than 0.3% relative error. These numerical identities are strong
lineage evidence, but the transformations are not documented in the paper.

| Family / parameter | Published | Historical (`Par.mat` or `Parameters.m`) | Concordance |
| --- | ---: | ---: | --- |
| R, T, F | 8.3144621, 310, 96485.3365 | exact; `RTF=26.7137302361` mV | **exact match**. |
| Ion valences | `z_Na=z_K=+1`, `z_Cl=-1` | Not stored; the signs are hard-coded consistently in Nernst and balance expressions | **implementation convention** equivalent to the published values. |
| Membrane capacitance | Appears as `C_m`, no table value | Absent because both historical solvers impose potential QSS | **equivalent after algebra** for QSS; full dynamic value **unresolved**. |
| Gₙₐᵗ, Gₖᵗ | 12.46, 0.9 nS | 357.8633121, 25.89317836; division by height gives 12.469105 and 0.902201 | Hidden geometry **implementation convention**. |
| G_Ae4 | 0.66 amol/um³ | `g4=18.95026072`; division by height gives 0.66028783 | Hidden geometry **implementation convention**; kinetics still differ. |
| G_Ae2 | 0.01807 fmol/s | `g2=0.3990848725`; multiplying by mean area and converting amol to fmol gives 0.018077015 | **equivalent after algebra** under recovered area convention. |
| G_Nhe1 | 0.0305 fmol/s | `g1=0.6747373081`; area conversion gives 0.030563014 | **equivalent after algebra** for activity; kinetic exponent differs. |
| G_CaCC | 71.3 nS | `gcl=2.046811983e7`; `gcl/(height*1e4)=71.317491` | Family-specific hidden scaling; **calibration choice**. |
| G_CaKC | 30.4 nS | `gk=8.749517823e7`; `gk/(height*1e5)=30.486125` | Different hidden power-of-ten scaling; **calibration choice**. |
| K_CaCC, eta1 | 0.26 microM, 1.46 | 100.26, 1.49 | **substantive model difference**. |
| K_CaKC, eta2 | 0.182 microM, 2.54 | 100.3, 1.7 | **substantive model difference**. |
| Pₐ, Pᵇ, Pₜ | `4.32e-12`, `5.15e-11`, `2.6e-13` L² mol⁻¹ s⁻¹ | `b1=5.79346122e-4`, `b2=5.45083717e-5`, `b3=7.31969563e-6` in undocumented units | No common conversion. Historical `b1/b2=10.6286`, while published `P_a/P_b=0.08388`. **calibration choice** / **unresolved** units. |
| alpha_Nkcc1 | 2.15 amol/um³ | `aNkcc1=0.0063812`; `aNkcc1*1e4/height=2.22341` | Near but not exact hidden scale; **calibration choice**. |
| Nkcc `a1,a3` | 157.5, 1.0306 | 157.55, 1.0306 | **exact match** to rounding. |
| Nkcc `a2,a4` | `2.0096e7`, `1.3852e6` labelled mM⁻⁴ | `2.0096e-5`, `1.3852e-6` with mM states | **publication typo likely**, corroborating the Palk corrigendum. |
| alpha_NaK | 4.84 amol/um³ | `aNaK=0.00138020584`; `aNaK*1e5/height=4.80908` | Hidden family-specific scale; **unresolved**. |
| NaK `r,alpha1` | `1.305e6`, 0.641 labelled mM-based | same numbers, applied to M concentrations via `s=1e-3` | **implementation convention**; published units remain **unresolved**. |
| G_Ae2, KCl, KB | 0.01807 fmol/s, 5.6 mM, `1e4` mM | area-equivalent G, 5.6 mM, `1e-4` mM | G/KCl agree; `K_B` is **publication typo likely**. |
| G_Nhe1, KH, KNa | 0.0305 fmol/s, `4.5e-4`, 15 mM | area-equivalent G, same constants | Numerical parameters agree; flux exponent is a **substantive model difference**. |
| Ae4 `kplus` | `1.92e-2` mM⁻⁴ s⁻¹ | `1.92e-2`, used with mM states | **exact match** and supports the literal mM concentration basis for the historical code. |
| Ae4 `kminus` | `1.3e-5` mM⁻⁴ s⁻¹ | effective `0.001593058814`; original base `1.341840734e-5` multiplied by 118.7219 | Hidden **calibration choice** and **substantive model difference**. |
| Buffer `k1,kminus1` | 11, `2.6e4` | first multiplied by 0.012 to 0.132 and 312, then by `gb=0.19965256` | Hidden **calibration choice**. |
| P_CO2 | `1.97e-13 s^-1` | 1970 in undocumented units | `10^16` numerical separation; **unresolved**. |

## Reported-observable and phenotype concordance

The paper reports time courses rather than a machine-readable output table, so
the comparison below separates its stated values and qualitative directions
from translated execution of the archived equations. The historical values are
from the seven-state BDF reproduction of the archived `0.05 -> 0.55` calcium
step at `t=100`, without retuning. Its endpoint is `t=200`; it is not assigned
the paper's minute scale because the code multiplies the RHS by `10^4`. The
fixed-lumen five-state variant is included only where it changes the knockout
conclusion.

| Published observable | Historical translated result | Assessment | Status |
| --- | --- | --- | --- |
| WT intracellular response: Cl, K, and HCO3 decrease; Na increases | At `t=200`, WT changes from `(Cl,K,HCO3,Na)=(50,120,10,25)` to `(33.4576,114.6933,6.7300,36.9342)` mM. | All four response directions agree, but HCO3 starts from the different historical baseline. | Qualitative agreement only; quantitative reproduction **unresolved**. |
| WT CO2 and pH change little | CO2 changes from `6.60012370` to `6.60012747` mM, but pH changes from `6.9100` to `6.83165`. | CO2 agrees qualitatively; the pH response does not. The historical pH is constrained by electroneutrality and a different acid/base subsystem. | pH: **substantive model difference**. |
| WT luminal Na, K, and Cl increase | `(Na_l,K_l,Cl_l)` changes from `(118.7,5.6,124.3)` to `(124.5930,6.9581,131.5511)` mM. | The three directions agree. Exact published curves are unavailable for a quantitative check. | Qualitative agreement only; magnitude **unresolved**. |
| Cell volume falls by `27.3%` at maximal stimulation | Seven-state WT endpoint volume is `0.9915659` pL from `1.3` pL, a `23.7257%` fall. The five-state endpoint is `1.0579597` pL, an `18.6185%` fall. | The seven-state number is similar but not equal and comes from a different input and unknown time mapping. | **unresolved**, not a numerical reproduction. |
| Apical voltage depolarizes from `-50.24` to about `-41` mV and later hyperpolarizes to `-55` mV; basolateral voltage hyperpolarizes by `14.7` mV | Persistent-step seven-state endpoint is `V_a=-38.5456` and `V_b=-81.0811` mV, changes of `+11.6944` and `-18.2811` mV. There is no archived agonist-removal phase. | Initial directions are related, but the archived protocol cannot produce the published biphasic/applied-then-removed trace. | **substantive model difference** in protocol and channel gates. |
| Total flow rises from `0.24e-08` to `3.2e-08` microlitre/min, about `13`-fold, then returns after agonist removal | Historical rest is `0.0005314155` in undocumented units; endpoint and largest one-unit post-step sample are `2.55936`- and `2.67155`-fold rest. Calcium is never returned to baseline. | Absolute units cannot be compared, and even the dimensionless fold response is far smaller. | Units **unresolved**; response/protocol a **substantive model difference**. |
| AE2 knockout has no significant effect on state or flow | Seven-state AE2-KO flow differs from WT by `-0.09477%` at the endpoint and `-0.26182%` at the first post-step maximum; five-state differences are `+0.000007%` and `+0.20229%`. | This qualitative phenotype survives both historical reductions under their own protocol. | Qualitative agreement; exact reproduction **unresolved**. |
| AE4 knockout lowers Cl, raises Na, K, and HCO3, leaves CO2/pH nearly unchanged | Relative to seven-state WT at `t=200`, AE4-KO has lower Cl (`32.3961` vs `33.4576`), higher K (`118.2900` vs `114.6933`) and HCO3 (`34.4224` vs `6.7300`), but **lower** Na (`33.3220` vs `36.9342`) and pH `7.19845` vs `6.83165`; CO2 remains nearly equal. | The Na and pH results contradict the published phenotype. The historical knockout starts away from its own KO steady state and encodes Na-only Ae4. | **substantive model difference**. |
| AE4 knockout lowers all luminal ions, hyperpolarizes apical voltage, and slightly depolarizes basolateral voltage | At `t=200` all three luminal values are slightly below WT and `V_a` is `0.855` mV more negative, but `V_b` is also `0.780` mV more negative rather than depolarized. | Lumen and apical directions agree; the basolateral direction does not. | Partial directional agreement; overall **substantive model difference**. |
| AE4 knockout lowers stimulated flow by about `24%` after an initially similar peak | Seven-state AE4-KO is `+10.1877%` at the first post-step maximum and `-0.14439%` at `t=200`. Five-state AE4-KO is `-9.6703%` at its first maximum and `+5.9019%` at `t=200`, where `max |dx/dt|=0.0722` after hidden scaling. | Neither archived dynamical variant approaches the central `-24%` steady phenotype; the five-state endpoint is not settled. | **substantive model difference** and failed output reproduction. |
| Knockout flux compensation: reciprocal Ae2/Ae4 increase, little Nkcc change, slight Nhe1 increase in AE2-KO and large Nhe1 decrease in AE4-KO | At seven-state endpoints, AE2-KO gives Ae4 `+47.70%`, Nkcc `+17.55%`, Nhe1 `-0.44%`; AE4-KO gives Ae2 `+4.38%`, Nkcc `+34.35%`, Nhe1 `-65.98%`. | Reciprocal exchanger and AE4-KO Nhe1 directions agree. The AE2-KO Nhe1 sign and claimed minimal Nkcc compensation do not. | Partial match with **substantive model difference**. |

Thus, matching resting values and several WT response directions do not extend
to the key quantitative observables. Most importantly, the archived equations
fail the published Ae4-knockout flow phenotype under the only reconstructable
archived protocol.

## Recomputed historical resting residual

Using the seven-state initial condition in `Salivary_ex.m`, the exact values in
`Par.mat`, and the archived formulas before the final global RHS multiplier:

| Quantity | Recomputed historical value |
| --- | ---: |
| Hᵢ | `0.0001230268770` mM (pH 6.9100000003) |
| CO₂ᵢ | `6.6001236986` mM |
| qₐ, qᵇ | `0.000463341313973`, `0.000463341313973` |
| qₜ, q_total | `0.000068074162788`, `0.000531415476761` |
| Vₐ, Vᵇ | `-50.2399999999999`, `-62.7999999999999` mV |
| CaCC `I/F`, tight Na `I/F`, tight K `I/F` | `-0.0660549437613`, `0.0630790170915`, `0.00297592666986` |
| Apical QSS residual | `-4.03e-17` historical flux units |
| CaKC `I/F`, NaK cycle flux | `0.0453088160576`, `0.0207461277037` |
| Basolateral QSS residual | `+4.03e-17` historical flux units |
| Maximum 12-row residual during original calibration | `1.39e-17` historical flux units |
| Maximum unscaled seven-state residual magnitude | `5.24e-12` in the mixed historical state units |

The original 12-row residual is at floating-point zero because
`Original/Saliva_Ae4.m` solves for the
conductances and activities from this same resting row. It validates internal
algebra and ancestry, not independent reproduction of the published Table 1
row. In particular, replacing the historical Na-only Ae4 allocation with the
published Na/K allocation immediately destroys the Na and K closures, and the
historical HCO₃, impermeant charge, bath bicarbonate, calcium, water
permeabilities, gates, and CO₂ closure are not the published values.

## Concordance conclusion

Historical evidence resolves or strongly suggests four conventions: the
conservation sign `d omega_i/dt=q_b-q_a`, mM-state conversion of the corrected
Nkcc1 coefficients, a consistent signed `I/F` convention for channel/tight
fluxes, and geometry mappings for Ae2/Nhe1/Ae4 and tight-junction fitted values.
It also identifies likely missing-minus-sign evidence for Ae2 `K_B`.

It does **not** supply the published 2018 executable model. The Na-only Ae4
stoichiometry, alternative channel gates, first-order Nhe1 proton law,
algebraic acid/base subsystem, different calcium step, altered bath values,
family-specific scaling factors, and conflicting time multipliers are
substantive. Those differences must not be silently repaired or promoted from
historical implementation evidence to published biology.
