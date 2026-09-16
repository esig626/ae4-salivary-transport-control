# Pre-Task52 algebraic gate

**Status:** FORMAL DEDUCTION + source-scale calculations from already published/frozen results.  
**Scientific-compute budget consumed:** zero.  
**Purpose:** exhaust the algebra before the single remaining Codex run. Task 52 should be numerical only.

This note starts from verified Task 51 closeout `75b50cdc056876ce7c8e83e0a1e2ba502a80c904`. It does not alter Task 51, the preserved Task 50 benchmark, or any scientific source/model file.

## 1. Exact conserved-coordinate equations with VRAC

Use the existing equal-routing notation:

- `N`: inward NKCC1 cycles, each `(+1 Na,+1 K,+2 Cl)`;
- `A`: positive AE4 cycle, source `(-A/2,-A/2,+A,-2A,-2A)` in `(Na,K,Cl,TIC,TA)`;
- `B`: inward NBC cycles, each `(+1 Na,+2 TIC,+2 TA)`;
- `H`: NHE1 inward-Na / outward-H cycles, source `(+1 Na,+1 TA)`;
- `E`: AE2 inward-Cl / outward-HCO3 cycles, source `(+1 Cl,-1 TIC,-1 TA)`;
- `P`: Na/K-pump cycles, source `(-3 Na,+2 K)`;
- `K`: total K-channel efflux;
- `C`: total CO2 influx to the cell;
- `J = J_TMEM16A + J_VRAC`: total apical chloride export.

The intracellular amount balances are therefore

`Na' = N + H + B - A/2 - 3P`

`K' = N - A/2 + 2P - K`

`Cl' = 2N + A + E - J`

`TA' = H + 2B - E - 2A`

`TIC' = C + 2B - E - 2A`.

These are exact for the declared Task40/48/51 architecture; adding VRAC changes only `J` and current closure.

Eliminating `E` from the Cl and TA equations gives

`J = 2N + H + 2B - A - TA' - Cl'`.  (1)

Eliminating `N+B-A/2` with the Na balance gives the equivalent exact transient identity

`J = 6P - H + 2Na' - TA' - Cl'`.  (2)

Also

`TIC' - TA' = C - H`.  (3)

At a true stationary state:

`C = H`,

`J = 2N + H + 2B - A`,  (4)

`J = 6P - H`.  (5)

**Consequence:** neither AE4, NKCC1 nor the partition between TMEM16A and VRAC appears explicitly in the stationary pump/NHE identity. A new apical Cl channel does not by itself create an AE4-specific sustained phenotype.

## 2. Exact knockout-versus-WT condition

For any two stationary states, define `Delta X = X_KO - X_WT`. Since `A_KO=0`, subtracting (4) gives

`Delta J = 2 Delta N + Delta H + 2 Delta B + A_WT`.  (6)

Equivalently from (5):

`Delta J = 6 Delta P - Delta H`.  (7)

Therefore a sustained AE4-KO deficit (`Delta J < 0`) requires

`2 Delta N + Delta H + 2 Delta B < -A_WT`.  (8)

This is an unusually demanding sign condition. The lost WT AE4 flux enters (6) with a **positive** sign in `KO-WT`: if the rest of the network were unchanged, removing AE4 would increase the export that the stationary balances require. Increased NKCC compensation makes the sign problem worse.

If the stimulated KO remains in the forward physiological regime `H_KO >= 0` and `B_KO >= 0`, then from (4)

`J_KO >= 2 N_KO`.  (9)

This gives a simple necessary NKCC suppression scale for a large sustained deficit.

## 3. What the frozen Task 40 fluxes imply

Task40 60–600 s integrated fluxes are retained only as a scale diagnostic, not as Task52 predictions.

WT:

- `N_WT = 78.4360885485 fmol` cycles;
- `A_WT = 41.0555019619 fmol`;
- `B_WT = 39.1102731990 fmol`;
- `H_WT = 4.05656455364 fmol`;
- `P_WT = 32.7214891030 fmol`;
- `J_WT = 207.327262871 fmol`.

AE4-null:

- `N_KO = 96.6045923152 fmol`;
- `A_KO = 0`;
- `B_KO = 5.39392724140 fmol`;
- `H_KO = 1.76217987392 fmol`;
- `P_KO = 33.3297787214 fmol`;
- `J_KO = 198.699313925 fmol`.

The broad-window loading fraction is

`2 N_WT / J_WT = 0.756640371`.

Under the favourable lower-bound assumption `H_KO,B_KO >= 0`, a quasi-sustained 35% deficit (`J_KO/J_WT=0.65`) requires approximately

`N_KO/N_WT <= 0.859060691`,

i.e. **NKCC must be at least 14.1% lower in KO than WT**.

For the Task50 null secretion ratio `0.697388494`, the analogous threshold is

`N_KO/N_WT <= 0.921690833`,

i.e. **at least 7.83% lower**.

The actual Task40 broad-window value is

`N_KO/N_WT = 1.231634495`,

or 23.16% higher. This is the opposite direction.

Using the 600 s WT endpoint as a quasi-stationary scale (`N_WT=0.1483651634`, `J_WT=0.3768556745 fmol/s`) makes the requirement stronger:

- 35% target: `N_KO/N_WT <= 0.82552`, at least **17.45% below WT**;
- Task50 ratio: `N_KO/N_WT <= 0.88570`, at least **11.43% below WT**.

These endpoint/broad-window thresholds are scale diagnostics, not global theorems, because finite storage remains nonzero.

## 4. Exact finite-window identity: small Task40 deficit is a storage effect

Integrating (1) over `[t0,t1]` gives exactly

`Integral J dt = 2 Integral N dt + Integral H dt + 2 Integral B dt - Integral A dt - Delta_state(TA+Cl)`.  (10)

For KO minus WT over the Task40 60–600 s window, the source-flux contribution is

`2 Delta N + Delta H + 2 Delta B + A_WT = +7.665432900 fmol`.

That term alone favours **more** KO export.

Observed integrated export difference is

`Delta J = -8.627948946 fmol`.

Therefore the differential storage term in `(TA+Cl)` is

`+16.293381846 fmol`,

large enough to reverse the source-flux sign temporarily.

The pump form gives the same result:

`6 Delta P - Delta H = +5.944122390 fmol`,

with a corresponding storage contribution of `-14.572071336 fmol` in (2).

Thus the small inherited Task40 deficit is not supported by a sustained source-flux disadvantage. It is primarily a finite-storage transient.

For scale only, keeping the same Task40 source-flux difference while demanding a 35% broad-window deficit would require about `80.23 fmol` of KO-minus-WT differential `(TA+Cl)` storage over the window. That is ~91% of the Task40 WT resting chloride pool (`~88.03 fmol`). This comparison is illustrative, not a hard bound, because TA also contributes and Task52 parameters/states may differ.

## 5. Palk NKCC law creates a strong compensation pressure

The active source-fixed NKCC core is

`N = alpha * M * f(x)`

with

`x = Na_i K_i Cl_i^2`

and

`f(x) = (A1 - A2 x)/(A3 + A4 x)`

where `A1=157.5`, `A2=2.0096e-5`, `A3=1.0306`, `A4=1.3852e-6`.

Its derivative is exactly

`f'(x) = -(A2*A3 + A4*A1)/(A3 + A4*x)^2 < 0`.  (11)

So for a common activity multiplier, inward NKCC is strictly decreasing in the intracellular product `Na_i K_i Cl_i^2`.

At the Task40 WT 600 s operating point, `x_WT=5.673826456e6`.

To suppress KO NKCC enough to satisfy the broad-window 35% lower-bound threshold (`N_KO/N_WT<=0.85906`) requires

`x_KO >= 5.910279e6`,

about **4.17% above WT**.

For the Task50-ratio threshold (`<=0.92169`), it requires

`x_KO >= 5.803089e6`,

about **2.28% above WT**.

But the actual Task40 KO endpoint has `x_KO=5.316152e6`, below WT, and therefore larger NKCC flux.

The experimentally measured resting chloride phenotype pushes even harder in that direction. If Na and K were unchanged, changing Cl from WT `50.1 mM` to AE4-KO `36.5 mM` reduces `x` by

`(36.5/50.1)^2 = 0.5308`.

Around the Task40 Palk operating point that alone would increase inward NKCC by roughly **3.8-fold**. To obtain the 4.17%-above-WT `x` required for the favourable 35% threshold despite that chloride ratio, the product `Na_i*K_i` would have to be about **1.96-fold WT**. This is a scale calculation, not a prediction of the true KO state.

**Interpretation:** the measured low-Cl KO state naturally increases the compensatory drive of the current Palk law unless other intracellular cations change dramatically.

## 6. Raising Na to suppress NKCC fights the desired secretion deficit twice

The Na/K-pump law contains

`P proportional to Na_i^3/(Na_i^3+K_Na^3)`

at fixed extracellular K, hence

`partial P / partial Na_i > 0`.  (12)

For the mechanistic Cha NHE1 law, at fixed pH and extracellular concentrations define `x=Na_i/K_Nai`. Its Na-dependent rate can be written

`H(x) = m * (aB - C d x) / [(a+d)(1+x) + B + Cx]`,

where `m,a,B,C,d` are positive constants independent of `x`. Differentiating gives

`dH/dx = -m [C d W + Z a B] / denominator^2 < 0`  (13)

for positive `W,Z` obtained from the denominator coefficients.

Thus increasing intracellular Na:

- suppresses NKCC through (11), **but**
- increases pump activity through (12), and
- decreases NHE flux through (13).

Because stationary total Cl export is `J=6P-H`, both of the latter effects increase `J`.

So one obvious route for suppressing compensatory NKCC — raising intracellular Na — simultaneously pushes the stationary export identity in the wrong direction. Suppressing NKCC through higher intracellular Cl is contrary to the measured low-Cl KO phenotype; suppressing it through higher K may be possible mathematically but requires the final numerical model to demonstrate it without violating other states/fluxes.

## 7. Osmotic sign of AE4 is opposite to the intuitive loading story

The water model counts tracked intracellular osmotic particles as `Na + K + Cl + TIC + fixed osmoles`.

Summing the exact amount balances gives the tracked mobile-osmole source

`M' = 4N + H + 3B - 2A - P - K - J + C`.  (14)

Per positive cycle/source:

- NKCC: `+4` osmotic particles;
- NHE: `+1`;
- NBC: `+3`;
- AE2: `0`;
- AE4: **`-2`**;
- pump: `-1`;
- K efflux: `-1`;
- apical Cl export: `-1`;
- CO2 influx: `+1`.

The AE4 sign follows directly from the current equal-routing source vector:

`-1/2 Na -1/2 K +1 Cl -2 TIC = -2 particles`.

Therefore, **within this model AE4 is a net osmotic sink**. At an otherwise identical state and protocol, deleting AE4 increases intracellular osmotic accumulation by `2*A_WT`, not decreases it.

If KO NKCC also compensates upward, the KO-minus-WT osmotic loading shift contains

`4 Delta N + 2 A_WT`,

again favouring more KO swelling.

Under a common-state acute comparison this sign is exact because AE4 is electroneutral and the other inherited source/current laws see the same state. With independently solved genotype rests the state differences can alter the result, but those chronic rests are already experimentally wrong and cannot be invoked as validated biology.

This is a major warning for a genotype-independent positive-swelling VRAC gate: the current equations do not naturally produce less gate activation in AE4 KO.

## 8. VRAC is negative feedback on swelling

Task51 uses

`G_V = g_scale * beta * s`,

`s = max(V/V_ref - 1,0)`.

An outward VRAC Cl flux removes one tracked intracellular osmotic particle per chloride, so outward VRAC is negative feedback on the very swelling that activates it.

Task51 numerical trials show the expected sign:

- `g_scale=0`: mean 300–600 s swelling `5.750804%`, endpoint `6.133046%`;
- `g_scale=1e-6 S`: mean swelling `2.436621%`, endpoint `2.421298%`.

Thus increasing a positive outward VRAC conductance cannot be expected to raise a 5.75% swelling response to the source 12.5% target. Swelling is primarily an upstream-loading constraint; electrophysiology is the appropriate source for the channel scale.

## 9. TMEM16A and VRAC are electrically parallel in the Task51 implementation

At a fixed state, both apical Cl pathways use the same chloride reversal potential `E_Cl`:

`I_TMEM16A = g_CaCC (V_a-E_Cl)`

`I_VRAC = g_V (V_a-E_Cl)`.

Therefore, after the common voltage is solved,

`I_VRAC / (I_TMEM16A + I_VRAC) = g_V/(g_CaCC+g_V)`  (15)

exactly whenever the denominator is nonzero.

The new channel supplies additional total apical Cl conductance but no new electrochemical reversal. Electrical closure consequently causes substantial channel substitution rather than additive export.

Task51 illustrates this strongly. Over 60–600 s:

- `g=0`: TMEM16A export `148.4198 fmol`;
- high-g trial: TMEM16A `18.9829 fmol`, VRAC `147.4373 fmol`, total `166.4202 fmol`.

About **87.8%** of added VRAC export replaced lost TMEM16A export, while total export rose only ~12.1%. This high-g trial is not a physiological calibration; it exposes the topology.

## 10. Newly recovered Catalan Supporting Information: electrophysiological scale

Task51 recorded the official SI as inaccessible. A web-accessible mirror of the Catalan Supporting Information has now been retrieved. This **supersedes only the access limitation**, not the historical Task51 record.

The SI states for the nystatin-perforated IPR recordings:

- internal: `58 mM NMDG-Cl`, `60 mM NMDG-glutamate`, `4.5 mM EGTA`, `10 mM HEPES`, pH 7.2;
- external: `90 mM NaCl`, `25 mM NaHCO3`, `4.3 mM KCl`, `1 mM CaCl2`, `1 mM MgCl2`, `5 mM glucose`, `10 mM HEPES`, pH 7.4;
- isotonic solution adjusted to `315 mOsm/kg`;
- IPR recordings use the nystatin-perforated configuration; Fig. S6 presents current density and the IPR-induced I-V relation.

Counting chloride from NaCl/KCl/CaCl2/MgCl2 gives approximately `Cl_o=98.3 mM` and pipette `Cl_i=58 mM` for an order-of-magnitude Nernst calculation. At room temperature:

`E_Cl approximately -13.5 mV`.

At a `-60 mV` clamp, the chloride driving-force magnitude is about `46.5 mV`.

Using the Task49 source extraction at that condition:

- total IPR-induced current magnitude `208.6 +/- 35.2 pA` -> ohmic-equivalent whole-cell conductance about **4.49 nS**;
- DCPIB-blocked magnitude `107.8 +/- 17.2 pA` -> about **2.32 nS**.

Dividing by the reported `12.5%` swelling to map onto the Task49/51 law gives approximate conductance-scale sensitivities

- total-current map: `g_scale approximately 3.6e-8 S`;
- DCPIB-component map: `g_scale approximately 1.9e-8 S`.

These are **not identified apical conductances**. Qualifications:

1. the patch current is whole-cell and the apical fraction is unknown;
2. the current is outwardly rectifying, so a single ohmic conductance is only a local equivalent;
3. pipette chloride is used as the available protocol concentration, not asserted to equal an exact physiological intracellular free chloride at all times;
4. current and swelling are not guaranteed to be simultaneous matched-cell measurements.

Nevertheless, they establish a useful scale. The Task51 numerical upper bracket `g_scale=1e-6 S` was about 30–50 times larger than these whole-cell source mappings.

## 11. Scale comparison with the working phenotype benchmark

The preserved Task50 WT combined-stimulus effective CaCC conductance is about `31.4 nS`.

At `12.5%` swelling, the source-mapped VRAC whole-cell conductance scale is only about `2.3–4.5 nS`, or roughly **7–14% of the WT combined-stimulus CaCC conductance**.

Even an unrealistically favourable comparison in which WT receives all of that extra conductance and KO receives none gives the bare conductance ratios

- using 4.49 nS: `31.4/(31.4+4.49) approximately 0.875`;
- using 2.32 nS: ratio approximately `0.931`.

By comparison, a conductance ratio equal to the Task50 null secretion ratio (`0.6974`) would require about `13.6 nS` of WT-only extra conductance, and a ratio of `0.65` would require about `16.9 nS` — several times the source-mapped whole-cell current scale.

This is **not a secretion-ratio theorem** because voltage, state and water feedback are nonlinear. It is a strong scale warning. More importantly, the actual literature mechanism does not give WT-only VRAC; equation (14) suggests KO can swell at least as strongly under shared-state conditions.

## 12. Algebraic go/no-go criteria for the final numerical shot

Task52 should not search mechanisms. It should numerically evaluate the already declared beta-NKCC/swelling/VRAC architecture against the following algebraic requirements.

### Gate A — independently constrained channel scale

Use only electrophysiological source mappings/sensitivities for `g_scale` (e.g. DCPIB-sensitive and total-current mappings). Do not tune `g_scale` to AE4 secretion.

### Gate B — upstream loading needed for mouse swelling

With `g_scale` fixed independently, one scalar effective beta-NKCC gain may be inferred from the independent mouse WT IPR swelling amplitude. Report whether the required gain is compatible with the 2.5–3-fold rat salivary evidence and the broader sixfold context. Do not use AE4 phenotype here.

A rough local scaling from Task51 (`M_beta=2.5` -> ~5.75% mean swelling with zero VRAC) suggests that reaching 12.5% could require `M_beta` of order 4 or more, and positive VRAC should increase the required loading. This is only a pre-numerical estimate.

### Gate C — NKCC sign/magnitude in KO

For a substantial sustained deficit, the final frozen model must overcome (6). In the forward H/NBC regime, a KO NKCC flux that remains equal to or above WT is already a strong failure signal. On Task40 scales, the required NKCC suppression is roughly 8–17% below WT depending whether the comparison is to the Task50 ratio or 35% phenotype and whether broad-window or endpoint scale is used.

If `N_KO >= N_WT` while H/B remain nonnegative and storage rates have relaxed, the literature mechanism cannot provide the desired sustained deficit in this chassis.

### Gate D — source flux versus storage

Any apparent finite-window deficit must be decomposed exactly using (10). If the deficit is produced mainly by temporary `(TA+Cl)` storage while the source-flux term favours KO, classify it as transient storage compensation rather than a sustained mechanistic explanation.

### Gate E — osmotic genotype direction

Report WT versus KO `4N + H +3B -2A -P-K-J+C` and VRAC gate activation. If the KO swells/activates VRAC at least as much as WT because AE4 is an osmotic sink and NKCC compensates, the proposed pathway has the wrong genotype direction for explaining the secretion loss.

### Gate F — pump/NHE identity

Report `6P-H` alongside total Cl export. A sustained 30–35% deficit requires the KO pump/NHE combination to move correspondingly downward. Increased Na-driven pump demand or decreased NHE makes this harder, not easier.

## 13. Pre-Task52 conclusion

The algebra supports three different levels of claim:

1. **Established:** beta-NKCC removes Task49's exact initiation deadlock and can generate swelling from AE4-null rest.
2. **Plausible physiology but structurally adverse for the phenotype:** swelling-gated VRAC is a credible cAMP secretory pathway, but in the current conservation model AE4 is an osmotic sink, NKCC tends to compensate upward after AE4 loss, and VRAC is electrically parallel with TMEM16A. Those facts do not naturally generate less KO secretion.
3. **Still numerical:** whether the independently constrained channel scale plus a mouse-swelling-constrained beta-NKCC gain can, through nonlinear state/current/water feedback, nevertheless drive `P`, `H`, `B`, NKCC and storage far enough to produce the observed KO deficit.

Therefore Task52 has one meaningful numerical job: test that narrow possibility and decompose it against the exact identities above. If it fails, the preserved Task50 construction becomes the quantitative demonstration of the **additional unidentified beta-conditioned network coupling required by the data**.

## References / provenance

- Catalan MA et al. *PNAS* 2015, DOI `10.1073/pnas.1415739112`, PMID `25646474`, PMCID `PMC4343136`.
- Catalan et al. Supporting Information mirror retrieved 2026-09-16: `https://esdocs.com/doc/477676/download-supporting-information--pdf-` (mirror used because official SI access failed during Task51; protocol text cross-checks the article's nystatin-perforated electrophysiology description).
- Paulais M, Turner RJ. *J Clin Invest* 1992, DOI `10.1172/JCI115695`, PMID `1313447`, PMCID `PMC442971`.
- Chaib N et al. 1998, PMID `9880083`, DOI `10.1016/S0196-9781(98)00134-X`.
- Tanimura A et al. 1995, DOI `10.1074/jbc.270.42.25252`, PMID `7559664`.
- Kurihara K et al. 1999, DOI `10.1152/ajpcell.1999.277.6.C1184`, PMID `10600770`.
- Kurihara K et al. 2002, DOI `10.1152/ajpcell.00352.2001`, PMID `11880270`.
- Pena-Munzenmayer G et al. *JBC* 2015, DOI `10.1074/jbc.M114.612895`, PMID `25745107`, PMCID `PMC4409235`.
- Pena-Munzenmayer G et al. *J Gen Physiol* 2016, DOI `10.1085/jgp.201611571`.
- Pena-Munzenmayer G et al. *AJP Gastrointest Liver Physiol* 2021, DOI `10.1152/ajpgi.00145.2021`, PMID `34585968`, PMCID `PMC8887885`.
- Repository equations/data: `src/modern_full_model/ae4_equal_cation_routing.py`, `nkcc1_palk2010.py`, `nhe1_cha2009.py`, `membranes.py`, `water.py`; Task40 frozen results; Task49/51 VRAC implementation and outputs; preserved Task50 benchmark.
