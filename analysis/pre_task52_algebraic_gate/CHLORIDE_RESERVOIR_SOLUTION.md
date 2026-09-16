# Pre-Task52 chloride-reservoir solution

**Status:** FORMAL DEDUCTION + source-scale calculation.  
**Scientific compute:** zero.  
**Task52 implication:** Palk/Benjamin NKCC is not to be used in the final candidate model.

This note supersedes the Palk-dependent Task52 recommendation in `PRE_TASK52_ALGEBRAIC_GATE.md`. The earlier analysis remains historically valid as a diagnosis of why Tasks 39–51 compensated AE4 loss, but it is not the architecture proposed for Task52.

## 1. Reverse the causal interpretation

The final candidate is not

`AE4 -> swelling -> VRAC`

and not

`AE4 expression -> apical-channel multiplier`.

Instead:

`AE4 maintains the intracellular chloride reservoir and supplies beta-stimulated chloride loading`;

while

`beta/IPR provides a shared additional apical chloride demand through the independently observed TMEM16A-independent anion conductance`.

The shared apical pathway need not be weaker in AE4 KO. Its role is to increase chloride demand. Genotype specificity comes from chloride availability: WT begins with more intracellular chloride and retains beta/PKA-activated AE4 loading; AE4 KO begins chloride-depleted and has no AE4 loading.

This combines four independent source facts:

1. AE4 KO acinar cells have resting `[Cl]_i = 36.50 +/- 1.60 mM` versus WT `50.10 +/- 1.50 mM` (Peña-Münzenmayer et al., JBC 2015).
2. The same study found no detectable genotype-specific loss or increase of isolated NKCC functional activity and attributed the secretion phenotype to AE4-dependent chloride uptake rather than altered NKCC/NHE activity.
3. AE4 is an electroneutral monovalent-cation-dependent Cl/HCO3 exchanger whose proposed physiological role is chloride accumulation in secretory acinar cells; the 2016 thermodynamic analysis explicitly supports AE4-mediated chloride accumulation above equilibrium in exchange for KHCO3 efflux.
4. Beta/cAMP/PKA activates AE4 (Peña-Münzenmayer et al. 2021), while Catalán et al. 2015 independently establish a beta/IPR TMEM16A-independent, DCPIB/NPPB-sensitive apical anion-secretory pathway associated with swelling.

No direct AE4-to-channel molecular interaction is required.

## 2. Palk/Benjamin NKCC is excluded from Task52

The Task39–40 Palk reduction creates strong automatic state-driven compensation after AE4 loss: approximately +23–33% NKCC loading in the later reconstructions. That compensation is the main reason AE4 deletion is almost erased numerically.

This is not a required conservation result. It is a kinetic-law result. The primary 2015 experiment did not detect a genotype difference in the isolated NKCC assay and explicitly argued against increased NKCC functional activity as compensation for AE4 loss.

Therefore Task52 must not use `src/modern_full_model/nkcc1_palk2010.py` as the NKCC mechanism. Historical Palk results remain valid negative diagnostics and are not deleted.

Task52 may represent NKCC only as a genotype-shared experimentally constrained supply term or with a non-Palk saturating/reversible sensitivity that does not introduce a phenotype-informed genotype multiplier. Exact equality of physiological WT/KO NKCC flux is an experimental idealisation, not a source fact; any central shared-flux representation must be labelled accordingly.

## 3. Exact chloride-reservoir balance

Let

- `x_g(t)` be intracellular chloride amount for genotype `g`;
- `N_g` be NKCC chloride loading divided by two in cycle units;
- `A_W(t) >= 0` be WT AE4 chloride loading and `A_K=0`;
- `E_g` be signed AE2 chloride loading;
- `J_g` be total apical chloride export, including TMEM16A and any beta/IPR auxiliary conductance.

Then exactly

`x'_g = 2 N_g + A_g + E_g - J_g`.  (R1)

For WT minus AE4 KO define `delta = x_W - x_K`. Then

`delta' = 2(N_W-N_K) + A_W + (E_W-E_K) - (J_W-J_K)`.  (R2)

If NKCC and AE2 are represented centrally as genotype-shared functional supply,

`N_W=N_K`, `E_W=E_K`,

so

`delta' = A_W - DeltaJ`,  where `DeltaJ = J_W-J_K`.  (R3)

Integrating over `[0,T]` gives the exact finite-window budget

`Integral_0^T DeltaJ dt = delta(0) + Integral_0^T A_W dt - delta(T)`.  (R4)

This is the key result. Cumulative WT export advantage can come from two independently measured/supported chloride sources:

1. the larger initial WT chloride reservoir `delta(0)`;
2. continuing WT AE4 loading during stimulation.

No AE4-dependent channel multiplier appears.

The initial reservoir is not automatically all secreted: `delta(T)` must be retained explicitly. Task52 numerics therefore have a single decisive question: how much of the available reservoir-plus-AE4 chloride advantage is discharged over 600 s under shared apical demand?

## 4. Source-scale chloride budget

Using the Task37 pre-Palk WT cell-volume scale `V0 = 1.448426083 pL`, the measured resting chloride difference gives

`deltaCl0 = (50.10 - 36.50) mM * V0 = 19.6986 fmol`.

Using only the reported Cl SEMs, the propagated concentration-measurement contribution is about `3.18 fmol` (cell-volume uncertainty not included).

The pre-Palk Task37 WT trajectory gives, on the common 60 s output grid over 0–600 s,

- integrated WT AE4 chloride loading approximately `41.25 fmol`;
- integrated WT apical CaCC chloride export approximately `225.70 fmol`.

Thus the maximum available chloride advantage before final storage is

`deltaCl0 + Integral A_W approximately 60.95 fmol`,

which is about `27.0%` of that WT apical chloride export scale.

The accepted Task41/Task50 construction produced the desired fluid phenotype with a `26.88%` reduction in null apical CaCC chloride export (`207.327 -> 151.594 fmol` over 60–600 s).

The near equality of these scales is not an identity and must not be called validation. It is nevertheless the strongest algebraic sufficiency clue in the project: **the experimentally measured missing chloride reservoir plus ordinary WT AE4 loading contains approximately the same amount of chloride as the artificial apical-export reduction Task50 had to impose.**

Task38/40 largely erased this available difference by allowing large state-driven NKCC compensation. For example, the pre-Palk Task38 null gained about `41.9 fmol` extra NKCC chloride loading over WT on 60–600 s, while Task40 gained about `36.3 fmol`; these terms consume most of the reservoir-plus-AE4 budget.

## 5. Chloride electrochemical driving-force effect

For any apical chloride conductance with lumen chloride fixed,

`E_Cl = V_T ln([Cl]_i/[Cl]_l)`.

Therefore the measured WT-to-KO resting chloride difference changes chloride reversal by

`Delta E_Cl = V_T ln(36.5/50.1) = -8.4646 mV` at 37 C.  (R5)

This shift is independent of the assumed lumen chloride because it is a ratio of intracellular concentrations.

In the frozen Task40 resting electrical context (`[Cl]_l approximately 145.924 mM`, `V_a approximately -31.429 mV`), zero apical chloride driving occurs at

`[Cl]_i* = [Cl]_l exp(V_a/V_T) approximately 45.0 mM`.

Thus, at that same model electrical/luminal state:

- measured WT `50.1 mM` lies on the outward-secretory side;
- measured AE4-KO `36.5 mM` lies on the opposite side of chloride reversal.

This is a model-context diagnostic, not a claim that the model lumen is measured physiology. The robust source-level statement is the `-8.46 mV` shift.

## 6. Voltage closure cannot remove the chloride-reservoir effect

Let

- `G_c` be total apical chloride conductance (TMEM16A plus auxiliary beta conductance);
- `G_ka` be apical K conductance;
- `G_p` be total paracellular electrical conductance;
- `D = G_kb + G_p + dI_NBC/dV_b > 0` be the positive linearised basolateral differential conductance.

Linearising the exact two-membrane current closure with respect to `E_Cl` gives

`dV_a/dE_Cl = G_c D / [(G_ka+G_c+G_p)D - G_p^2]`.  (R6)

Hence the outward chloride driving force `D_Cl = E_Cl - V_a` obeys

`dD_Cl/dE_Cl = alpha`,

with

`alpha = [(G_ka+G_p)D - G_p^2] / [(G_ka+G_c+G_p)D - G_p^2] > 0`.  (R7)

So membrane-voltage relaxation can attenuate the effect of lower intracellular chloride, but it cannot reverse it. Lower `[Cl]_i` necessarily lowers outward apical chloride driving in the stable positive-conductance regime.

On the frozen pre-Palk stimulated scale, the WT apical chloride driving force is only a few millivolts. A source-imposed `8.46 mV` shift in `E_Cl` therefore has ample leverage; a rough linearised circuit evaluation gives an `alpha` of order `0.1–0.2`, leaving an effective driving-force loss of order `1 mV`, comparable with a large fraction of the WT stimulated chloride driving force. This scale estimate must be replaced by exact Task52 numerical closure, not treated as a fitted parameter.

## 7. Comparison theorem: the sign is correct without genotype-specific channels

Suppose the two genotypes have common NKCC/AE2 supply and the same apical-channel law, and write total apical export as `J(x,t)` where the exact electrical closure satisfies

`partial J / partial x > 0`

through (R7).

Then from (R3):

`delta' = A_W(t) - [J(x_W,t)-J(x_K,t)]`,

with measured `delta(0)>0` and `A_W(t)>=0`.

At any hypothetical first crossing `delta=0`, the export difference is zero and

`delta' = A_W >= 0`.

Therefore `delta(t)` cannot cross from positive to negative. Consequently

`x_W(t) >= x_K(t)` and `J_W(t) >= J_K(t)`

for the reduced common-supply system.  (R8)

This gives the desired phenotype direction **without any AE4-dependent channel conductance**.

## 8. Why a shared beta/IPR auxiliary conductance helps rather than hurts

Locally write

`J_W-J_K approximately kappa(t) delta`, with `kappa(t)>0` proportional to the effective apical chloride demand after voltage closure.

Then

`delta' = A_W - kappa delta`.  (R9)

Increasing a shared beta/IPR apical conductance increases `kappa`. It therefore discharges the WT-versus-KO chloride-reservoir difference faster and increases cumulative export separation through (R4), even though the conductance itself is identical in both genotypes.

This reverses the Task49/51 interpretation:

- VRAC does **not** need to encode genotype specificity;
- VRAC/shared beta secretory demand exposes genotype specificity already present in chloride supply.

In this interpretation, Catalán's pathway is a demand amplifier, while AE4 is the supply/reservoir pathway.

For approximately constant `kappa` and AE4 loading `A`, (R9) has

`delta(t) = A/kappa + [delta(0)-A/kappa] exp(-kappa t)`.

The cumulative export advantage remains exactly

`delta(0)+A*T-delta(T)`.

Thus the final numerical question is not whether a mysterious coupling exists. It is whether the shared apical conductance is sufficient to discharge most of the approximately `61 fmol` available chloride advantage within 600 s.

## 9. Task52 numerical architecture implied by the algebra

The last Codex shot should test this one candidate, not search.

### Frozen source constraints

1. **No Palk/Benjamin NKCC.** Do not import or call `nkcc1_palk2010.py` in the candidate model.
2. Initialise WT and AE4 KO with the measured genotype resting chloride concentrations `50.10` and `36.50 mM`, with pH `6.91` and `6.89` respectively. Unmeasured coordinates must be handled transparently; do not refit them to secretion.
3. Retain the independently established beta/PKA activation of AE4 in WT; AE4 is exactly zero in KO.
4. Use the same apical secretory conductances in WT and KO. No AE4-expression multiplier on TMEM16A, VRAC or water.
5. Treat Catalán beta/IPR swelling/anion conductance as a **shared demand pathway**. Prefer the independently source-mapped electrophysiological scale; do not tune it to the AE4 phenotype.
6. Constrain NKCC to have no genotype-specific capacity/regulatory multiplier, consistent with the 2015 isolated NKCC assay. A central genotype-shared effective supply representation is permitted as an explicit experimental idealisation. If a state-dependent non-Palk sensitivity is included, it must be predeclared and cannot be selected from the AE4 output.
7. Freeze all choices before evaluating the AE4 secretion phenotype or Task50 benchmark.

### Required outputs

Task52 must report

- cumulative WT-minus-KO apical chloride export;
- decomposition `delta(0) + Integral A_W - delta(T)` and any deviations caused by non-shared NKCC/AE2 terms;
- WT/KO `E_Cl`, `V_a`, and chloride driving force through time;
- remaining chloride-reservoir difference `delta(T)`;
- NKCC WT/KO difference explicitly;
- cumulative water/secretion ratio;
- comparison with Task50's `26.88%` apical-export reduction and `30.26%` fluid deficit;
- CCh-only and beta-containing specificity if available without adding a new calibration.

### Go/no-go

The reservoir mechanism is numerically successful if independently constrained shared apical demand plus measured chloride states and WT AE4 loading yield a substantial persistent AE4-null secretion deficit without a genotype-specific channel multiplier or post-reveal tuning.

If the model only fails because a chosen NKCC kinetic law generates large automatic KO compensation inconsistent with the experimental NKCC constraint, that is a failure of that NKCC representation, not a proof against the chloride-reservoir mechanism.

## 10. References

- Peña-Münzenmayer G et al. *J Biol Chem* 2015. DOI `10.1074/jbc.M114.612895`; PMID `25745107`; PMCID `PMC4409235`. Resting Cl/pH phenotype, secretion phenotype, isolated NKCC/NHE assays and beta-associated uptake contrast.
- Peña-Münzenmayer G et al. *J Gen Physiol* 2016. DOI `10.1085/jgp.201611571`; PMID `27114614`; PMCID `PMC4845690`. AE4 electroneutral cation-dependent Cl/HCO3 mechanism and physiological interpretation as a chloride-accumulating pathway.
- Peña-Münzenmayer G et al. *AJP Gastrointest Liver Physiol* 2021. DOI `10.1152/ajpgi.00145.2021`; PMID `34585968`; PMCID `PMC8887885`. Beta/cAMP/PKA activation of AE4.
- Catalán MA et al. *PNAS* 2015. DOI `10.1073/pnas.1415739112`; PMID `25646474`; PMCID `PMC4343136`. TMEM16A-independent beta/IPR secretion, swelling and DCPIB/NPPB-sensitive anion conductance.
- Repository: Task37/38 pre-Palk outputs; Task40/41 algebraic balances; verified Task51 closeout `75b50cdc056876ce7c8e83e0a1e2ba502a80c904`; preserved Task50 benchmark `b16c30094b95f61a73d8f1977cd58e79c7bb50f6`.
