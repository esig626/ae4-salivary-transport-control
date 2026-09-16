# Mandatory research ledger continuation — R51I

**Binding status:** append-only continuation of `docs/MANDATORY_RESEARCH_LEDGER.md`, `docs/MANDATORY_RESEARCH_LEDGER_PRE_TASK52_APPENDIX.md` and `docs/MANDATORY_RESEARCH_LEDGER_R51H.md`.

**Budget:** Task51 is complete. Exactly **one** funded/scientific Codex shot remains: Task52. The numerical diagnostic recorded here was executed outside Codex and does not consume that shot.

Full diagnostic: `analysis/pre_task52_algebraic_gate/REDUCED_RESERVOIR_SIMULATION.md`.

## R51I — reduced chloride-reservoir simulation

**Classification:** REDUCED-MODEL NUMERICAL DIAGNOSTIC. Not a full production-model prediction and not independent validation.

A reduced two-genotype chloride-amount model was simulated using frozen pre-Palk Task37 WT forcing, measured genotype resting chloride, genotype-shared non-AE4 chloride supply, WT-only AE4 loading and an identical apical-channel law in WT and KO.

No Palk/Benjamin NKCC law, Task50 multiplier, genotype-specific apical conductance, phenotype fit, optimiser, stationary solve or parameter search was used.

### Central construction

The frozen Task37 WT chloride balance is used to reconstruct shared non-AE4 supply:

`S(t)=d x_ref/dt + J_ref(t) - A_ref(t)`.

WT receives `S+A`; AE4 KO receives `S`. Initial chloride is imposed from the 2015 source at `50.10 mM` WT and `36.50 mM` KO.

A common chloride-export response is linearised from the exact current-closure theorem in R51H:

`delta(E_Cl-V_a)=alpha*delta E_Cl`, `alpha>0`.

`alpha` is not fitted and is profiled from `0.08` to `0.20` around the pre-Task52 circuit scale estimate.

### Result

Across `alpha=0.08–0.20`, cumulative AE4-null apical chloride-export deficit is approximately **21.16–25.64%**.

Central `alpha=0.12`:

- WT export `212.5538 fmol`;
- KO export `162.5725 fmol`;
- export deficit **23.5147%**;
- initial WT chloride advantage `19.6986 fmol`;
- integrated WT AE4 loading `41.7168 fmol`;
- final WT-KO reservoir difference `11.4340 fmol`;
- cumulative WT-minus-KO export `49.9813 fmol`.

The exact reduced balance closes:

`19.6986 + 41.7168 - 11.4340 = 49.9813 fmol`.

Approximately **81% of the available reservoir-plus-AE4 chloride advantage is discharged as extra WT apical chloride export** in the central reduced diagnostic.

Task50's successful target-calibrated construction reduced null CaCC export by `26.88%` and fluid secretion by `30.26%`. The reduced reservoir mechanism therefore produces a chloride-export effect of the same order without an AE4-dependent channel multiplier.

### Sensitivities

Using one reported SEM in opposite directions for the resting chloride measurements gives:

- at `alpha=0.12`, approximately `21.67–25.32%` export deficit;
- at `alpha=0.20`, approximately `23.77–27.47%`.

A diagnostic KO-specific increase in the reconstructed non-AE4 chloride supply progressively erases the effect. At `alpha=0.12`:

- 0% extra KO supply -> `23.51%` deficit;
- +2% -> `22.27%`;
- +5% -> `20.40%`;
- +10% -> `17.32%`;
- +20% -> `11.24%`.

This is consistent with the historical observation that the later +20–30% automatic NKCC compensation masked AE4 loss. It does **not** convert the isolated NKCC nonsignificant genotype comparison into exact physiological flux equality.

### Decision

R51I strengthens the R51H go decision. Task52 is now justified as one full-model numerical test of the chloride-reservoir/supply mechanism.

Task52 must:

1. exclude Palk/Benjamin from the candidate model;
2. use measured genotype resting chloride;
3. retain WT beta-activated AE4 and zero AE4 in KO;
4. use identical apical conductance laws across genotypes;
5. constrain non-AE4 chloride supply without phenotype-selected KO compensation;
6. preserve the exact chloride budget and report `delta(0)`, integrated AE4 loading, non-AE4 supply difference, `delta(T)` and cumulative apical export difference;
7. compute full water/secretion dynamics and compare to the preserved Task50 benchmark only after the model is frozen;
8. perform no mechanism search or post-reveal retuning.

If the full conserved model does not retain a substantial secretion deficit under these constraints, Task50 remains the proof that an additional unidentified beta-conditioned network coupling is required.
