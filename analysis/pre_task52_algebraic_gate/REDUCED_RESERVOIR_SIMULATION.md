# Reduced chloride-reservoir simulation

**Classification:** REDUCED-MODEL NUMERICAL DIAGNOSTIC. This is not the full repository ODE model and does not consume the final Task52/Codex shot.

## Purpose

Test whether the chloride-reservoir mechanism identified in `CHLORIDE_RESERVOIR_SOLUTION.md` is numerically capable of generating a substantial WT-versus-AE4-KO apical chloride-export difference before spending the final full-model numerical shot.

The reduced simulation deliberately excludes the Palk/Benjamin NKCC law and any AE4-expression-dependent apical-channel multiplier.

## Construction

The reference time-dependent quantities are the frozen **pre-Palk Task37 WT CCh+IPR** outputs sampled every 60 s:

- cell volume `V_ref(t)`;
- intracellular chloride `Cl_ref(t)`;
- apical voltage `V_a,ref(t)`;
- WT AE4 chloride loading `A_ref(t)`;
- apical chloride export `J_ref(t)`.

Piecewise monotone cubic interpolation is used only to define smooth forcing between saved points.

The non-AE4 chloride supply is reconstructed from the frozen WT chloride balance:

`S(t) = d[Cl amount]_ref/dt + J_ref(t) - A_ref(t)`.

The central reduced hypothesis imposes the same `S(t)` on WT and AE4 KO. WT additionally receives `A_ref(t)`; KO receives zero AE4 loading.

Initial intracellular chloride is imposed from the experiment:

- WT `50.10 mM`;
- AE4 KO `36.50 mM`.

The Task37 reference volume at onset, `1.448426083 pL`, converts this to an initial WT-minus-KO chloride-amount advantage of `19.6986 fmol`.

The common apical chloride current response is linearised from the exact two-membrane voltage-closure result derived in R51H:

`delta(E_Cl - V_a) = alpha * delta E_Cl`, with `alpha > 0`.

At each time the reduced export law is

`J(C,t) = J_ref(t) + [J_ref(t)/D_ref(t)] * alpha * V_T * ln(C/Cl_ref(t))`,

where `D_ref = E_Cl,ref - V_a,ref`. This keeps the same apical-channel law in WT and KO and varies only chloride availability. `alpha` is **not fitted**. It is profiled over the predeclared algebraic scale `0.08–0.20` around the rough circuit estimate `~0.1–0.2`.

Cell volume is held to the frozen WT reference trajectory in this reduced diagnostic. Full voltage, water, Na/K, TIC/TA and lumen feedback are deferred to Task52.

## Main result

| alpha | WT Cl export 0–600 s (fmol) | KO Cl export 0–600 s (fmol) | KO export deficit | final WT-KO reservoir (fmol) | WT Cl at 600 s (mM) | KO Cl at 600 s (mM) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.08 | 213.6253 | 168.4155 | 21.16% | 16.2056 | 52.911 | 41.322 |
| 0.10 | 212.9492 | 164.9787 | 22.53% | 13.4448 | 53.394 | 43.779 |
| 0.12 | 212.5538 | 162.5725 | **23.51%** | 11.4340 | 53.677 | 45.500 |
| 0.15 | 212.2449 | 160.1240 | 24.56% | 9.2945 | 53.898 | 47.251 |
| 0.20 | 212.0586 | 157.6933 | **25.64%** | 7.0500 | 54.031 | 48.990 |

Thus the reduced reservoir mechanism produces a robust **~21–26% cumulative apical chloride-export deficit** across the unfitted voltage-closure sensitivity interval.

The accepted Task41/50 construction required a `26.88%` reduction in null apical CaCC chloride export to generate its `30.26%` cumulative fluid deficit. The reduced reservoir result is therefore of the same order before any genotype-specific channel multiplier is introduced.

## Central alpha=0.12 budget audit

For `alpha=0.12`:

- measured initial reservoir advantage: `19.6986 fmol`;
- integrated WT AE4 loading: `41.7168 fmol`;
- total available reservoir + AE4 advantage: `61.4154 fmol`;
- final residual reservoir advantage: `11.4340 fmol`;
- exact cumulative WT-minus-KO export difference: `49.9813 fmol`.

The balance closes numerically:

`19.6986 + 41.7168 - 11.4340 = 49.9813 fmol`.

So about **81% of the available chloride advantage is discharged into extra WT apical chloride export** in the central reduced simulation.

The resulting cumulative export deficit is `23.51%`.

## Measurement sensitivity

Using one reported SEM in the conservative direction for both resting chloride measurements (`WT=48.6 mM`, `KO=38.1 mM`) versus the opposite direction (`WT=51.6`, `KO=34.9`):

- at `alpha=0.12`, export deficit spans approximately `21.67–25.32%`;
- at `alpha=0.20`, approximately `23.77–27.47%`.

This is a measurement-scale sensitivity only; the SEMs are not treated as a joint confidence interval.

## NKCC-compensation sensitivity

Because the central model uses genotype-shared non-AE4 supply, a diagnostic KO supply increase was applied without fitting.

At `alpha=0.12`:

| KO non-AE4 supply increase | export deficit |
| ---: | ---: |
| 0% | 23.51% |
| +2% | 22.27% |
| +5% | 20.40% |
| +10% | 17.32% |
| +20% | 11.24% |

This demonstrates why the later automatic +20–30% NKCC compensation erased the AE4 phenotype. It also identifies the decisive full-model uncertainty: whether a source-compatible non-Palk NKCC representation keeps genotype-dependent supply changes modest enough for the measured chloride reservoir to matter.

The primary isolated NKCC experiment did not detect a genotype difference, but that is not converted here into a claim of exact physiological flux equality.

## Interpretation

This reduced simulation **does not prove the full model will reproduce the saliva phenotype**. It omits the coupled Na/K, acid-base, lumen, water and exact voltage dynamics.

It does establish a stronger go decision for Task52:

1. the measured resting chloride deficit plus WT AE4 loading has enough chloride mass to matter;
2. after electrochemical attenuation, a large fraction of that mass is discharged over 600 s under shared secretory demand;
3. the resulting `21–26%` chloride-export deficit is close to the `26.88%` artificial export reduction in Task50;
4. the mechanism fails progressively as large KO-specific non-AE4 chloride compensation is reintroduced.

Therefore Task52 should implement exactly this chloride-reservoir architecture in the full conserved model, without Palk/Benjamin and without the Task50 multiplier, and determine the actual fluid-secretion result.
