# Task 39 final report

**PALK_NKCC1_FULL_VALIDATION_COMPLETE**

The source-fixed Palk/Benjamin NKCC1 law passes source reproduction, exact
Task 31 REST nesting, and all three 600-second trajectory gates. It does
**not** suppress NKCC1 compensation in the frozen Task 38 architecture.
Compensation increases to **29.6543% at AE4=0.05** and **33.1013% at null**.
The final cumulative secretion deficits are only **0.482162%** and
**0.627320%**, respectively. Endpoint flow exceeds WT in both perturbations.
The proposed mechanistic/phenotype hypothesis is therefore unsupported by
this single permitted replacement.

Work is confined to `codex/task-39-palk-nkcc1-full-validation`, based on
prepared head `4aa1a1454dfb042d55523d103f10f27d3de77132`. GitHub comparison
confirmed that only Task 39 instructions were added after frozen Task 38
commit `d4bce14c5caed4e05587e237f002814f0f27e226`.

## Source, nesting, and scientific isolation

[Source/equation verification](source_equation_verification.md) records the
primary-paper and historical-code checks, exact fourth-order M-to-mM
conversion, reversal, and inward `1 Na : 1 K : 2 Cl` source signature.

The independent resting calculation gives
`X_rest=4941065.334966328 mM^4`, `F_rest=7.391062769440939`, and
`alpha_eff=0.017334746894096052 fmol/s`. The scale was hashed before dynamics.
It reproduces the accepted `0.1281222023866353 fmol/s` resting NKCC1 cycle flux.
No parameter was inferred from a stimulated or genotype trajectory.

All seven source tests pass. The one permitted literal coding correction
replaced an inappropriate bitwise conversion assertion with exact Decimal
and numerical-precision checks; it changed no coefficient or scientific
code. The initial and corrected logs are retained.

The complete 13-row REST RHS matches bit for bit. Differences in all amount
and volume RHS rows, all three voltages, and water outflow are exactly zero.
Regulatory RHS and conservation gates pass. The accepted Task 31 R09 core
state hash remains `d77c1907e3e7124061721d606148a733f858cce2ef591c4e9f16e5cb7fc2469d`.
Every production initial vector, including the regulatory coordinate, has
hash `879a650761e9044691ab6b2a9b1bfc3be8f951ece176925117596cc2baa31abe`.

The existing N1 multiplier is **1.0 at REST and 1.75 under stimulation**.
It multiplies the complete Palk cycle flux, preserves reversal, and is
verified at every monitored trajectory state. The generic law remains
selectable but is inactive in every Task 39 production layer. The new
selector is separate from inherited parameter records, so whole-cell,
AE4, NBC, regulation, and stimulus hashes are unchanged. The four existing
source edits only evaluate/select/forward the NKCC1 core law; the other 113
prepared files remained byte-identical in the executed snapshot.

Cha NHE1 amount remains `2.339370005697548e-05 fmol`; NBC remains electrogenic
`1 Na : 2 HCO3` at `G_NBC=0.11570913197464398 fmol/s`; the NHE1 multiplier
remains 1.0/2.3. All other mechanisms, geometry, bath, and stimulus are frozen.
No stationary/resting solve was performed. The inherited acid-base and NBC
algebraic closure routines are unchanged.

## Gated trajectory results

WT, 5%, and null were run in that order from the same accepted WT REST
vector, with CCh 0.3 uM, isoproterenol 5 uM, Ca 0.25 uM, and beta occupancy 1.
The inherited onset convention is REST at exactly zero and active thereafter;
continuous Radau integration starts at 1e-6 seconds with the unchanged state.

All use production Radau: `rtol=1e-7`, amount `atol=1e-10 fmol`, volume
`atol=1e-12 pL`, regulatory `atol=1e-10`, maximum step 2 seconds. Each run
completed with 333 accepted steps, 934 checked samples, 2336 RHS evaluations,
and no scientific gate failure. Checks cover accepted endpoints plus dense
output at integer seconds, not an assertion of a continuous-time bound
between samples.

| Metric | WT | AE4=0.05 | AE4=0.00 |
| --- | ---: | ---: | ---: |
| Cumulative secretion, 0–600 s (pL) | 0.958698900 | 0.954076419 | 0.952684789 |
| Cumulative ratio to Task 39 WT | 1 | 0.995178381 | 0.993726799 |
| Cumulative deficit (%) | 0 | **0.482162** | **0.627320** |
| Mean flow, 60–600 s (pL/s) | 0.00157517140 | 0.00156870451 | 0.00156627566 |
| Mean-flow ratio to WT | 1 | 0.995894487 | 0.994352529 |
| Flow at 600 s (pL/s) | 0.00151477902 | 0.00155405248 | 0.00155587387 |
| Endpoint-flow ratio to WT | 1 | 1.025926860 | 1.027129271 |
| Endpoint deficit (%) | 0 | **−2.592686** | **−2.712927** |

Negative endpoint deficits denote increased flow. WT passes all three
secretion activation gates: mean and cumulative flow exceed 1.10 times
their resting comparators, and endpoint flow exceeds REST. The perturbations
also satisfy these thresholds, although their required decision gates are
physiology and conservation, as in Task 38.

Task 39 WT cumulative secretion is **0.971331006 times** the frozen Task 37
value `0.986995055625111 pL`: a 2.86690% decrease caused by this replacement.

All tracked amounts and volumes are positive and finite. Monitored ranges:

| Quantity | WT | AE4=0.05 | AE4=0.00 |
| --- | --- | --- | --- |
| Na_i (mM) | 11.6361–20.6007 | 11.6361–17.7339 | 11.6361–17.5681 |
| K_i (mM) | 107.037–116.763 | 112.562–116.763 | 112.770–116.763 |
| Cl_i (mM) | 51.2065–60.3050 | 51.7763–60.3050 | 51.7344–60.3050 |
| pH_i | 6.91000–7.03559 | 6.91000–7.20212 | 6.91000–7.22239 |
| HCO3_i (mM) | 4.90937–7.75517 | 4.90937–14.4270 | 4.90937–15.5029 |
| TIC_i (mM) | 5.67175–8.76746 | 5.67175–15.5789 | 5.67175–16.6854 |
| Cell volume (pL) | 1.35372–1.46310 | 1.44835–1.49867 | 1.44835–1.52000 |
| Largest conservation residual / tolerance | 0.00612843 | 0.00592237 | 0.00606803 |

Every inherited charge, current, carbon, water, buffer, and speciation
tolerance passes. Across all runs, the largest basolateral current residual
is `2.48801e-23 A` versus `1e-20 A`; the largest water-accounting residual is
`2.16841e-19 pL/s` versus `1e-12 pL/s`. Per-diagnostic maxima are preserved in
the verification JSON files.

## Integrated transport and compensation

All integrals below use 60–600 s. Positive loading integrates `max(J_Cl,0)`;
signed AE2 is also retained. Quadrature uses the inherited one-second grid,
with t=0 and 1e-6 for cumulative secretion. These are flux integrals in the
current one-cell model, not activity assays.

| Quantity (fmol unless stated) | WT | AE4=0.05 | AE4=0.00 |
| --- | ---: | ---: | ---: |
| Positive NKCC1 chloride loading | 145.291036 | 188.376129 | 193.384277 |
| NKCC1 cycles | 72.645518 | 94.188064 | 96.692138 |
| Positive AE4 chloride loading | 40.479639 | 4.298243 | 0 |
| Signed AE2 chloride flux | −0.404498 | 0.435009 | 0.511099 |
| Positive AE2 chloride loading | 0 | 0.435484 | 0.511247 |
| NBC cycles | 38.153270 | 8.618984 | 5.404307 |
| NBC HCO3 influx | 76.306540 | 17.237969 | 10.808614 |
| NBC integrated outward current (C) | 3.68123e-9 | 8.31606e-10 | 5.21436e-10 |
| NHE1 influx | 4.019203 | 1.908558 | 1.768121 |
| Total Na/K pump cycles | 34.492511 | 33.475694 | 33.389813 |
| CaCC chloride export | 200.392343 | 198.816014 | 198.452935 |
| Paracellular chloride return | 74.968719 | 74.023846 | 73.871380 |

| Positive chloride-loading share (%) | WT | AE4=0.05 | AE4=0.00 |
| --- | ---: | ---: | ---: |
| NKCC1 | 78.2099 | 97.5487 | 99.7363 |
| AE4 | 21.7901 | 2.2258 | 0 |
| AE2 | 0 | 0.2255 | 0.2637 |

WT NKCC1 loading is above the requested 65–75% experimental-context band.
This contextual comparison was not a fitting or trajectory-failure gate.
AE2 remains a negligible positive loader: below 0.27% in both perturbations.
Its positive-loading ratio to WT is undefined because WT has zero positive
AE2 loading; the signed integral changes direction and is reported separately.

| Perturbation | Task 38 NKCC1 compensation | Task 39 compensation | Change |
| --- | ---: | ---: | ---: |
| AE4=0.05 | +24.1302% | **+29.6543%** | +5.5242 percentage points |
| AE4=0.00 | +26.8944% | **+33.1013%** | +6.2070 percentage points |

Compensation is relative to each law's own WT positive NKCC1 loading integral.
The Palk law therefore increases relative compensation. This statement does
not imply that every absolute Palk flux exceeds the corresponding Task 38 flux.

Relative to Task 39 WT, integrated NBC cycles fall to 0.225904 and 0.141647;
NHE1 to 0.474860 and 0.439918; pump cycles to 0.970521 and 0.968031; CaCC
export to 0.992134 and 0.990322; and paracellular return to 0.987396 and
0.985363. Residual AE4 loading at 5% expression is 0.106183 times WT because
the concentration-dependent flux per expressed carrier changes. No capacity
or regulation parameter was altered to obtain these responses.

## Mechanistic signs and sustained secretion

At 600 s:

| Intracellular quantity | WT | AE4=0.05 | AE4=0.00 |
| --- | ---: | ---: | ---: |
| Na_i (mM) | 20.600656 | 17.731120 | 17.489924 |
| K_i (mM) | 107.037165 | 112.879455 | 113.488003 |
| Cl_i (mM) | 51.206535 | 51.776290 | 51.751216 |
| pH_i | 7.035590 | 7.202125 | 7.222392 |
| X_i (mM^4) | 5781844.748 | 5365533.712 | 5315926.854 |

AE4 loss raises K by **5.84229/6.45084 mM**, but lowers Na by
**2.86954/3.11073 mM**. Thus K has the direction described by the
[2018 mechanism](https://pmc.ncbi.nlm.nih.gov/articles/PMC5792321/), while Na
does not. Both Na perturbations are slightly above WT at 60 s and below WT
by 120 s; the sustained Na response is a decrease. pH rises by
0.166535/0.186803, rather than remaining near WT.

For the prescribed law,

\[
F'(X)=-\frac{a_2a_3+a_4a_1}{(a_3+a_4X)^2}<0.
\]

The endpoint product X is 7.20032%/8.05829% below WT despite the rise in K.
The reduced product increases inward NKCC1 flux at the identical activity
multiplier. This algebra and the observed cation responses explain why
installing the source core alone does not remove the compensation in this
architecture. They do not identify a unique cause for all downstream responses.

Chloride initially falls relative to WT: the most negative published
60-second-grid differences are **−2.21683/−2.42484 mM at 120 s**. This is a
transient change, not a sustained chloride deficit. At 600 s chloride is
**0.569755/0.544681 mM above WT**.

The cumulative secretion deficit first grows and then collapses:

| Time (s) | Cumulative deficit, 5% (%) | Cumulative deficit, null (%) | Flow deficit, 5% (%) | Flow deficit, null (%) |
| --- | ---: | ---: | ---: | ---: |
| 0 | undefined (0/0) | undefined (0/0) | 0 | 0 |
| 60 | 1.0456 | 1.1197 | 2.2918 | 2.4709 |
| 120 | 1.8695 | 2.0262 | 2.9417 | 3.2348 |
| 180 | 2.1540 | 2.3594 | 2.4706 | 2.7836 |
| 240 | 2.1335 | 2.3618 | 1.6371 | 1.9197 |
| 300 | 1.9532 | 2.1867 | 0.7407 | 0.9664 |
| 360 | 1.6927 | 1.9196 | −0.1062 | 0.0501 |
| 420 | 1.3958 | 1.6081 | −0.8653 | −0.7823 |
| 480 | 1.0870 | 1.2796 | −1.5288 | −1.5181 |
| 540 | 0.7799 | 0.9497 | −2.1015 | −2.1590 |
| 600 | 0.4822 | 0.6273 | −2.5927 | −2.7129 |

These tabulated extrema are not claims about exact continuous peak times.
The complete cumulative difference in pL at each published time is in
`flow_and_cumulative_comparison.csv`.

## Held-out phenotype and limits

The null prediction has a **0.627320% cumulative deficit** and a **2.712927%
endpoint-flow surplus**. It does not reproduce the task's held-out
approximately 35% experimental deficit or the original 2018 model's
approximately 24% flow deficit. The latter is a flow-rate comparison, not
an identical 0–600 s cumulative protocol. Both present observables are
reported to keep that distinction visible. Neither held-out value was used
to select or adjust any parameter.

The Task 38 cumulative deficits were 1.568465% and 1.839548%; the Palk
replacement makes the already small cumulative phenotypes smaller.
This result concerns the published NKCC1 concentration-response core inside
the current Task 38 architecture, including its inherited N1 activity
multiplier. It is not a reproduction or refutation of the entire 2018 model.

## Evidence and numerical budget

Results are in `results/39_palk_nkcc1_full_validation/`. The three
`*_timeseries.csv` files contain the complete requested panel at
0, 60, ..., 600 s: flow, Na/K/Cl, pH, HCO3, TIC, cell volume, NKCC1 cycles and
chloride loading, AE4 and signed AE2, NBC cycles/HCO3/current, NHE1, total pump
cycles, CaCC export, paracellular return, all voltages, lumen osmolarity,
and X. The corresponding `*_verification.json` and `*_integrated_fluxes.csv`
files retain all gates, solver evidence, and signed/positive integrals.
`integrated_flux_comparison.csv`, `intracellular_comparison.csv`,
`nkcc1_compensation_comparison.csv`, and `comparison_summary.json` provide
the ratios, differences, and mechanism comparisons.

Exactly **three intended integrations** were executed, with **zero numerical
retries**, **zero stationary solves**, **zero optimisation calls**, and
**zero sweeps**. One scientific worker and one thread in each loaded BLAS
pool were verified. Total recorded numerical execution, including source
tests and nesting, was **6.507844 seconds**, below the 900-second limit.
Subsequent work only summarized saved text results and prepared publication.
No further scientific execution is authorized by this completed task.

Publication note: during execution the remote branch was rebased to
`8c10ac041c36a9feff1b538930daf7ea32e730f0` (repository lineage consolidation).
An exact Git-blob comparison verified that all scientific inputs and task
instructions are unchanged from the prepared head. Among the 117 prepared
files, only README, the ignore file, and the archive index differ or are
absent in that rebased repository tree. Publication preserves that updated
tree and history and applies only this Task 39 patch. The recorded numerical
freeze remains the original executed snapshot; no calculations were rerun.
