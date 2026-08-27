# Historical MATLAB execution and reproduction

## Outcome

The archive contains an exactly constructed resting-state calibration and two later reduced dynamical functions. The calibration arithmetic in `Original/Saliva_Ae4.m` closes all 12 displayed balances to a maximum absolute residual of `1.39e-17`, and exact transformations of that calibration reproduce all 45 fields in `Par.mat` to a maximum relative discrepancy of `2.51e-15`. This is strong implementation-lineage evidence.

It is not a reproduction of the central 2018 result. The only archived full driver is internally incomplete, the later driver does not run a WT/AE2-KO/AE4-KO comparison, the archived calcium protocol differs from the published one, and translated runs of the available equations do not produce the reported approximately 24% AE4-knockout flow reduction. No retuning was performed.

Evidence levels in this note are:

- **Historical implementation evidence:** source expressions, dependencies, fitted baseline, and archived protocol.
- **Translated execution evidence:** direct evaluation of those expressions and SciPy BDF trajectories because no MATLAB-compatible runtime is installed.
- **Published comparison:** targets taken from the 2018 article and the completed Task 10 audit.

## Evidence integrity and runtime

No file under `archive/` was written or executed in place. The seven MATLAB/MAT inputs were copied byte-for-byte to `tmp/forensics/historical_execution_copy/`; the helper checks each source/copy SHA-256 pair before use.

| Archived input | SHA-256 |
| --- | --- |
| `Original/Parameters.m` | `fa19bfabd63f0798bbbf72e269a4d469cf00b2b41434bc7d32f614163db20eef` |
| `Original/Saliva_Ae4.m` | `7d72a055ab240f658287643520558a592a618845945c01bbdc4cb454f84fdcf3` |
| `Par.m` | `06ae901e1c045e234ec08d808bf909539b5d6d8cf009ebd98797791b7d313eff` |
| `Par.mat` | `fea3bdbe49dd794d43204bb234bd2b31e6215d784d920f46d3b4300f7068e1ca` |
| `Salivary.m` | `a6524879c8d948d8be41cdc081c2423b1bda8510943ddd3ce76fab41dc87ed8e` |
| `Salivary2.m` | `94e8ddc12e26d4cac0437f3e3d946014c91a1dcffa8fac6072aed56ac010b923` |
| `Salivary_ex.m` | `5f3ecaa2473a3bb409a1ff3dac980dcf356897f42923c50d092ae847f1c6f760` |

Neither `matlab` nor `octave` is available on `PATH`. Literal source execution was therefore impossible. Even with a native runtime, `Original/Saliva_Ae4.m:135` calls `Salivary(t,x,50e-3,par)` with four arguments, while the only archived `Salivary.m` declares six arguments, `Salivary(t,x,ca,par,g4,g2)`. There is no `Original/Salivary.m`. The original calibration block can be evaluated, but its final ODE call cannot be completed from the archived dependency set as found.

The forensic helper is `analysis/11_forensic_reconstruction/tools/reproduce_historical_matlab.py`. It directly transcribes the immutable MATLAB arithmetic and uses `scipy.integrate.solve_ivp(method="BDF")` as a documented substitute for `ode15s`. The archived tolerances, `RelTol=1e-6` and `AbsTol=1e-6`, are retained. The discontinuous calcium step is integrated as two segments so it is not stepped across silently.

Reproduction commands from the repository root are:

```bash
command -v matlab
command -v octave
python analysis/11_forensic_reconstruction/tools/reproduce_historical_matlab.py --check
python -m unittest discover -s tests/forensics -p 'test_historical_reproduction.py' -v
```

The translated BDF trajectories are solver-specific forensic evidence. Exact baseline arithmetic and parameter transformations are separately tested and do not depend on the BDF trajectories.

## Dependency combinations and entry points

| Component | Supported dependency combination | Finding |
| --- | --- | --- |
| `Original/Parameters.m` | Script workspace for `Original/Saliva_Ae4.m` | Defines geometry, bath/resting values, raw M-based NKCC coefficients, water coefficients, and preliminary exchanger constants. It is not a parameter function. |
| `Original/Saliva_Ae4.m` calibration block | `Original/Parameters.m`; also loads outer `Par.mat` before the ODE section | Constructs conductances and densities algebraically so the selected resting row is a steady state, then saves an unarchived `Conductances.mat`. The calibration portion is self-consistent. |
| `Original/Saliva_Ae4.m` ODE section | Requires a missing four-argument `Salivary` variant | Not runnable from the archive as found. Arbitrarily supplying the later six-argument function would change the dependency contract. |
| `Par.m` | Standalone generator for `Par.mat` | Manually records the Original calibration and exact coefficient transformations. It calls `clear` and `save('Par.mat')`; it was not run because doing so in `archive/` would overwrite evidence. |
| `Salivary_ex.m` first solve | `Par.mat` + six-argument `Salivary.m` | Runnable dependency set in MATLAB/Octave. Despite the comment “Full Model,” it explicitly sets `g2=par.g2*0`, so this solve is AE2 knockout, not WT. |
| `Salivary_ex.m` second solve | `Par.mat` + `Salivary2.m` | Uses the unmodified `par.g2` and `par.g4`, hence is WT for the fixed-lumen reduction. It is not the matching comparator implied by the first block. |
| `Salivary.m` | `Par.mat` struct plus explicit `g4,g2` arguments | Seven states: luminal Na/K, cell height, and intracellular Na/K/Cl/HCO3. H is eliminated by electroneutrality, CO2 algebraically, luminal Cl as Na+K, and both potentials by QSS formulas. |
| `Salivary2.m` | `Par.mat` struct | Five-state fixed-lumen reduction. Activities are embedded as `par.g2` and `par.g4`; five-state KO runs below therefore use an in-memory forensic copy of `par` with one field set to zero. The archive and its driver are unchanged. |

No archived MATLAB file generates the continuation, determinant, or sensitivity figures in the historical TeX. Those figures refer to a further three-state nondimensional system and MATCONT workflow whose implementation is absent from the seven archived MATLAB files.

## Exact baseline reconstruction

### Original calibration

Direct evaluation of `Original/Parameters.m` followed by `Original/Saliva_Ae4.m:13-123` gives:

- `CO2_i = 6.600123698596884` mM and `pH_i = 6.91`;
- `q_a = q_b` to `2e-19` in code units;
- `q_total = 0.000531415476760618` in undocumented code units; and
- a maximum absolute residual of `1.39e-17` across the 12 balance, volume, and QSS rows.

This exact closure is by construction. The script first evaluates each unscaled transporter/channel law, then solves algebraically for every fitted conductance or activity needed to cancel the selected resting fluxes. It is therefore evidence of internal implementation consistency, not independent validation of the baseline.

Every `Par.mat` field is recovered from this calibration. Important exact transformations include:

- M-based NKCC values `a2=2.0096e7` and `a4=1.3852e6` become `2.0096e-5` and `1.3852e-6` for later mM-state evaluation;
- the Original reverse AE4 constant is multiplied by `118.7219`, giving stored `k2=0.00159305881398359`;
- geometric `delta` is multiplied by `1e-15`, so `delta*Hi` is a volume in litres; and
- the fitted conductances/activities are copied into the later `par` names.

The maximum absolute field difference is `1.49e-8`, occurring only because a roughly `8.75e7` conductance was rounded at the final displayed digits; the maximum relative field difference is `2.51e-15`.

### Later outer-model baseline

At the archived initial state and calcium value `0.05`, both later functions recover the fitted baseline. For the seven-state function the maximum unscaled RHS magnitude is `5.24e-12`; after its hidden `10000` multiplier it is `5.24e-8`. The five-state function has the same maximum unscaled residual and a `1000` multiplier.

| Quantity | Historical evaluation | Published resting target | Status |
| --- | ---: | ---: | --- |
| intracellular Na, K, Cl | `25`, `120`, `50` mM | `25`, `120`, `50` mM | Exact selected values |
| intracellular HCO3 | `10` mM | `12.1` mM | Substantive baseline difference |
| intracellular pH | `6.9100000003` | `6.91` | Exact selected value |
| intracellular CO2 | `6.6001236986` mM | `6.6` mM | Exact to displayed precision |
| luminal Na, K, Cl | `118.7`, `5.6`, `124.3` mM | same | Cl is algebraically Na+K |
| apical, basolateral potentials | `-50.24`, `-62.8` mV | same | Exact fitted QSS values |
| cell volume | `1.3000000000` pL | `1.3` pL | From `delta*Hi` |
| total flow | `0.0005314154768` code units | reported in microlitre/min | Not comparable: code unit/normalization is undocumented |

Thus the historical convention resolves the Task 10 baseline electrical residual only by using its fitted current/flux normalization, not by reproducing the published Table 2 conductances literally.

## Conventions exposed by executable arithmetic

Several published blockers have identifiable historical causes:

- **Water:** `Parameters.m` first makes `b2=7*b1`, then multiplies only `b1` by `74.4`. The stored ratio is therefore `b1/b2=10.6285714286`. At the selected baseline the two code osmotic brackets are `0.7997659713` and `8.5003697524`, yielding exact `q_a=q_b`. This is an implementation-specific calibration, not the published permeability pair.
- **NKCC:** the Original M-valued coefficients plus explicit `s=1e-3` concentration conversion are algebraically identical to the later mM coefficients in `Par.mat`. This independently confirms the Palk-corrigendum conversion used in Task 10.
- **Electrical QSS:** the calibration solves `GtNa`, `GtK`, `GCaCC`, and `GCaKC` from the target voltages and resting flux. The later QSS formulas recover both potentials exactly. The numerical conductances and implicit units differ drastically from the published nS table, so this is fitted historical normalization rather than a transparent paper/code unit match.
- **Acid/base:** the outer dynamics eliminate H by electroneutrality and CO2 by the algebraic equality of the code's aggregate CO2 term and a scaled buffer term. This yields `CO2_i=6.6001237` mM without imposing the paper's separate equilibrium formula that gave `3.52` mM in Task 10. The code retains the printed aggregate CO2 sign, uses `PCO2=1970`, and introduces the fitted multiplier `gb=0.1996525589`; it does not resolve the physical sign or published units.
- **AE2:** the code stores `KB=1e-4`, rather than the `10^4` value transcribed from the published table. This makes the selected AE2 flux positive and is consistent across Original and later code.
- **AE4 and cations:** both historical dynamical functions put the entire AE4 cycle into the Na balance and no AE4 term into the K balance. They do not implement the published Na/K fractional partition. Their stored `g4`, effective reverse constant, and units also differ from Table 7. This is a substantive model difference, not merely a hidden scale conversion.
- **NHE1:** the historical proton factor is first power, whereas the published reconstructed expression is squared. This is another substantive equation difference.
- **Time:** `Salivary.m` multiplies every RHS by `10000`; `Salivary2.m` uses `1000`. Both switch calcium from `0.05` to `0.55` when `t>100` and label plots in seconds. These hidden time scalings and the step input do not match the published qualitative minute-scale calcium trace.

## Archived-protocol trajectories

The translated harness starts every scenario from the archived WT initial state, integrates calcium `0.05` over `t=0..100`, then `0.55` over `t=100..200`, and makes no parameter adjustment other than setting the named exchanger activity to zero. Seven-state knockouts use the function's explicit activity arguments. Five-state knockouts require an in-memory copy of `par`, because `Salivary2.m` does not expose activity arguments and the archived driver invokes it only with WT fields.

The table reports the archived-protocol endpoint at `t=200`, not an assumed steady state. “First post-step sample” is the maximum sampled `q_total` on the one-unit output grid; it occurs at `t=101` in all runs.

| Function | Scenario | Endpoint flow / WT | First post-step sample / WT | Max endpoint RHS after hidden scale |
| --- | --- | ---: | ---: | ---: |
| seven-state `Salivary` | WT | `1.000000` | `1.000000` | `3.49e-5` |
| seven-state `Salivary` | AE2 KO | `0.999052` | `0.997382` | `5.58e-5` |
| seven-state `Salivary` | AE4 KO | `0.998556` | `1.101877` | `5.71e-5` |
| five-state `Salivary2` | WT | `1.000000` | `1.000000` | `3.52e-5` |
| five-state `Salivary2` | AE2 KO | `1.000000` | `1.002023` | `5.58e-5` |
| five-state `Salivary2` | AE4 KO | `1.059019` | `0.903297` | `7.22e-2` |

The seven-state endpoint is close to stationary in its own raw residual scale, but its AE4 knockout changes endpoint flow by only `-0.144%` and initially raises, rather than lowers, the sampled stimulated flow. The five-state AE4 knockout initially lowers the sampled flow by `9.67%`, but the `t=200` endpoint is still moving (`max |dx|=0.0722` after the hidden scale) and its flow is `5.90%` above WT. AE2 knockout has a negligible flow effect in both versions.

None of these archived-protocol outputs reproduces the approximately 24% AE4-knockout flow loss reported in the 2018 article or the approximately 25% continuation result described in the historical TeX. The continuation result must have come from the separate three-state nondimensional/MATCONT workflow, whose code is not in the archive.

Full state, potential, flux, pH, volume, and flow trajectories for all six translated runs are in `results/11_forensic_reconstruction/reproduction_trajectories.csv`. Baseline residuals are in `reproduction_residuals.csv`; solver metadata, every `Par.mat` value, calibration values, endpoint states, and comparison metrics are in `reproduction_summary.json`.

## Reproduction verdict

1. **Baseline:** reproduced exactly as a historical fitted construction. It matches many displayed 2018 resting values but uses HCO3 `10` mM, different equations, and hidden calibration/scaling conventions.
2. **WT dynamics:** translated successfully for the archived step protocol; no machine-readable published trajectory exists for exact numeric comparison, and the input protocols differ.
3. **AE2 phenotype:** negligible flow change is reproduced qualitatively.
4. **AE4 phenotype:** not reproduced. Neither available dynamical variant yields the reported central decrease under its archived protocol.
5. **Implementation identity:** the files are directly lineaged, model-related calibration/reduction evidence, but this execution does not support treating either later function as the executable implementation that generated the 2018 figures.
6. **Physiological-map gate:** failed. A fitted baseline plus nonmatching knockout dynamics, substantive equation differences, undocumented units, and missing native dependencies are insufficient to validate a physiological parameter-to-observation map.
