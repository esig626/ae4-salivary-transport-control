# Identity evidence and conditional gate decision

This is the final identity-evidence synthesis based on the archive audit, equation lineage, public-literature cross-check, and no-retuning historical reproduction. [final_assessment.md](final_assessment.md) is the controlling project decision if its integrated classification differs from any supporting classification here.

## Decision summary

The archive supports four different propositions at different grades:

| Proposition | Grade | Reason |
|---|---|---|
| The historical materials belong to the same AE4 salivary-secretion project/model family as the 2018 article | **CONFIRMED** | Explicit manuscript attribution, matching topology and variables, exact target states, distinctive parameter fingerprints, and the same AE2/AE4 question |
| `Original/Parameters.m` and `Original/Saliva_Ae4.m` are the direct calibration source of `Par.m`/`Par.mat` and the reduced MATLAB functions | **CONFIRMED** | All 45 stored parameters are reproduced from the Original calculation or exact unit transformations to floating-point precision |
| The archived snapshot is an immediate analytical or development branch of the 2018 model | **PLAUSIBLE BUT UNPROVEN** | Strong structural and numerical ancestry, but no complete prepublication dependency set or generating-script chain survives |
| The exact archived MATLAB snapshot generated the published 2018 simulations and figures | **UNLIKELY** | Substantive equation differences, a postpublication parameter bundle, broken dependencies, no matching published figures, and failure to reproduce the central AE4-knockout flow reduction |

For the prompt's operative question—“is the historical code actually the 2018 implementation?”—the recommended grade is therefore **`UNLIKELY` for the exact archived snapshot**. This does not deny direct ancestry. It rejects the stronger claim that these files, as preserved, are the publication-generating implementation.

The conditional clean-reconstruction gate is **closed**, and the Task 10 identifiability project should **not** be reopened on this historical branch.

## Meaning of “identity”

This report treats “the 2018 implementation” as the executable equation/parameter/protocol set that generated the article's central baseline, WT, AE2-knockout, AE4-knockout, and figure results. Sharing authors, names, scientific intent, or even most equations establishes lineage, not exact identity. A postpublication reduction may be a genuine descendant while still being the wrong object for reproducing or inferring the 2018 model.

## Evidence matrix

The weight column concerns exact implementation identity. “For” and “against” are kept separate so a strong lineage fingerprint is not mistaken for output authentication.

| Dimension | Evidence for identity or direct ancestry | Evidence against exact snapshot identity | Weight | Dimension verdict |
|---|---|---|---|---|
| Equation identity | Same cell/lumen/interstitium topology; NKCC1, Na/K pump, CaCC, CaKC, AE2, AE4, NHE1, buffer, water and tight-junction terms; same seven dynamic quantities after algebraic elimination; cell-height balance corresponding to \(d\omega_i/dt=q_b-q_a\); luminal electroneutrality \(Cl_l=Na_l+K_l\); quasi-steady membrane potentials | Historical dynamics assigns the **entire** AE4 cation flux to Na and none to K, whereas published Eqs. 2–3 partition it between Na and K; historical NHE1 uses first-power H gating whereas the paper uses a square; external bicarbonate is hard-coded as 21 mM in exchangers; channel gates and water coefficients differ; Ca and time protocols differ | **Very high** | Strong common skeleton; **substantive model differences** prevent exact equation identity |
| Parameter identity | Exact or near-exact transformed fingerprints for AE2, AE4, NHE1, tight-junction conductances, NKCC, target voltages, geometry, and \(\Psi_l\); `Par.mat` matches `Par.m` field-for-field | Resting bicarbonate and bath values differ; water ordering is reversed by a hidden multiplier; channel half-activation values are offset by approximately 100; conductances are recalibrated to compensate; multiple normalizations are undocumented | **High** | Strong direct lineage, mixed exact identity |
| Output reproduction | Historical prescribed WT row closes its own 12 residuals to \(1.39\times10^{-17}\), and the seven-state WT row is a numerical steady state under its own equations | Seven-state AE4 knockout gives only \(-0.144\%\) endpoint flow change, not the published approximately \(-24\%\); the five-state fork gives \(+5.90\%\); no archived path reproduces the published WT/KO figure set | **Decisive** | **Against exact identity** |
| Figure correspondence | The Golubitsky manuscript plots AE2/AE4 continuation and sensitivity results for a reduction explicitly attributed to the Sigüenza model | The archived figures are determinant, normalized-flow, state, and sensitivity plots, not duplicates of 2018 Figs. 1–10; the archived MATLAB driver does not save them; the MATCONT/continuation source is absent | **High** | Analysis-of relationship, not generating-figure identity |
| Chronology | EPS assets dated June 2017 and the Golubitsky manuscript build dated September 2017 are contemporary with the paper's online-first period and call the source model “Sigüenza et al. (2017)” | `Par.mat` was created 10 October 2018, after the article's 5 December 2017 online publication and 2018 issue; no embedded date ties the preserved MATLAB snapshot to the submitted simulations | **Medium-high** | Supports contemporary ancestry; weighs against `Par.mat` as prepublication proof |
| Dependency structure | `Original` calibration quantities transform exactly into the outer parameter bundle; `Salivary.m` is a mathematically recognizable elimination/reduction of the calibration ledger | `Original/Saliva_Ae4.m` calls `Salivary` with too few arguments and writes a missing `Conductances.mat`; `Salivary_ex.m` calls `Salivary2` with too few arguments; the preserved driver runs AE2 knockout rather than WT; no complete figure-generation chain exists | **High** | Direct internal lineage, but **not a self-contained release** |
| External provenance | 2019, 2020, 2021, and 2022 primary papers explicitly reuse descendants of the 2018 secretion model; later literature independently confirms the conservation sign and corrected NKCC scale | No 2018 correction, supplement, or public code was located; later models change geometry, activities, water parameters, channel placement, calcium, and AE4 equations | **Medium** | Confirms model-family continuity, not this snapshot's identity |

## Quantitative parameter fingerprints

The historical model uses cell height \(H_{i0}=28.7\ \mu\mathrm{m}\) and mean membrane area \(A=45.296167\ \mu\mathrm{m}^2\). The following transformations are too specific to be accidental and are the strongest numerical evidence for a direct common calibration lineage.

| Quantity | Published value | Normalization suggested by historical geometry | Expected historical-scale value | Archived value | Relative difference |
|---|---:|---|---:|---:|---:|
| NHE1 activity | \(0.0305\) fmol/s | \(1000G/A\), converting fmol to amol per mean area | 0.673346 | `g1 = 0.674737` | +0.207% |
| AE2 activity | \(0.01807\) fmol/s | \(1000G/A\) | 0.398930 | `g2 = 0.399085` | +0.0388% |
| AE4 activity | \(0.66\) amol/\(\mu\mathrm{m}^3\) | archived per-area coefficient divided by height | `g4 / 28.7 = 0.660288` | `g4 = 18.950261` | +0.0436% after division |
| Tight Na conductance | 12.46 nS | archived coefficient divided by height | `gtna / 28.7 = 12.469105` | `gtna = 357.863312` | +0.0731% after division |
| Tight K conductance | 0.9 nS | archived coefficient divided by height | `gtk / 28.7 = 0.902201` | `gtk = 25.893178` | +0.245% after division |
| Luminal impermeant osmolyte | 48.8 mM | none | 48.8 | `xl = 48.8001357` | +0.000278% |
| NKCC \(a_2,a_4\) | corrigendum: \(2.0096\times10^7,1.3852\times10^6\) with M concentrations | multiply fourth-order coefficients by \(10^{-12}\) for mM | \(2.0096\times10^{-5},1.3852\times10^{-6}\) | exact in `Par.mat` | floating-point exact |

These fingerprints establish derivation from the same target model. They do **not** establish which side of a later fork generated the paper: several are explicitly obtained by solving the historical baseline equations, and published values may themselves be rounded outputs of that calibration.

## Material mismatches

| Feature | 2018 publication | Historical snapshot | Identity effect |
|---|---|---|---|
| AE4 cation balance | \(-[Na_i/(Na_i+K_i)]J_{Ae4}\) in Na and \(-[K_i/(Na_i+K_i)]J_{Ae4}\) in K | `-j4` in Na; no `j4` term in K | Different dynamical model and mechanism |
| NHE1 H gate | Squared intracellular-H saturation factor | First power | Different nonlinearity and sensitivity |
| Intracellular bicarbonate baseline | 12.1 mM | 10 mM | Published baseline not reproduced |
| External bicarbonate | 42.9 mM in the table | 40 mM assigned, 21 mM hard-coded in AE2/AE4 | Multiple incompatible bath conventions |
| Channel half-activation | \(K_{CaCC}=0.26\ \mu\mathrm{M}\), \(K_{CaKC}=0.182\ \mu\mathrm{M}\) | 100.26 and 100.3 in code-native units, with enormous fitted conductances | Compensating calibration, not a transparent unit conversion |
| Water coefficients | Published/Palk ordering has basolateral much larger than apical | `b1` is multiplied by 74.4 after setting `b2=7*b1` and `b3=0.94*b1`, leaving \(b_2/b_1=0.09409\) and \(b_3/b_1=0.01263\) | Resolves the historical row but reverses the public ordering |
| Calcium | Published multi-minute stimulation trace | Fixed 0.05, abrupt step to 0.55 at code time 100 | Different input protocol |
| Time | Article reports physical minutes | Hidden whole-RHS multiplier 10,000 in seven-state code and 1,000 in five-state code | No common physical time without extra evidence |
| State representation | Full paper DAE with explicit balances/constraints | Seven-state and five-state reductions with H, CO2, luminal Cl, and voltages eliminated | Reduction could be valid, but only if all constitutive terms were otherwise identical; they are not |

## Reproduction evidence

Neither MATLAB nor Octave is installed. The reproduction therefore evaluates a direct transcription of the immutable MATLAB arithmetic and uses SciPy BDF in place of `ode15s`. The preserved `Original` driver could not have been run literally against the archived outer function in any case because its call omits the two required exchanger arguments. Full diagnostics are in [reproduction_summary.json](../../results/11_forensic_reconstruction/reproduction_summary.json), [reproduction_residuals.csv](../../results/11_forensic_reconstruction/reproduction_residuals.csv), and the final reproduction report.

### Baseline

- The 12-equation `Original` calibration row has maximum absolute unscaled residual \(1.39\times10^{-17}\). This verifies the algebra and the `Original` to `Par.mat` lineage, but it is expected because the driver solves activities and conductances from that same row.
- The seven-state WT initial row has a maximum scaled residual of about \(5.24\times10^{-8}\), so it is also a numerical steady state under the historical conventions.
- That row has \(HCO_{3,i}=10\) mM, not the paper's 12.1 mM. Internal closure is therefore not reproduction of the complete published baseline.

### Knockout phenotype without retuning

| Historical run | Endpoint flow ratio to historical WT | Endpoint percent change | Peak-sample flow ratio | Published comparison |
|---|---:|---:|---:|---|
| Seven-state WT | 1.000000 | 0.000% | 1.000000 | Control |
| Seven-state AE2 knockout | 0.999052 | -0.0948% | 0.997382 | Qualitatively consistent with “no significant effect” |
| Seven-state AE4 knockout | 0.998556 | **-0.144%** | 1.101877 | **Inconsistent** with approximately -24% |
| Five-state AE4 knockout | 1.059019 | **+5.90%** | 0.903297 | Wrong endpoint direction |

The central discriminator is not reproduced. Agreement of AE2 knockout with a negligible effect is insufficient because many related branches can preserve that qualitative result. The large AE4-knockout discrepancy is strong evidence against this snapshot as the published generating implementation.

## Integrated interpretation

The most coherent history is:

1. A common full secretion-model calibration generated the distinctive target row and parameter family.
2. A sodium-only or otherwise earlier AE4 balance was used in the preserved calibration/reduction branch.
3. The Golubitsky project reduced and analysed that branch in 2017–2018.
4. The published article used a different or later branch for its Na/K-partitioned AE4 equations and central knockout result.
5. `Par.mat` and the outer functions preserve a postpublication state of the analysis branch, with signature drift and exploratory drivers.

This is an inference from combined evidence, not a recovered commit history. An alternative is that the paper's printed Na/K allocation differs from its private generating code. The unchanged historical phenotype argues against that alternative: the preserved sodium-only branch does not generate the paper's approximately 24% AE4-knockout reduction either.

## Conditional clean-reconstruction gates

The prompt permits a clean implementation only if identity reaches at least `VERY LIKELY` **and** the historical implementation reproduces the central published baseline/phenotype sufficiently well.

| Gate | Required evidence | Current result | Decision |
|---|---|---|---|
| Identity | At least `VERY LIKELY` | Exact snapshot graded `UNLIKELY` | **FAIL** |
| Published baseline | Resting states, potentials, balances, and flow under one coherent convention | Historical row closes internally but uses different bicarbonate, water, channel, and exchanger conventions | **FAIL** |
| Central phenotype | WT, AE2-KO, and approximately 24% AE4-KO flow reduction without tuning | AE4-KO is -0.144% in seven-state and +5.90% in five-state | **FAIL** |
| Dependency integrity | Complete executable historical version or an unambiguous dependency combination | Broken arities, missing output dependency, no publication figure generator | **FAIL** |
| Convention provenance | Each repair traceable as published, historical evidence, or new decision | NKCC and volume sign are resolved; AE4, water, acid/base, electrical, Ca, and time conventions remain mixed | **FAIL** |
| No arbitrary retuning | Central reproduction before calibration changes | Central reproduction already fails as found | **FAIL** |

**Decision:** do not create `src/reconstructed_2018/` from this snapshot. Doing so would require choosing among substantive equations and calibrations, not merely transcribing a recovered implementation.

## Conditions that could change the decision

Reconsideration would require materially new evidence, such as:

1. a complete, independently dated prepublication code bundle or author-supplied generating archive;
2. explicit resolution of the AE4 Na/K allocation, NHE1 exponent, bicarbonate baths, channel gating, water calibration, calcium input, and time normalization;
3. unchanged reproduction of the published resting state, WT outputs, negligible AE2-KO effect, and approximately 24% AE4-KO flow reduction;
4. a traceable link from the executable version to at least one central published figure or underlying numeric dataset; and
5. regression tests separating published facts, recovered historical conventions, and genuinely new modelling decisions.

A parameter retune that merely forces the knockout percentage would not satisfy these conditions.

## Identifiability and discrimination recommendation

The physiological steady-state map required by Task 10 has not been recovered. The historical branch's AE4 activity-to-flow response is demonstrably not the published phenotype, so computing \(F_u\), \(F_\theta\), observation geometry, or mechanism overlap on this branch would answer a different-model question.

The appropriate project classification remains **`STOP`** for reopening identifiability/discrimination. The analytical Task 10 results remain valid only as conditional mathematics under their stated assumptions; they must not be promoted to a physiological result until a validated model map passes the gates above.

## Uncertainties and limits

- Literal MATLAB/`ode15s` execution was unavailable. The direct arithmetic transcription is highly diagnostic but is not bit-for-bit MATLAB execution.
- The full Palk corrigendum body was blocked by the publisher during this audit; its indexed correction and coefficients were cross-checked against the later publisher supplement.
- The absence of indexed 2018 code or a correction is negative search evidence, not proof of nonexistence.
- File metadata establishes the `Par.mat` creation date but not the date at which every source line was first written.
- A missing private ancestor could explain the shared fingerprints and published outputs. It would be new evidence; it is not present in this repository.
