# Task52D central production summaries

**52D COMPLETE: all three authorised cases reached 600 s with Radau and passed all frozen sampled physical, conservation and reservoir gates.**

Classification: NUMERICAL PREDICTION under the frozen matched-supply and onset-projection MODEL IDEALISATIONS. This is the central52D production boundary only; no52E sensitivities/protocol controls or52F comparison/final classification were performed.

The scientific inputs are exactly52C at `36f7c3d0bb2353bf55ac0ed9e22582720097c829`, freeze SHA `9cebb8fcf85557643f6c7838ff57f9fb9603f5a7479cbe6fb63ba8c40b417ba2`. Authorisation is `f2b56f903d33a0ac5f5c292db92cead8b6157f1f`. All50 immutable scientific/input hashes remain unchanged. The failure at `9c0fde42b5eb16746413801fae555c84cb32e201` remains preserved.

## Secretion and total apical chloride export

Headline integrals use the frozen one-second trapezoidal rule, including zero/epsilon onset samples. Deficits are100×(1−KO/WT).

| Case | G_aux (S) | WT fluid (pL) | KO fluid (pL) | Fluid deficit (%) | WT Cl export (fmol) | KO Cl export (fmol) | Cl deficit (%) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| case_01 | 0 | 0.9729547286 | 0.7779801601 | 20.039429 | 228.486672438 | 186.884020174 | 18.207912 |
| case_02 | 2.32e-09 | 0.9786120576 | 0.7821483215 | 20.075753 | 229.632497270 | 187.766308003 | 18.231822 |
| case_03 | 4.49e-09 | 0.9827438302 | 0.7851967712 | 20.101582 | 230.468167539 | 188.410807218 | 18.248663 |

| Case | WT mean flow60–600 (pL/s) | KO mean flow60–600 (pL/s) | Broad deficit (%) | WT flow600 (pL/s) | KO flow600 (pL/s) | Endpoint deficit (%) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| case_01 | 0.001614478343 | 0.001295413977 | 19.762691 | 0.001562785633 | 0.001235633930 | 20.933882 |
| case_02 | 0.001623209257 | 0.001302003487 | 19.788316 | 0.001569295997 | 0.001241322247 | 20.899419 |
| case_03 | 0.001629574940 | 0.001306819642 | 19.806103 | 0.001574014570 | 0.001245475263 | 20.872698 |

## Independent chloride-reservoir audit

`delta(t)=nCl_WT(t)−nCl_KO(t)` is stored at every observation. Independently evaluated8-node Gauss flux quadrature tests `D_J=delta(0)+2D_N+D_A+D_E−delta(600)`, with4-node cross-check and separately recorded epsilon contribution. All quantities below arefmol; N is in cycle units before the displayed factor2.

| Case | delta(0) | 2D_N | D_A | D_E | delta(600) | D_J | D_J−RHS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| case_01 | 19.698594732 | 0.000000000 | 46.636713087 | 0.000000000 | 24.735189577 | 41.600118311 | 6.875e-08 |
| case_02 | 19.698594732 | 0.000000000 | 46.652417676 | 0.000000000 | 24.487456436 | 41.863556044 | 7.225e-08 |
| case_03 | 19.698594732 | 0.000000000 | 46.663679031 | 0.000000000 | 24.307622169 | 42.054651669 | 7.488e-08 |

All reservoir errors are below1e-5fmol, every individual/differential8-versus4 discrepancy is below1e-6fmol, and the onset bounds pass. D_N and D_E are exactly zero in these matched-supply cases. The coarser headline-trapezoid budget errors are separately recorded in each JSON and are not used as pointwise conservation errors.

| Case | Genotype | Total Cl Gauss8 (fmol) | TMEM16A (fmol) | Auxiliary (fmol) |
| --- | --- | ---: | ---: | ---: |
| case_01 | WT | 228.484671477 | 228.484671477 | 0.000000000 |
| case_01 | KO | 186.884553166 | 186.884553166 | 0.000000000 |
| case_02 | WT | 229.630337913 | 199.310260419 | 30.320077494 |
| case_02 | KO | 187.766781869 | 162.974311376 | 24.792470493 |
| case_03 | WT | 230.465886816 | 178.046452361 | 52.419434455 |
| case_03 | KO | 188.411235147 | 145.557125466 | 42.854109681 |

## Output coverage and software repair

Complete compact trajectories and per-case summaries are under `output/reporter_recovery_52D/cases/`. Each NPZ preserves26 paired states, all original named diagnostic values, integer-second and accepted-step records, delta(t), raw accepted states, cumulative Gauss8 flux integrals and an explicit presence mask. Each union has156 columns; exactly four absent entries are two optional current fields at time zero for two genotypes. Missing values are NaN with a false mask, never invented zero; JSON missing endpoints use null. All present diagnostics and states are finite. The original physical/conservation gates remain unchanged.

Per-case JSONs preserve complete reservoir audits, onset terms,4/8 checks, physical residual maxima, all ten cumulative observation windows, endpoint electrochemical/flow states, extrema and separately integrated transporter/current components. CSV and `output/summary_52D.json` retain machine-readable stage summaries.

Independent file-only review checked all three definitions/onsets, summary/trajectory hashes, full integer grids, masks, delta, matched supply/shared conductances, KO AE4 zero, physical limits, reservoir and trapezoidal summaries. No model was evaluated by this review. A pre-model output-path setup error was preserved separately and corrected only in authorised bookkeeping; it created no case or numerical evaluation.

The successful fresh execution used one scientific process and one thread per BLAS pool; three Radau attempts, no BDF fallback, no fits/searches or extra cases. Counts are3 parent factories,5 constructor core evaluations,21,583 paired evaluations and43,166 local evaluations. The prior failed process’s unrecorded count remains unknown and is not added to these fresh counts.

No scientific parameter, transport law, diagnostic definition, onset, solver tolerance or quadrature changed. No Palk/Benjamin orTask41/Task50 coupling was loaded. Stop after independently verified52D publication;52E and52F remain unexecuted.
