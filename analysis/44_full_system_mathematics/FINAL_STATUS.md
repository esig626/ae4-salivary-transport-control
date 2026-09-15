# Task 44 final status

Complete full system mathematical analysis, standalone report and verification on the separate branch `analysis/task-44-full-system-mathematics`. This status records the published report checkpoint below. The containing status commit is a subsequent administrative checkpoint; its final branch head is reported separately, avoiding a self-referential commit SHA.

## Exact checkpoints

- Initial remote Task 44 checkpoint: `b7b775a0277d93263c4de339745806e5658ec9c2`.
- Independently replayed numerical checkpoint: `f783785a863440df459e9c5530beba4c7eb59110`.
- Published complete report and verification checkpoint: `fac497b5d09f219f4d197bec665c1dbd700e7610`.
- Final companion Task 43 head checked here: `547f113d9123ab1976774639a69483faf40cef34`. Its crosswalk uses the exact numerical checkpoint above, with seven source file hashes.
- Frozen scientific base: `c4d4f207f702d8eee09ab94d2d90ae90552dd641`.
- Unchanged main head: `1aba5ee753035613bcdfbd39b1e5a9d8cdcf61ea`.

## Completed results and their status

| Status | Result |
| --- | --- |
| Proved or exact | All 13 production states retained; two exact charge constraints leave 11 independent dynamics. Full nondimensional equations retain carbon, alkalinity, changing cell and lumen volumes, finite lumen composition, algebraic pH and current closure, and AE4 regulation. |
| Proved or exact | Carbonate and electrical uniqueness under explicit positivity/slope assumptions, smooth dependence within strict executable brackets, local well posedness, exact pH derivatives, alkalinity storage budget and conditional capacity obstruction. |
| Proved or exact | Equal routing gives `J_Cl,a = 6P - J_H + 2 dn_Na,i/dt - dn_A,i/dt - dn_Cl,i/dt`. The stationary explicit AE4 term cancels; state mediated sensitivity remains. The general coefficient is `nu_Cl - 2 nu_Na + nu_A`; the five joint source projections have rank five. |
| Local numerical | All 26 class/expression cases yield nonsingular stable local roots; 22 pass physiological gates. Null source classes coincide, so these are not 26 distinct roots at one parameter point. WT and null slow decay times are 345.82 and 459.92 s. Modal input/output projection and nonnormality are quantified. |
| Local numerical | Independently checked IFT expression and parameter derivatives; stationary null secretion deficit 0.94617%; WT local NKCC chloride compensation 92.901%. Sixteen uncertain parameter sensitivities cover secretion, pH, the null effect, alkalinity shares/capacity, compensation and slow modes. |
| Local numerical | A physiological NBC absent equilibrium exists: pH 6.96833154 and secretion 0.00153598288 pL/s, 3.035% below WT. Its AE4 cycle is 0.00632651 fmol/s, compared with WT 0.07507245. NBC supplies 95.48% of the WT AE4 alkalinity demand; the conditional no NBC capacity bound covers 57.96% of that specified demand. |
| Finite time numerical | Frozen Task 40 observables reproduce. Over 60–600 s the relative integrated NKCC increase is 23.16345%, while the fraction of missing AE4 chloride replaced is 88.50704%. These differ from the local derivative compensation and the 0–600 s cumulative secretion observable. |
| Finite time numerical | Within the one target selected Task 41 family, the first numerical crossings require 89.50727% (adaptive) or 89.50806% (original sampled observable) AE4 dependent stimulated CaCC recruitment. All inherited 600 s gates pass. Eight uncertain parameters have independently checked sensitivities of both crossings. |
| Finite time numerical | Separate continued stimulus null cases cross the upper pH gate at 1486.31 s (selected case) and 1484.75 s (adaptive crossing). Their late mathematical roots have pH above 7.3. This does not change the frozen 600 s results. |
| Unresolved | Global uniqueness, global stability/positivity, a justified singular reduction, global exclusion of unsampled earlier inverse crossings, finite range robustness without defensible uncertainty intervals, native pathway identity/capacity, molecular stoichiometry and compatible kinetics, and independent validation of the CaCC coupling. |

## Failures, caveats and corrected earlier claims

- The four physiological exclusions remain visible: C1 WT fails intracellular sodium; C2 WT fails sodium and pH; C4a and C4b WT fail pH. Solver success and physiological acceptance are separate.
- Thermodynamic opposition is retained at C0 WT, C3b WT/5% and C4b WT/5% under the inherited scalar kinetics and substituted source vectors. No class was rescued or refitted.
- The inherited checkpoint's four unsuccessful continuation attempts are preserved as an aggregate record. It contains no recoverable per-case seeds or errors. Documented nearby seeds now recover all 26 case roots; the historical failures are not evidence of nonexistence.
- Analysis differentiation across the production regulatory clamp produced a spurious 60 s regulatory mode. An inward derivative gives the correct 30 s mode. The 346/460 s chemical modes remain. Production source was not changed.
- Equal routing cancellation does not prove zero total AE4 sensitivity or a small phenotype for every parameterisation. The stationary null effect is not the 600 s cumulative 3.8575% effect.
- Universal NBC necessity for secretion is not supported: the valid NBC absent root selects much lower AE4 loading. The exact obstruction applies only to a declared loading demand and admissible capacity bounds. No specific NBC isoform is identified.
- The original selected Task 41 point has a sampled deficit of 30.26115%, slightly below the requested 30.3% comparison. Refined adaptive and sampled crossings remain distinct. The comparison is experimental mean minus one reported standard error, not a confidence bound or parameter uncertainty interval.
- The first numerical inverse crossing is not a proved global lower bound. Task 41 remains target selected. Later pH failure limits sustained physiological interpretation.
- Strong modelled NKCC compensation remains in tension with the isolated 2015 NKCC assay. A shared inherited WT initial state is not equivalent to established knockout animals. A slow eigenvalue alone does not establish the experimental delay.

## Verification

The numerical checkpoint was independently cloned and checked out detached. Its generated outputs and inverse caches were removed only inside that clean checkout, then all 18 analysis commands were rerun. All 47 authoritative JSON/CSV files were compared (46 regenerated results and one unchanged historical source record). All 331,157 numeric values matched bit for bit. Only runtime and ancestor-validated verification commit metadata were excluded. All 532 inverse trial cache files were regenerated without reuse. All 1,751 tracked files outside Task 44 remained unchanged.

The independent full RHS reconstruction verifies 84 random valid states across all seven source classes. Separate dimensional/dimensionless trajectories agree within 7.524e-10 in scaled state, and cumulative dimensional water recovery agrees within 1.443e-15 pL. Independent nearby roots, two derivative steps and matched eigenvalues verify local derivatives and spectra. Independent 13-state integration and adaptive quadrature verify the inverse observables.

The supplementary verifier independently solves 128 perturbed WT/null parameter roots and 256 compensation expression roots with every inherited Task 40 gate, including its equal routing architecture checks. All 384 pass. Its clean replay matches all 7,567 numeric values exactly with no exclusions. The five NBC continuation roots also pass an independent complete gate review.

The 26-page report compiles without unresolved citations/references, overfull boxes or TeX warnings. All pages and all five figures were visually reviewed; all 17 required scientific topics were independently checked. Caption, heading and rounding issues were corrected. Two separate builds produce the identical final PDF, and its 26 page texts and rendered pixels match the reviewed version exactly after metadata controls were added. `verify_report.py` also verifies that report preparation did not change the 47 numerical payloads or numerical analysis scripts.

## Parameters and experiments needing constraint

The largest local sensitivities prioritise functional pump capacity and apical distribution; effective NKCC capacity and stimulation gain; NBC capacity, saturation width and its assumed loading allocation; NHE abundance and stimulation; the finite cell buffer pool; and CaCC conductance/recruitment. No defensible experimental uncertainty intervals are available for these effective settings, so no finite range robustness claim is made.

Most urgent are matched stimulated CaCC currents and surface recruitment during acute AE4 loss; simultaneous sodium/bicarbonate flux, current and pH measurements for alkalinity support; ion conditioned NKCC/NHE and pump measurements; and TIC, buffering, pH, chloride and volume time courses. Source stoichiometry requires combined cation, chloride, carbon and current/reversal measurements. Within the prescribed inverse family the null retains about 10.49% of stimulated CaCC recruitment; this is a directly testable family-specific requirement.

## Deliverables and preservation

- `analysis/44_full_system_mathematics/report.tex`
- `analysis/44_full_system_mathematics/AE4_full_system_mathematics.pdf`
- `analysis/44_full_system_mathematics/output/`: source audit, scales, root/spectrum/IFT/compensation/source-vector/sensitivity results, figures, clean replay and report verification.
- `analysis/44_full_system_mathematics/inverse_detail/`: ordered scan, both refined crossings, original-point reproduction, long diagnostics, all parameter sensitivities and independent verification. Optional trial caches are excluded; authoritative results are retained.
- `README.md`, `RECOVERY.md`, analysis/verification scripts, TeX fragments, bibliography, deterministic `build.sh` and normalised `build.log`.
- `output/final_repository_audit.json`: exact reviewed refs, permitted path differences, preservation and branch separation checks.

The final report checkpoint and artifact hashes identify the reviewed content. The final PDF SHA256 is `3a86965653ff2cc4834ee99772f48edc4dbbba8b6b6b0960ff3110c91b63c60d`; remaining report/figure hashes are in `output/report_verification.json`.

`main`, the manuscript, Tasks 1–42, production model source and frozen outputs were untouched. No model fit was performed. The old Na/K/Cl-only reduction was not reinstated. Neither analysis branch was merged, and no force push was used.
