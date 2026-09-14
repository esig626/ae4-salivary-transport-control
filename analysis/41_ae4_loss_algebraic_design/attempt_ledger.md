Task 41 records every attempted design, including failures.

The latest user request authorized target-directed algebra and amendment as
new work after Task 40. The original no-new-mechanism/held-out rules still
describe Task 40; its files and results were not changed. The Task 41 scope
and each candidate's intent were written before the corresponding numerical
work. Earlier self-imposed candidate stop ledgers remain locked; subsequent
candidates are separate attempts under the new user task.

| Candidate | Algebraic design | Result and disposition | Production |
|---|---|---|---|
| 01 | Shared finite NKCC1 ceiling; initial guide 27.5% sustained loss | C=0.11825333 fmol/s converged but limiting null Cl=27.3225 mM and pH=7.34098 fail. One amended WT rest passed. | WT and 5% complete; 5% cumulative loss 13.8617%. Null stops at 577.752874 s when Cl=29.997678 mM. |
| 01 supporting solves | Null stationary pH=7.29; optional AE4 gain | Initial root, root retry, and least-squares attempts did not solve the constrained equations. Exact water elimination gave C=0.17100089 and only 8.5915% sustained loss. Two AE4-gain inverse attempts were not accepted. | None beyond the selected cap cohort. |
| 02 | Shared ceiling + CaCC conductance, 33% sustained loss and null pH=7.28 | Initial root stepped outside the positive-volume manifold; counter not saved. Damped retry failed to converge (scaled residual 0.07455). | None. |
| 03 | Shared ceiling + coordinated AE4/NBC gain; null Cl=40 mM and 33% loss | Did not converge (residual 0.08134). | None. |
| 03 variant | First solve null Cl=34 mM for C, then WT AE4-only gain for 35% loss | Null root converged at C=0.13880567, pH=7.32360. WT inverse did not converge (residual 0.36026). | None. |
| 04 | Shared ceiling + pump Na half-saturation; pH=7.28 and 33% loss | Converged at C=0.08759126, KNa=16.50607 mM, but null Cl=20.27302 mM. | None. |
| 05 | Shared ceiling + pump Na half-saturation + CaCC factor; Cl=40, pH=7.28, 33% loss | Did not converge (residual 0.06702). | None. |
| 06 | pH-dependent CaCC conductance, fixed width 0.02 pH units; 30% loss | Converged at half-pH=7.27067134, but limiting null pH=7.319139. | None. |
| 07 | AE4-dependent CaCC recruitment; 30% sustained-loss guide | Two stationary roots converged; residual null recruitment b=0.08140757. Limiting null pH=7.319139 remains outside the box. | Three 600 s cases pass all gates. Cumulative losses 26.1848% / 35.6915%; null just above range. |
| 08 selected | Same recruitment law; one closed-form correction using two existing null trajectories | b=0.10511872843288, targeting approximately 30% cumulative loss. No new algebraic root or trajectory during parameter selection. | Three 600 s cases pass. Cumulative losses 23.1634% / 30.2612%. |

There were 17 algebraic equation-solver calls, including two
least-squares equation solves, one new WT resting solve (candidate 01 only),
and nine production integration attempts: eight complete and one scientific
stop. No genotype-specific resting solve, routing-fraction search, parameter
grid, trajectory optimizer, or parallel scientific worker was used. Parameter
selection and the final correction did use the target; this is calibration
of a hypothesis, not independent validation. There were 51151 recorded
algebraic residual evaluations, plus the unavailable count from candidate
02's initial exception. About 59.92 seconds of numerical
execution are recorded; that exception's timing was not saved.

The early `capacity_design.json` stores the parent NKCC1 field inside its
`observables` block. Its explicit `limited_nkcc1_fmol_s` and substituted
`raw_rhs` are the candidate values; production was tested against those.
Do not interpret parent-field values or any nonconverged iterate as an
accepted amended-model result. Candidate-03 AE4-only rows explicitly label
their observables as being before the source substitution.

All original execution inputs and candidate 01/07/08 frozen hashes were
verified again during final reporting. Empty initial candidate-02 stdout is
retained; its exception is recorded in the candidate's design result.
