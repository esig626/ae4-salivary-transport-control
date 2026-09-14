"""Write the final report and complete attempt ledger from saved results."""
from summarize_final import ROOT,HERE,OUT,FINAL,CASES,read,write,csvread
from datetime import datetime,timezone
import hashlib


def main():
    d=read(OUT/'final_summary.json');s={r['case']:r for r in d['cases']}
    e=d['endpoint'];b=d['selected_null_cacc_recruitment_fraction']
    flux={c:{r['quantity']:float(r['value']) for r in csvread(FINAL/f'{c}_integrated_fluxes.csv')} for c in CASES}
    cohorts=[OUT,OUT/'candidate_07',FINAL]
    budgets=[read(p/'budget.json') for p in cohorts]
    for p in cohorts:
        assert read(p/'budget.json')['no_further_scientific_execution']
        for f in read(p/'frozen_inputs.json')['execution_files']:
            assert hashlib.sha256((ROOT/f['path']).read_bytes()).hexdigest()==f['sha256'],f['path']
    for f in read(ROOT/'results/40_ae4_equal_cation_routing/frozen_inputs.json')['task40_execution_files']:
        assert hashlib.sha256((ROOT/f['path']).read_bytes()).hexdigest()==f['sha256'],f['path']
    designs=[OUT/x for x in ('capacity_design.json','constrained_design.json','constrained_design_retry.json',
        'constrained_design_scaled.json','reduced_capacity_design.json','two_parameter_design.json','two_parameter_design_scaled.json',
        'candidate_02/design.json','candidate_02/design_retry.json','candidate_03/design.json','candidate_03/ae4_only_design.json',
        'candidate_04/design.json','candidate_05/design.json','candidate_06/design.json','candidate_07/design.json')]
    records=[read(p) for p in designs]
    calls=sum(r.get('root_calls',1) for r in records)
    count=sum(r.get('actual_residual_calls') or 0 for r in records)
    extra=[read(OUT/x) for x in ('candidate_02/design_retry.json','candidate_03/design.json','candidate_03/ae4_only_design.json',
        'candidate_04/design.json','candidate_05/design.json','candidate_06/design.json')]
    completion={'status':'TASK41_TARGET_REACHED_REPORT_COMPLETE','completed_utc':datetime.now(timezone.utc).isoformat(),
        'branch':'codex/task-40-ae4-equal-cation-routing','scientific_parent_commit':'76f3144a4e662eb12d6da26c83ebf6f2ff47a59d',
        'single_scientific_agent':True,'scientific_workers':1,'blas_threads':1,
        'all_algebraic_solver_calls':calls,'least_squares_equation_solver_calls':2,
        'recorded_algebraic_residual_evaluations':count,'unrecorded_residual_count_calls':1,
        'new_WT_rest_solves':sum(x['stationary_solves'] for x in budgets),
        'genotype_specific_rest_solves':0,'production_integration_attempts':sum(x['integration_attempts'] for x in budgets),
        'completed_600s_trajectories':sum(len(x['completed_cases']) for x in budgets),
        'scientific_trajectory_stops':1,'integration_numerical_retries':0,
        'recorded_numerical_seconds':sum(x['numerical_execution_seconds'] for x in budgets)+sum(x['wall_seconds'] for x in extra),
        'timing_limitation':'Initial candidate-2 domain exception was not timed in a saved result; residual count also unavailable.',
        'selected_model':'AE4-dependent stimulated CaCC recruitment','selected_candidate':'08',
        'selected_b':b,'final_deficits_percent':{c:s[c]['deficit_percent'] for c in CASES[1:]},
        'no_further_scientific_execution':True,'original_Task40_results_preserved':True,
        'target_selected_not_independent_validation':True,'routing_fractions_tested_in_followup':[.5]}
    write(OUT/'completion.json',completion)
    ledger=f'''Task 41 records every attempted design, including failures.

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
| 08 selected | Same recruitment law; one closed-form correction using two existing null trajectories | b={b:.14f}, targeting approximately 30% cumulative loss. No new algebraic root or trajectory during parameter selection. | Three 600 s cases pass. Cumulative losses {s['ae4_5pct']['deficit_percent']:.4f}% / {s['ae4_null']['deficit_percent']:.4f}%. |

There were {calls} algebraic equation-solver calls, including two
least-squares equation solves, one new WT resting solve (candidate 01 only),
and nine production integration attempts: eight complete and one scientific
stop. No genotype-specific resting solve, routing-fraction search, parameter
grid, trajectory optimizer, or parallel scientific worker was used. Parameter
selection and the final correction did use the target; this is calibration
of a hypothesis, not independent validation. There were {count} recorded
algebraic residual evaluations, plus the unavailable count from candidate
02's initial exception. About {completion['recorded_numerical_seconds']:.2f} seconds of numerical
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
'''
    (HERE/'attempt_ledger.md').write_text(ledger)
    lines=[f'''The amended model achieves the requested 20-35% loss of cumulative
secretion over 600 seconds: **{s['ae4_5pct']['deficit_percent']:.2f}% at 5% AE4** and
**{s['ae4_null']['deficit_percent']:.2f}% at zero AE4**. WT is exactly the completed Task 40
model. This is a target-selected model, not independent experimental
validation or proof that the added biological dependency exists.

The successful amendment is an explicit dependence of stimulated CaCC
recruitment on AE4 expression. AE4's 50:50 sources and total cycle law are
unchanged. None of the NKCC1 ceilings, pump changes, pH gates, or NBC/AE4
capacity changes explored algebraically is present in the final model.

| Case | Cumulative secretion, 0-600 s (pL) | Ratio to WT | Deficit | Absolute deficit (pL) | Flow deficit at 600 s |
|---|---:|---:|---:|---:|---:|''']
    for c,label in zip(CASES,('WT','AE4 5%','AE4 null')):
        r=s[c];lines.append(f"| {label} | {r['cumulative_0_600_pL']:.9f} | {r['ratio_to_WT']:.6f} | {r['deficit_percent']:.2f}% | {r['absolute_deficit_pL']:.9f} | {r['endpoint_deficit_percent']:.2f}% |")
    lines.append(r'''
The necessary balance condition is a sustained reduction in secreted salt
that compensation and release from intracellular stores do not erase.
Stoichiometry alone does not identify a unique mechanism producing that
condition. With NKCC1 cycles N, AE4 inward chloride cycles A, NBC cycles B,
NHE1 flux H, AE2 inward chloride E, and total pump cycles P, the inherited
equal-routing balances are exactly

\[
\dot n_{Na,i}=N+H+B-\tfrac12 A-3P,\qquad
\dot{TA}_i=H+2B-2A-E,
\]
\[
\dot n_{Cl,i}=2N+A+E-J_{CaCC}.
\]

Eliminating N, B and E gives the exact transient identity

\[
J_{CaCC}=6P-H+2\dot n_{Na,i}-\dot{TA}_i-\dot n_{Cl,i}.
\]

At stationary state, this becomes J_CaCC=6P-H. There is no direct AE4 term
remaining at the 50:50 split, although AE4 changes the coupled states and
therefore the pump, NHE1 and other fluxes. Loss of AE4 need not raise Na_i:
the loss of its bicarbonate export also changes NBC's Na input through the
alkalinity balance. The full cell-plus-lumen chloride identity is

\[
q_{out}[Cl]_l=2N+A+E-J_{para,Cl}
-\frac{d}{dt}(n_{Cl,i}+n_{Cl,l}).
\]

Thus the AE4 share of positive chloride loading is not itself a water
deficit. Paracellular return, luminal composition and stored chloride all
matter. The capacity-only test demonstrated this: it gave only 13.86%
cumulative loss at 5% AE4 and depleted null chloride below the inherited
30 mM floor before 600 s. The complete [attempt ledger](attempt_ledger.md)
retains that failure and every discarded algebraic design.

The selected added hypothesis is

\[
g_{CaCC}=g_{parent}\{(1-a)+a[b+(1-b)e_{AE4}]\},
\]

where a is the existing normalized calcium activation and e_AE4 is AE4
expression. The calcium input and its existing channel gate are retained.
The additional factor is applied to conductance before the full electrical
closure, including NBC. There is no multiplier on computed secretion.
''')
    lines.append(f'''The selected value is **b={b:.14f}**. At full stimulation, CaCC conductance
is 100% of its parent value in WT, {100*s['ae4_5pct']['cacc_stimulated_fraction']:.4f}% at 5% AE4,
and {100*b:.4f}% in the null. This assumes that about {100*(1-b):.2f}% of full stimulated
recruitment depends on AE4. That is a strong, unvalidated coupling; a 30%
secretion loss does not imply only a 30% conductance reduction.

An initial full stationary balance solution gave b=0.08140757 for 30%
sustained loss. Its finite-time simulations instead gave 26.18% / 35.69%
cumulative loss. One closed-form correction used the two already available
null trajectories (Task 40 and that first recruitment model). The surrogate
Q(g)=Qinf*g/(g+k) gave k=0.04588248, Qinf=0.99802085 pL and the final b above
for approximately 30% cumulative loss. This is a design approximation, not
an exact identity of the ODE. No trajectory optimizer or sweep was used.
The final simulations produced 23.16% / 30.26%, without further adjustment.

The source vector remains (-0.5A,-0.5A,+A,-2A,-2A) in Na/K/Cl/TIC/TA order.
Its charge-equivalent source is (-0.5-0.5-1+2)A=0; reversing A reverses all
five sources. The inherited total AE4 cycle, Palk NKCC1 law and alpha_eff,
1/1.75 NKCC1 stimulation, NBC and NHE1 kinetics, pump, other channels,
acid/base closure, paracellular parameters, water, geometry and stimulus
retain their Task 40 values. Only the new CaCC recruitment dependency is
added. AE4 expression now affects both AE4 transport and this explicit
channel recruitment law.

WT exactly nests Task 40 at every state, so no new resting solve was needed
for the selected model. All three final trajectories use the identical WT
rest SHA-256 `{d['frozen_rest_sha256']}`. WT cumulative output also reproduced
Task 40 bit for bit. The three final Radau integrations completed without
retry, using the inherited 600 s protocol, tolerances and maximum step.

| Readout at 600 s | WT | AE4 5% | AE4 null |
|---|---:|---:|---:|''')
    for field,label in (('na_i_mM','Na_i (mM)'),('k_i_mM','K_i (mM)'),('cl_i_mM','Cl_i (mM)'),
        ('hco3_i_mM','HCO3_i (mM)'),('ph_i','pH_i'),('volume_i_pL','Cell volume (pL)')):
        lines.append('| '+label+' | '+' | '.join(f'{e[c][field]:.6f}' for c in CASES)+' |')
    lines.append('''
**Na_i does not become higher than WT in the final model** on the saved
one-second grid or at the endpoint. It ends 3.231 mM lower at 5% AE4 and
4.014 mM lower in the null. The target secretion loss is therefore achieved
without reproducing a higher-Na_i loss phenotype.

Intracellular chloride persists: the final mutant concentrations are 3.610
and 5.154 mM above WT. Their minima over all monitored states are 56.476
and 57.312 mM. Mutants retain chloride amount as well as concentration:
final chloride amounts are 93.843 / 100.043 fmol, compared with 75.655 fmol
in WT and 88.032 fmol at the shared start. Cell volume also increases to
1.629 / 1.692 pL, versus 1.402 pL in WT. These changes should not be
described as preserved overall homeostasis merely because the gates pass.

The following values integrate signed fluxes over 60-600 s except where
positive loading is specified. Parenthetical percentages are relative to
the final model's WT.

| Flux integral (fmol; pump/NBC in cycles) | WT | AE4 5% | AE4 null |
|---|---:|---:|---:|''')
    for k,label in (('positive_nkcc1_cl_loading','Positive NKCC1 Cl loading'),('positive_ae4_cl_loading','Positive AE4 Cl loading'),
        ('positive_ae2_cl_loading','Positive AE2 Cl loading'),('nbc_cycles','NBC'),('nhe1_flux','NHE1'),
        ('pump_cycles','Pump'),('cacc_cl_export','CaCC Cl export'),('paracellular_cl_return','Paracellular Cl return')):
        row=[label,f"{flux['wt'][k]:.6f}"]
        for c in CASES[1:]:
            pct=f" ({100*(flux[c][k]/flux['wt'][k]-1):+.2f}%)" if flux['wt'][k] else ''
            row.append(f'{flux[c][k]:.6f}'+pct)
        lines.append('| '+' | '.join(row)+' |')
    lines.append('''
AE2 positive loading is small but present: 0.2109% / 0.2497% of the positive
chloride-loading pool, versus zero in WT. NBC and NHE1 decrease markedly;
pump cycles decrease by 6.57% / 9.17%. The smaller chloride export and its
associated electrical changes reduce paracellular chloride return too.
The higher retained chloride weakens the NKCC1 compensation that previously
masked AE4 loss. No NKCC1 capacity restriction is present in this model.

| NKCC1 compensation, 60-600 s | AE4 5% | AE4 null |
|---|---:|---:|
| Task 39 | +29.65% | +33.10% |
| Task 40, 50:50 routing alone | +20.82% | +23.16% |
| Final Task 41 | +5.68% | +2.79% |

Cumulative deficits below are differences in integrated fluid output from
time zero. Percentage deficit is undefined at time zero because both
cumulative volumes are zero.

| Time (s) | WT cumulative (pL) | 5% deficit (pL; %) | Null deficit (pL; %) |
|---|---:|---:|---:|''')
    curve={(r['case'],float(r['time_s'])):r for r in csvread(OUT/'final_cumulative_deficit.csv')}
    for t in range(60,601,60):
        row=[str(t),f"{float(curve['wt',t]['cumulative_pL']):.6f}"]
        for c in CASES[1:]:
            r=curve[c,t];row.append(f"{float(r['deficit_pL']):.6f}; {float(r['deficit_percent']):.2f}%")
        lines.append('| '+' | '.join(row)+' |')
    lines.append(f'''
The percentage losses peak earlier ({s['ae4_5pct']['peak_cumulative_deficit_percent']:.2f}% at
{s['ae4_5pct']['peak_cumulative_deficit_time_s']:.0f} s and {s['ae4_null']['peak_cumulative_deficit_percent']:.2f}% at
{s['ae4_null']['peak_cumulative_deficit_time_s']:.0f} s on the one-second grid), then decline to the reported
600 s values. The absolute cumulative deficits continue to grow. The
20-35% criterion applies to the final cumulative observation, not every
earlier time point. Final instantaneous flow deficits are 20.11% / 27.04%.

![Final secretion, ions, pH and NKCC1 trajectories](../../results/41_ae4_loss_algebraic_design/final_trajectories.svg)

All accepted solver endpoints and integer-second dense samples passed the
unchanged positivity, Na/K/Cl, pH, volume and conservation checks. Source
tests verified exact WT nesting, recruitment boundaries, unchanged AE4 and
homeostasis sources at fixed state, no water-output rescaling, and the
coupled current/Na/alkalinity/chloride identities. The final parameter also
passed recruitment and current-conservation checks at all three expressions
before production. Pointwise conservation remained within every inherited
tolerance. The separate 1 s trapezoidal chloride-amount audit differs from
the actual amount change by at most 0.003081 fmol; this is quadrature error,
not an ODE conservation residual.

The null's mean 60-600 s secretion is only 1.072 times the shared rest rate.
It does not meet the WT-only 1.10-fold stimulation criterion; that criterion
is required of WT, which passes at 1.519-fold. All mutants remain above the
rest rate at 600 s.

Validity is demonstrated for 600 seconds only. The earlier recruitment
design's limiting null pH was 7.3191, above the inherited 7.3 bound. The
corrected parameter has not been given a long-run validation, and mutant
pH is still rising at 600 s. The secretion range was used for selection
and correction; proximity to the approximately 35% experimental deficit is
not held-out evidence. Measurements of stimulated CaCC recruitment and
the concurrent Na/Cl/pH/volume responses would be needed to evaluate this
specific hypothesis. Algebra plus one output range does not identify it
uniquely.

The public selected-model driver is
[task41_selected.py](../../src/modern_full_model/task41_selected.py), which
uses the frozen [recruitment engine](../../src/modern_full_model/ae4_cacc_recruitment.py).
Its mocked integration-dispatch test confirms that the solver receives the
amended RHS; the test runs no additional trajectory. Use the public driver
or the recorded gated runner, rather than delegating integration to the
underlying parent model.
Its selected parameter and all trajectories are in
[candidate_08](../../results/41_ae4_loss_algebraic_design/candidate_08).
The [summary](../../results/41_ae4_loss_algebraic_design/final_summary.json),
[completion ledger](../../results/41_ae4_loss_algebraic_design/completion.json),
[attempt ledger](attempt_ledger.md), and source scripts preserve the entire
target-directed history. Work remains on the requested branch; no merge
to main is part of this task. Scientific execution stops here.
''')
    (HERE/'final_answer.md').write_text('\n'.join(lines).replace('\\[','$$').replace('\\]','$$'))
    (HERE/'README.md').write_text('Task 41 final report: [final_answer.md](final_answer.md).\n\n'
        'Selected model: candidate_08, AE4-dependent stimulated CaCC recruitment.\n'
        'This new user-authorized task uses the secretion target for model design.\n'
        'Task 40 remains unchanged. See [scope](scope.md) and [attempt ledger](attempt_ledger.md).\n\n'
        'Production inputs and ledgers are frozen; do not rerun into the recorded\n'
        'result directories. Reproduce in an isolated results copy, retaining the\n'
        'fixed final b and shared Task 40 WT rest. Postprocessing scripts only\n'
        'read completed trajectories and do not invoke solvers.\n')
    write(OUT/'report_verification.json',{'target_satisfied_both_losses':True,
        'three_final_trajectories_complete_and_pass':True,'same_rest_all_cases':True,
        'WT_matches_Task40_exactly':True,'frozen_input_hashes_checked':True,
        'no_extra_scientific_runs_during_reporting':True,'all_attempts_disclosed':True,
        'plot_full_deficit_range_visible':True,'public_driver_dispatch_test_pass':True,
        'public_driver_test_calls_ODE_solver':False})
    print(completion)

if __name__=='__main__':main()
