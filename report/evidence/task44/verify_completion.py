"""Final cross-file checks for Task 45, in addition to independent recomputation."""
from core import *
import subprocess
P=Path(__file__).resolve().parent
O=P/'output'
def read(path): return json.loads((P/path).read_text())
def main():
    reference()
    eq=read('output/equilibrium_review.json')
    assert len(eq['equilibrium_rows'])==26
    assert sum(r['full_inherited_physiological'] for r in eq['equilibrium_rows'])==22
    assert all(r['maximum_conservation_tolerance_ratio']<1 for r in eq['equilibrium_rows'])
    assert all(r['matched_spectrum_max_relative_error']<1e-5 for r in eq['equilibrium_rows'])
    assert all(r['max_real_eigenvalue_per_s']<0 for r in eq['equilibrium_rows'])
    assert len(eq['production_reference'])==3
    assert all(r['cumulative_absolute_error']<1e-8 for r in eq['production_reference'])
    ex=read('output/exact_verification.json');assert ex['status']=='passed' and ex['states']==84
    sens=read('output/sensitivity_extended.json');v=read('output/sensitivity_verification.json')
    assert len(sens['parameters'])==16 and v['status']=='passed'
    assert v['max_IFT_relative_error']<1e-4 and v['max_eigenvalue_method_relative_difference']<1e-5
    assert v['independent_nearby_WT_null_roots']==128
    base=sens['baseline'];assert .009<base['null_stationary_deficit']<.01
    assert .90<base['NBC_alkalinity_share']<1 and base['noNBC_capacity_gap']>0
    threshold=read('inverse_detail/threshold_sensitivities.json')
    assert len(threshold)==16
    assert all(r['maximum_validation_relative_error']<2e-3 for r in threshold)
    assert all(all(x['perturbed_gates_pass'] for x in r['estimates']) for r in threshold)
    inverse=read('output/inverse_threshold.json')
    assert inverse['tight_check']['physiological'] and abs(inverse['tight_check']['adaptive_deficit']-.303)<1e-6
    long=read('output/inverse_long_time_summary.json')
    assert 600<long['ph_gate_crossing_s']<3600 and long['end_ph']>7.3
    scale=read('output/rescaling_verification.json')
    assert scale['dimensionless_max_scaled_error']<1e-7
    original=read('output/test_summary.json');assert original['success']
    changed=subprocess.check_output(['git','diff','--name-only','b7b775a0277d93263c4de339745806e5658ec9c2','--'],cwd=ROOT,text=True).splitlines()
    assert all(x.startswith('analysis/44_full_system_mathematics/') for x in changed),changed
    required=['core.py','run_analysis.py','verify.py','verify_exact.py','exact_results.tex','equilibrium_review.py','equilibrium_results.tex','sensitivity_analysis.py','sensitivity_results.tex','inverse_analysis.py','inverse_verify.py','inverse_results.tex','literature_comparison.tex','sources.bib']
    assert all((P/f).is_file() for f in required)
    # Verify both the amount scaling and the augmented water observable explicitly.
    _,stim,y0,_,_=reference();m=clone(stim);sc=Scaling(m,y0)
    dim=integrate(m,y0);nd=integrate(m,y0,scaled=sc)
    water_error=abs(float(dim.y[-1,-1])-float(nd.y[-1,-1])*sc.V0)
    assert water_error<1e-8
    summary=dict(status='passed',analysed_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
      equilibrium_roots=26,physiological_roots=22,random_valid_exact_states=84,uncertain_equilibrium_parameters=16,
      independently_solved_nearby_roots=128,inverse_parameter_observable_pairs=16,
      dimensional_recovered_cumulative_water_error_pL=water_error,
      protected_path_differences=[],required_analysis_files=required,
      scope='Cross-file and fresh dimensional recovery checks. Full clean-checkout numerical replay is recorded separately.')
    write_json(O/'completion_verification.json',summary)
    print('COMPLETION_CHECKS passed; all full-system result families present.')
if __name__=='__main__':main()
