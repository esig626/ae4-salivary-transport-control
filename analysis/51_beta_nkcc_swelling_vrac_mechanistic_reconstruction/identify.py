"""One prospective scalar solve against independent WT swelling only."""
import json
from scipy.optimize import brentq
from checkpoint import HERE,OUT,verified,dump,sha,integrity
from trajectory import run_case

def main():
    verified('51B');integrity()
    target=OUT/'parameter_freeze.json'
    if target.exists():raise RuntimeError('Identification already frozen')
    attempts=[];cache={}
    def residual(g):
        if g in cache:return cache[g]['residual']
        case=run_case(f'calibration_{len(attempts):02d}',track='track2',genotype_name='WT',arm='IPR_ONLY',
            beta_multiplier=2.5,conductance_S=float(g))
        entry=dict(g_V_S=float(g),case=case['name'],status=case['status'],
            physiologically_admissible=case.get('physiologically_admissible'),
            residual=case.get('swelling_mean_300_600',0)-.125 if case['status']=='COMPLETE' else None)
        attempts.append(entry);cache[g]=entry;dump(OUT/'identification_attempts.json',attempts)
        if case['status']!='COMPLETE':raise RuntimeError('Numerical calibration failure, preserved case')
        return entry['residual']
    freeze=dict(route='PREDECLARED_SWELLING_FALLBACK',M_beta=2.5,M_beta_sensitivity=3.,
        g_V_S=None,parameter_status='UNIDENTIFIED',target_relative_swelling=.125,
        observation_window_s=[300,600],observation='Mean relative cell volume minus one',
        source_mapping_status='ASSUMPTION/SENSITIVITY; source narrative magnitude, not measured averaging window',
        cross_species_status='Rat submandibular central, rat parotid sensitivity; not mouse confidence bounds',
        current_mapping='UNRESOLVED; SI inaccessible and apical/current/volume mapping unidentified',
        model_sha256=sha(HERE/'task51_model.py'),runner_sha256=sha(HERE/'trajectory.py'),
        identification_code_sha256=sha(HERE/'identify.py'),heldout_data_used=False,
        predeclared_scalar_bracket_S=[0.,1e-6],xtol_S=1e-13,maxiter=40)
    try:
        low=residual(0.);high=residual(1e-6)
        if low*high>0:
            freeze.update(outcome='FALLBACK_NOT_BRACKETED',endpoint_residuals=[low,high],
                qualification='No accepted scalar conductance. Two endpoints do not prove absence of interior or other roots; no search or target rescue is authorised.')
        else:
            root,info=brentq(residual,0.,1e-6,xtol=1e-13,rtol=1e-8,maxiter=40,full_output=True)
            selected=cache[root]
            freeze.update(scalar_iterations=info.iterations,scalar_converged=bool(info.converged),selected_case=selected['case'])
            if info.converged and selected['physiologically_admissible'] and root>0:
                freeze.update(g_V_S=float(root),parameter_status='INDEPENDENT_SWELLING_CALIBRATED_CONDITIONAL',outcome='IDENTIFIED_CONDITIONALLY')
                run_case('validation_tmem16a_off',track='track2',genotype_name='WT',arm='IPR_ONLY',
                    beta_multiplier=2.5,conductance_S=float(root),tmem16a_off=True,purpose='INDEPENDENT_BETA_VALIDATION',parameter_status=freeze['parameter_status'])
                freeze['validation_cases']={'WT_IPR':selected['case'],'VRAC_off':'calibration_00','TMEM16A_off':'validation_tmem16a_off'}
            else:freeze.update(outcome='FALLBACK_INADMISSIBLE',rejected_root_S=float(root))
    except Exception as exc:
        freeze.update(outcome='FALLBACK_EXECUTION_FAILURE',failure=str(exc))
    freeze['attempts']=attempts
    dump(target,freeze)
    print(json.dumps(freeze),flush=True)

if __name__=='__main__':main()
