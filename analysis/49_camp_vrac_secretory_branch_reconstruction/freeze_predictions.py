"""49D prospective record after the independent identification gate fails."""
import json,hashlib
from pathlib import Path
from checkpoint import HERE,ROOT,OUT,dump,verified

def main():
    dep=verified('49C')
    freeze=json.loads((OUT/'parameter_freeze.json').read_text())
    assert freeze['calibrated_conductance_S'] is None
    assert hashlib.sha256((HERE/'vrac_model.py').read_bytes()).hexdigest()==freeze['model_source_sha256']
    p48=ROOT/'analysis/48_joint_experimental_constraint_reconstruction/output'
    predictions=[]
    for group,genotype in [('AE4_WT','WT'),('AE4_KO','AE4_KO'),('AE2_control','WT'),('AE2_KO','AE2_KO')]:
        rest=p48/'rests'/f'{genotype}.json'
        cch=p48/'protocol_predictions'/f'{genotype}__CCH_ONLY.json'
        for arm in ['REST','CCH_ONLY','IPR_ONLY','CCH_IPR']:
            row=dict(cohort=group,genotype=genotype,protocol=arm,calibrated_numerical_prediction=None)
            if arm in ['REST','CCH_ONLY']:
                path=rest if arm=='REST' else cch
                row.update(parameter_independent_status='exact nesting; reuse published parent',
                    inherited_artifact=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            elif genotype=='AE4_KO' and arm=='IPR_ONLY':
                row.update(parameter_independent_status='exact-rest invariant for every finite nonnegative conductance',
                    chemical_state='constant exact genotype equilibrium',vrac_current_A=0,
                    chloride_concentration_slope_mM_s=0,SPQ_slope_per_s=0,
                    water_flow='constant resting outflow; no beta-induced increment',
                    condition='exact own-genotype equilibrium; common time-independent optical map')
            else:row.update(parameter_independent_status='unavailable without independently identified conductance',
                reason='49C identification gate failed; no arbitrary scale promoted')
            predictions.append(row)
    dump(OUT/'pre_AE4_prediction.json',dict(checkpoint='49D',verified_49C_sha=dep['sha'],
        calibrated_prediction_generated=False,status='PRE_AE4_NUMERICAL_PREDICTION_UNAVAILABLE',
        new_trajectories=0,heldout_evaluation_performed=False,parameter_freeze=freeze,
        family='one predeclared law, all finite g>=0; not a candidate-family search',
        predictions=predictions,proof='At exact KO rest, inherited core is beta-independent and positive-swelling gate is zero. The constant core with evolving disconnected AE4 activation solves the augmented ODE. Local uniqueness gives the invariant for every finite g.',
        numerical_rest_residual_rule='Numerical equilibrium error is not an authorised source of activating swelling.'))
    print('49D numerical unavailability and parameter-independent predictions frozen; no phenotype comparison')

if __name__=='__main__':main()
