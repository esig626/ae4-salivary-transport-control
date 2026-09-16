"""49E: reveal only after remote 49D; no integration, inference or retuning."""
import json
from checkpoint import HERE, ROOT, OUT, dump, sha, verified

dependency = verified('49D')
frozen = json.loads((OUT/'pre_AE4_prediction.json').read_text())
assert frozen['parameter_freeze']['calibrated_conductance_S'] is None
ledger_path = ROOT/'analysis/48_joint_experimental_constraint_reconstruction/constraints.json'
ledger = json.loads(ledger_path.read_text())
row = next(r for r in ledger['ae4_table1'] if r['observable']=='initial_Cl_uptake_IPR')
mean, se, n = row['AE4_KO']
result = dict(
    checkpoint='49E', verified_49D_sha=dependency['sha'],
    frozen_prediction_sha256=sha(OUT/'pre_AE4_prediction.json'),
    heldout_ledger_sha256=sha(ledger_path),
    status='PREDECLARED_SWELLING_GATE_FAILS_KO_IPR_INITIATION; FULL_NUMERICAL_RECONSTRUCTION_UNAVAILABLE',
    calibration_after_reveal=False, new_trajectories=0, parameter_updates=0,
    exact_prediction_comparison=dict(
        observable='AE4_KO initial_Cl_uptake_IPR', units=row['units'],
        predicted_mean=0.0, observed_mean=mean, observed_SE=se, n=n,
        model_minus_data_in_reported_SE=-mean/se,
        statistical_limit='Descriptive discrepancy in reported SE units, not a calibrated p-value or joint likelihood.',
        observation_map='A constant chemical state gives zero SPQ slope for any fixed time-independent optical calibration.',
        scope='Every finite nonnegative conductance in the single frozen positive-swelling law at exact own-genotype rest.'),
    questions=[
        dict(id=1, question='KO CCh versus combined invariance removed?', status='not quantitatively established',
             reason='The extra current can depend on beta at swollen states, but no conductance was identified. KO IPR-only still has an exact no-initiation invariant. This does not prove combined-stimulus equality for all states.'),
        dict(id=2, question='CCh preserved with selective beta impairment?', status='CCh exactly nests; beta prediction unavailable'),
        dict(id=3, question='KO IPR uptake in correct direction?', status='fails: exact zero instead of positive 0.20 +/- 0.03'),
        dict(id=4, question='Combined early preservation and later loss?', status='unavailable'),
        dict(id=5, question='Ten-minute cumulative deficit?', status='unavailable; 35 percent not fitted'),
        dict(id=6, question='AE2 secretion neutrality?', status='unavailable for beta-containing protocols'),
        dict(id=7, question='Task41 immediate-deficit and AE4-TMEM16A coupling avoided?', status='no direct AE4-expression factor; timing unavailable'),
        dict(id=8, question='Resting mismatch?', status='unchanged: beta=0 exactly nests parent; StageII follows remote49E')],
    interpretation='Reject this sole predeclared gate as a sufficient reconstruction, not the experimental existence of a cAMP/VRAC pathway.',
    deferred_numeric_observations='All other beta-containing AE4/AE2 endpoints remain unavailable, including combined uptake and secretion. No old-parent output is relabelled as a new VRAC prediction.')
dump(OUT/'heldout_comparison.json', result)
print('49E: exact KO IPR slope 0 versus 0.20 +/- 0.03; no fit or new trajectory.')
