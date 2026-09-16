"""49C: consume published WT data; no trajectory, optimiser or AE4 phenotype."""
import csv,json,hashlib
from pathlib import Path
from checkpoint import HERE,ROOT,OUT,dump,verified

def main():
    verified('49B')
    source=ROOT/'analysis/48_joint_experimental_constraint_reconstruction/output/protocol_predictions/WT__IPR_ONLY.csv'
    rows=list(csv.DictReader(source.open()))
    v=[float(r['volume_pL']) for r in rows]
    sw=[x/v[0]-1 for x in v]
    gate=[max(x,0.) for x in sw]
    # These are cached samples, not a continuous enclosure or a rerun.
    audit=dict(source=str(source.relative_to(ROOT)),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        samples=len(rows),time_s=[float(r['time_s']) for r in rows],
        volume_pL=v,fractional_swelling=sw,positive_swelling_gate=gate,
        initial_volume_pL=v[0],final_volume_pL=v[-1],maximum_saved_gate=max(gate),
        final_swelling_percent=100*sw[-1],new_trajectories=0,
        interpretation='No activation at published samples; no assertion about unsaved transient crossings.',
        source_narrative_swelling_percent=dict(mean=12.5,SEM=.2),
        source_curve_discrepancy='Fig5C digitised endpoint approximately16%; retained separately from narrative12.5%.',
        formal_assay_standardised_residual=None,
        observation_qualification='Model finite lumen and source isolated-acinar volume experiment are not silently equated.')
    dump(OUT/'independent_activation_audit.json',audit)
    frozen=dict(status='BLOCKED_INDEPENDENT_CALIBRATION',law='g_V beta max(V_i/V_rest_genotype - 1,0)',
        calibrated_conductance_S=None,finite_identified_envelope_S=None,
        diagnostic_fixture_conductance_is_parameter=False,production_eligible=False,
        new_parameters_inferred=0,optimisation_calls=0,production_trajectories=0,
        beta_NKCC_changed=False,AE4_heldout_data_used=False,post_reveal_retuning=False,
        blockers=[
          'No verified matched blocked apical current, swelling and chloride driving-force map.',
          'Perforated whole-cell current is not a measured apical conductance; total and blocked components differ.',
          'WT current-parent IPR samples shrink instead of activating the predeclared swelling gate.',
          'Independent gland shapes/ratios have not supplied a demonstrated finite identifiable scale under these activation/observation restrictions.'
        ],
        limits='This is an unresolved/failed identification gate, not proof that the unavailable supplement lacks data or that every possible observation map is impossible.',
        source_manifest_sha256=hashlib.sha256((HERE/'input/pnas2015/repository_source_manifest.json').read_bytes()).hexdigest(),
        frozen_inherited_parameter_sha256=hashlib.sha256((ROOT/'analysis/47_dynamic_potassium_recycling_reconstruction/input/active_parameters.json').read_bytes()).hexdigest(),
        model_source_sha256=hashlib.sha256((HERE/'vrac_model.py').read_bytes()).hexdigest(),
        next='Publish immutable unavailable-prediction record at49D, then49E limitations, thenStageII after remote49E.')
    dump(OUT/'parameter_freeze.json',frozen)
    dump(OUT/'independent_constraint_status.json',[
      dict(observation='TMEM16A-independent IPR output',implementation='separate current verified',numerical_prediction=None,reason='No calibrated conductance; ideal acute deletion differs from chronic source knockout'),
      dict(observation='DCPIB/NPPB inhibition',implementation='zero-current mask nests parent',numerical_prediction=None,reason='No drug potency or calibrated conductance; source blocked and total currents retained separately'),
      dict(observation='IPR swelling and delay',implementation='existing volume dynamics only',numerical_prediction='cached parent sample audit',reason='No positive gate at published times; source swelling positive'),
      dict(observation='chloride/bicarbonate dependence',implementation='blocked assay domain',numerical_prediction=None,reason='Inherited substrate/bath/carbonate restrictions'),
      dict(observation='CCh-only unchanged branch',implementation='exact RHS nesting tested',numerical_prediction='inherited published CCh trajectory, by exact nesting',reason='beta=0 for every CCh-only state')])
    print(json.dumps(dict(status=frozen['status'],maximum_saved_gate=max(gate),final_swelling_percent=100*sw[-1],calibrated_conductance_S=None)))

if __name__=='__main__':main()
