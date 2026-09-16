"""Post-51D file-only comparison. No model, fit or trajectory imports."""
import json
from checkpoint import HERE, ROOT, OUT, verified, integrity, dump, sha, git


def main():
    boundary = verified('51D')
    integrity()
    destination = OUT / 'heldout_comparison.json'
    if destination.exists():
        raise RuntimeError('Held-out comparison already exists')
    receipt_path = str((OUT / 'checkpoints/51D.json').relative_to(ROOT))
    receipt = json.loads(git('show', boundary['sha'] + ':' + receipt_path))
    frozen = ['task51_model.py', 'trajectory.py', 'identify.py',
              'SOURCE_ARCHITECTURE_FREEZE.md', 'output/parameter_freeze.json',
              'output/pre_AE4_predictions.json']
    frozen_hashes = {}
    for relative in frozen:
        path = HERE / relative
        assert sha(path) == receipt['artifact_sha256'][str(path.relative_to(ROOT))]
        frozen_hashes[relative] = sha(path)
    prediction = json.loads((OUT / 'pre_AE4_predictions.json').read_text())
    assert prediction['heldout_comparison_performed'] is False
    assert json.loads((OUT / 'parameter_freeze.json').read_text())['g_V_S'] is None
    benchmark_path = ROOT / 'analysis/50_minimal_beta_conditioned_effective_coupling/output/inherited_phenotype.json'
    benchmark = json.loads(benchmark_path.read_text())
    benchmark_rows = [dict(case=row['case'], cumulative_pL=row['cumulative_0_600_pL'],
                           deficit_percent=row['deficit_percent'],
                           endpoint_deficit_percent=row['endpoint_deficit_percent'],
                           task51_cumulative_pL=None, task51_deficit_percent=None,
                           task51_minus_benchmark_percent=None)
                      for row in benchmark['cases']]
    rests = []
    observations = {'WT': (50.10, 1.50, 6.91, .07), 'AE4_KO': (36.50, 1.60, 6.89, .02)}
    for genotype, (cl, cl_sem, ph, ph_sem) in observations.items():
        path = ROOT / 'analysis/48_joint_experimental_constraint_reconstruction/output/rests' / (genotype + '.json')
        cached = json.loads(path.read_text())['observables']
        rests.append(dict(genotype=genotype, source='Unchanged Task48 cache by exact beta-zero nesting',
                          file=str(path.relative_to(ROOT)), sha256=sha(path),
                          model_cl_mM=cached['cl_i_mM'], observed_cl_mM=cl, observed_cl_sem_mM=cl_sem,
                          cl_difference_mM=cached['cl_i_mM']-cl,
                          model_ph=cached['ph_i'], observed_ph=ph, observed_ph_sem=ph_sem,
                          ph_difference=cached['ph_i']-ph))
    data = dict(
        classification='UNRESOLVED_NO_INDEPENDENTLY_IDENTIFIED_FULL_MODEL',
        immutable_reveal_boundary=boundary,
        post_reveal_frozen_hashes=frozen_hashes,
        heldout_source=dict(doi='10.1074/jbc.M114.612895', pmcid='PMC4409235',
                            provenance='Task48 source-checked ledger sections1.1-1.3; no new extraction'),
        secretion=dict(observed_AE4_null_deficit_percent=35., observed_error_percent=4.7,
                       task51_deficit_percent=None, task51_minus_experiment_percent=None,
                       persistence='UNAVAILABLE_NO_IDENTIFIED_G'),
        uptake=dict(units='1e-3 s^-1, source SPQ recovery metric',
                    observed_AE4_KO_CCH=dict(mean=2.30, sem=.10),
                    observed_AE4_KO_CCH_IPR=dict(mean=.90, sem=.09),
                    observed_AE4_KO_IPR=dict(mean=.20, sem=.03),
                    task51_SPQ_CCH=None, task51_SPQ_CCH_IPR=None, task51_SPQ_IPR=None,
                    initial_loading_sign='POSITIVE_FORMAL_LOCAL_PREDICTION',
                    qualification='Initial chemical loading is not the fluorescence recovery slope; no numerical agreement claimed'),
        task50=dict(classification=benchmark['classification'], file=str(benchmark_path.relative_to(ROOT)),
                    sha256=sha(benchmark_path), cases=benchmark_rows,
                    archive_branch='archive/task-50-working-effective-coupling-benchmark',
                    pinned_sha='b16c30094b95f61a73d8f1977cd58e79c7bb50f6',
                    imported_model=False),
        AE2_specificity=dict(observed='Approximately secretion-neutral in source whole-gland experiment',
                             model='UNAVAILABLE_FOR_BETA_PROTOCOLS; CCh-only exact parent IVP retained',
                             control='Shared WT model alias; experimental strain means not fitted'),
        compensation=dict(excess_NKCC_compensation=None, beta_demand_exposes_AE4_dependence=None,
                          explanation='No accepted matched WT/AE4 beta trajectories'),
        physical_constraints=dict(calibration_trials='Both passed inherited sampled physical/conservation checks',
                                  accepted_model='Unavailable', independent_swelling='Bracket endpoints below target',
                                  chronic_rest='Inherited experimental Cl/pH mismatch persists'),
        unchanged_chronic_rests=rests,
        conclusion='Local beta-blind initiation deadlock removed conditionally; held-out quantitative mechanism remains unpredicted',
        new_model_evaluations=0, new_trajectories=0, new_fits=0, new_rest_solves=0, retuned=False)
    dump(destination, data)
    print(json.dumps(dict(boundary=boundary['sha'], classification=data['classification'],
                          frozen_hashes_unchanged=len(frozen_hashes), new_model_evaluations=0)))


if __name__ == '__main__':
    main()
