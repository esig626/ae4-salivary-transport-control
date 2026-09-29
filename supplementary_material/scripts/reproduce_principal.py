"""Rerun the three archived principal cases without changing their equations."""
from pathlib import Path
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'model'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=('case_01', 'case_02', 'case_03', 'all'), default='case_02')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    protected = (ROOT / 'data', ROOT / 'model', ROOT / 'parameters', ROOT / 'provenance')
    if any(output == p or p in output.parents for p in protected):
        parser.error('Choose a new output directory outside the archived inputs.')
    if output.exists():
        parser.error('The output directory must not already exist.')

    # The original runner sets numerical thread limits before importing NumPy.
    import run_frozen_cases as runner
    from run_reporter_recovery_52D import encode_observations, patched_attempt_source
    from threadpoolctl import threadpool_limits
    from verify_files import verify
    verify()
    freeze = json.loads((ROOT / 'parameters/parameter_and_case_freeze.json').read_text())
    cases = [c for c in freeze['case_matrix']['cases']
             if c['id'] in ('case_01', 'case_02', 'case_03')
             and (args.case == 'all' or c['id'] == args.case)]
    onsets = freeze['projected_onsets']['central']
    y0 = runner.np.concatenate([onsets[g]['state'] for g in ('WT', 'KO')])
    settings = freeze['case_matrix']['solver']

    # Reuse the archived numerical routine and its recorded serializer repair.
    # Only its one-time branch/HEAD launcher is replaced by this portable entry.
    _, source = patched_attempt_source()
    namespace = dict(runner.__dict__, encode_observations=encode_observations)
    exec(compile(source, 'archived_run_attempt', 'exec'), namespace)
    run_attempt = namespace['run_attempt']
    output.mkdir(parents=True)
    receipts = []
    with threadpool_limits(limits=1):
        for case in cases:
            directory = output / case['id']
            directory.mkdir()
            model = runner.PairedModel(protocol=case['protocol'], G_aux_S=case['G_aux_S'],
                                       ko_supply_multiplier=case['ko_supply_multiplier'])
            if not runner.dependency_audit()['pass']:
                raise RuntimeError('A model dependency was loaded from outside the supplied snapshot.')
            result = run_attempt(case, model, y0.copy(), 'Radau', directory, settings)
            if result['status'] == 'solver_status_failure':
                result = run_attempt(case, model, y0.copy(), 'BDF', directory, settings)
            receipts.append({'case': case['id'], 'status': result['status'],
                             'headline_0_600': result['headline_0_600']})
            print(case['id'], result['status'])
    (output / 'reproduction_summary.json').write_text(json.dumps(receipts, indent=2, allow_nan=False) + '\n')
    if any(row['status'] != 'completed' for row in receipts):
        raise SystemExit('At least one case did not complete its original numerical or physiological checks.')


if __name__ == '__main__':
    main()
