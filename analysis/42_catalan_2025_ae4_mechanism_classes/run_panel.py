"""One worker; each invocation stops at a required publication boundary."""
from validation_common import *
import argparse
import subprocess
from datetime import datetime, timezone


def require_publication(stage):
    records = json.loads((HERE/'publication_history.json').read_text())['checkpoints']
    assert stage in records, f'Publish {stage} before the next execution stage'
    assert len(records[stage]['sha']) == 40
    return records[stage]['sha']


def invoke(name, script, args, log_name):
    process = subprocess.run([sys.executable, str(HERE/script), *args],
        env={**os.environ, 'AE4_TASK42_CLASS': name}, capture_output=True, text=True)
    (RESULTS/name/log_name).write_text(process.stdout+process.stderr)
    if process.returncode:
        raise RuntimeError(f'{name}/{script}: inspect the saved execution log')


def run_rests():
    source_sha = require_publication('source_checks')
    lock = RESULTS/'panel_started.json'
    assert not lock.exists(), 'Stationary panel already started; inspect its ledger'
    preflight = json.loads((RESULTS/'preflight.json').read_text())
    assert preflight['status'] == 'PREFLIGHT_PASS'
    for entry in preflight['execution_files']:
        assert hashlib.sha256((ROOT/entry['path']).read_bytes()).hexdigest() == entry['sha256']
    write_json(lock, {'started_utc': datetime.now(timezone.utc).isoformat(),
        'preflight_sha256': hashlib.sha256((RESULTS/'preflight.json').read_bytes()).hexdigest(),
        'source_checks_remote_sha': source_sha, 'class_order': list(CLASSES), 'workers': 1})
    for name in CLASSES:
        invoke(name, 'run_rest.py', [], 'rest_execution.log')
        record = json.loads((RESULTS/name/'wt_rest.json').read_text())
        print(f"{name}: WT rest {record['status']}", flush=True)


def run_class(name):
    require_publication('wt_rests')
    for preceding in list(CLASSES)[:list(CLASSES).index(name)]:
        require_publication(preceding)
    record = json.loads((RESULTS/name/'wt_rest.json').read_text())
    if not record['admissible']:
        print(f'{name}: WT rest failed; no genotype trajectories permitted', flush=True)
        return
    for case in ('wt', 'ae4_5pct', 'ae4_null'):
        invoke(name, 'run_case.py', [case], case+'_console.log')
        record = json.loads((RESULTS/name/(case+'_verification.json')).read_text())
        print(f"{name}/{case}: {record['status']}", flush=True)


def complete():
    for name in CLASSES:
        require_publication(name)
    assert not (RESULTS/'panel_complete.json').exists(), 'Panel already locked complete'
    budgets = [json.loads((RESULTS/name/'budget.json').read_text()) for name in CLASSES]
    expected_integrations = 3*sum(b.get('wt_rest_pass', False) for b in budgets)
    assert all(b['no_further_scientific_execution'] for b in budgets)
    summary = {'status': 'ALL_PREDECLARED_CLASSES_RECORDED',
        'completed_utc': datetime.now(timezone.utc).isoformat(),
        **{key: sum(b[key] for b in budgets) for key in (
            'stationary_solves', 'integration_attempts', 'numerical_retries',
            'stationary_numerical_retries', 'numerical_execution_seconds',
            'optimisation_calls', 'parameter_sweeps', 'genotype_rest_solves')},
        'intended_production_integrations': expected_integrations,
        'workers': 1, 'blas_threads': 1, 'classes': budgets,
        'interpretation_unlocked': False,
        'next_action': 'Publish the complete prediction checkpoint and stop before phenotype reveal'}
    assert summary['stationary_solves'] == 7
    assert expected_integrations <= summary['integration_attempts'] <= expected_integrations+7
    write_json(RESULTS/'panel_complete.json', summary)
    print('All seven class records complete; phenotype interpretation remains locked.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=('rests', 'class', 'complete'))
    parser.add_argument('class_id', nargs='?', choices=tuple(CLASSES))
    args = parser.parse_args()
    if args.stage == 'rests': run_rests()
    elif args.stage == 'class':
        assert args.class_id is not None
        run_class(args.class_id)
    else: complete()
