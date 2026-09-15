"""Import exact committed companion results; perform no new sensitivity runs."""
from pathlib import Path
import argparse
import csv
import hashlib
import io
import json
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE_COMMIT = 'f783785a863440df459e9c5530beba4c7eb59110'
PREFIX = 'analysis/44_full_system_mathematics/'
FILES = {
    'stationary': 'output/sensitivity_extended.json',
    'stationary_verification': 'output/sensitivity_verification.json',
    'inverse': 'inverse_detail/threshold_sensitivities.json',
    'inverse_verification': 'inverse_detail/independent_verification.json',
    'crossings': 'inverse_detail/crossings.json',
    'long_time': 'inverse_detail/long_time_summary.json',
    'nbc_continuation': 'output/nbc_continuation.csv',
}


def inventory_name(name):
    if name.startswith('parameters.'):
        return name[len('parameters.'):]
    return {
        'nkcc1.alpha_eff_fmol_s': 'nkcc.alpha_eff_fmol_s',
        'regulation.ae4_regulatory_model.tau_activation_s': 'regulation.tau_ae4_s',
        'regulation.ae4_regulatory_model.gain.fully_activated_increment': 'regulation.fully_activated_increment',
        'regulation.nkcc1_regulatory_model.fully_activated_multiplier': 'nkcc_activation.fully_activated_multiplier',
    }.get(name, name)


def products(source_commit=SOURCE_COMMIT):
    commit = subprocess.check_output(
        ['git', 'rev-parse', '--verify', source_commit + '^{commit}'], cwd=ROOT,
        text=True).strip()
    assert commit == source_commit, 'Use the full immutable source commit SHA.'
    raw = {key: subprocess.check_output(
        ['git', 'show', commit + ':' + PREFIX + path], cwd=ROOT)
        for key, path in FILES.items()}
    data = {key: json.loads(value) for key, value in raw.items()
            if key != 'nbc_continuation'}
    assert data['stationary_verification']['status'] == 'passed'
    assert data['inverse_verification']['status'] == 'PASS'
    inventory = {r['name']: r for r in json.loads(
        (HERE / 'output/parameter_inventory.json').read_text())}
    rows = []
    for p in data['stationary']['parameters']:
        name = inventory_name(p['inventory_path'])
        assert name in inventory and inventory[name]['value'] == p['value'], name
        rows.append({
            'inventory_parameter': name, 'parameter': p['parameter'],
            'active_value': p['value'], 'classification': inventory[name]['classification'],
            'uncertainty_status': inventory[name]['uncertainty_status'],
            'constraint_experiment': p['experiment'], 'results': p['results'],
            'inverse_family_sensitivities': [
                {k: v for k, v in q.items() if k != 'estimates'}
                for q in data['inverse'] if q['parameter'] == p['parameter']],
        })
    assert len(rows) == 16
    assert sum(len(r['inverse_family_sensitivities']) for r in rows) == len(data['inverse']) == 16
    nbc = list(csv.DictReader(io.StringIO(raw['nbc_continuation'].decode())))
    nbc_full = next(r for r in nbc if float(r['NBC_fraction']) == 1.)
    nbc_null = next(r for r in nbc if float(r['NBC_fraction']) == 0.)
    assert nbc_null['physiological'] == 'True'
    no_nbc = {
        'physiological': True,
        'intracellular_pH': float(nbc_null['ph_i']),
        'secretion_pL_s': float(nbc_null['q_pL_s']),
        'secretion_decrease_fraction': 1 - float(nbc_null['q_pL_s']) / float(nbc_full['q_pL_s']),
        'AE4_cycle_fmol_s': float(nbc_null['J4']),
        'full_NBC_AE4_cycle_fmol_s': float(nbc_full['J4']),
        'interpretation': 'A physiological no-NBC equilibrium survives at much lower AE4 loading; the capacity obstruction concerns a specified AE4 demand, not universal NBC necessity for secretion.',
    }
    result = {
        'source_branch': 'analysis/task-44-full-system-mathematics',
        'source_commit': commit,
        'imported_source_files': [
            {'path': PREFIX + path, 'sha256': hashlib.sha256(raw[key]).hexdigest()}
            for key, path in FILES.items()],
        'scope': 'Read-only import of committed Task 44 stationary and prescribed inverse-family results; no new sensitivity calculation on Task 43.',
        'baseline': data['stationary']['baseline'], 'parameters': rows,
        'limitations': data['stationary']['limitations'],
        'stationary_verification': data['stationary_verification'],
        'inverse_verification': data['inverse_verification'],
        'no_NBC_equilibrium': no_nbc,
        'inverse_crossings': {
            key: {k: value[k] for k in ['rho', 'b', 'ae4_dependent_recruitment_percent']}
            for key, value in data['crossings'].items()
            if isinstance(value, dict) and 'rho' in value},
        'uncertainty': 'Numerical difference steps are verification choices, not uncertainty intervals. No finite parameter change that overturns a claim is established by these local results.',
    }
    flat = [{k: r[k] for k in ['inventory_parameter', 'parameter', 'active_value', 'classification']} | r['results'] for r in rows]
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, list(dict.fromkeys(k for row in flat for k in row)), lineterminator='\n')
    writer.writeheader()
    writer.writerows(flat)
    return {
        'task44_sensitivity_crosswalk.json': json.dumps(result, indent=2, sort_keys=True) + '\n',
        'task44_sensitivity_crosswalk.csv': stream.getvalue(),
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-commit', default=SOURCE_COMMIT)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    for name, contents in products(args.source_commit).items():
        path = HERE / 'output' / name
        if args.check:
            assert path.read_text() == contents, 'Saved crosswalk differs from exact committed source: ' + name
        else:
            path.write_text(contents)
    print(('Verified' if args.check else 'Imported') + ' 16 stationary and 16 inverse parameter/observable records from ' + args.source_commit)
