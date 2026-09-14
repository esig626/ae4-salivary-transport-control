"""Record the finite panel and its algebra before any stationary solve."""
from validation_common import *
from datetime import datetime, timezone
import csv
from modern_full_model.ae4_catalan2025 import forward_affinity


def main():
    assert not (RESULTS / 'preflight.json').exists(), 'Preflight is already frozen'
    log = (HERE / 'focused_tests.log').read_text()
    assert 'Ran 7 tests' in log and log.rstrip().endswith('OK')
    rest, stim, y40, frozen = parent_models_and_state()
    e = rest.evaluate(0., y40)
    o, p = e.diagnostics.observables, rest.parameters
    inside = dict(o.cell_concentrations_mM)
    outside = {'na': p.bath.na_mM, 'k': p.bath.k_mM, 'cl': p.bath.cl_mM,
               'hco3': o.bath_acid_base.hco3_mM, 'co3': o.bath_acid_base.co3_mM}
    voltage = e.diagnostics.membranes.v_basolateral_V
    panel = []
    for name, c in CLASSES.items():
        assert c.exact_charge == 0
        a = forward_affinity(c, inside, outside, voltage_V=voltage,
                             thermal_voltage_V=p.constants.thermal_voltage_V)
        a.update(class_id=name, source_coefficients_na_k_cl_tic_ta=c.conserved_coefficients,
                 exact_charge_equivalents=str(c.exact_charge),
                 affinity_kJ_mol=a['affinity_over_rt']*p.constants.gas_constant_J_mol_K*p.constants.temperature_K/1000,
                 inherited_j4_fmol_s=e.diagnostics.ae4.cl_cell_fmol_s,
                 scalar_direction_agrees_at_reference=a['affinity_over_rt']*e.diagnostics.ae4.cl_cell_fmol_s >= 0)
        panel.append(a)
    execution_files = [ROOT / 'src/modern_full_model/ae4_catalan2025.py',
                       ROOT / 'tests/test_task42_catalan2025.py', *HERE.glob('*.py')]
    data = {'status': 'PREFLIGHT_PASS', 'frozen_utc': datetime.now(timezone.utc).isoformat(),
            'prepared_head': PREPARED_HEAD, 'branch': BRANCH,
            'task40_reference_sha256': frozen['state_sha256'],
            'reference_state_vector': y40.tolist(), 'reference_observables': asdict(o),
            'concentrations_i_mM': inside, 'concentrations_o_mM': outside,
            'v_basolateral_V': voltage, 'constants': asdict(p.constants),
            'acid_base_parameters': asdict(p.acid_base), 'classes': panel,
            'source_tests_pass': True, 'numerical_solves_so_far': 0,
            'execution_files': [{'path': str(f.relative_to(ROOT)),
                                 'sha256': hashlib.sha256(f.read_bytes()).hexdigest()} for f in sorted(execution_files)],
            'firewall_note': 'The required repository instructions and Task 40 report contain a phenotype reference. No held out experimental source or Task 41 content is inspected in this execution. No target value enters the implementation, gates, execution or reporting. All seven classes and all parameters are fixed by the supplied contract. Stop before phenotype reveal.',
            'execution_note': 'Fresh execution from verified remote parent inputs; previous Task 42 scripts were reviewed and adapted for incremental publication; no previous Task 42 numerical outputs were reused.'}
    RESULTS.mkdir(parents=True, exist_ok=True)
    write_json(RESULTS / 'preflight.json', data)
    with (RESULTS / 'thermodynamic_affinities.csv').open('w', newline='') as f:
        keys = ('class_id', 'exact_charge_equivalents', 'chemical_delta_g_over_rt',
                'electrical_delta_g_over_rt', 'affinity_over_rt', 'affinity_kJ_mol',
                'forward_direction', 'inherited_j4_fmol_s', 'scalar_direction_agrees_at_reference')
        w = csv.DictWriter(f, fieldnames=keys);w.writeheader()
        w.writerows({key: a[key] for key in keys} for a in panel)
    for c in CLASSES:
        directory = RESULTS / c;directory.mkdir()
        write_json(directory / 'budget.json', {'class_id': c, 'source_tests_pass': True,
            'stationary_solves': 0, 'stationary_residual_evaluations': 0,
            'started_cases': [], 'completed_cases': [], 'integration_attempts': 0,
            'numerical_retries': 0, 'stationary_numerical_retries': 0,
            'numerical_execution_seconds': 0., 'no_further_scientific_execution': False,
            'status': 'PREFLIGHT_FROZEN', 'optimisation_calls': 0, 'parameter_sweeps': 0,
            'genotype_rest_solves': 0, 'workers': 1, 'blas_threads': 1})
    print('PREFLIGHT_PASS: seven source, charge and thermodynamic audits recorded; no stationary or production solves.')


if __name__ == '__main__':
    main()
