"""Independent Task 46 milestone 01 inspection; no integrations or root solves.

Reevaluate saved dense states, independently sum trapezoids, and compare with
pinned Task 40/44 references. Preserve original replay outputs, including its
omitted volume-extrema comparison; this verifier explicitly covers that field.
"""
from __future__ import annotations
import os
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
             'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
import csv
from dataclasses import replace
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
PINS = {'task43': '547f113d9123ab1976774639a69483faf40cef34',
        'task44': 'e5fa9840bf96be3147c6118daee4942468c78f8b'}


def read(path):
    return json.loads(path.read_text())


def rows(path):
    with path.open() as stream:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]


def git(path, *args):
    return subprocess.check_output(['git', '-C', str(path), *args], text=True).strip()


def trapezoid(data, field, start=0.):
    selected = [r for r in data if r['time_s'] >= start]
    return math.fsum((b['time_s']-a['time_s'])*(a[field]+b[field])/2
                     for a, b in zip(selected, selected[1:]))


class Audit:
    def __init__(self):
        self.numeric = 0
        self.identical = 0
        self.gates = 0
        self.failures = []
        self.differences = []
        self.categories = {}

    def number(self, label, expected, actual, *, atol=1e-12, rtol=1e-10):
        expected, actual = float(expected), float(actual)
        self.numeric += 1
        category = label.split('.')[0]
        self.categories[category] = self.categories.get(category, 0)+1
        finite = math.isfinite(expected) and math.isfinite(actual)
        allowed = atol+rtol*max(abs(expected), abs(actual))
        delta = abs(expected-actual)
        if finite and expected == actual:
            self.identical += 1
        else:
            record = dict(label=label, expected=expected, actual=actual,
                          absolute_error=delta, allowed_error=allowed)
            if not finite or delta > allowed:
                self.failures.append(record)
            elif len(self.differences) < 30:
                self.differences.append(record)

    def gate(self, label, condition, detail=None):
        self.gates += 1
        if not condition:
            self.failures.append(dict(label=label, detail=detail))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source43', required=True, type=Path)
    parser.add_argument('--source44', required=True, type=Path)
    parser.add_argument('--replay', type=Path, default=HERE/'output/checkpoint_01_replay')
    parser.add_argument('--output', type=Path,
                        default=HERE/'output/checkpoint_01_independent_verification.json')
    args = parser.parse_args()
    assert not args.output.exists(), 'Independent verification records are immutable.'
    audit = Audit()
    pins = read(HERE/'source_pins.json')
    audit.gate('source.pins', pins['sources'] == PINS)
    sources = {'task43': args.source43.resolve(), 'task44': args.source44.resolve()}
    source_status = {}
    for name, path in sources.items():
        head, status = git(path, 'rev-parse', 'HEAD'), git(path, 'status', '--porcelain')
        audit.gate('source.'+name+'.head', head == PINS[name], head)
        audit.gate('source.'+name+'.clean_before', not status, status)
        source_status[name] = dict(head=head, clean_before=not status)
    source, replay = sources['task44'], args.replay.resolve()
    input_hashes = {str(p.relative_to(replay)): hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in replay.glob('*') if p.is_file()}
    verification = read(replay/'verification.json')
    audit.gate('replay.status', verification['status'] == 'passed')
    audit.gate('replay.source', verification['source_sha'] == PINS['task44'])
    audit.gate('replay.no_failed_comparisons', not verification['failed_comparisons'])
    audit.gate('replay.reported_count', verification['numerical_comparisons'] == len(verification['comparisons']))
    audit.gate('replay.reported_identical', verification['bit_identical'] ==
               sum(c['expected'] == c['actual'] for c in verification['comparisons']))
    for c in verification['comparisons']:
        audit.gate('replay.record.'+c['label'], c['passed'] and
                   c['bit_identical'] == (c['expected'] == c['actual']) and
                   c['absolute_error'] == abs(c['expected']-c['actual']) and
                   abs(c['expected']-c['actual']) <= c['allowed_error'])
    audit.gate('replay.case_count', verification['production_cases'] == 3)
    audit.gate('replay.root_count', verification['continuation_roots'] == 8)
    audit.gate('replay.nearby_root_count', verification['independent_nearby_roots'] == 4)
    for name, digest in verification['source_python_sha256'].items():
        audit.gate('source.hash.'+name, hashlib.sha256((source/name).read_bytes()).hexdigest() == digest)
    sys.path.insert(0, str(source/'analysis/44_full_system_mathematics'))
    import core
    import numpy as np
    rest, stim, y0, task, _ = core.reference()
    frozen = source/'results/40_ae4_equal_cation_routing'
    reference_rest = read(frozen/'wt_rest.json')
    rest_evaluation = rest.evaluate(0., y0)
    rest_row, bad, ratios = task.diagnose(rest, 0., y0, rest_evaluation)
    audit.gate('rest.full_inherited_gates', not bad, bad)
    for key, value in reference_rest['observables'].items():
        audit.number('rest.'+key, value, rest_row[key], atol=1e-20 if key.endswith('_A') else 1e-12)
    audit.gate('rest.stationary_amount_rhs', max(abs(rest_evaluation.rhs[k]) for k in
               [0,1,2,3,4,6,7,8,9,10]) < 1e-9)
    audit.gate('rest.stationary_volume_rhs', max(abs(rest_evaluation.rhs[k]) for k in [5,11]) < 1e-10)
    summaries, trajectories, totals, integrals = {}, {}, {}, {}
    sample_gates, maximum_conservation_ratio = 0, 0.
    flux_fields = {'nkcc1_cycles':'nkcc1_cycle_inward_fmol_s',
                   'ae4_signed_cl_flux':'ae4_cl_inward_fmol_s'}
    for case, expression in [('wt',1.), ('ae4_5pct',.05), ('ae4_null',0.)]:
        summary = read(replay/(case+'_summary.json'))
        reference = read(frozen/(case+'_verification.json'))
        data, states = rows(replay/(case+'_trajectory.csv')), rows(replay/(case+'_states.csv'))
        summaries[case], trajectories[case] = summary, data
        expected_grid = [0., 1e-6, *map(float, range(1,601))]
        audit.gate('trajectory.'+case+'.grid', [r['time_s'] for r in data] == expected_grid)
        audit.gate('trajectory.'+case+'.state_grid', [r['time_s'] for r in states] == expected_grid)
        audit.gate('trajectory.'+case+'.expression', summary['expression'] == expression)
        audit.gate('trajectory.'+case+'.counter', summary['checked_states'] == reference['checked_state_count'] == 934)
        audit.gate('trajectory.'+case+'.steps', summary['accepted_steps'] == reference['accepted_solver_steps'])
        audit.gate('trajectory.'+case+'.nfev', summary['nfev'] == reference['integration_attempts'][0]['nfev'])
        audit.gate('trajectory.'+case+'.solver', summary['solver'] == reference['integration_attempts'][0]['solver'])
        audit.gate('trajectory.'+case+'.no_failed_gates', summary['failed_gates'] == [])
        audit.gate('trajectory.'+case+'.initial_state', summary['initial_state'] == reference_rest['state_vector'])
        for field in reference['minima']:
            audit.number('extrema.'+case+'.min.'+field, reference['minima'][field], summary['minima'][field], atol=1e-20 if field.endswith('_A') else 1e-12)
            audit.number('extrema.'+case+'.max.'+field, reference['maxima'][field], summary['maxima'][field], atol=1e-20 if field.endswith('_A') else 1e-12)
        for field, ratio in summary['maximum_conservation_ratios'].items():
            audit.gate('trajectory.'+case+'.conservation.'+field, math.isfinite(ratio) and ratio <= 1)
            audit.number('ratios.'+case+'.'+field, reference['max_conservation_ratios'][field], ratio)
        genotype = replace(core.WT, ae4_expression=expression)
        for index, (row, state) in enumerate(zip(data, states)):
            y = np.array([state[k] for k in stim.state_names])
            ev = stim.evaluate(row['time_s'], y, genotype=genotype)
            evaluated, bad, ratios = task.diagnose(stim, row['time_s'], y, ev)
            audit.gate('state.'+case+'.'+str(index)+'.inherited_gates', not bad, bad)
            sample_gates += 1
            maximum_conservation_ratio = max(maximum_conservation_ratio, *ratios.values())
            for field, value in evaluated.items():
                audit.number('diagnostics.'+case+'.'+str(index)+'.'+field, row[field], value,
                             atol=1e-20 if field.endswith('_A') else 1e-12)
        for i, name in enumerate(stim.state_names):
            audit.number('endpoint.'+case+'.'+name, reference['last_checked_state_vector'][i], states[-1][name])
        for row in rows(frozen/(case+'_timeseries.csv')):
            actual = next(r for r in data if r['time_s'] == row['time_s'])
            for field, value in row.items():
                audit.number('frozen_diagnostics.'+case+'.'+str(row['time_s'])+'.'+field, value, actual[field],
                             atol=1e-20 if field.endswith('_A') else 1e-12)
        cumulative = 0.
        audit.number('quadrature.'+case+'.initial', 0., data[0]['cumulative_outflow_0_t_pL'])
        for a, b in zip(data, data[1:]):
            cumulative += (b['time_s']-a['time_s'])*(a['q_out_pL_s']+b['q_out_pL_s'])/2
            audit.number('quadrature.'+case+'.t'+str(b['time_s']), b['cumulative_outflow_0_t_pL'], cumulative)
        totals[case] = trapezoid(data, 'q_out_pL_s')
        mean = trapezoid(data, 'q_out_pL_s', 60.)/540.
        audit.number('secretion.'+case+'.total', reference['secretion']['cumulative_0_600_pL'], totals[case])
        audit.number('secretion.'+case+'.summary_total', summary['cumulative_0_600_pL'], totals[case])
        audit.number('secretion.'+case+'.mean', reference['secretion']['mean_60_600_pL_s'], mean)
        activation = {'mean':mean>1.1*rest_row['q_out_pL_s'],
                      'cumulative':totals[case]>1.1*600*rest_row['q_out_pL_s'],
                      'endpoint':data[-1]['q_out_pL_s']>=rest_row['q_out_pL_s']}
        audit.gate('activation.'+case+'.record', activation == summary['activation_checks'])
        audit.gate('activation.'+case+'.scope', summary['activation_required'] == (case == 'wt'))
        audit.gate('activation.'+case+'.wt_required', case != 'wt' or all(activation.values()))
        integrals[case] = {name:trapezoid(data, field, 60.) for name, field in flux_fields.items()}
        for name, actual in integrals[case].items():
            audit.number('integrals.'+case+'.'+name, summary['integrals_60_600'][name], actual)
    comparison = read(frozen/'comparison_summary.json')
    deficits = {case:1-total/totals['wt'] for case,total in totals.items()}
    for case in ['wt','ae4_5pct','ae4_null']:
        audit.number('deficit.'+case+'.replay', verification['summary'][case]['cumulative_deficit'], deficits[case])
        if case != 'wt':
            audit.number('deficit.'+case+'.reference', comparison['perturbations'][case]['cumulative_deficit_percent']/100, deficits[case])
    roots = read(replay/'equilibria.json')['roots']
    saved_roots = {r['ae4_expression']:r for r in read(source/'analysis/44_full_system_mathematics/output/equilibria.json')['roots'] if r['mechanism'] == 'C0'}
    audit.gate('roots.expressions', [r['ae4_expression'] for r in roots] == [1.,.9,.75,.5,.25,.1,.05,0.])
    by_expression = {}
    for record in roots:
        expression = record['ae4_expression']
        by_expression[expression] = record
        reference = saved_roots[expression]
        audit.gate('root.'+str(expression)+'.physiology', record['full_inherited_gates_passed'] and record['physiological'] and not record['failed_gates'])
        audit.gate('root.'+str(expression)+'.success', record['root_success'])
        audit.gate('root.'+str(expression)+'.residual', record['residual_max'] < 1e-7)
        audit.gate('root.'+str(expression)+'.jacobian', record['jacobian_relative_step_error'] < 1e-4)
        audit.gate('root.'+str(expression)+'.conservation', max(record['conservation_ratios'].values()) < 1)
        eig = record['eigenvalues_per_s']
        audit.gate('root.'+str(expression)+'.spectrum_count', len(eig) == 11)
        audit.gate('root.'+str(expression)+'.stable', all(real < 0 for real, imag in eig))
        audit.gate('root.'+str(expression)+'.regulatory_mode', min(abs(complex(*v)+1/30) for v in eig) < 1e-9)
        for k, (expected, actual) in enumerate(zip(reference['eigenvalues_per_s'], eig)):
            for component in [0,1]:
                audit.number('eigenvalues.'+str(expression)+'.'+str(k)+'.'+str(component), expected[component], actual[component])
        audit.number('slowmode.'+str(expression)+'.definition', -1/max(v[0] for v in eig), record['slowest_decay_s'])
        audit.number('slowmode.'+str(expression)+'.reference', reference['slowest_decay_s'], record['slowest_decay_s'])
        for field, actual in record['observables'].items():
            audit.number('root_observables.'+str(expression)+'.'+field, reference['observables'][field], actual,
                         atol=1e-20 if field.endswith('_A') else 1e-12)
        for k, actual in enumerate(record['state']):
            audit.number('root_states.'+str(expression)+'.'+str(k), reference['state'][k], actual)
    implicit = read(replay/'implicit_verification.json')
    review = read(source/'analysis/44_full_system_mathematics/output/equilibrium_review.json')
    ref_implicit = next(r for r in review['implicit_sensitivity_two_step_validation'] if r['expression'] == 1.)
    audit.gate('ift.reference_identity', implicit == ref_implicit)
    audit.gate('ift.steps', [r['step'] for r in implicit['trials']] == [1e-4,5e-5])
    for trial in implicit['trials']:
        audit.gate('ift.component_gate.'+str(trial['step']), max(trial['component_relative_errors']) < 1e-3)
        audit.gate('ift.state_gate.'+str(trial['step']), trial['state_derivative_relative_error'] < 1e-3)
        audit.gate('ift.root_gate.'+str(trial['step']), trial['nearby_root_max_residual'] < 1e-7)
    extra = 2*(integrals['ae4_null']['nkcc1_cycles']-integrals['wt']['nkcc1_cycles'])
    lost = integrals['wt']['ae4_signed_cl_flux']-integrals['ae4_null']['ae4_signed_cl_flux']
    metrics = dict(extra_NKCC_chloride_fmol=extra,lost_AE4_chloride_fmol=lost,
                   signed_integrated_compensation=extra/lost,
                   relative_NKCC_increase=extra/(2*integrals['wt']['nkcc1_cycles']),
                   steady_null_deficit=1-by_expression[0.]['observables']['q_pL_s']/by_expression[1.]['observables']['q_pL_s'],
                   local_gain=-2*implicit['implicit_derivatives'][1]/implicit['implicit_derivatives'][2])
    replay_metrics = read(replay/'compensation.json')
    reference_metrics = read(source/'analysis/44_full_system_mathematics/output/compensation_comparison.json')
    audit.gate('compensation.integration_window', replay_metrics['integration_window_s'] == [60.,600.])
    for key, actual in metrics.items():
        audit.number('compensation.'+key+'.replay', replay_metrics[key], actual)
        expected = ref_implicit['chloride_compensation_gain'] if key == 'local_gain' else reference_metrics[key]
        audit.number('compensation.'+key+'.reference', expected, actual)
    for name, digest in input_hashes.items():
        audit.gate('output.unchanged.'+name, hashlib.sha256((replay/name).read_bytes()).hexdigest() == digest)
    for name, path in sources.items():
        clean = not git(path, 'status', '--porcelain')
        audit.gate('source.'+name+'.clean_after', clean)
        source_status[name]['clean_after'] = clean
    result = dict(status='passed' if not audit.failures else 'failed',
                  verification_scope='Read saved products and freshly evaluate saved states; no integration, root solve, parameter search or output overwrite.',
                  source_pins=source_status, numeric_comparisons=audit.numeric,
                  bit_identical_numeric_values=audit.identical,
                  tolerated_numeric_values=audit.numeric-audit.identical-len([f for f in audit.failures if 'absolute_error' in f]),
                  logical_or_gate_checks=audit.gates,numeric_categories=audit.categories,
                  independent_dense_state_evaluations=sample_gates,
                  maximum_dense_conservation_tolerance_ratio=maximum_conservation_ratio,
                  production_checked_states_reported={case:s['checked_states'] for case,s in summaries.items()},
                  continuation_roots_checked=len(roots),eigenvalues_checked=sum(len(r['eigenvalues_per_s']) for r in roots),
                  deficits=deficits,compensation=metrics,
                  verification_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  replay_input_sha256=input_hashes, failures=audit.failures,
                  tolerated_difference_examples=audit.differences,
                  issues_preserved=[dict(issue='Initial replay comparison loop used absent field cell_volume_pL and skipped volume extrema.',
                                         resolution='This independent verifier compares all 34 minimum and maximum diagnostic fields, including volume_i_pL, for every case; existing replay output is preserved.'),
                                    dict(issue='Accepted solver endpoint states between integer seconds are represented by replay summary extrema and gate maxima, not separate state CSV rows.',
                                         resolution='Compared those summary extrema, conservation maxima and counters against frozen production reference; independently reevaluated all 602 saved grid states for each case.')])
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
    print(json.dumps({key:result[key] for key in ['status','numeric_comparisons','bit_identical_numeric_values','tolerated_numeric_values','logical_or_gate_checks','independent_dense_state_evaluations','eigenvalues_checked','failures']}))
    if audit.failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
