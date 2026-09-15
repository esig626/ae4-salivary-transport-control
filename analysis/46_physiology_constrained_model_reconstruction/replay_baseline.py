"""Task 46 milestone 01: replay pinned Task 40/44 without modifying inputs.

Production stepping and quadrature follow Task 40; equilibrium continuation and
IFT calculations call the unchanged Task 44 functions. No fitting is performed.
"""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
            'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import csv
from dataclasses import asdict, replace
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

TASK44_SHA = 'e5fa9840bf96be3147c6118daee4942468c78f8b'
TASK46 = Path(__file__).resolve().parent

def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False,
                               default=lambda value: value.item()) + '\n')

def git(root, *args):
    return subprocess.check_output(['git', '-C', str(root), *args], text=True).strip()

def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=TASK46/'output/checkpoint_01_replay')
    args = parser.parse_args()
    source, out = args.source.resolve(), args.output.resolve()
    assert git(source, 'rev-parse', 'HEAD') == TASK44_SHA
    assert not git(source, 'status', '--porcelain'), 'Pinned input checkout must be clean.'
    assert not out.exists(), 'Checkpoint outputs are immutable; choose a new scratch destination for replay.'
    out.mkdir(parents=True)
    manifests = {str(p.relative_to(source)): hashlib.sha256(p.read_bytes()).hexdigest()
                 for folder in [source/'src', source/'analysis/44_full_system_mathematics']
                 for p in folder.rglob('*.py')}
    sys.path.insert(0, str(source/'analysis/44_full_system_mathematics'))
    import core
    import equilibrium_review as review
    import numpy as np
    import scipy
    from scipy.integrate import Radau, cumulative_trapezoid
    from scipy.optimize import linear_sum_assignment
    from modern_full_model.validation import PRODUCTION_RADAU
    from modern_full_model.model import WT, Genotype
    rest, stim, y0, task, source_audit = core.reference()
    # Import constants/diagnostic field definitions only; do not call the old
    # task driver, change its frozen completion ledger, or weaken its guards.
    sys.path.insert(0, str(source/'analysis/40_ae4_equal_cation_routing'))
    historical = load_module('task40_readonly_driver', source/'analysis/40_ae4_equal_cation_routing/run_case.py')
    comparisons = []
    def compare(label, expected, actual, rtol=1e-8, atol=1e-10):
        expected, actual = float(expected), float(actual)
        delta = abs(expected-actual)
        allowed = atol + rtol*max(abs(expected), abs(actual))
        comparisons.append(dict(label=label, expected=expected, actual=actual,
                                absolute_error=delta, allowed_error=allowed,
                                bit_identical=expected == actual, passed=delta <= allowed))

    rest_ev = rest.evaluate(0., y0)
    rest_row, rest_fail, rest_ratios = task.diagnose(rest, 0., y0, rest_ev)
    assert not rest_fail
    frozen_rest = json.loads((source/'results/40_ae4_equal_cation_routing/wt_rest.json').read_text())
    compare('rest.flow', frozen_rest['observables']['q_out_pL_s'], rest_row['q_out_pL_s'])
    records = {}
    for case, expression in [('wt',1.), ('ae4_5pct',.05), ('ae4_null',0.)]:
        genotype = replace(WT, ae4_expression=expression)
        grid, states, minima, maxima, maximum_ratios = [], [], {}, {}, {}
        checked, accepted = 0, 0
        failures = []
        def inspect(t, y):
            nonlocal checked
            ev = stim.evaluate(float(t), y, genotype=genotype)
            row, bad, ratios = task.diagnose(stim, float(t), y, ev)
            checked += 1
            if bad:
                failures.extend(dict(time_s=float(t), **b) for b in bad)
                raise RuntimeError(f'{case}: physiological gate failed at {t}: {bad}')
            for k, v in row.items():
                if k != 'time_s':
                    minima[k] = min(minima.get(k,v),v)
                    maxima[k] = max(maxima.get(k,v),v)
            for k, v in ratios.items():
                maximum_ratios[k] = max(maximum_ratios.get(k,0.),float(v))
            if t == 0. or t == 1e-6 or float(t).is_integer():
                if not grid or grid[-1]['time_s'] != t:
                    grid.append(row)
                    states.append([float(t), *map(float,y)])
        inspect(0.,y0)
        inspect(1e-6,y0)
        spec = PRODUCTION_RADAU
        solver = Radau(lambda t,y:stim.rhs(t,y,genotype=genotype),1e-6,y0.copy(),600.,
                       rtol=spec.rtol,atol=spec.atol_vector(stim.state_names),max_step=spec.max_step_s)
        next_time = 1.
        while solver.status == 'running':
            solver.step()
            if solver.status == 'failed':
                raise RuntimeError(f'{case}: original production Radau failed')
            accepted += 1
            dense = solver.dense_output()
            while next_time <= solver.t:
                inspect(next_time,dense(next_time))
                next_time += 1.
            if solver.t != next_time-1.:
                inspect(float(solver.t),solver.y)
        assert solver.t == 600.
        ts = np.array([r['time_s'] for r in grid])
        assert np.array_equal(ts,np.r_[0.,1e-6,np.arange(1.,601.)])
        q = np.array([r['q_out_pL_s'] for r in grid])
        cumulative = cumulative_trapezoid(q,ts,initial=0.)
        for row,value in zip(grid,cumulative):
            row['cumulative_outflow_0_t_pL'] = float(value)
        total = float(np.trapezoid(q,ts))
        window = ts >= 60.
        mean = float(np.trapezoid(q[window],ts[window])/540.)
        integrals = {}
        for name,field,positive,unit in historical.FLUXES:
            values = np.array([r[field] for r in grid])[window]
            if positive: values = np.maximum(values,0.)
            integrals[name] = float(np.trapezoid(values,ts[window]))
        frozen = json.loads((source/f'results/40_ae4_equal_cation_routing/{case}_verification.json').read_text())
        compare(case+'.cumulative',frozen['secretion']['cumulative_0_600_pL'],total,atol=1e-8)
        compare(case+'.mean',frozen['secretion']['mean_60_600_pL_s'],mean)
        compare(case+'.endpoint',frozen['secretion']['endpoint_pL_s'],q[-1])
        original = list(csv.DictReader((source/f'results/40_ae4_equal_cation_routing/{case}_timeseries.csv').open()))
        lookup = {r['time_s']:r for r in grid}
        for row in original:
            t = float(row['time_s'])
            for field,expected in row.items():
                compare(f'{case}.t{t}.{field}',expected,lookup[t][field],
                        rtol=2e-5,atol=1e-20 if field.endswith('_A') else 1e-10)
        for row in csv.DictReader((source/f'results/40_ae4_equal_cation_routing/{case}_integrated_fluxes.csv').open()):
            compare(case+'.integral.'+row['quantity'],row['value'],integrals[row['quantity']],atol=1e-9)
        # The separately recorded independent verifier covers volume_i_pL.
        for field in ['na_i_mM','k_i_mM','cl_i_mM','ph_i']:
            if field in frozen['minima']:
                compare(case+'.minimum.'+field,frozen['minima'][field],minima[field])
                compare(case+'.maximum.'+field,frozen['maxima'][field],maxima[field])
        activation = {'mean':mean>1.1*rest_row['q_out_pL_s'],
                      'cumulative':total>1.1*600*rest_row['q_out_pL_s'],
                      'endpoint':q[-1]>=rest_row['q_out_pL_s']}
        assert case != 'wt' or all(activation.values())
        record = dict(expression=expression,solver=asdict(spec),accepted_steps=accepted,
                      checked_states=checked,nfev=solver.nfev,failed_gates=failures,
                      maximum_conservation_ratios=maximum_ratios,minima=minima,maxima=maxima,
                      cumulative_0_600_pL=total,mean_60_600_pL_s=mean,endpoint=grid[-1],
                      initial_state=y0.tolist(),final_state=solver.y.tolist(),integrals_60_600=integrals,
                      activation_checks=activation,activation_required=case=='wt')
        records[case] = record
        core.write_csv(out/(case+'_trajectory.csv'),grid)
        core.write_csv(out/(case+'_states.csv'),[dict(zip(['time_s',*stim.state_names],row)) for row in states])
        dump(out/(case+'_summary.json'),record)
        print('PRODUCTION_PASS',case,total,checked,flush=True)

    # Exact source functions and original continuation order. The seed is a
    # freshly integrated constant-input WT state, never a saved solved root.
    constant = core.clone(stim)
    scale = core.Scaling(constant,y0)
    transient = core.integrate(constant,y0)
    seed = transient.y[:-1,-1]
    roots = {}
    saved = json.loads((source/'analysis/44_full_system_mathematics/output/equilibria.json').read_text())
    originals = {r['ae4_expression']:r for r in saved['roots'] if r['mechanism']=='C0'}
    for expression in [1.,.9,.75,.5,.25,.1,.05,0.]:
        genotype = Genotype('AE4',ae4_expression=expression)
        record = core.stationary(constant,scale,seed,genotype)
        record.update(mechanism='C0',ae4_expression=expression)
        seed = np.array(record['state'])
        ev = constant.evaluate(1.,seed,genotype=genotype)
        _,bad,ratios = task.diagnose(constant,1.,seed,ev)
        assert not bad and max(ratios.values())<1
        record['full_inherited_gates_passed'] = True
        record['conservation_ratios'] = ratios
        ref = originals[expression]
        for field,value in ref['observables'].items():
            compare(f'equilibrium.{expression}.{field}',value,record['observables'][field],
                    atol=1e-19 if field.endswith('_A') else 1e-9)
        for i,(expected,actual) in enumerate(zip(ref['state'],record['state'])):
            compare(f'equilibrium.{expression}.state{i}',expected,actual)
        a = np.array([complex(*v) for v in ref['eigenvalues_per_s']])
        b = np.array([complex(*v) for v in record['eigenvalues_per_s']])
        ii,jj = linear_sum_assignment(abs(a[:,None]-b[None,:]))
        for i,j in zip(ii,jj):
            compare(f'equilibrium.{expression}.eigenvalue{i}.real',a[i].real,b[j].real,rtol=1e-3,atol=1e-10)
            compare(f'equilibrium.{expression}.eigenvalue{i}.imag',a[i].imag,b[j].imag,rtol=1e-3,atol=1e-10)
        compare(f'equilibrium.{expression}.slow_mode',ref['slowest_decay_s'],record['slowest_decay_s'],rtol=1e-4)
        roots[expression] = record
    # Two-step independent nearby-root verification is the original Task44
    # verifier. This is baseline replay, not a new causal/sensitivity study.
    implicit = review.ift_validation(constant,scale,roots[1.])
    assert max(max(t['component_relative_errors']) for t in implicit['trials'])<1e-3
    ref_review = json.loads((source/'analysis/44_full_system_mathematics/output/equilibrium_review.json').read_text())
    ref_ift = next(r for r in ref_review['implicit_sensitivity_two_step_validation'] if r['expression']==1.)
    compare('local_compensation',ref_ift['chloride_compensation_gain'],implicit['chloride_compensation_gain'],rtol=1e-5)
    wt_i = records['wt']['integrals_60_600']
    null_i = records['ae4_null']['integrals_60_600']
    extra = 2*(null_i['nkcc1_cycles']-wt_i['nkcc1_cycles'])
    lost = wt_i['ae4_signed_cl_flux']-null_i['ae4_signed_cl_flux']
    compensation = dict(extra_NKCC_chloride_fmol=extra,lost_AE4_chloride_fmol=lost,
                        signed_integrated_compensation=extra/lost,
                        relative_NKCC_increase=extra/(2*wt_i['nkcc1_cycles']),
                        steady_null_deficit=1-roots[0.]['observables']['q_pL_s']/roots[1.]['observables']['q_pL_s'])
    ref_comp = json.loads((source/'analysis/44_full_system_mathematics/output/compensation_comparison.json').read_text())
    for field,value in compensation.items():
        compare('compensation.'+field,ref_comp[field],value,rtol=1e-7,atol=1e-10)
    for case in ['wt','ae4_5pct','ae4_null']:
        records[case]['cumulative_deficit'] = 1-records[case]['cumulative_0_600_pL']/records['wt']['cumulative_0_600_pL']
    compensation.update(integration_window_s=[60.,600.],local_gain=implicit['chloride_compensation_gain'])
    dump(out/'equilibria.json',{'roots':list(roots.values()),'seed_source':'fresh constant-stimulus integration'})
    dump(out/'implicit_verification.json',implicit)
    dump(out/'compensation.json',compensation)
    for name,digest in manifests.items():
        assert hashlib.sha256((source/name).read_bytes()).hexdigest()==digest,name
    assert not git(source,'status','--porcelain'), 'Pinned source/result bytes were changed.'
    failed = [r for r in comparisons if not r['passed']]
    result = dict(status='passed' if not failed else 'failed',source_sha=TASK44_SHA,
                  task46_parent_sha=git(TASK46,'rev-parse','HEAD'),source_python_sha256=manifests,
                  numpy=np.__version__,scipy=scipy.__version__,python=sys.version,
                  numerical_comparisons=len(comparisons),bit_identical=sum(r['bit_identical'] for r in comparisons),
                  failed_comparisons=failed,comparisons=comparisons,production_cases=3,
                  continuation_roots=8,independent_nearby_roots=4,
                  all_production_and_root_gates_passed=True,source_checkout_clean=True,
                  scientific_equations_or_parameters_changed=False,
                  summary={k:dict(cumulative_0_600_pL=v['cumulative_0_600_pL'],
                                 cumulative_deficit=v['cumulative_deficit'],
                                 checked_states=v['checked_states']) for k,v in records.items()},
                  compensation=compensation)
    dump(out/'verification.json',result)
    print('BASELINE_COMPARISON',result['status'],result['numerical_comparisons'],
          'bit_identical',result['bit_identical'],'failures',len(failed),flush=True)
    if failed: raise SystemExit(1)

if __name__ == '__main__':
    main()
