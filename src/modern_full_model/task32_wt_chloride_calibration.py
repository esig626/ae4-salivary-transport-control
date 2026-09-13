"""Task 32: bounded sequential WT-only calibration of two abundance scalars.

No transport equation is changed. The AE4 dimensionless activity multiplies
the existing carrier_amount_fmol, independently of genotype expression.
The only multivariate solve is stationary state closure, never parameter fit.
Run with OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=1.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import hashlib
import json
import math
import os
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np
from scipy.optimize import least_squares

from . import task31_nhe1_repair as prior
from . import task30_nhe1_repair as inherited
from .task30_nhe1_repair import (
    load_background, state_from_coordinates, COORDINATE_NAMES,
    COORDINATE_LOWER, COORDINATE_UPPER, INDEPENDENT_ROWS,
    INDEPENDENT_RHS_SCALES, OMITTED_CHARGE_ROWS, AMOUNT_ROWS, VOLUME_ROWS,
    ROOT_SCALED_TOLERANCE, OMITTED_ROW_TOLERANCE, CHARGE_TOLERANCE_FMOL,
    CURRENT_TOLERANCE_A, REGULATORY_TOLERANCE_S_INV,
    BOUNDARY_RELATIVE_TOLERANCE, JACOBIAN_RANK_RELATIVE_TOLERANCE,
    _conservation_audit, _rest_observables, _write_csv,
)
from .model import WT, AE4_NULL
from .nhe1_cha2009 import published_cha_kinetics
from .validation import sha256_object

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / 'results/32_wt_chloride_calibration'
ANALYSIS = ROOT / 'analysis/32_wt_chloride_calibration'
PREPARED_HEAD = '6027f68dfcf8113c13d726b91d0b38cd0349ef47'
BRANCH = 'codex/task-32-wt-chloride-calibration'
LIMITS = dict(stationary_solver_calls=18, stationary_residual_evaluations=30000,
              ae4_scalar_evaluations=6, nhe1_scalar_evaluations=4,
              stimulated_integrations=2)
AE4_BOUNDS = (1e-6, 1e6)
NHE1_BOUNDS = (1e-8, 0.1)
PER_SOLVE_RESIDUAL_CEILING = 2400
WT_PH_INTERVAL = (6.84, 6.98)
WT_CL_INTERVAL = (48.6, 51.6)


class BoundedStop(RuntimeError):
    pass


class LocalSolveStop(RuntimeError):
    pass


def write_json(path, value):
    """Atomic UTF-8 checkpoint; no nonfinite JSON numbers."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True,
                                   allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


@dataclass
class Budget:
    stationary_solver_calls: int = 0
    stationary_residual_evaluations: int = 0
    ae4_scalar_evaluations: int = 0
    nhe1_scalar_evaluations: int = 0
    stimulated_integrations: int = 0
    numerical_wall_seconds: float = 0.0

    def __post_init__(self):
        self._clock = time.perf_counter()

    def update_time(self):
        now = time.perf_counter()
        self.numerical_wall_seconds += now - self._clock
        self._clock = now

    def check(self):
        self.update_time()
        if self.numerical_wall_seconds >= 1800:
            raise BoundedStop('1800 s numerical execution limit')
        for name, limit in LIMITS.items():
            if getattr(self, name) >= limit:
                raise BoundedStop(name + ' limit reached; no further numerical case')

    def count(self, name):
        self.check()
        setattr(self, name, getattr(self, name) + 1)

    def count_residual(self):
        # A scalar/solver count is reserved at case start. Reaching that count
        # permits this last case to finish, but never another numerical case.
        self.update_time()
        if self.numerical_wall_seconds >= 1800:
            raise BoundedStop('1800 s numerical execution limit')
        if self.stationary_residual_evaluations >= LIMITS['stationary_residual_evaluations']:
            raise BoundedStop('30000 residual evaluations exhausted')
        self.stationary_residual_evaluations += 1
        if self.stationary_residual_evaluations % 100 == 0:
            self.save()

    def save(self):
        self.update_time()
        write_json(RESULTS / 'budget.json', {**asdict(self), 'limits': LIMITS,
            'scientific_workers': 1, 'blas_threads': 1,
            'per_solve_residual_ceiling': PER_SOLVE_RESIDUAL_CEILING,
            'accepted_identical_cases_rerun': 0,
            'limits_respected': all(getattr(self, k) <= v for k, v in LIMITS.items())
                and self.numerical_wall_seconds <= 1800})


def build_model(background, *, carrier_amount_fmol, ae4_activity, stimulated=False):
    """Scale existing AE4 capacity, preserving kinetics, routing and genotype."""
    if not math.isfinite(ae4_activity) or ae4_activity <= 0:
        raise ValueError('AE4 activity must be finite and positive')
    if not math.isfinite(carrier_amount_fmol) or carrier_amount_fmol <= 0:
        raise ValueError('NHE1 amount must be finite and positive')
    scaled = replace(background, ae4_parameters=replace(background.ae4_parameters,
        carrier_amount_fmol=background.ae4_parameters.carrier_amount_fmol * ae4_activity))
    return prior.build_model(scaled, carrier_amount_fmol=carrier_amount_fmol,
                             stimulated=stimulated)


def flux_row(evaluation, model):
    row = inherited._flux_row(evaluation, model)
    row['paracellular_cl_lumen_to_bath_fmol_s'] = inherited.current_to_fmol_s(
        evaluation.diagnostics.membranes.currents_A['para_cl'], valence=-1,
        faraday_C_mol=model.parameters.constants.faraday_C_mol)
    return row


def acceptance(solution, *, require_wt_intervals=False):
    """Null measurements never enter this stationary admissibility gate."""
    o = solution['observables']
    core = (solution['numerical_pass'] and solution['finite_capacity']['pass']
            and solution['exact_ae4_zero'] and o['na_i_mM'] <= 30
            and o['volume_i_pL'] < 3)
    return bool(core and (not require_wt_intervals or (
        WT_PH_INTERVAL[0] <= o['ph_i'] <= WT_PH_INTERVAL[1]
        and WT_CL_INTERVAL[0] <= o['cl_i_mM'] <= WT_CL_INTERVAL[1])))


def allocation(solution):
    f = solution['fluxes']
    ae4, nkcc, ae2 = (f[k] for k in ('ae4_cl_cell_fmol_s',
        'nkcc1_cl_in_fmol_s', 'ae2_cl_in_hco3_out_fmol_s'))
    a, n, b = max(ae4, 0), max(nkcc, 0), max(ae2, 0)
    positive, net = a+n+b, ae4+nkcc+ae2
    return {'label': solution['label'], 'background': solution['background'],
        'genotype': solution['genotype'], 'admissible': solution['admissible'],
        'ae4_signed_cl_fmol_s': ae4, 'nkcc1_signed_cl_fmol_s': nkcc,
        'ae2_signed_cl_fmol_s': ae2, 'ae4_positive_cl_fmol_s': a,
        'nkcc1_positive_cl_fmol_s': n, 'ae2_positive_cl_fmol_s': b,
        'total_positive_basolateral_cl_fmol_s': positive,
        'net_basolateral_cl_fmol_s': net,
        'ae4_over_positive_nkcc1': a/n if n else None,
        'ae4_over_positive_ae2': a/b if b else None,
        'ae4_fraction_positive_basolateral': a/positive if positive else None,
        'ae4_over_net_basolateral': a/net if net else None}


def solve_stationary(background, *, carrier_amount_fmol, ae4_activity,
                     genotype, start_coordinates, label, budget, scalar=None):
    """Task 31 bounded state-closure method and gates, without phenotype gates."""
    model = build_model(background, carrier_amount_fmol=carrier_amount_fmol,
                        ae4_activity=ae4_activity)
    key = sha256_object({'background': background.alias, 'parameters': model.parameters,
        'ae4_parameters': model.ae4_parameters, 'genotype': asdict(genotype),
        'start_coordinates': list(start_coordinates), 'kinetics': published_cha_kinetics()})
    path = RESULTS / 'cases.json'
    cases = json.loads(path.read_text()) if path.exists() else {}
    if key in cases:
        return cases[key]
    # A successful case is reusable independently of the proposed warm start.
    for case in cases.values():
        if (case['background'] == background.alias and case['ae4_expression'] == genotype.ae4_expression
            and case['carrier_amount_fmol'] == carrier_amount_fmol
            and case['ae4_activity'] == ae4_activity and case['numerical_pass']):
            return case
    budget.check()
    if scalar:
        setattr(budget, scalar+'_scalar_evaluations', getattr(budget, scalar+'_scalar_evaluations') + 1)
    budget.stationary_solver_calls += 1
    write_json(RESULTS / 'in_progress.json', {'key': key, 'label': label,
        'carrier_amount_fmol': carrier_amount_fmol, 'ae4_activity': ae4_activity,
        'ae4_expression': genotype.ae4_expression})
    budget.save()
    start = np.asarray(start_coordinates, dtype=float)
    if np.any(start <= COORDINATE_LOWER) or np.any(start >= COORDINATE_UPPER):
        raise ValueError('start must be strictly inside inherited coordinate bounds')
    local = 0
    best_x, best_cost = start.copy(), math.inf

    def counted(candidate, *, audit=False):
        nonlocal local, best_x, best_cost
        ceiling = PER_SOLVE_RESIDUAL_CEILING if audit else PER_SOLVE_RESIDUAL_CEILING - 13
        if local >= ceiling:
            raise LocalSolveStop('per-solve actual residual ceiling')
        budget.count_residual()
        local += 1
        try:
            state = state_from_coordinates(model, candidate)
            raw = model.rhs(0, state, genotype=genotype)
            values = raw[INDEPENDENT_ROWS] / INDEPENDENT_RHS_SCALES
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
            values = np.full(10, 1e6)
        cost = float(values@values)/2
        if cost < best_cost:
            best_x, best_cost = np.array(candidate), cost
        return values

    try:
        fit = least_squares(counted, start, bounds=(COORDINATE_LOWER, COORDINATE_UPPER),
            max_nfev=2400, xtol=1e-12, ftol=1e-12, gtol=1e-12, x_scale='jac')
    except LocalSolveStop as exc:
        final = best_x.copy()
        f0 = counted(final, audit=True)
        jac = np.empty((10, 10))
        for j in range(10):
            dx = math.sqrt(np.finfo(float).eps)*max(abs(final[j]), 1)
            if final[j]+dx >= COORDINATE_UPPER[j]:
                dx = -dx
            xj = final.copy()
            xj[j] += dx
            jac[:, j] = (counted(xj, audit=True)-f0)/dx
        fit = SimpleNamespace(x=final, jac=jac, success=False, message=str(exc),
            nfev=None, cost=float(f0@f0)/2, optimality=float(np.max(np.abs(jac.T@f0))))
    x = np.asarray(fit.x)
    state = state_from_coordinates(model, x)
    budget.count_residual()
    local += 1
    e = model.evaluate(0, state, genotype=genotype)
    budget.count_residual()
    local += 1
    later = model.evaluate(600, state, genotype=genotype)
    raw = np.asarray(e.rhs)
    scaled = raw[INDEPENDENT_ROWS]/INDEPENDENT_RHS_SCALES
    span = COORDINATE_UPPER - COORDINATE_LOWER
    singular = np.linalg.svd(np.asarray(fit.jac)*span[None, :], compute_uv=False)
    rank = int(np.sum(singular > JACOBIAN_RANK_RELATIVE_TOLERANCE*singular[0]))
    distances = np.minimum(x-COORDINATE_LOWER, COORDINATE_UPPER-x)/span
    boundary_hits = [n for n, d in zip(COORDINATE_NAMES, distances)
                     if d <= BOUNDARY_RELATIVE_TOLERANCE]
    charge = max(abs(float(v)) for v in e.diagnostics.state_charge_fmol.values())
    current = max(abs(float(v)) for v in e.diagnostics.membranes.current_residuals_A.values())
    regulatory = float(np.max(np.abs(raw[12:]))) if raw.size > 12 else 0.
    conservation, ratio, conservation_pass = _conservation_audit(e)
    capacity = prior._finite_capacity_audit(e, model)
    obs = _rest_observables(e, state)
    full_rhs = bool(np.max(np.abs(scaled)) <= ROOT_SCALED_TOLERANCE
        and np.max(np.abs(raw[OMITTED_CHARGE_ROWS])) <= OMITTED_ROW_TOLERANCE
        and regulatory <= REGULATORY_TOLERANCE_S_INV)
    positive = bool(np.all(np.isfinite(state)) and np.all(state[:12] > 0)
        and all(math.isfinite(v) for v in obs.values())
        and all(v > 0 for k, v in obs.items() if k.endswith('_mM')))
    numerical = bool(fit.success and full_rhs and charge <= CHARGE_TOLERANCE_FMOL
        and current <= CURRENT_TOLERANCE_A and conservation_pass and rank == 10
        and not boundary_hits and positive and np.array_equal(raw, later.rhs))
    exact_zero = all(float(v) == 0 for v in (e.diagnostics.ae4.na_cell_fmol_s,
        e.diagnostics.ae4.k_cell_fmol_s, e.diagnostics.ae4.cl_cell_fmol_s,
        e.diagnostics.ae4.hco3_cell_fmol_s, e.diagnostics.ae4.transported_charge_fmol_s)
        ) if genotype.ae4_expression == 0 else True
    solution = {'background': background.alias, 'root_id': background.root_id,
        'label': label, 'genotype': genotype.name, 'ae4_expression': genotype.ae4_expression,
        'carrier_amount_fmol': carrier_amount_fmol, 'ae4_activity': ae4_activity,
        'ae4_carrier_amount_fmol': model.ae4_parameters.carrier_amount_fmol,
        'optimizer_success': bool(fit.success), 'optimizer_message': str(fit.message),
        'optimizer_nfev': fit.nfev, 'actual_residual_evaluations': local,
        'cost': float(fit.cost), 'optimality': float(fit.optimality),
        'coordinates': x.tolist(), 'core_state': state[:12].tolist(),
        'core_state_sha256': sha256_object(state[:12].tolist()), 'raw_rhs': raw.tolist(),
        'scaled_independent_rhs': scaled.tolist(),
        'max_abs_scaled_independent_rhs': float(np.max(np.abs(scaled))),
        'max_abs_amount_rhs_fmol_s': float(np.max(np.abs(raw[AMOUNT_ROWS]))),
        'max_abs_volume_rhs_pL_s': float(np.max(np.abs(raw[VOLUME_ROWS]))),
        'max_abs_omitted_rhs_fmol_s': float(np.max(np.abs(raw[OMITTED_CHARGE_ROWS]))),
        'max_abs_regulatory_rhs_s_inv': regulatory, 'max_abs_bulk_charge_fmol': charge,
        'max_abs_current_residual_A': current, 'max_dimensionless_conservation_ratio': ratio,
        'conservation_residuals': conservation, 'normalized_jacobian_singular_values': singular.tolist(),
        'independent_jacobian_rank': rank, 'boundary_hits': boundary_hits,
        'rest_protocol_time_invariant_rhs': bool(np.array_equal(raw, later.rhs)),
        'positivity_pass': positive, 'full_rhs_pass': full_rhs, 'numerical_pass': numerical,
        'conservation_pass': conservation_pass, 'finite_capacity': capacity,
        'exact_ae4_zero': exact_zero, 'observables': obs, 'fluxes': flux_row(e, model),
        'active_parameter_sha256': sha256_object(model.parameters),
        'ae4_parameter_sha256': sha256_object(model.ae4_parameters)}
    solution['admissible'] = acceptance(solution)
    solution['wt_calibration_pass'] = acceptance(solution, require_wt_intervals=True) if genotype.ae4_expression == 1 else None
    cases[key] = solution
    write_json(path, cases)
    budget.save()
    (RESULTS / 'in_progress.json').unlink()
    print(json.dumps({k: solution[k] for k in ('label', 'carrier_amount_fmol',
        'ae4_activity', 'ae4_expression', 'actual_residual_evaluations',
        'max_abs_scaled_independent_rhs', 'numerical_pass', 'admissible', 'observables')}), flush=True)
    return solution


def bracketed_stage(anchor, evaluate, *, scalar, target, interval, bounds, record):
    """One sign-bracketed scalar solve in log activity; never a parameter grid.

    Bracket discovery takes one factor-ten directional step, then extends only
    in the direction inferred from the two closed WT endpoints. Interpolation
    is permitted only after a sign bracket exists. Stop at the measured interval,
    using its centre as the interpolation target. An unclosed endpoint has no sign.
    """
    field = 'ae4_activity' if scalar == 'ae4' else 'carrier_amount_fmol'
    observable = 'cl_i_mM' if scalar == 'ae4' else 'ph_i'
    history = [anchor]
    record.update({'scalar': scalar, 'observable': observable, 'target': target,
        'acceptance_interval': interval, 'bounds': bounds, 'anchor': anchor['label'],
        'evaluations': [], 'sign_brackets': []})
    first_direction = -1 if anchor['observables'][observable] > target else 1
    while True:
        last = history[-1]
        if not last['numerical_pass']:
            raise BoundedStop('unclosed WT scalar endpoint; no valid residual sign for '+record['stage'])
        value = last['observables'][observable]
        if interval[0] <= value <= interval[1]:
            record['accepted_label'] = last['label']
            return last
        points = sorted((math.log(r[field]), r['observables'][observable]-target, r)
                        for r in history)
        brackets = [(a, b) for a, b in zip(points, points[1:]) if a[1]*b[1] < 0]
        if brackets:
            a, b = min(brackets, key=lambda pair: pair[1][0]-pair[0][0])
            record['sign_brackets'].append([a[2]['label'], b[2]['label']])
            u = (a[0]*b[1]-b[0]*a[1])/(b[1]-a[1])
            if not a[0] + 0.01*(b[0]-a[0]) < u < b[0] - 0.01*(b[0]-a[0]):
                u = (a[0]+b[0])/2
        elif len(history) == 1:
            u = points[0][0] + first_direction*math.log(10)
        else:
            a, b = points[0], points[-1]
            slope = (b[1]-a[1])/(b[0]-a[0])
            if abs(slope) < 1e-12:
                raise BoundedStop('no resolvable WT scalar sensitivity for '+record['stage'])
            direction = -math.copysign(1, last['observables'][observable]-target)*math.copysign(1, slope)
            endpoint = a[0] if direction < 0 else b[0]
            u = endpoint + direction*math.log(10)
        trial = float(np.clip(math.exp(u), *bounds))
        if any(math.isclose(trial, r[field], rel_tol=1e-13, abs_tol=0) for r in history):
            raise BoundedStop('finite positive scalar bound exhausted without a bracket for '+record['stage'])
        seed = min(history, key=lambda r: abs(math.log(r[field]/trial)))
        result = evaluate(trial, seed)
        history.append(result)
        record['evaluations'].append(result['label'])


def cached_wt(background):
    directory = ROOT / 'results/31_nhe1_mechanistic_repair'
    freeze = json.loads((directory/'density_freeze.json').read_text())
    calibration = json.loads((directory/'calibration.json').read_text())
    selected = [r for r in calibration['evaluations'] if r['background'] == 'R09'
        and r['ae4_expression'] == 1 and r['carrier_amount_fmol'] == freeze['carrier_amount_fmol']]
    if len(selected) != 1 or not selected[0]['admissible']:
        raise ValueError('missing unique accepted Task 31 WT cache')
    r = dict(selected[0])
    model = build_model(background, carrier_amount_fmol=freeze['carrier_amount_fmol'], ae4_activity=1.)
    if (r['active_parameter_sha256'] != sha256_object(model.parameters)
        or r['ae4_parameter_sha256'] != sha256_object(model.ae4_parameters)
        or r['core_state_sha256'] != sha256_object(r['core_state'])
        or freeze['kinetics_sha256'] != sha256_object(published_cha_kinetics())):
        raise ValueError('Task 31 WT cache provenance mismatch')
    r.update(label='Task31_cached_R09_WT', ae4_activity=1.,
        ae4_carrier_amount_fmol=model.ae4_parameters.carrier_amount_fmol,
        reused_cached_state=True, task32_residual_evaluations=0,
        task31_residual_evaluations=r['actual_residual_evaluations'])
    r['actual_residual_evaluations'] = 0
    r['admissible'] = acceptance(r)
    r['wt_calibration_pass'] = acceptance(r, require_wt_intervals=True)
    return r


def calibrate(background, budget, report):
    current = cached_wt(background)
    report['cached_start'] = current
    report['stages'] = []

    def stage(name, scalar, anchor, target, interval, bounds):
        record = {'stage': name}
        report['stages'].append(record)
        def evaluate(value, seed):
            amount = value if scalar == 'nhe1' else anchor['carrier_amount_fmol']
            activity = value if scalar == 'ae4' else anchor['ae4_activity']
            result = solve_stationary(background, carrier_amount_fmol=amount,
                ae4_activity=activity, genotype=WT, start_coordinates=seed['coordinates'],
                label=f'{name}_{len(record["evaluations"])+1}', budget=budget, scalar=scalar)
            report['last_candidate'] = result
            return result
        try:
            return bracketed_stage(anchor, evaluate, scalar=scalar, target=target,
                interval=interval, bounds=bounds, record=record)
        finally:
            write_json(RESULTS/'wt_calibration.json', report)

    current = stage('AE4_initial', 'ae4', current, 50.10, WT_CL_INTERVAL, AE4_BOUNDS)
    if not WT_PH_INTERVAL[0] <= current['observables']['ph_i'] <= WT_PH_INTERVAL[1]:
        current = stage('NHE1_once', 'nhe1', current, 6.91, (6.90999, 6.91001), NHE1_BOUNDS)
        # The same stationary state already supplies chloride: no duplicate solve.
        report['post_nhe1_chloride_check'] = current['observables']['cl_i_mM']
        if not WT_CL_INTERVAL[0] <= current['observables']['cl_i_mM'] <= WT_CL_INTERVAL[1]:
            current = stage('AE4_final_once', 'ae4', current, 50.10, WT_CL_INTERVAL, AE4_BOUNDS)
    else:
        report['nhe1_recalibration'] = 'skipped; pH already in WT interval'
    report['final_candidate'] = current
    report['accepted'] = acceptance(current, require_wt_intervals=True)
    if not report['accepted']:
        raise BoundedStop('permitted sequential procedure ended outside simultaneous WT acceptance')
    return current


def freeze_parameters(background, wt, budget):
    if not acceptance(wt, require_wt_intervals=True):
        raise BoundedStop('cannot freeze an unaccepted WT calibration')
    model = build_model(background, carrier_amount_fmol=wt['carrier_amount_fmol'],
                        ae4_activity=wt['ae4_activity'])
    record = {'status': 'accepted_WT_frozen_before_genotype_access',
        'prepared_head': PREPARED_HEAD, 'calibration_background': 'R09',
        'carrier_amount_fmol': wt['carrier_amount_fmol'], 'ae4_activity': wt['ae4_activity'],
        'ae4_carrier_amount_fmol': model.ae4_parameters.carrier_amount_fmol,
        'whole_cell_parameters': asdict(model.parameters),
        'ae4_parameters': asdict(model.ae4_parameters),
        'kinetics': asdict(published_cha_kinetics()),
        'accepted_wt_state_sha256': wt['core_state_sha256'], 'budget_at_freeze': asdict(budget)}
    path = RESULTS/'frozen_parameters.json'
    if path.exists() and json.loads(path.read_text()) != record:
        raise BoundedStop('immutable frozen parameter record already exists')
    write_json(path, record)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    (RESULTS/'frozen_parameters.sha256').write_text(digest+'  frozen_parameters.json\n')
    return record, digest


def export_tables(rows):
    rows = sorted(rows, key=lambda r: (not r.get('reused_cached_state',False), r['label']))
    _write_csv(RESULTS/'resting_states.csv', [{
        **{k:r.get(k) for k in ('label','background','genotype','ae4_expression',
            'carrier_amount_fmol','ae4_activity','ae4_carrier_amount_fmol',
            'numerical_pass','admissible','wt_calibration_pass','actual_residual_evaluations',
            'max_abs_scaled_independent_rhs','max_abs_amount_rhs_fmol_s',
            'max_abs_volume_rhs_pL_s','max_abs_current_residual_A',
            'max_abs_bulk_charge_fmol','max_dimensionless_conservation_ratio',
            'independent_jacobian_rank')}, **r['observables']} for r in rows])
    _write_csv(RESULTS/'flux_ledger.csv', [{'label':r['label'], 'background':r['background'],
        'genotype':r['genotype'], 'numerical_pass':r['numerical_pass'],
        'admissible':r['admissible'], **r['fluxes']} for r in rows])
    _write_csv(RESULTS/'chloride_allocation.csv', [allocation(r) for r in rows])


def contract():
    return {'repository':'esig626/ae4-salivary-transport-control', 'branch':BRANCH,
        'prepared_head':PREPARED_HEAD, 'scientific_source':'Task 31 frozen Cha eight-state/Mod2',
        'workers':1, 'blas_threads':1, 'limits':LIMITS, 'numerical_seconds_limit':1800,
        'calibration_data':{'background':'R09','genotype':'WT',
            'ph_centre':6.91,'ph_interval':WT_PH_INTERVAL,'cl_centre_mM':50.10,
            'cl_interval_mM':WT_CL_INTERVAL},
        'scalars':{'nhe1':'effective carrier amount (fmol)',
            'ae4':'dimensionless multiplier of inherited AE4 carrier_amount_fmol'},
        'genotype_expression':{'WT':1.,'exact_null':0.},
        'sequence':['AE4 with frozen Task31 NHE1','NHE1 only if pH outside interval',
            'one chloride check from same stationary state','at most one final AE4 correction',
            'simultaneous WT acceptance','write and hash freeze','exact-null evaluation'],
        'scalar_method':'directional factor-ten bracket discovery; safeguarded log false-position only inside a sign bracket',
        'ae4_bounds':AE4_BOUNDS,'nhe1_bounds_fmol':NHE1_BOUNDS,
        'ae4_stop':'measured WT chloride interval; centre is root target',
        'nhe1_stop':'6.91 +/- 1e-5 when a correction is required',
        'unclosed_scalar_endpoint':'stop; no phenotype residual sign assigned',
        'cached_anchor':'accepted Task31 R09 WT; no new solve or scalar count',
        'stationary_method':'inherited bounded least_squares in 10 physical state coordinates; no parameter optimiser',
        'freeze_gate':'WT intervals, full closure, numerical/conservation/current/rank/positive gates; Na<=30 mM, volume<3 pL',
        'null_gate':'same closure gates and ordinary Na/volume; no experimental pH or chloride gate',
        'null_continuation':'direct expression=0 first; if unclosed, bounded 0.5,0.25,0.125,0 from last closed state',
        'validation_only':['null chloride','null pH','secretion'],
        'frozen_laws':['Cha kinetics','AE4 topology and stoichiometry and donor routing',
            'AE2','NKCC1','pump','K','CaCC','CO2','paracellular','water','bath','hydraulic','regulatory'],
        'r10':'only after an accepted R09 matched pair; frozen R09 scalars; no recalibration',
        'dynamics':'only admissible matched R09 WT/exact-null; at most two 0-600 s CCh+IPR, Ca=0.25 uM integrations'}


def run_calibration():
    for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
        if os.environ.get(name) != '1':
            raise RuntimeError(name+' must equal 1 before Python starts')
    if (RESULTS/'execution_summary.json').exists():
        return json.loads((RESULTS/'execution_summary.json').read_text())
    if (RESULTS/'in_progress.json').exists():
        raise BoundedStop('interrupted case checkpoint exists; counters must not be reset or case silently rerun')
    old = json.loads((RESULTS/'budget.json').read_text()) if (RESULTS/'budget.json').exists() else {}
    budget = Budget(**{k:old.get(k,0) for k in Budget.__dataclass_fields__})
    write_json(RESULTS/'contract.json', contract())
    report = {'accepted':False, 'null_data_used':False, 'secretion_data_used':False}
    stop, wt, freeze, digest = None, None, None, None
    try:
        background = load_background('R09')
        wt = calibrate(background, budget, report)
        freeze, digest = freeze_parameters(background, wt, budget)
    except BoundedStop as exc:
        stop = str(exc)
        report['stop_reason'] = stop
        print('BOUNDED_STOP: '+stop, flush=True)
    finally:
        report['accepted'] = wt is not None
        write_json(RESULTS/'wt_calibration.json', report)
        if freeze is None:
            write_json(RESULTS/'frozen_parameters.json', {'status':'not_frozen_calibration_failed',
                'accepted_parameters':None,'exact_null_evaluated':False,'reason':stop})
        cases = json.loads((RESULTS/'cases.json').read_text()) if (RESULTS/'cases.json').exists() else {}
        rows = ([report['cached_start']] if 'cached_start' in report else [])+list(cases.values())
        export_tables(rows)
        budget.save()
        summary = {'calibration_accepted':wt is not None,'stop_reason':stop,
            'accepted_wt':wt,'frozen_parameters_sha256':digest,
            'validation_status':'pending_frozen_validation' if wt else 'blocked_by_failed_WT_calibration',
            'exact_null_evaluated':False,'stimulated_integrations':0,'r10_evaluated':False}
        write_json(RESULTS/'execution_summary.json', summary)
    return summary


if __name__ == '__main__':
    run_calibration()
