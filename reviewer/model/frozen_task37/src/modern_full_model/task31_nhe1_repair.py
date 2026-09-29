"""Task 31 continuation: one source-fixed Cha model and one WT density root.

No historical calculations, alternate starts, parameter grids or transporter
changes. Completed cases are immutable checkpoints keyed by all active inputs.
Run with OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=1.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import json
import math
import os
from pathlib import Path
import time
from types import SimpleNamespace
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import least_squares

from . import task30_nhe1_repair as inherited
from .task30_nhe1_repair import (
    Background, load_background, state_from_coordinates,
    COORDINATE_NAMES, COORDINATE_LOWER, COORDINATE_UPPER,
    INDEPENDENT_ROWS, INDEPENDENT_RHS_SCALES, OMITTED_CHARGE_ROWS,
    AMOUNT_ROWS, VOLUME_ROWS, ROOT_SCALED_TOLERANCE, OMITTED_ROW_TOLERANCE,
    CHARGE_TOLERANCE_FMOL, CURRENT_TOLERANCE_A, REGULATORY_TOLERANCE_S_INV,
    BOUNDARY_RELATIVE_TOLERANCE, JACOBIAN_RANK_RELATIVE_TOLERANCE,
    ORDINARY_NA_UPPER_MM, CELL_VOLUME_UPPER_PL, WT_TARGET_PH,
    EXACT_NULL_DIAGNOSTIC_PH, _conservation_audit, _rest_observables,
    _flux_row, _write_json, _write_csv,
)
from .nhe1_cha2009 import NHE1_CHA_2009, published_cha_kinetics, AVOGADRO_PER_MOL
from .model import AE4_NULL, WT, Genotype, ModernFullModel
from .camp_pka import AE4Construct, R1EffectiveActivation, RegulatoryGain
from .ae4_routing_only import routing_only_adapter
from .nkcc_stimulation import N1AlgebraicNkcc1, attach_stimulated_nkcc1
from .validation import sha256_object, PRODUCTION_RADAU

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results/31_nhe1_mechanistic_repair"
PER_SOLVE_RESIDUAL_CEILING = 2400
DENSITY_BOUNDS_FMOL = (1e-8, 0.1)
INITIAL_DENSITY_FMOL = 1e-4  # Numerical starting value, not source N.
CALIBRATION_PH_TOLERANCE = 1e-5


class BudgetStop(RuntimeError):
    pass


class LocalSolveStop(RuntimeError):
    pass


@dataclass
class Budget:
    stationary_solver_calls: int = 0
    stationary_residual_evaluations: int = 0
    scalar_density_evaluations: int = 0
    stimulated_integrations: int = 0
    numerical_wall_seconds: float = 0.0

    def __post_init__(self):
        self._clock = time.perf_counter()

    def check_time(self):
        now = time.perf_counter()
        self.numerical_wall_seconds += now - self._clock
        self._clock = now
        if self.numerical_wall_seconds >= 1800:
            raise BudgetStop("1800 s numerical execution ceiling")

    def count(self, field, limit):
        self.check_time()
        if getattr(self, field) >= limit:
            raise BudgetStop(field + " ceiling")
        setattr(self, field, getattr(self, field) + 1)

    def count_stationary_residual(self):
        self.count("stationary_residual_evaluations", 24000)

    def save(self):
        # Persist the stop even when the wall-time guard itself was reached.
        now = time.perf_counter()
        self.numerical_wall_seconds += now - self._clock
        self._clock = now
        _write_json(RESULTS / "budget.json", {
            **asdict(self), "scientific_workers": 1, "blas_threads": 1,
            "per_solve_residual_ceiling": PER_SOLVE_RESIDUAL_CEILING,
            "accepted_identical_cases_rerun": 0,
            "limits_respected": (self.stationary_solver_calls <= 16 and
                self.stationary_residual_evaluations <= 24000 and
                self.scalar_density_evaluations <= 6 and
                self.stimulated_integrations <= 4 and self.numerical_wall_seconds <= 1800),
        })


def task31_parameters(background, *, carrier_amount_fmol):
    base = background.inherited_parameters
    return replace(base, homeostasis=replace(base.homeostasis,
        nhe1_model=NHE1_CHA_2009,
        nhe1_cha_carrier_amount_fmol=float(carrier_amount_fmol)))


def build_model(background, *, carrier_amount_fmol, stimulated=False):
    regulation = R1EffectiveActivation(tau_activation_s=30.0,
        gain=RegulatoryGain(basal_capacity_multiplier=1.0,
            fully_activated_increment=0.25, coupling_scale=1.0),
        construct=AE4Construct.WT)
    base = ModernFullModel(task31_parameters(background,
        carrier_amount_fmol=carrier_amount_fmol),
        stimulus=inherited._protocol(stimulated=stimulated),
        regulatory_model=regulation, ae4_parameters=background.ae4_parameters,
        ae4_evaluator=routing_only_adapter)
    return attach_stimulated_nkcc1(base, N1AlgebraicNkcc1(
        fully_activated_multiplier=1.75, resting_calcium_uM=0.058,
        stimulated_calcium_uM=0.10))  # Unchanged inherited NKCC activation law.


def _finite_capacity_audit(evaluation, model):
    audit = inherited._finite_capacity_audit(evaluation, model)
    capacity = (model.parameters.homeostasis.nhe1_cha_carrier_amount_fmol
                * published_cha_kinetics().turnover_bound_per_s)
    audit["capacities"]["nhe1_fmol_s"] = capacity
    realised = audit["realised_absolute_fluxes"]["nhe1_fmol_s"]
    audit["utilisation"]["nhe1_fmol_s"] = realised/capacity if capacity else (0 if realised == 0 else math.inf)
    audit["pass"] = max(audit["utilisation"].values()) <= 1 + 1e-12
    return audit


def solve_stationary(model, *, genotype, start_coordinates, background,
                     label, budget):
    """One direct bounded solve; actual finite-difference calls are counted."""
    key = sha256_object({"background": background.alias,
        "parameters": model.parameters, "ae4_parameters": background.ae4_parameters,
        "genotype": asdict(genotype), "kinetics": published_cha_kinetics()})
    path = RESULTS / "rest_checkpoints" / (key + ".json")
    if path.exists():
        return json.loads(path.read_text())
    budget.count("stationary_solver_calls", 16)
    start = np.asarray(start_coordinates, dtype=float)
    if np.any(start <= COORDINATE_LOWER) or np.any(start >= COORDINATE_UPPER):
        raise ValueError("start must be strictly inside inherited coordinate bounds")
    local_residual_evaluations = 0
    best_x, best_cost = start.copy(), math.inf

    def counted(candidate, *, audit=False):
        nonlocal local_residual_evaluations, best_x, best_cost
        # Reserve the final forward-difference Jacobian and two gate RHS calls.
        ceiling = PER_SOLVE_RESIDUAL_CEILING if audit else PER_SOLVE_RESIDUAL_CEILING - 13
        if local_residual_evaluations >= ceiling:
            raise LocalSolveStop("per-solve actual residual ceiling")
        budget.count_stationary_residual()
        local_residual_evaluations += 1
        try:
            state = state_from_coordinates(model, candidate)
            raw = model.rhs(0.0, state, genotype=genotype)
            values = raw[INDEPENDENT_ROWS] / INDEPENDENT_RHS_SCALES
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
            values = np.full(10, 1e6)
        cost = float(values @ values)/2
        if cost < best_cost:
            best_x, best_cost = np.array(candidate), cost
        return values

    try:
        fit = least_squares(counted, start,
            bounds=(COORDINATE_LOWER, COORDINATE_UPPER), max_nfev=2400,
            xtol=1e-12, ftol=1e-12, gtol=1e-12, x_scale="jac")
    except LocalSolveStop as exc:
        final_x = best_x.copy()
        f0 = counted(final_x, audit=True)
        jac = np.empty((10, 10))
        for j in range(10):
            dx = math.sqrt(np.finfo(float).eps)*max(abs(final_x[j]), 1.0)
            if final_x[j]+dx >= COORDINATE_UPPER[j]:
                dx = -dx
            xj = final_x.copy()
            xj[j] += dx
            jac[:, j] = (counted(xj, audit=True)-f0)/dx
        fit = SimpleNamespace(x=final_x, jac=jac, success=False,
            message=str(exc), nfev=None, cost=float(f0@f0)/2,
            optimality=float(np.max(np.abs(jac.T@f0))))
    coordinates = np.asarray(fit.x, dtype=float)
    state = state_from_coordinates(model, coordinates)
    budget.count_stationary_residual()
    local_residual_evaluations += 1
    evaluation = model.evaluate(0.0, state, genotype=genotype)
    budget.count_stationary_residual()
    local_residual_evaluations += 1
    later = model.evaluate(600.0, state, genotype=genotype)
    raw = np.asarray(evaluation.rhs, dtype=float)
    scaled = raw[INDEPENDENT_ROWS] / INDEPENDENT_RHS_SCALES
    span = COORDINATE_UPPER - COORDINATE_LOWER
    jacobian = np.asarray(fit.jac, dtype=float) * span[None, :]
    singular = np.linalg.svd(jacobian, compute_uv=False)
    threshold = JACOBIAN_RANK_RELATIVE_TOLERANCE * singular[0] if singular.size else math.inf
    rank = int(np.sum(singular > threshold))
    relative_boundary_distance = np.minimum(
        coordinates - COORDINATE_LOWER, COORDINATE_UPPER - coordinates
    ) / span
    boundary_hits = [
        name
        for name, distance in zip(COORDINATE_NAMES, relative_boundary_distance)
        if distance <= BOUNDARY_RELATIVE_TOLERANCE
    ]
    charge_max = max(abs(float(value)) for value in evaluation.diagnostics.state_charge_fmol.values())
    current_max = max(abs(float(value)) for value in evaluation.diagnostics.membranes.current_residuals_A.values())
    regulatory_max = float(np.max(np.abs(raw[12:]))) if raw.size > 12 else 0.0
    conservation, conservation_ratio, conservation_pass = _conservation_audit(evaluation)
    finite_capacity = _finite_capacity_audit(evaluation, model)
    observables = _rest_observables(evaluation, state)
    full_rhs_pass = bool(
        np.max(np.abs(scaled)) <= ROOT_SCALED_TOLERANCE
        and np.max(np.abs(raw[OMITTED_CHARGE_ROWS])) <= OMITTED_ROW_TOLERANCE
        and regulatory_max <= REGULATORY_TOLERANCE_S_INV
    )
    numerical_pass = bool(
        fit.success
        and full_rhs_pass
        and charge_max <= CHARGE_TOLERANCE_FMOL
        and current_max <= CURRENT_TOLERANCE_A
        and conservation_pass
        and rank == 10
        and not boundary_hits
        and np.all(state[:12] > 0.0)
        and np.array_equal(raw, later.rhs)
    )
    ordinary_na_pass = bool(observables["na_i_mM"] <= ORDINARY_NA_UPPER_MM)
    ph_diagnostic_centre = (
        EXACT_NULL_DIAGNOSTIC_PH
        if genotype.ae4_expression == 0.0
        else WT_TARGET_PH
    )
    ph_diagnostic_pass = bool(
        abs(observables["ph_i"] - ph_diagnostic_centre)
        <= (0.02 if genotype.ae4_expression == 0.0 else 0.07)
    )
    volume_pass = bool(observables["volume_i_pL"] <= CELL_VOLUME_UPPER_PL)
    potassium_pass = bool(60.0 <= observables["k_i_mM"] <= 210.0)
    # The inherited 100 mM TIC boundary is a broad diagnostic, not a fit target.
    carbon_sanity_pass = bool(observables["tic_i_mM"] < 100.0 and observables["hco3_i_mM"] < 100.0)
    physiological_pass = bool(ordinary_na_pass and ph_diagnostic_pass and volume_pass
                              and potassium_pass and carbon_sanity_pass)
    exact_ae4_zero = all(
        float(value) == 0.0
        for value in (
            evaluation.diagnostics.ae4.na_cell_fmol_s,
            evaluation.diagnostics.ae4.k_cell_fmol_s,
            evaluation.diagnostics.ae4.cl_cell_fmol_s,
            evaluation.diagnostics.ae4.hco3_cell_fmol_s,
            evaluation.diagnostics.ae4.transported_charge_fmol_s,
        )
    ) if genotype.ae4_expression == 0.0 else True
    admissible = bool(
        numerical_pass
        and physiological_pass
        and finite_capacity["pass"]
        and exact_ae4_zero
    )
    solution = {
        "background": background.alias,
        "root_id": background.root_id,
        "routing_family": background.routing_family,
        "genotype": genotype.name,
        "ae4_expression": float(genotype.ae4_expression),
        "label": label,
        "carrier_amount_fmol": float(model.parameters.homeostasis.nhe1_cha_carrier_amount_fmol),
        "nhe1_model": model.parameters.homeostasis.nhe1_model,
        "optimizer_success": bool(fit.success),
        "optimizer_message": str(fit.message),
        "optimizer_nfev": int(fit.nfev) if fit.nfev is not None else None,
        "actual_residual_evaluations": int(local_residual_evaluations),
        "cost": float(fit.cost),
        "optimality": float(fit.optimality),
        "coordinates": coordinates.tolist(),
        "core_state": state[:12].tolist(),
        "core_state_sha256": sha256_object(state[:12].tolist()),
        "raw_rhs": raw.tolist(),
        "scaled_independent_rhs": scaled.tolist(),
        "max_abs_scaled_independent_rhs": float(np.max(np.abs(scaled))),
        "max_abs_amount_rhs_fmol_s": float(np.max(np.abs(raw[AMOUNT_ROWS]))),
        "max_abs_volume_rhs_pL_s": float(np.max(np.abs(raw[VOLUME_ROWS]))),
        "max_abs_omitted_rhs_fmol_s": float(np.max(np.abs(raw[OMITTED_CHARGE_ROWS]))),
        "max_abs_regulatory_rhs_s_inv": regulatory_max,
        "max_abs_bulk_charge_fmol": float(charge_max),
        "max_abs_current_residual_A": float(current_max),
        "max_dimensionless_conservation_ratio": conservation_ratio,
        "conservation_residuals": conservation,
        "normalized_jacobian_singular_values": singular.tolist(),
        "independent_jacobian_rank": rank,
        "independent_jacobian_nullity": 10 - rank,
        "boundary_hits": boundary_hits,
        "rest_protocol_time_invariant_rhs": bool(np.array_equal(raw, later.rhs)),
        "positivity_pass": bool(np.all(state[:12] > 0.0)),
        "full_rhs_pass": full_rhs_pass,
        "conservation_pass": conservation_pass,
        "numerical_pass": numerical_pass,
        "ordinary_na_pass": ordinary_na_pass,
        "ph_diagnostic_centre": ph_diagnostic_centre,
        "ph_diagnostic_pass": ph_diagnostic_pass,
        "volume_pass": volume_pass,
        "potassium_pass": potassium_pass,
        "carbon_sanity_pass": carbon_sanity_pass,
        "chloride_context_pass": abs(observables["cl_i_mM"] - (36.50 if genotype.ae4_expression == 0.0 else 50.10)) <= (1.60 if genotype.ae4_expression == 0.0 else 1.50),
        "finite_capacity": finite_capacity,
        "exact_ae4_zero": exact_ae4_zero,
        "physiological_pass": physiological_pass,
        "admissible": admissible,
        "observables": observables,
        "fluxes": _flux_row(evaluation, model),
        "active_parameter_sha256": sha256_object(model.parameters),
        "inherited_parameter_sha256": background.inherited_whole_cell_sha256,
        "ae4_parameter_sha256": background.inherited_ae4_sha256,
    }

    _write_json(path, solution)
    budget.save()
    print(json.dumps({k: solution[k] for k in ("background", "genotype",
        "label", "carrier_amount_fmol", "actual_residual_evaluations",
        "max_abs_scaled_independent_rhs", "numerical_pass", "admissible",
        "observables")}), flush=True)
    return solution


def calibrate_r09(background, budget):
    path = RESULTS / 'calibration.json'
    history = json.loads(path.read_text())['evaluations'] if path.exists() else []
    while len(history) < 6:
        if not history:
            amount, start = INITIAL_DENSITY_FMOL, background.saved_coordinates
        else:
            last = history[-1]
            if not last['numerical_pass']:
                raise BudgetStop('R09 density evaluation did not produce a numerical REST root')
            if abs(last['observables']['ph_i'] - WT_TARGET_PH) <= CALIBRATION_PH_TOLERANCE:
                return last, history
            start = last['coordinates']
            u1 = math.log(last['carrier_amount_fmol'])
            f1 = last['observables']['ph_i'] - WT_TARGET_PH
            if len(history) == 1:
                u = u1 + math.copysign(math.log(1.25), -f1)
            else:
                prev = history[-2]
                u0 = math.log(prev['carrier_amount_fmol'])
                f0 = prev['observables']['ph_i'] - WT_TARGET_PH
                if abs(f1-f0) < 1e-12:
                    raise BudgetStop('WT scalar secant has no resolvable pH sensitivity')
                u = u1 - f1*(u1-u0)/(f1-f0)
                # Safeguarded secant; if a sign bracket exists stay inside it.
                bracket = [(math.log(r['carrier_amount_fmol']), r['observables']['ph_i']-WT_TARGET_PH)
                           for r in history]
                opposite = [v for v, f in bracket if f*f1 < 0]
                if opposite:
                    other = min(opposite, key=lambda v: abs(v-u1))
                    lo, hi = sorted((u1, other))
                    if not lo < u < hi:
                        u = (lo+hi)/2
                else:
                    u = float(np.clip(u, u1-math.log(4), u1+math.log(4)))
            amount = float(np.clip(math.exp(u), *DENSITY_BOUNDS_FMOL))
            if any(amount == r['carrier_amount_fmol'] for r in history):
                raise BudgetStop('WT scalar iteration would repeat an evaluated density')
        budget.count('scalar_density_evaluations', 6)
        solution = solve_stationary(build_model(background, carrier_amount_fmol=amount),
            genotype=WT, start_coordinates=start, background=background,
            label=f'R09_WT_density_{len(history)+1}', budget=budget)
        history.append(solution)
        _write_json(path, {'target_ph': WT_TARGET_PH, 'ph_tolerance': CALIBRATION_PH_TOLERANCE,
            'method': 'bounded safeguarded secant in log carrier amount',
            'bounds_fmol': DENSITY_BOUNDS_FMOL,
            'initial_value_note': 'numerical starting amount; no cardiac density imported',
            'evaluations': history})
        if solution['numerical_pass'] and abs(solution['observables']['ph_i']-WT_TARGET_PH) <= CALIBRATION_PH_TOLERANCE:
            return solution, history
    raise BudgetStop('six WT-only scalar evaluations exhausted without the pH root')


def stimulate(background, rest, genotype, budget):
    from scipy.integrate import solve_ivp
    amount = rest['carrier_amount_fmol']
    path = RESULTS / 'dynamic_checkpoints' / f'{background.alias}_{genotype.name}.json'
    if path.exists():
        return json.loads(path.read_text())
    model = build_model(background, carrier_amount_fmol=amount, stimulated=True)
    rest_model = build_model(background, carrier_amount_fmol=amount)
    y0 = np.r_[rest['core_state'], model.initial_state()[12:]]
    onset = model.evaluate(0, y0, genotype=genotype)
    onset_match = np.array_equal(onset.rhs, rest_model.rhs(0, y0, genotype=genotype))
    if not onset_match:
        raise BudgetStop('stimulus onset does not match REST')
    budget.count('stimulated_integrations', 4)
    def rhs(t, y):
        budget.check_time()
        return model.rhs(t, y, genotype=genotype)
    grid = np.r_[1e-6, np.arange(1., 601.)]
    sol = solve_ivp(rhs, (1e-6, 600), y0, method=PRODUCTION_RADAU.method,
        rtol=PRODUCTION_RADAU.rtol, atol=PRODUCTION_RADAU.atol_vector(model.state_names),
        t_eval=grid, max_step=PRODUCTION_RADAU.max_step_s)
    times = np.r_[0., sol.t]
    states = np.column_stack((y0, sol.y))
    rows, ratios, capacities = [], [], []
    for t, state in zip(times, states.T):
        budget.check_time()
        e = model.evaluate(float(t), state, genotype=genotype)
        rows.append({'time_s': float(t), **_rest_observables(e, state), **_flux_row(e, model)})
        ratios.append(_conservation_audit(e)[1])
        capacities.append(_finite_capacity_audit(e, model)['pass'])
    integrated = {key: float(np.trapezoid([r[key] for r in rows], times)) for key in (
        'water_lumen_outflow_pL_s', 'nhe1_na_in_h_out_fmol_s', 'ae2_cl_in_hco3_out_fmol_s',
        'ae4_cl_cell_fmol_s', 'nkcc1_cl_in_fmol_s', 'cacc_cl_cell_to_lumen_fmol_s')}
    summary = {'background': background.alias, 'genotype': genotype.name,
        'carrier_amount_fmol': amount, 'calcium_uM': 0.25, 'protocol': '600_s_CCH_IPR',
        'solver': PRODUCTION_RADAU.label, 'solver_success': bool(sol.success),
        'solver_message': sol.message, 'onset_matches_rest_rhs': bool(onset_match),
        'complete': bool(sol.t[-1] == 600), 'integrated_fluxes': integrated,
        'integrated_flux_units': 'pL for water; fmol for ions',
        'cumulative_secretion_pL': integrated['water_lumen_outflow_pL_s'],
        'max_dimensionless_conservation_ratio': max(ratios),
        'capacity_pass': all(capacities), 'positive_core': bool(np.all(states[:12] > 0)),
        'finite': bool(np.all(np.isfinite(states))),
        'ranges': {key: {'min': min(r[key] for r in rows), 'max': max(r[key] for r in rows)}
            for key in rows[0] if key not in ('nhe1_model', 'time_s')},
        'endpoint': rows[-1]}
    summary['valid'] = bool(summary['solver_success'] and summary['complete'] and
        summary['positive_core'] and summary['finite'] and summary['capacity_pass'] and
        max(ratios) <= 1 and min(r['water_lumen_outflow_pL_s'] for r in rows) >= 0)
    # Full 1-second samples are used for integration/ranges; compact 10-second
    # text trajectories retain all requested observables and fluxes.
    compact = [r for r in rows if r['time_s'] == 0 or r['time_s'] % 10 == 0 or r is rows[-1]]
    _write_csv(RESULTS / f'trajectory_{background.alias}_{genotype.name}.csv', compact)
    _write_json(path, summary)
    budget.save()
    return summary


def export_rest(rest):
    _write_csv(RESULTS / 'resting_states.csv', [
        {**inherited._flatten_rest(r), 'carrier_amount_fmol': r['carrier_amount_fmol'],
         'chloride_context_pass': r['chloride_context_pass'],
         'potassium_pass': r['potassium_pass'], 'carbon_sanity_pass': r['carbon_sanity_pass']}
        for r in rest])
    _write_csv(RESULTS / 'flux_ledger.csv', [inherited._flatten_flux(r) for r in rest])


def main():
    for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
        if os.environ.get(name) != '1':
            raise RuntimeError(name + ' must equal 1 before Python starts')
    old = json.loads((RESULTS / 'budget.json').read_text())
    budget = Budget(**{k: old.get(k, 0) for k in Budget.__dataclass_fields__})
    rest, dynamics, stop = [], [], None
    try:
        backgrounds = {alias: load_background(alias) for alias in ('R09', 'R10')}
        wt09, history = calibrate_r09(backgrounds['R09'], budget)
        amount = wt09['carrier_amount_fmol']
        rest.append(wt09)
        freeze_path = RESULTS / 'density_freeze.json'
        freeze = {'carrier_amount_fmol': amount,
            'carrier_count': amount*AVOGADRO_PER_MOL/1e15,
            'relative_to_source_N_489900': amount*AVOGADRO_PER_MOL/1e15/489900,
            'kinetics': asdict(published_cha_kinetics()),
            'kinetics_sha256': sha256_object(published_cha_kinetics()),
            'calibration_background': 'R09', 'genotype': 'WT',
            'calibration_ph': wt09['observables']['ph_i'],
            'frozen_before_R10_and_exact_null': True}
        if freeze_path.exists() and json.loads(freeze_path.read_text()) != freeze:
            raise BudgetStop('attempt to change the frozen density')
        _write_json(freeze_path, freeze)
        export_rest(rest)
        wt10 = solve_stationary(build_model(backgrounds['R10'], carrier_amount_fmol=amount),
            genotype=WT, start_coordinates=backgrounds['R10'].saved_coordinates,
            background=backgrounds['R10'], label='R10_WT_frozen_density', budget=budget)
        rest.append(wt10)
        export_rest(rest)
        # The two predeclared null cases are solved exactly once from matched WT.
        for wt in (wt09, wt10):
            b = backgrounds[wt['background']]
            null = solve_stationary(build_model(b, carrier_amount_fmol=amount),
                genotype=AE4_NULL, start_coordinates=wt['coordinates'],
                background=b, label=f'{b.alias}_exact_AE4_null', budget=budget)
            rest.append(null)
            export_rest(rest)
        for alias, b in backgrounds.items():
            pair = [r for r in rest if r['background'] == alias]
            if len(pair) == 2 and all(r['admissible'] for r in pair):
                for r in pair:
                    genotype = WT if r['ae4_expression'] == 1 else AE4_NULL
                    dynamics.append(stimulate(b, r, genotype, budget))
        if dynamics:
            flat = [{k: v for k, v in d.items() if not isinstance(v, dict)} for d in dynamics]
            for alias in backgrounds:
                pair = [d for d in flat if d['background'] == alias]
                if len(pair) == 2:
                    wt = next(d for d in pair if d['genotype'] == WT.name)
                    null = next(d for d in pair if d['genotype'] == AE4_NULL.name)
                    null['null_to_wt_secretion_ratio'] = null['cumulative_secretion_pL']/wt['cumulative_secretion_pL']
            _write_csv(RESULTS / 'trajectory_summary.csv', flat)
    except BudgetStop as exc:
        stop = str(exc)
        print('BOUNDED_STOP: ' + stop, flush=True)
    finally:
        # A publication problem never calls this numerical runner.
        budget.save()
        _write_json(RESULTS / 'execution_summary.json', {'stop_reason': stop,
            'selected_rest': rest, 'dynamics': dynamics,
            'matched_rest_admissible': {a: len([r for r in rest if r['background'] == a]) == 2
                and all(r['admissible'] for r in rest if r['background'] == a) for a in ('R09', 'R10')}})


if __name__ == '__main__':
    main()
