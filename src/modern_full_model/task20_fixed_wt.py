"""WT-only Task 20 preparation, structural feasibility, inverse solve, freeze.

No genotype runner or result is imported/read. An affinity sign certificate
precedes numerical minimization. An empty feasible set has no canonical
capacity optimum; no signed/unphysical capacity is passed to production.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from decimal import Decimal, localcontext
import hashlib
import io
import math
from pathlib import Path
import subprocess
import unittest

import numpy as np

from . import task18_fixed_wt as inverse
from .pooled_ae4 import MODEL_ID, PooledAE4Environment, evaluate_pooled_ae4, select_ae4_mechanism
from .run_calcium_fast_screen import load_freeze, write_json, write_rows
from .validation import sha256_file, sha256_object


REPO = Path(__file__).resolve().parents[2]
OUT = REPO / 'results/20_ae4_pooled_cation'
ANALYSIS = REPO / 'analysis/20_ae4_pooled_cation'
BRANCH = 'codex/task-20-ae4-pooled-cation-no-slip'
INTAKE = 'e45144492d682b0d2c3b581ef3b1b51c73e7397e'
INHERITANCE = 'e98789b276607807f7c9c88ff53605008a19cb66'
CONDITIONS = ('P0', 'P10')
CALCIUM = (0.10, 0.25, 0.50)
INFEASIBLE = 'STRUCTURALLY_INFEASIBLE_POOLED_AFFINITY'


def read(path):
    return inverse.read(path)


def relative(path):
    return str(Path(path).relative_to(REPO))


def git(*args):
    return subprocess.check_output(['git', *args], cwd=REPO, text=True).strip()


def ensure_unfrozen():
    if (OUT / 'frozen_wt_manifest.json').exists():
        raise FileExistsError('Task 20 WT freeze is immutable')
    if git('branch', '--show-current') != BRANCH:
        raise AssertionError('Task 20 may run only on its requested branch')


def environment_from_evaluation(model, evaluation):
    ci = evaluation.diagnostics.observables.cell_concentrations_mM
    bath = model.parameters.bath
    return PooledAE4Environment(ci['na'], ci['k'], ci['cl'], ci['hco3'],
                                bath.na_mM, bath.k_mM, bath.cl_mM,
                                evaluation.diagnostics.observables.bath_acid_base.hco3_mM)


def independent_thermodynamics(environment):
    """50-digit mass-action calculation independent of the implemented log sum."""
    with localcontext() as context:
        context.prec = 50
        d = {k: Decimal.from_float(float(v)) for k, v in asdict(environment).items()}
        ci = d['na_i_mM'] + d['k_i_mM']
        co = d['na_o_mM'] + d['k_o_mM']
        ratio = (d['cl_o_mM'] * ci * d['hco3_i_mM']**2 /
                 (d['cl_i_mM'] * co * d['hco3_o_mM']**2))
        hco3_reversal = d['hco3_o_mM'] * (d['cl_i_mM']*co/(d['cl_o_mM']*ci)).sqrt()
        return {'activity_ratio_decimal50': str(ratio),
                'affinity_decimal50': str(ratio.ln()),
                'independent_affinity': float(ratio.ln()),
                'hco3_i_at_zero_affinity_mM': float(hco3_reversal),
                'cl_i_at_zero_affinity_mM': float(d['cl_i_mM'] * ratio)}


def assay_tests():
    ensure_unfrozen()
    stream = io.StringIO()
    suite = unittest.defaultTestLoader.discover(str(REPO/'tests'), pattern='test_task20_pooled_ae4.py')
    cases = [case.id() for group in suite for subgroup in group for case in subgroup]
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/'assay_test_results.txt').write_text(stream.getvalue())
    write_json(OUT/'assay_tests.json', {
        'run_utc': inverse.now(), 'tests_run': result.testsRun,
        'passed': result.wasSuccessful(), 'failures': len(result.failures),
        'errors': len(result.errors), 'skipped': len(result.skipped), 'test_ids': cases,
        'phenotype_independent': True,
        'test_source_sha256': sha256_file(REPO/'tests/test_task20_pooled_ae4.py'),
        'mechanism_source_sha256': sha256_file(REPO/'src/modern_full_model/pooled_ae4.py'),
    })
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    print(stream.getvalue())


def prepare():
    ensure_unfrozen()
    if (OUT/'contract.json').exists():
        raise FileExistsError('WT contract already exists')
    assay = read(OUT/'assay_tests.json')
    if not assay['passed'] or assay['tests_run'] < 12:
        raise AssertionError('Required assays must pass before whole-cell work')
    for key, path in (('test_source_sha256', REPO/'tests/test_task20_pooled_ae4.py'),
                      ('mechanism_source_sha256', REPO/'src/modern_full_model/pooled_ae4.py')):
        if sha256_file(path) != assay[key]:
            raise AssertionError('Rerun assays after a source edit')
    manifest, verification = load_freeze()
    inputs = dict(verification['verified_input_sha256'])
    for root in manifest['roots']:
        path = REPO/'results/17_fixed_wt_flux_repartition/fixed_wt_payloads'/f'{inverse.slug(root)}__inherited_baseline.json'
        inputs[relative(path)] = sha256_file(path)
    for path in (REPO/'src/modern_full_model/task18_fixed_wt.py',
                 REPO/'analysis/13_state_resolved_ae4/evidence_freeze.md',
                 REPO/'analysis/12_ae4_mechanism_reconstruction/experimental_evidence.md'):
        inputs[relative(path)] = sha256_file(path)
    write_json(OUT/'contract.json', {
        'task': '20_AE4_POOLED_CATION_NO_SLIP', 'created_utc': inverse.now(),
        'branch': BRANCH, 'intake_commit': INTAKE, 'inheritance_commit': INHERITANCE,
        'model_id': MODEL_ID, 'conditions': CONDITIONS, 'calcium_uM': CALCIUM,
        'root_ids': sorted(manifest['roots']), 'weights': {'Na': 1., 'K': 1.},
        'state_rule': 'Exact inherited conserved amounts, both volumes and regulatory coordinates',
        'loading_rules': {'P0': 'Match inherited signed AE4 chloride source',
                          'P10': '0.10 of inherited positive basolateral chloride-loading pool'},
        'other_positive_loaders': 'Proportional scaling; inherited positive pool preserved',
        'canonical_objective': ['minimize largest absolute natural log capacity fold',
                                'then minimize sum of squared natural log capacity folds'],
        'objective_capacities': inverse.OBJECTIVE_CAPACITIES,
        'ae4_reference_gauge': 'attempt rate 1/s, inherited carrier as reference; unmeasured convention',
        'hard_capacity_rule': 'nonnegative, no inherited capacity boxes imposed',
        'feasibility_precheck': 'J_target*A must be nonnegative before optimization',
        'infeasibility_scope': 'Sign certificate covers every nonnegative AE4 capacity and every other capacity at the fixed state',
        'genotype_information_in_objective': False, 'genotype_information_read_before_freeze': False,
        'genotype_expression': {'AE4': 0.05, 'AE2': 0.0}, 'exact_ae4_zero_attempts': 0,
        'genotype_policy': 'Only frozen feasible WT candidates; production Radau, 600 s, no refitting',
        'no_feasible_candidates_policy': 'Explicit NOT_RUN rows, null outcomes, no substitute WT states',
        'prefreeze_input_hashes': inputs, 'intake_verification': verification,
    })


def solve():
    ensure_unfrozen()
    contract = read(OUT/'contract.json')
    for name, digest in contract['prefreeze_input_hashes'].items():
        if sha256_file(REPO/name) != digest:
            raise AssertionError(f'Inherited input changed: {name}')
    manifest, _ = load_freeze()
    rows, reference_rows = [], []
    for root_index, root in enumerate(sorted(manifest['roots'])):
        legacy, state, before = inverse.baseline(manifest, root)
        exact_bytes = state.tobytes()
        pooled = select_ae4_mechanism(legacy, MODEL_ID)
        after = pooled.evaluate(0., state)
        e = environment_from_evaluation(legacy, before)
        thermo = independent_thermodynamics(e)
        f = evaluate_pooled_ae4(e, pooled.ae4_parameters)
        if abs(f.affinity - thermo['independent_affinity']) > 1e-13:
            raise AssertionError('Independent thermodynamic check failed')
        if after.diagnostics.ae4.cl_cell_fmol_s != f.j_ae4_fmol_s:
            raise AssertionError('Production adapter changed pooled cycle flux')
        a4 = after.diagnostics.ae4
        if a4.na_cell_fmol_s * a4.k_cell_fmol_s < 0:
            raise AssertionError('Pooled model generated opposed cation sources')
        inherited_allocation = inverse.allocation(before, None)
        lplus = inherited_allocation['inherited_positive_loading_pool_fmol_s']
        ref = inverse.reference_capacities(legacy)
        reference_row = {
            'root_id': root, 'routing_family': inverse.routing(manifest, root),
            **asdict(e), **thermo, 'pooled_affinity': f.affinity,
            'legacy_na_source_fmol_s': before.diagnostics.ae4.na_cell_fmol_s,
            'legacy_k_source_fmol_s': before.diagnostics.ae4.k_cell_fmol_s,
            'legacy_cl_source_fmol_s': before.diagnostics.ae4.cl_cell_fmol_s,
            'pooled_na_source_at_reference_capacity_fmol_s': a4.na_cell_fmol_s,
            'pooled_k_source_at_reference_capacity_fmol_s': a4.k_cell_fmol_s,
            'pooled_cl_source_at_reference_capacity_fmol_s': a4.cl_cell_fmol_s,
            'pooled_reference_state_is_accepted_wt': False,
            'inherited_positive_loading_pool_fmol_s': lplus,
            'pooled_no_slip': True,
        }
        reference_rows.append(reference_row)
        for condition_index, condition in enumerate(CONDITIONS):
            target = before.diagnostics.ae4.cl_cell_fmol_s if condition == 'P0' else .10*lplus
            target_share = target/lplus
            rho = (lplus-target)/inherited_allocation['signed_basolateral_cl_fmol_s']['NKCC1']
            formal_capacity = target/(f.j_ae4_fmol_s/ref['AE4']) if f.j_ae4_fmol_s else None
            incompatible = target != 0 and (target*f.affinity < 0 or f.affinity == 0)
            row = {
                'candidate_id': inverse.slug(root)+'__'+condition, 'root_id': root,
                'routing_family': inverse.routing(manifest, root), 'condition': condition,
                'target_ae4_cl_fmol_s': target, 'target_ae4_share': target_share,
                'inherited_positive_loading_pool_fmol_s': lplus,
                'required_other_positive_loader_multiplier': rho,
                'pooled_affinity': f.affinity, 'target_times_affinity': target*f.affinity,
                'formal_signed_carrier_required_fmol': formal_capacity,
                'formal_signed_capacity_is_applied': False,
                'fixed_wt_feasible': False, 'status': INFEASIBLE if incompatible else 'PENDING_INVERSE',
                'exact_inherited_state_bytes_retained': state.tobytes() == exact_bytes,
                'state_coordinates_changed': 0,
                'complete_state_bytes_sha256': inverse.state_bytes_sha256(state),
                'pooled_no_slip': True, 'capacity_optimizer_run': False,
                'canonical_objective_has_feasible_point': False,
                'largest_abs_log_fold': None, 'largest_fold_multiplier': None,
            }
            payload = {
                'row': row, 'contract_sha256': sha256_file(OUT/'contract.json'),
                'core_state': manifest['roots'][root]['core_state'],
                'core_state_sha256': manifest['roots'][root]['core_state_sha256'],
                'complete_state': state.tolist(),
                'complete_state_bytes_sha256': inverse.state_bytes_sha256(state),
                'reference_capacities': ref, 'baseline_whole_cell_parameters': asdict(legacy.parameters),
                'baseline_ae4_parameters': asdict(legacy.ae4_parameters),
                'pooled_reference_parameters': asdict(pooled.ae4_parameters),
                'accepted_capacities': None, 'capacity_multipliers': None,
                'environment': asdict(e), 'independent_thermodynamics': thermo,
                'pooled_reference_flux': asdict(f),
                'baseline_residuals': inverse.residual_summary(before),
                'unrebalanced_reference_rhs_diagnostic': inverse.residual_summary(after),
                'accepted_wt_closure': None, 'phenotype_information_used': False,
                'feasibility_certificate': {
                    'proof': 'kf/kr=exp(A)<1 implies J=E*g*(kf-kr)<=0 for E,g>=0; P0 and P10 demand J>0. Other capacities and voltages do not enter A at a fixed conserved state.',
                    'capacity_bound_independent': True, 'barrier_gauge_independent': True,
                    'electrical_reclosure_cannot_change_affinity': True,
                    'nonlinear_optimizer_failure': False,
                    'canonical_minimax_and_tiebreak': 'undefined over empty feasible set' if incompatible else 'pending',
                },
            }
            if not incompatible:
                # General feasible path reuses the inherited method. P0 gets
                # its explicit share so it cannot enter Task 18's baseline bypass.
                allocation = {**inherited_allocation, 'target_share': target_share,
                    'ae4_capacity_multiplier': target/f.j_ae4_fmol_s,
                    'nkcc1_capacity_multiplier': rho, 'other_positive_loader_common_multiplier': rho}
                solution = inverse.solve_canonical(pooled, after, allocation,
                                                   seed=200000+101*root_index+condition_index)
                accepted, closure = inverse.verify_solution(pooled, state, before, solution, allocation)
                logs = inverse.log_changes(solution['capacities'], ref)
                row.update(fixed_wt_feasible=True, status='FIXED_WT_EXACT_REBALANCE_FEASIBLE',
                    capacity_optimizer_run=True, canonical_objective_has_feasible_point=True,
                    largest_abs_log_fold=float(np.max(np.abs(logs))),
                    largest_fold_multiplier=float(math.exp(np.max(np.abs(logs)))))
                payload.update(accepted_capacities=solution['capacities'], accepted_wt_closure=closure,
                    accepted_whole_cell_parameters=asdict(accepted.parameters),
                    accepted_ae4_parameters=asdict(accepted.ae4_parameters),
                    capacity_multipliers={k:solution['capacities'][k]/ref[k] for k in ref},
                    allocation=allocation, stage1_optimisation=solution['stage1'],
                    stage2_optimisation=solution['stage2'], search_diagnostics=solution['search'])
            if state.tobytes() != exact_bytes:
                raise AssertionError('Inherited WT state was modified')
            write_json(OUT/'wt_parameter_payloads'/f"{row['candidate_id']}.json", payload)
            rows.append(row)
            print(row['candidate_id'], row['status'], f'A={f.affinity:.10g}', flush=True)
    write_rows(OUT/'legacy_vs_pooled_reference.csv', reference_rows)
    write_rows(OUT/'wt_inverse_solutions.csv', rows)


def freeze():
    ensure_unfrozen()
    contract = read(OUT/'contract.json')
    paths = sorted((OUT/'wt_parameter_payloads').glob('*.json'))
    payloads = [read(path) for path in paths]
    expected = {(r,c) for r in contract['root_ids'] for c in CONDITIONS}
    actual = {(p['row']['root_id'],p['row']['condition']) for p in payloads}
    if len(payloads) != 20 or actual != expected:
        raise AssertionError('All 20 feasibility decisions must be present')
    manifest, _ = load_freeze()
    for p in payloads:
        model, state, evaluation = inverse.baseline(manifest, p['row']['root_id'])
        if p['complete_state'] != state.tolist() or p['complete_state_bytes_sha256'] != inverse.state_bytes_sha256(state):
            raise AssertionError('Fixed state mismatch')
        if not p['row']['fixed_wt_feasible']:
            if p['row']['status'] != INFEASIBLE or p['row']['target_times_affinity'] >= 0:
                raise AssertionError('Missing infeasibility certificate')
            if p['accepted_capacities'] is not None:
                raise AssertionError('Infeasible decision cannot have accepted capacities')
        elif p['accepted_wt_closure']['max_abs_scaled_independent_rhs'] > inverse.SPEC.root_scaled_tolerance:
            raise AssertionError('Accepted WT closure failed')
    files = [p for p in OUT.rglob('*') if p.is_file()]
    files += [REPO/'src/modern_full_model/pooled_ae4.py', Path(__file__),
              REPO/'tests/test_task20_pooled_ae4.py', REPO/'tests/test_task20_fixed_wt.py',
              ANALYSIS/'evidence_boundary.md', ANALYSIS/'mechanism.md', ANALYSIS/'wt_rebalancing.md']
    hashes = {relative(p):sha256_file(p) for p in files}
    write_json(OUT/'frozen_wt_manifest.json', {
        'manifest_id': 'TASK20_FIXED_WT_PRE_GENOTYPE_V1', 'created_utc': inverse.now(),
        'branch': BRANCH, 'intake_commit': INTAKE, 'inheritance_commit': INHERITANCE,
        'contract_sha256': sha256_file(OUT/'contract.json'), 'hashes': hashes,
        'all_decisions': [p['row'] for p in payloads],
        'feasible_parameterisations': [
            {'candidate_id':p['row']['candidate_id'], 'payload_path':relative(path),
             'payload_sha256':sha256_file(path)} for path,p in zip(paths,payloads)
            if p['row']['fixed_wt_feasible']],
        'all_ten_roots_retained_in_decisions': True, 'both_legacy_routing_families_retained': True,
        'state_solves_performed': 0, 'genotype_outcomes_used': False, 'genotype_result_files_read': 0,
        'pushed_checkpoint_required_before_genotype': True,
        'feasibility_decisions_complete': True,
        'canonical_selection_complete_for_feasible_candidates': True,
        'no_feasible_candidates': not any(p['row']['fixed_wt_feasible'] for p in payloads),
    })
    print('WT decisions frozen. Push and remotely verify this checkpoint before post-freeze work.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=('assays','prepare','solve','freeze'))
    args = parser.parse_args()
    {'assays':assay_tests, 'prepare':prepare, 'solve':solve, 'freeze':freeze}[args.stage]()


if __name__ == '__main__':
    main()
