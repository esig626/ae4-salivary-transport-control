"""Independent feasibility certificates and immutable WT decision checks."""
from dataclasses import replace
from decimal import Decimal
import math

import numpy as np

from src.modern_full_model import task20_fixed_wt as task
from src.modern_full_model import task18_fixed_wt as inverse
from src.modern_full_model.pooled_ae4 import MODEL_ID, select_ae4_mechanism
from src.modern_full_model.run_calcium_fast_screen import load_freeze
from src.modern_full_model.validation import sha256_file


def payloads():
    return [task.read(p) for p in sorted((task.OUT/'wt_parameter_payloads').glob('*.json'))]


def test_every_root_and_condition_has_an_exact_state_and_a_sign_certificate():
    manifest, _ = load_freeze()
    items = payloads()
    assert len(items) == 20
    assert {(p['row']['root_id'],p['row']['condition']) for p in items} == {
        (r,c) for r in manifest['roots'] for c in task.CONDITIONS}
    for p in items:
        row = p['row']
        model, state, evaluation = inverse.baseline(manifest, row['root_id'])
        assert p['core_state'] == manifest['roots'][row['root_id']]['core_state']
        assert p['complete_state'] == state.tolist()
        assert p['complete_state_bytes_sha256'] == inverse.state_bytes_sha256(state)
        # Independent raw-concentration mass-action inequality, without using
        # either production pooled affinity or the certificate helper.
        e = p['environment']
        numerator = e['cl_o_mM']*(e['na_i_mM']+e['k_i_mM'])*e['hco3_i_mM']**2
        denominator = e['cl_i_mM']*(e['na_o_mM']+e['k_o_mM'])*e['hco3_o_mM']**2
        assert numerator < denominator
        assert Decimal(p['independent_thermodynamics']['activity_ratio_decimal50']) < 1
        assert row['target_ae4_cl_fmol_s'] > 0
        assert row['target_times_affinity'] < 0
        assert row['formal_signed_carrier_required_fmol'] < 0
        assert row['formal_signed_capacity_is_applied'] is False
        assert row['fixed_wt_feasible'] is False
        assert row['status'] == task.INFEASIBLE
        assert p['accepted_capacities'] is None
        assert p['accepted_wt_closure'] is None
        assert row['largest_abs_log_fold'] is None
        assert p['baseline_residuals']['max_abs_scaled_independent_rhs'] <= inverse.SPEC.root_scaled_tolerance
        assert p['baseline_residuals']['max_conservation_tolerance_ratio'] <= 1


def test_only_p0_p10_with_correct_loading_and_no_thirty_percent_condition():
    manifest, _ = load_freeze()
    for p in payloads():
        row = p['row']
        model, state, evaluation = inverse.baseline(manifest, row['root_id'])
        loading = inverse.allocation(evaluation, None)
        pool = loading['inherited_positive_loading_pool_fmol_s']
        target = evaluation.diagnostics.ae4.cl_cell_fmol_s if row['condition'] == 'P0' else .1*pool
        assert row['target_ae4_cl_fmol_s'] == target
        assert row['target_ae4_share'] == target/pool
        rho = row['required_other_positive_loader_multiplier']
        assert math.isclose(target + rho*loading['signed_basolateral_cl_fmol_s']['NKCC1'], pool, rel_tol=2e-15)
    assert task.read(task.OUT/'contract.json')['conditions'] == ['P0','P10']


def test_capacity_changes_and_voltage_reclosure_cannot_reverse_fixed_state_affinity():
    manifest, _ = load_freeze()
    for root in manifest['roots']:
        legacy, state, baseline = inverse.baseline(manifest, root)
        model = select_ae4_mechanism(legacy, MODEL_ID)
        reference = inverse.reference_capacities(model)
        affinity = model.evaluate(0., state).diagnostics.ae4.diagnostics['affinity']
        assert affinity < 0
        for factor in (1e-3, 1e3):
            capacities = {name:value*factor for name,value in reference.items()}
            altered = inverse.model_with_capacities(model, capacities)
            result = altered.evaluate(0., state)
            ae4 = result.diagnostics.ae4
            assert ae4.diagnostics['affinity'] == affinity
            assert ae4.cl_cell_fmol_s < 0
            assert ae4.na_cell_fmol_s >= 0 and ae4.k_cell_fmol_s >= 0
            assert inverse.residual_summary(result)['max_conservation_tolerance_ratio'] <= 1
        zero = select_ae4_mechanism(model, MODEL_ID,
                                  parameters=replace(model.ae4_parameters, carrier_amount_fmol=0.))
        assert zero.evaluate(0., state).diagnostics.ae4.cl_cell_fmol_s == 0
        # Zero capacity is an assay feasibility boundary, not an AE4 genotype run.


def test_wt_preparation_has_no_genotype_execution_or_outcome_inputs():
    import ast
    tree = ast.parse((task.REPO/'src/modern_full_model/task20_fixed_wt.py').read_text())
    imports = [node.module or '' for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)]
    assert not any('genotype' in name or 'postfreeze' in name for name in imports)
    contract = task.read(task.OUT/'contract.json')
    assert not contract['genotype_information_in_objective']
    assert not contract['genotype_information_read_before_freeze']
    assert not any('ae4_5pct' in name or 'ae2_results' in name for name in contract['prefreeze_input_hashes'])
    for p in payloads():
        assert not p['phenotype_information_used']
        assert not p['row']['capacity_optimizer_run']
        assert p['feasibility_certificate']['capacity_bound_independent']
        assert not p['feasibility_certificate']['nonlinear_optimizer_failure']


def test_freeze_preserves_all_decisions_even_when_no_candidate_is_feasible():
    path = task.OUT/'frozen_wt_manifest.json'
    if not path.exists():
        # Prefreeze invocation; the same test checks hashes after checkpoint.
        return
    frozen = task.read(path)
    assert len(frozen['all_decisions']) == 20
    assert frozen['feasible_parameterisations'] == []
    assert frozen['no_feasible_candidates'] is True
    assert frozen['genotype_outcomes_used'] is False
    assert frozen['genotype_result_files_read'] == 0
    for name,digest in frozen['hashes'].items():
        assert sha256_file(task.REPO/name) == digest, name
