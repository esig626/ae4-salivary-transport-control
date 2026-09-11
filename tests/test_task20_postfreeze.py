"""No imputed pooled outcomes, complete accounting, verified legacy reuse."""
from unittest.mock import patch

from src.modern_full_model import task20_fixed_wt as task
from src.modern_full_model import task20_postfreeze as post


def test_checkpoint_gate_verifies_remote_receipt_and_all_frozen_files():
    frozen = post.checked_checkpoint()
    assert len(frozen['all_decisions']) == 20
    assert frozen['feasible_parameterisations'] == []
    original_read = task.read

    def corrupt_receipt(path):
        data = original_read(path)
        if path.name == 'pushed_wt_checkpoint.json':
            data = {**data, 'remote_ref_verified':False}
        return data

    with patch.object(task,'read',side_effect=corrupt_receipt):
        try:
            post.checked_checkpoint()
        except AssertionError:
            pass
        else:
            raise AssertionError('Unverified remote checkpoint was accepted')


def test_each_required_pooled_dynamic_table_retains_all_sixty_unavailable_rows():
    expected = {(r,c,ca) for r in task.read(task.OUT/'contract.json')['root_ids']
                for c in task.CONDITIONS for ca in task.CALCIUM}
    for name in ('wt_dynamic_results.csv','ae4_5pct_results.csv','ae2_results.csv'):
        rows = post.csv_rows(task.OUT/name)
        assert len(rows) == 60
        assert {(r['root_id'],r['condition'],float(r['calcium_uM'])) for r in rows} == expected
        for row in rows:
            assert row['status'] == post.UNAVAILABLE
            assert row['reason'] == task.INFEASIBLE
            assert row['simulation_attempted'] == 'False'
            assert row['fixed_wt_feasible'] == 'False'
            for field in ('total_0_600_pL','secretion_ratio','solver_success','numerical_gate_pass'):
                assert row[field] == ''


def test_four_conceptual_conditions_are_complete_with_hash_checked_legacy_values():
    rows = post.csv_rows(task.OUT/'final_comparison.csv')
    assert len(rows) == 120
    assert {r['condition'] for r in rows} == {'legacy_inherited','legacy_P10','P0','P10'}
    tables, verified = post.verify_legacy()
    expected = post.legacy_comparison(tables)
    legacy = [r for r in rows if r['architecture']=='legacy_SR2']
    assert len(legacy) == len(expected) == 60
    for actual, source in zip(legacy,expected):
        for name in ('R_AE4','R_AE2','wt_total_0_600_pL','ae4_5pct_total_0_600_pL','ae2_total_0_600_pL'):
            assert float(actual[name]) == source[name]
        assert actual['source_sha256'] == source['source_sha256']
    pooled = [r for r in rows if r['architecture']==post.MODEL_ID]
    assert len(pooled) == 60
    assert all(r['R_AE4']==r['R_AE2']=='' for r in pooled)
    assert all(r['comparison_to_pooled_available']=='False' for r in rows)
    assert task.read(task.OUT/'legacy_reuse_manifest.json')['frozen_inherited_inputs_unchanged'] == 187


def test_native_state_diagnostics_distinguish_targets_from_accepted_genotype_states():
    rows = post.csv_rows(task.OUT/'state_flux_summary.csv')
    assert len(rows) == 60
    assert {r['genotype'] for r in rows} == {'WT','AE4_5pct','AE2_loss'}
    for row in rows:
        assert row['status'] == post.UNAVAILABLE
        assert row['reference_pooled_is_accepted_wt'] == 'False'
        assert float(row['reference_pooled_na_source_fmol_s']) > 0
        assert float(row['reference_pooled_k_source_fmol_s']) > 0
        assert float(row['reference_pooled_cl_source_fmol_s']) < 0
        for field in ('accepted_resting_cl_i_mM','accepted_resting_ph_i','accepted_resting_na_i_mM',
                      'accepted_resting_k_i_mM','delta_cl_i_mM','delta_ph_i','delta_na_i_mM','delta_k_i_mM',
                      'NKCC1_cycles_fmol_s','NHE1_fmol_s','AE2_cl_fmol_s','pump_apical_cycles_fmol_s',
                      'pump_basolateral_cycles_fmol_s','lower_cl_after_ae4_loss'):
            assert row[field] == ''


def test_report_does_not_conflate_infeasibility_with_a_genotype_result_or_optimizer_failure():
    summary = task.read(task.OUT/'postfreeze_summary.json')
    assert summary['primary_classification'] == post.CLASSIFICATION
    assert summary['scientific_status'] == task.INFEASIBLE
    assert 'not numerical nonconvergence' in summary['classification_limitation']
    for name in ('new_wt_trajectories','new_ae4_trajectories','new_ae2_trajectories',
                 'genotype_continuations','exact_ae4_zero_attempts','capacity_refits_after_freeze'):
        assert summary[name] == 0
    assert summary['ae4_direction_P0'] is None
    assert summary['ae4_direction_P10'] is None
    assert summary['pooled_ae2_neutrality'] is None
    assert summary['ae4_loss_resting_cl_direction'] is None
    assert summary['experimental_deficit_used_as_target'] is False
