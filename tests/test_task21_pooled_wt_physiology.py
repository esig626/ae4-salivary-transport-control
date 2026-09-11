"""Task 21 regression gates use production equations and saved WT-only evidence."""
from copy import deepcopy
from dataclasses import asdict
import hashlib
import json
import numpy as np
import pytest

from modern_full_model import task21_calibration as task
from modern_full_model.task21_wt_checkpoint import candidate_model, equivalent
from modern_full_model import task21_review as review
from modern_full_model.task21_flux_search import coordinates, linear_system
from modern_full_model.validation import CONSERVATION_RESIDUAL_TOLERANCES, PRODUCTION_RADAU


def candidates():
    return [task.read(p) for p in sorted((task.OUT/'wt_candidates').glob('*.json'))]


@pytest.mark.parametrize('field',list(task.PHYSIOLOGY))
def test_physiology_contract_enforces_each_boundary(field):
    values={k:(b['lower']+b['upper'])/2 for k,b in task.PHYSIOLOGY.items()}
    band=task.PHYSIOLOGY[field]
    for edge in ('lower','upper'):
        values[field]=band[edge]
        assert not task.physiology_failures(values)
    for invalid in (band['lower']-1e-6,band['upper']+1e-6,float('nan'),float('inf')):
        values[field]=invalid
        assert field in task.physiology_failures(values)


def test_all_rest_valid_witnesses_have_positive_no_slip_production_closure():
    saved=candidates();assert saved
    for c in saved:
        model=candidate_model(c,.10);ev=model.evaluate(0.,c['complete_state'])
        residual=task.inverse.residual_summary(ev);ae4=ev.diagnostics.ae4
        assert ae4.diagnostics['affinity']>0 and ae4.cl_cell_fmol_s>0
        assert ae4.na_cell_fmol_s<=0 and ae4.k_cell_fmol_s<=0
        assert ae4.na_cell_fmol_s+ae4.k_cell_fmol_s == pytest.approx(-ae4.cl_cell_fmol_s,rel=1e-13)
        assert ae4.tic_cell_fmol_s==pytest.approx(-2*ae4.cl_cell_fmol_s,rel=1e-13)
        passes=all(residual[k]<=t for k,t in (
            ('max_abs_scaled_independent_rhs',task.SPEC.root_scaled_tolerance),
            ('max_abs_omitted_rhs',task.SPEC.omitted_row_raw_tolerance),
            ('max_abs_current_A',task.SPEC.current_tolerance_A),
            ('max_abs_state_charge_fmol',task.SPEC.charge_tolerance_fmol),
            ('max_conservation_tolerance_ratio',1.)))
        assert passes==c['resting_gate_pass']
        assert residual==c['closure']
        assert not task.physiology_failures(c['observables'])
        assert asdict(model.parameters)==c['whole_cell_parameters']
        assert asdict(model.ae4_parameters)==c['ae4_parameters']


def test_deduplication_requires_state_and_parameter_agreement():
    a=candidates()[0];b=deepcopy(a);assert equivalent(a,b)
    b['complete_state'][0]*=1+5e-7;assert equivalent(a,b)
    b['complete_state'][0]*=1+1e-4;assert not equivalent(a,b)
    b=deepcopy(a);b['whole_cell_parameters']['homeostasis']['ae2_capacity_fmol_s']*=1.001
    assert not equivalent(a,b)


def test_all_ten_seeds_and_complete_dynamic_panels_are_preserved():
    manifest,_=task.load_freeze();inventory=task.read(task.OUT/'wt_candidate_inventory.json')
    for folder in ('seed_attempts','refinements','mean_flux_search','profile_calibration','joint_tie_break'):
        outputs=[task.read(p) for p in (task.OUT/folder).glob('*.json') if '.progress.' not in p.name]
        assert {r['root_id'] for r in outputs}==set(manifest['roots'])
        assert all(r['genotype_results_used'] is False for r in outputs)
    rows=task.read(task.OUT/'wt_dynamic_results.json')
    assert len(rows)==3*len(inventory['candidate_ids'])
    for cid in inventory['candidate_ids']:
        assert {r['calcium_uM'] for r in rows if r['solution_id']==cid}=={.10,.25,.50}
    for path in (task.OUT/'wt_trajectories').glob('*.json'):
        payload=task.read(path);assert payload['solver']==asdict(PRODUCTION_RADAU)
        if payload['trajectory']:
            worst=max(abs(v)/CONSERVATION_RESIDUAL_TOLERANCES[k]
                      for k,v in payload['trajectory']['max_abs_conservation_residuals'].items())
            assert worst==pytest.approx(payload['summary']['max_conservation_tolerance_ratio'])
            if worst>1:assert not payload['summary']['dynamic_gate_pass']


def test_mechanism_sources_and_archive_are_unchanged():
    contract=task.read(task.OUT/'precalibration_contract.json')
    for name,digest in contract['source_hashes'].items():
        assert task.sha256_file(task.REPO/name)==digest
    assert task.git('rev-parse','HEAD:archive')==contract['archive_tree_sha']
    assert not task.git('diff','--name-only','--','archive')
    raw=(task.OUT/'precalibration_contract.json').read_bytes()
    assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()==task.STAGE0_BLOB


def test_checkpoint_blocks_genotypes_without_receipt_or_admissible_wt(monkeypatch,tmp_path):
    monkeypatch.setattr(review,'OUT',tmp_path)
    with pytest.raises(FileNotFoundError):review.authorize_genotypes()
    monkeypatch.setattr(review,'checked_checkpoint',lambda:({'retained_solution_ids':[]},{}))
    with pytest.raises(RuntimeError,match='STOP_NO_ADMISSIBLE_WT'):review.authorize_genotypes()


def test_matched_ratio_arithmetic_rejects_invalid_pairs():
    wt=dict(solution_id='synthetic',calcium_uM=.1,dynamic_gate_pass=True,Q_0_600_pL=12.)
    loss={**wt,'Q_0_600_pL':9.}
    assert review.matched_ratio(wt,loss)==.75
    assert review.matched_ratio(wt,{**loss,'dynamic_gate_pass':False}) is None
    with pytest.raises(ValueError):review.matched_ratio(wt,{**loss,'solution_id':'different'})
    with pytest.raises(ValueError):review.matched_ratio(wt,{**loss,'calcium_uM':.5})
    with pytest.raises(ValueError):review.matched_ratio({**wt,'Q_0_600_pL':0.},loss)


def test_actual_frozen_checkpoint_and_no_genotype_refit():
    if not (task.OUT/'wt_checkpoint_receipt.json').exists():pytest.skip('WT checkpoint is not yet pushed')
    manifest,receipt=review.checked_checkpoint()
    assert manifest['genotype_evaluations']==0
    if not manifest['retained_solution_ids']:
        with pytest.raises(RuntimeError,match='STOP_NO_ADMISSIBLE_WT'):review.authorize_genotypes()
    if (task.OUT/'final_summary.json').exists():
        final=task.read(task.OUT/'final_summary.json')
        assert final['genotype_evaluations']==final['genotype_specific_refits']==0
        assert not final['experimental_loss_magnitude_compared']
