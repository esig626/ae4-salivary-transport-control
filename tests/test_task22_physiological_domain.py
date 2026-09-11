"""WT-only Task 22 domain, proof, and checkpoint ordering validation."""
from copy import deepcopy
from dataclasses import asdict
from decimal import Decimal
import json
import subprocess

import numpy as np
import pytest

from modern_full_model import task22_contract as c
from modern_full_model import task22_search as s
from modern_full_model.task22_certificate import certificate_for
from modern_full_model.task21_calibration import WTProblem, with_capacities
from modern_full_model.run_calcium_fast_screen import load_freeze
from modern_full_model.acid_base import _carbon_fractions, total_alkalinity_mM
from modern_full_model.validation import sha256_file


@pytest.mark.parametrize('name',list(c.PHYSIOLOGY))
def test_every_intracellular_and_luminal_boundary_is_hard(name):
    values={k:(b['lower']+b['upper'])/2 for k,b in c.PHYSIOLOGY.items()}
    band=c.PHYSIOLOGY[name]
    for edge in ('lower','upper'):
        values[name]=band[edge]
        assert not s.physiology_failures(values)
    for v in (band['lower']-1e-8,band['upper']+1e-8,float('nan'),float('inf'),-float('inf')):
        values[name]=v
        assert name in s.physiology_failures(values)


@pytest.mark.parametrize('side',['cell','lumen'])
def test_osmotic_ratios_are_relative_to_production_bath(side):
    osm=dict(cell=300.,lumen=300.,bath=300.)
    for v in (240.,360.):
        osm[side]=v
        assert not s.osmotic_failures(osm)
    for v in (239.999,360.001,float('nan'),float('inf')):
        osm[side]=v
        assert side in s.osmotic_failures(osm)
    assert s.osmotic_failures(dict(cell=300.,lumen=300.,bath=0.))==['bath']


@pytest.mark.parametrize('name',c.CAPACITIES)
def test_every_capacity_bound_zero_reference_and_fraction(name):
    reference={k:1. for k in c.CAPACITIES};values=reference.copy()
    for v in (.01,100.):
        values[name]=v
        assert not s.capacity_failures(values,reference,[0.,1.])
    for v in (.009999,100.001,-1.,float('nan'),float('inf')):
        values[name]=v
        assert name in s.capacity_failures(values,reference)
    assert not s.capacity_failures({'absent':0.},{'absent':0.})
    assert 'absent' in s.capacity_failures({'absent':1e-30},{'absent':0.})
    for fraction in (-1e-15,1+1e-15,float('nan'),float('inf')):
        assert 'fractions' in s.capacity_failures(reference,reference,[fraction])


def test_signed_pool_uses_positive_contributors_before_cancellation():
    signed=dict(AE4=2.,NKCC1=3.,AE2=-4.)
    result=s.anti_cancellation(dict(one=500.,two=-500.),signed)
    assert result['positive_loading_pool_fmol_s']==5
    assert result['passes'] and result['maximum_ratio']==100
    assert not s.anti_cancellation(dict(one=500.0001,two=-500.0001),signed)['passes']
    assert s.anti_cancellation(dict(one=700.),dict(AE4=2.,NKCC1=3.,AE2=2.))['passes']
    assert not s.anti_cancellation(dict(one=0.),dict(AE4=0.,NKCC1=-1.,AE2=0.))['passes']
    assert not s.anti_cancellation(dict(one=float('nan')),signed)['passes']


@pytest.fixture(scope='module')
def problems():
    manifest,_=load_freeze()
    return [WTProblem(manifest,r) for r in sorted(manifest['roots'])]


def test_production_ledger_preserves_individual_pump_and_exchanger_sources(problems):
    p=problems[0];ev=p.model.evaluate(0.,p.seed);ledger=s.flux_ledger(ev,p.model)
    flux=ledger['solute_fluxes_fmol_s']
    assert flux['pump_apical_Na']==-3*flux['pump_apical_cycles']
    assert flux['pump_basolateral_K']==2*flux['pump_basolateral_cycles']
    assert flux['AE4_tic']==-2*flux['AE4_cycles']
    assert flux['NKCC1_Cl']==2*flux['NKCC1_cycles']
    assert ledger['water']==asdict(ev.diagnostics.water)
    assert set(ledger['ratios'])==set(flux)


def test_search_transcription_and_osmolarity_match_production(problems):
    for p in problems:
        z=p.initial.copy();ev=p.model.evaluate(0.,p.fast(z)['state'])
        z[25:]=[100*ev.diagnostics.membranes.v_apical_V,100*ev.diagnostics.membranes.v_basolateral_V]
        fast=p.fast(z)
        assert np.max(abs(fast['rhs']-ev.rhs[:12]))<1e-11
        assert s.production_osmolarities(p,fast)==pytest.approx(ev.diagnostics.observables.osmolarities_mOsm,abs=1e-11)
        got=s.production_screen(p,z)
        assert not got['resting_gate_pass']
        assert 'thermodynamics:nonpositive_affinity' in got['resting_gate_failures']


def test_production_carbon_minus_alkalinity_identities(problems):
    # These identities are used in the proof, and hold away from equilibrium.
    for p in problems:
        for kind in ('inherited_root','WT_only_high_carbon'):
            z=s.starting_coordinates(p,dict(kind=kind));lo,hi=s.bounds(p);z=np.clip(z,lo,hi)
            fast=p.fast(z);model=with_capacities(p.model,fast['capacities']);ev=model.evaluate(0.,fast['state'])
            d=ev.diagnostics;r=ev.rhs;co2=d.co2_fluxes_fmol_s
            assert r[3]-r[4]==pytest.approx(co2['bath_to_cell']+co2['lumen_to_cell']-d.homeostasis.nhe1_inward_fmol_s,abs=1e-12)
            li=d.observables.lumen_concentrations_mM
            assert co2['lumen_to_cell']==pytest.approx(d.water.lumen_outflow_pL_s*(li['alkalinity']-li['tic'])+r[10]-r[9],abs=1e-12)
            assert d.membranes.cell_sources_fmol_s['tic']==d.membranes.cell_sources_fmol_s['alkalinity']==0.
            assert d.membranes.lumen_sources_fmol_s['tic']==d.membranes.lumen_sources_fmol_s['alkalinity']


def test_analytic_certificate_covers_all_backgrounds_and_tolerances(problems):
    for p in problems:
        cert=certificate_for(p)
        assert cert['proved']
        assert cert['classification']=='PROVEN_STRUCTURAL_CONTRADICTION'
        assert Decimal(cert['four_balance_residual_budget_fmol_s'])==Decimal('6e-9')
        assert Decimal(cert['bath_to_cell_CO2_lower_fmol_s'])>Decimal('.00005118')
        assert Decimal(cert['pooled_affinity_upper'])<Decimal('-.0764')
        assert Decimal(cert['lumen_TA_minus_TIC_upper_mM'])<Decimal('-.00645')
        assert not cert['optimization_results_used']


def test_certificate_fails_closed_when_its_ph_domain_assumption_changes(problems,monkeypatch):
    bounds=deepcopy(c.PHYSIOLOGY);bounds['ph_l']['upper']=9.
    import modern_full_model.task22_certificate as cert
    monkeypatch.setattr(cert,'PHYSIOLOGY',bounds)
    assert not cert.certificate_for(problems[0])['proved']


def test_carbon_fraction_monotonic_bound_against_production_speciation(problems):
    p=problems[0];a=p.p.acid_base
    fmax=_carbon_fractions(8.,a)
    for ph in np.linspace(6.8,8.,61):
        f=_carbon_fractions(float(ph),a)
        assert f[2]-f[0]<=fmax[2]-fmax[0]+1e-15
        for tic in (1.,10.,80.):
            ta=total_alkalinity_mM(ph=float(ph),total_carbon_mM=tic,
                buffer_total_mM=p.p.geometry.lumen_buffer_total_fmol/.02,
                buffer_pka=a.lumen_buffer_pka,parameters=a)
            assert ta-tic<-.00645


def test_contract_remote_bytes_and_all_frozen_sources():
    path=c.OUT/'precalibration_contract.json';receipt=c.read(c.OUT/'stage0_checkpoint.json')
    c.check_receipt(receipt,path)
    contract=c.read(path)
    assert len(contract['seed_root_ids'])==10 and len(contract['task21_witness_files'])==7
    for name,digest in contract['source_hashes'].items():
        assert sha256_file(c.REPO/name)==digest
    assert c.git('rev-parse','HEAD:archive')==contract['archive_tree_sha']
    assert not c.git('diff','--name-only','--','archive')


@pytest.mark.parametrize('mutation',[{'remote_ref_verified':False},{'remote_blob_sha':'wrong'},
    {'remote_branch':'main'},{'checkpoint_sha':c.BASE},{'remote_verified_utc':'2000-01-01T00:00:00Z'}])
def test_checkpoint_cannot_authorize_wrong_ref_bytes_commit_or_order(mutation):
    receipt={**c.read(c.OUT/'stage0_checkpoint.json'),**mutation}
    with pytest.raises((AssertionError,subprocess.CalledProcessError)):
        c.check_receipt(receipt,c.OUT/'precalibration_contract.json')


def test_search_rejects_missing_stage0(tmp_path,monkeypatch):
    monkeypatch.setattr(c,'OUT',tmp_path)
    with pytest.raises(FileNotFoundError):c.require_stage0()


def test_search_rejects_already_frozen_ensemble(tmp_path,monkeypatch):
    import shutil
    for name in ('precalibration_contract.json','stage0_checkpoint.json'):
        shutil.copy(c.OUT/name,tmp_path/name)
    (tmp_path/'wt_manifest.json').write_text('{}')
    monkeypatch.setattr(c,'OUT',tmp_path)
    monkeypatch.setattr(c,'check_receipt',lambda *a,**kw:None)
    with pytest.raises(AssertionError,match='frozen'):c.require_stage0()


def test_genotype_authorization_requires_complete_nonempty_validated_ensemble(tmp_path,monkeypatch):
    monkeypatch.setattr(c,'OUT',tmp_path)
    monkeypatch.setattr(c,'require_stage0',lambda **kw:None)
    monkeypatch.setattr(c,'check_receipt',lambda *a,**kw:None)
    with pytest.raises(FileNotFoundError):c.authorize_genotypes()
    (tmp_path/'wt_checkpoint_receipt.json').write_text('{}')
    def manifest(complete,ids):
        (tmp_path/'wt_manifest.json').write_text(json.dumps(dict(search_complete=complete,
            file_sha256={},retained_solution_ids=ids)))
    manifest(False,[])
    with pytest.raises(AssertionError,match='incomplete'):c.authorize_genotypes()
    manifest(True,[])
    with pytest.raises(RuntimeError,match='STOP_NO_ADMISSIBLE_WT'):c.authorize_genotypes()
    manifest(True,['test'])
    (tmp_path/'wt_solutions').mkdir()
    (tmp_path/'wt_solutions/test.json').write_text('{"resting_gate_pass":true}')
    (tmp_path/'wt_dynamic_results.json').write_text('[]')
    with pytest.raises(AssertionError,match='Incomplete WT dynamic'):c.authorize_genotypes()
    rows=[dict(solution_id='test',calcium_uM=x,dynamic_gate_pass=x!=.50) for x in c.CALCIUM]
    (tmp_path/'wt_dynamic_results.json').write_text(json.dumps(rows))
    with pytest.raises(AssertionError,match='WT dynamics failed'):c.authorize_genotypes()


def test_all_seed_attempts_preserved_with_hard_capacity_bounds_and_zero_genotypes():
    plan=c.read(c.OUT/'search_plan.json');jobs=plan['jobs']
    assert len(jobs)==27
    contract=c.read(c.OUT/'precalibration_contract.json')
    assert {j['root_id'] for j in jobs if j['kind']=='inherited_root'}==set(contract['seed_root_ids'])
    assert {j['source'] for j in jobs if j['kind']=='task21_resting_witness'}==set(contract['task21_witness_files'])
    for job in jobs:
        result=c.read(c.OUT/'attempts'/(job['attempt_id']+'.json'))
        assert result['genotype_evaluations']==0
        assert result['stage0_checkpoint_sha']==c.read(c.OUT/'stage0_checkpoint.json')['checkpoint_sha']
        assert result['started_utc']>=c.read(c.OUT/'stage0_checkpoint.json')['remote_verified_utc']
        candidate=result['candidate'];assert candidate is not None
        assert not s.capacity_failures(candidate['calibrated_capacities'],candidate['reference_capacities'])
        assert not candidate['resting_gate_pass']


def test_all_saved_wt_endpoints_replay_through_production(problems):
    by_root={p.root:p for p in problems}
    for path in sorted((c.OUT/'attempts').glob('*.json')):
        result=c.read(path);candidate=result['candidate'];problem=by_root[result['root_id']]
        model=with_capacities(problem.model,candidate['calibrated_capacities'])
        ev=model.evaluate(0.,candidate['complete_state'])
        np.testing.assert_array_equal(ev.rhs,candidate['complete_rhs'])
        assert s.inverse.residual_summary(ev)==candidate['closure']
        assert asdict(model.parameters)==candidate['whole_cell_parameters']
        assert s.flux_ledger(ev,model)==candidate['flux_ledger']
