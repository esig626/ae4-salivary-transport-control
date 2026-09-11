"""Task 17: fixed-state source algebra; no state optimizer or ODE calls.

The electrical source matrix is evaluated at the inherited concentrations,
speciation and voltages. Current equations are retained as explicit rows.
The native 14 nS total K source constraint is a separate equality row.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
from datetime import datetime, timezone
from itertools import combinations
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.optimize import linprog

from .calibration import CORE_INDEPENDENT_ROWS, WTCalibrationSpec
from .model import ModernFullModel
from .nkcc_stimulation import StimulatedNkcc1Model
from .run_calcium_fast_screen import load_freeze, write_json, write_rows
from .states import CORE_STATE_NAMES
from .task14_blind import build_model, routing, slug
from .task16_wt_allocation import check_inputs
from .validation import CONSERVATION_RESIDUAL_TOLERANCES, attach_basal_regulation, sha256_file, sha256_object

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / 'results/17_fixed_wt_flux_repartition'
ANALYSIS = REPO / 'analysis/17_fixed_wt_flux_repartition'
BRANCH = 'codex/task-17-fixed-wt-flux-repartition'
BASE = '4d4403f279fc36aec940c9e4759486d02f07e3f6'
INTAKE = '70a342487100edf8a8f407d02944ba4ec652d786'
CONDITIONS = ('inherited_baseline', 'share_0.10', 'share_0.30')
TARGETS = {'share_0.10': .10, 'share_0.30': .30}
AXES = ('NHE1', 'pump_apical', 'pump_basolateral', 'K_apical',
        'K_basolateral', 'CO2_basolateral', 'CO2_apical',
        'paracellular_Na', 'paracellular_K', 'paracellular_HCO3')
ROWS = (*CORE_STATE_NAMES, 'apical_current_fmol_charge_s',
        'basolateral_current_fmol_charge_s', 'total_K_capacity_fraction')
SPEC = WTCalibrationSpec()


def read(path):
    return json.loads(Path(path).read_text())


def git(*args):
    return subprocess.check_output(['git', *args], cwd=REPO, text=True).strip()


def now():
    return datetime.now(timezone.utc).isoformat()


def relative(path):
    return str(Path(path).relative_to(REPO))


def baseline(manifest, root_id):
    model = build_model(manifest, root_id, .10)
    # Copy the frozen amounts verbatim. No encode/decode round trip is used.
    state = attach_basal_regulation(model, manifest['roots'][root_id]['core_state'])
    state.flags.writeable = False
    return model, state, model.evaluate(0., state)


def source_ledger(model, evaluation):
    """Independently assemble named pathway vectors in the production layout."""
    d = evaluation.diagnostics
    h, a, e = d.homeostasis, d.ae4, d.membranes
    convert = 1e15 / model.parameters.constants.faraday_C_mol
    rows = {}

    def add(name, entries, ia=0., ib=0., include=True):
        v = np.zeros(14)
        for index, value in entries.items():
            v[index] = value
        v[12:14] = ia, ib
        rows[name] = (v, include)

    jna, jk = a.diagnostics['j_na_fmol_s'], a.diagnostics['j_k_fmol_s']
    add('AE4_Na_branch', {0: -jna, 2: jna, 3: -2*jna, 4: -2*jna})
    add('AE4_K_branch', {1: -jk, 2: jk, 3: -2*jk, 4: -2*jk})
    add('AE4_total', {0: -jna, 1: -jk, 2: jna+jk, 3: -2*(jna+jk), 4: -2*(jna+jk)}, include=False)
    n, q, z = h.nkcc1_inward_fmol_s, h.nhe1_inward_fmol_s, h.ae2_inward_fmol_s
    add('NKCC1', {0:n, 1:n, 2:2*n})
    add('AE2', {2:z, 3:-z, 4:-z})
    add('NHE1', {0:q, 4:q})
    pa, pb = e.pump_apical_fmol_s, e.pump_basolateral_fmol_s
    add('pump_apical', {0:-3*pa, 1:2*pa, 6:3*pa, 7:-2*pa}, ia=pa)
    add('pump_basolateral', {0:-3*pb, 1:2*pb}, ib=pb)
    ka, kb = [e.currents_A[k]*convert for k in ('k_apical','k_basolateral')]
    add('K_apical', {1:-ka, 7:ka}, ia=ka)
    add('K_basolateral', {1:-kb}, ib=kb)
    cl = -e.currents_A['cl_apical']*convert
    add('CaCC', {2:-cl, 8:cl}, ia=-cl)
    for name, key, idx, valence in (
        ('paracellular_Na','para_na',6,1), ('paracellular_K','para_k',7,1),
        ('paracellular_Cl','para_cl',8,-1), ('paracellular_HCO3','para_hco3',9,-1)):
        current = e.currents_A[key]*convert
        entries = {idx:-current/valence}
        if idx == 9:
            entries[10] = -current/valence
        add(name, entries, ia=-current, ib=current)
    cb, ca = d.co2_fluxes_fmol_s['bath_to_cell'], d.co2_fluxes_fmol_s['lumen_to_cell']
    add('CO2_basolateral', {3:cb})
    add('CO2_apical', {3:ca, 9:-ca})
    w = d.water
    add('water', {5:w.cell_volume_source_pL_s, 11:w.lumen_volume_source_pL_s})
    add('lumen_outflow', {6+i:-d.outflow_sources_fmol_s[k]
        for i,k in enumerate(('na','k','cl','tic','alkalinity'))})
    total = sum((v for v, include in rows.values() if include), np.zeros(14))
    np.testing.assert_allclose(total[:12], evaluation.rhs[:12], atol=4e-15, rtol=1e-10)
    np.testing.assert_allclose(total[12:], np.array([
        e.current_residuals_A['apical'],e.current_residuals_A['basolateral']])*convert,
        atol=4e-15, rtol=1e-10)
    return rows


def source_matrix(model, ledger):
    matrix = np.zeros((15,len(AXES)))
    matrix[:14] = np.column_stack([ledger[k][0] for k in AXES])
    fk = model.parameters.membranes.apical_k_fraction
    matrix[14,AXES.index('K_apical')] = fk
    matrix[14,AXES.index('K_basolateral')] = 1-fk
    return matrix


def rank_info(matrix):
    # Normalize columns before rank decisions; rows are rates except for the
    # explicit dimensionless K-source constraint, whose scaling is recorded.
    scale = np.maximum(np.max(np.abs(matrix),axis=1),1e-14)
    normalized = matrix/scale[:,None]
    normalized /= np.maximum(np.linalg.norm(normalized,axis=0),1e-14)
    u, s, vh = np.linalg.svd(normalized, full_matrices=True)
    rank = int(np.sum(s > SPEC.jacobian_rank_relative_tolerance*s[0]))
    # Lexicographic smallest full-column-rank basis of the complete source span.
    basis = None
    for subset in combinations(range(len(AXES)), rank):
        values = np.linalg.svd(normalized[:,subset], compute_uv=False)
        if values[-1] > SPEC.jacobian_rank_relative_tolerance*values[0]:
            basis = [AXES[i] for i in subset]
            break
    return {'rank':rank, 'nullity':len(AXES)-rank, 'row_scale':scale.tolist(),
            'singular_values':s.tolist(), 'smallest_full_rank_basis':basis}


def allocation(ledger, target):
    a = ledger['AE4_total'][0][2]
    positives = {k:ledger[k][0][2] for k in ('NKCC1','AE2') if ledger[k][0][2]>0}
    other = sum(positives.values())
    assert a > 0 and other > 0
    if target is not None and target not in TARGETS.values():
        raise ValueError('Only the two predeclared target shares are allowed')
    lam = 1. if target is None else target*(a+other)/a
    rho = 1. if target is None else (1-target)*(a+other)/other
    delta = (lam-1)*ledger['AE4_total'][0]
    for key in positives:
        delta += (rho-1)*ledger[key][0]
    return {'ae4_capacity_multiplier':lam, 'rho':rho, 'positive_other_loaders':positives,
            'A0_fmol_s':a, 'O0_fmol_s':other, 'Lplus_fmol_s':a+other}, delta


def allocated_only(model, values):
    hp = model.parameters.homeostasis
    changes = {}
    for name in values['positive_other_loaders']:
        key = {'NKCC1':'nkcc1_capacity_fmol_s','AE2':'ae2_capacity_fmol_s'}[name]
        changes[key] = getattr(hp,key)*values['rho']
    parameters = replace(model.parameters, homeostasis=replace(hp,**changes))
    ae4 = replace(model.ae4_parameters,carrier_amount_fmol=(
        model.ae4_parameters.carrier_amount_fmol*values['ae4_capacity_multiplier']))
    return StimulatedNkcc1Model(ModernFullModel(parameters=parameters,
        ae4_parameters=ae4, ae4_evaluator=model.ae4_evaluator,
        regulatory_model=model.regulatory_model, stimulus=model.stimulus))


def with_effective_multipliers(model, multipliers):
    """Map effective capacity coordinates back to existing model parameters.

    This utility does not certify eligibility. In particular callers must
    separately enforce the native total-K equality before accepting a model.
    """
    m=dict(zip(AXES,map(float,multipliers)))
    if not all(np.isfinite(v) and v>=0 for v in m.values()):
        raise ValueError('Effective capacities must be finite and nonnegative')
    hp,mp=model.parameters.homeostasis,model.parameters.membranes
    pa=mp.nak_capacity_fmol_s*mp.apical_pump_fraction*m['pump_apical']
    pb=mp.nak_capacity_fmol_s*(1-mp.apical_pump_fraction)*m['pump_basolateral']
    ka=mp.g_k_total_S*mp.apical_k_fraction*m['K_apical']
    kb=mp.g_k_total_S*(1-mp.apical_k_fraction)*m['K_basolateral']
    if pa+pb<=0 or ka+kb<=0:
        raise ValueError('Nonzero pump/K totals required by this coordinate map')
    hp=replace(hp,nhe1_capacity_fmol_s=hp.nhe1_capacity_fmol_s*m['NHE1'],
        co2_basolateral_permeability_fmol_s_mM=hp.co2_basolateral_permeability_fmol_s_mM*m['CO2_basolateral'],
        co2_apical_permeability_fmol_s_mM=hp.co2_apical_permeability_fmol_s_mM*m['CO2_apical'])
    mp=replace(mp,nak_capacity_fmol_s=pa+pb,apical_pump_fraction=pa/(pa+pb),
        g_k_total_S=ka+kb,apical_k_fraction=ka/(ka+kb),
        g_para_na_S=mp.g_para_na_S*m['paracellular_Na'],
        g_para_k_S=mp.g_para_k_S*m['paracellular_K'],
        g_para_hco3_S=mp.g_para_hco3_S*m['paracellular_HCO3'])
    return StimulatedNkcc1Model(ModernFullModel(parameters=replace(model.parameters,
        homeostasis=hp,membranes=mp),ae4_parameters=model.ae4_parameters,
        ae4_evaluator=model.ae4_evaluator,regulatory_model=model.regulatory_model,
        stimulus=model.stimulus))


def residual_summary(ev):
    return {'max_abs_rhs':float(np.max(np.abs(ev.rhs))),
        'max_abs_scaled_independent_rhs':float(np.max(np.abs(
            ev.rhs[list(CORE_INDEPENDENT_ROWS)]/SPEC.independent_rhs_scales))),
        'max_abs_omitted_rhs':float(np.max(np.abs(ev.rhs[[4,10]]))),
        'conservation_residuals':dict(ev.diagnostics.conservation_residuals),
        'max_conservation_tolerance_ratio':max(abs(v)/CONSERVATION_RESIDUAL_TOLERANCES[k]
            for k,v in ev.diagnostics.conservation_residuals.items()),
        'max_abs_state_charge_fmol':max(abs(v) for v in ev.diagnostics.state_charge_fmol.values()),
        'max_abs_current_A':max(abs(v) for v in ev.diagnostics.membranes.current_residuals_A.values())}


def prepare():
    assert git('branch','--show-current') == BRANCH
    assert not (OUT/'contract.json').exists(), 'Contract is immutable'
    inputs = check_inputs()
    manifest, verification = load_freeze()
    rows, matrices, ranks = [], [], {}
    for root_id in sorted(manifest['roots']):
        model, state, ev = baseline(manifest,root_id)
        ledger = source_ledger(model,ev)
        matrix = source_matrix(model,ledger)
        ranks[root_id] = rank_info(matrix)
        ranks[root_id]['without_native_K_equality'] = rank_info(matrix[:14])
        ranks[root_id]['baseline_residuals'] = residual_summary(ev)
        for name,(vector,include) in ledger.items():
            rows.append({'root_id':root_id,'pathway':name,'include_in_rhs_sum':include,
                **dict(zip(ROWS[:14],vector.tolist()))})
        for i,name in enumerate(ROWS):
            matrices.append({'root_id':root_id,'equation':name,
                **dict(zip(AXES,matrix[i].tolist()))})
    write_rows(OUT/'baseline_flux_ledger.csv',rows)
    write_rows(OUT/'fixed_state_source_matrix.csv',matrices)
    write_json(OUT/'fixed_state_rank.json',ranks)
    write_json(OUT/'baseline_input_hashes.json',inputs)
    contract = {'task':'17_FIXED_WT', 'created_utc':now(), 'branch':BRANCH,
        'scientific_baseline':BASE,'intake_commit':INTAKE,
        'conditions':CONDITIONS,'target_shares':[.10,.30], 'calcium_uM':[.10,.25,.50],
        'root_ids':sorted(manifest['roots']), 'state_coordinates_optimized':0,
        'source_matrix_constructed_before_target_solves':True,
        'source_matrix_sha256':sha256_file(OUT/'fixed_state_source_matrix.csv'),
        'rank_sha256':sha256_file(OUT/'fixed_state_rank.json'),
        'capacity_axes':AXES,
        'axis_coordinates':'Multipliers on existing effective capacities at each membrane; reconstructed total capacities and fractions retain the existing topology.',
        'selection':'Smallest full-rank source basis is recorded. Preserve all null directions in feasibility testing, so basis choice cannot discard feasible solutions.',
        'underdetermined_rule':'If feasible, minimize squared Euclidean distance of effective-capacity multipliers from one subject to all equalities/bounds; no genotype objective.',
        'feasibility_algorithm':'Exact affine source equations, normalized rank, HiGHS primal and dual-simplex feasibility, and independent analytical K-flux certificate.',
        'equality_rows':ROWS,
        'native_K_total_S':14e-9,
        'native_K_source_record':['src/modern_full_model/native_source_panel.py',
            'analysis/13B_modern_full_model/cation_topology.md',
            'analysis/13B_modern_full_model/evidence_ledger.md'],
        'K_rule':'Total K maximum is a literal source constraint, not an unmeasured WT capacity. Apical/basolateral allocation is allowed to vary; no new backgrounds.',
        'allowed_outer_relaxation':'For impossibility certification only, all pump/K fractions may range over [0,1], pump and paracellular capacities have no upper bound, and the NHE1 upper bound is removed. Infeasibility in this superset proves infeasibility within every inherited bracket.',
        'historical_bounds':'Early common-conductance and AE4-carrier fitting boxes belong to superseded calibrations; do not reject by importing them into the native source panel. NHE1 calibration bound [1e-6,0.2] fmol/s is checked if a physical solution exists.',
        'chloride_rules':'AE4 whole carrier scales; all other positive basolateral Cl loaders scale by one rho; signed counterfluxes, CaCC, paracellular Cl stay fixed.',
        'electrical_rule':'Fixed chloride sources with fixed CaCC and paracellular Cl conductances fix Va and Vt, hence Vb. Retain both current-closure equations explicitly.',
        'other_axes':'NHE1, neutral CO2 exchange, existing pump capacities/partition and non-Cl paracellular capacities are the maximal existing unmeasured balancing set; negative AE2 counterflux cannot change without breaking the fixed Cl equation.',
        'fixed_fields':['all WT amounts and volumes','acid-base chemistry','AE4 routing and rates',
            'calcium and regulation','CaCC and paracellular Cl','bath','geometry','water/outflow'],
        'wt_tolerances':asdict(SPEC),'conservation_tolerances':CONSERVATION_RESIDUAL_TOLERANCES,
        'input_count':len(inputs['input_files']),'intake_verification':verification,
        'phenotype_firewall':'No genotype payload is read in preparation/solve/freeze; frozen input byte hashes are verified without outcome access.',
        'postfreeze':'Push all decisions and feasible payloads first. Reuse hash-valid baseline; any feasible modified case gets new connected AE4 0.05 and AE2 0 roots and matched 600s WT/genotype trajectories at all three calcium values.'}
    write_json(OUT/'contract.json',contract)
    print(json.dumps({'source_ranks':{k:v['rank'] for k,v in ranks.items()},'target_solves':0}))


def k_certificate(model, ledger, delta):
    """Necessary K-flux interval, independent of all pump/para partitions."""
    # TA enforces dNHE=-delta_TA. Na then fixes the TOTAL pump change.
    dnhe = -delta[4]
    dpump = (delta[0]+dnhe)/3
    ka, kb = -ledger['K_apical'][0][1], -ledger['K_basolateral'][0][1]
    required_k = ka+kb+delta[1]+2*dpump
    fk = model.parameters.membranes.apical_k_fraction
    all_apical, all_basal = ka/fk, kb/(1-fk)
    low, high = min(all_apical,all_basal), max(all_apical,all_basal)
    gap = low-required_k
    return {'required_NHE1_flux_change':float(dnhe),'required_total_pump_flux_change':float(dpump),
        'inherited_total_K_outward_fmol_s':float(ka+kb),
        'required_total_K_outward_fmol_s':float(required_k),
        'minimum_K_outward_at_14nS_fmol_s':float(low),
        'maximum_K_outward_at_14nS_fmol_s':float(high),
        'K_lower_bound_gap_fmol_s':float(gap),
        'required_K_negative':bool(required_k<0),
        'native_source_bound_incompatible':bool(required_k<low-1e-9 or required_k>high+1e-9),
        'required_apical_K_fraction_signed':float((required_k-all_basal)/(all_apical-all_basal)),
        'certificate_domain':'All K partitions in [0,1]; arbitrary nonnegative existing pump/NHE/CO2/non-Cl paracellular capacities; fixed native total K and chloride sources.',
        'proof':'dNHE=-dTA; dP=(dNa+dNHE)/3; Kreq=K0+dK+2dP; Ktotal=Gtotal*[fK*driveA+(1-fK)*driveB] lies between its two endpoint values.'}


def solve():
    assert not (OUT/'frozen_fixed_wt_manifest.json').exists()
    c=read(OUT/'contract.json')
    assert sha256_file(OUT/'fixed_state_source_matrix.csv')==c['source_matrix_sha256']
    manifest,_=load_freeze()
    decisions=[]
    for root_id in sorted(manifest['roots']):
        model,state,ev=baseline(manifest,root_id)
        ledger=source_ledger(model,ev); matrix=source_matrix(model,ledger)
        norm=np.maximum(np.max(np.abs(matrix),axis=1),1e-12)
        for condition in CONDITIONS:
            target=TARGETS.get(condition)
            values,delta=allocation(ledger,target)
            modified=allocated_only(model,values); allocated=modified.evaluate(0.,state)
            np.testing.assert_allclose(allocated.rhs[:12]-ev.rhs[:12],delta[:12],atol=5e-15,rtol=1e-10)
            certificate=k_certificate(model,ledger,delta)
            b=matrix@np.ones(len(AXES))-np.r_[delta,0.]
            solvers=[]
            for method in ('highs-ds','highs-ipm'):
                result=linprog(np.zeros(len(AXES)),A_eq=matrix/norm[:,None],
                    b_eq=b/norm,bounds=[(0,None)]*len(AXES),method=method,
                    options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
                solvers.append({'method':method,'status':int(result.status),'message':result.message,
                    'success':bool(result.success),'multipliers':result.x.tolist() if result.success else None})
            change=np.linalg.lstsq(matrix/norm[:,None],-np.r_[delta,0.]/norm,rcond=1e-10)[0]
            signed=1+change
            signed_residual=matrix@change+np.r_[delta,0.]
            assert max(abs(signed_residual))<1e-10, 'Structural rank inconsistency requires separate handling'
            if target is None:
                feasible=True; accepted_multipliers=np.ones(len(AXES))
                assert all(s['success'] for s in solvers)
                accepted=residual_summary(allocated)
                assert accepted['max_abs_scaled_independent_rhs']<=SPEC.root_scaled_tolerance
                assert accepted['max_abs_omitted_rhs']<=SPEC.omitted_row_raw_tolerance
                assert accepted['max_conservation_tolerance_ratio']<=1
            else:
                # Do not replace a physical solution with a failure label if
                # the certificate and independent LP methods ever disagree.
                assert certificate['native_source_bound_incompatible']
                assert all(s['status']==2 for s in solvers), solvers
                feasible=False; accepted_multipliers=None; accepted={}
            candidate_id=slug(root_id)+'__'+condition
            a=allocated.diagnostics.ae4
            turnover=abs(a.diagnostics['j_na_fmol_s'])+abs(a.diagnostics['j_k_fmol_s'])
            actual_pool=a.cl_cell_fmol_s+sum(ledger[k][0][2]*values['rho'] for k in values['positive_other_loaders'])
            row={'candidate_id':candidate_id,'root_id':root_id,'routing_family':routing(manifest,root_id),
                'condition':condition,'target_share':target,**{k:v for k,v in values.items() if not isinstance(v,dict)},
                'other_positive_loaders':';'.join(values['positive_other_loaders']),
                'realized_ae4_share':a.cl_cell_fmol_s/actual_pool,
                'positive_loading_pool_change_fmol_s':actual_pool-values['Lplus_fmol_s'],
                'state_coordinates_changed':0,'exact_inherited_state_bytes_retained':True,
                'fixed_wt_feasible':feasible,
                'status':'FIXED_WT_FEASIBLE' if feasible else ('INFEASIBLE_NEGATIVE_REQUIRED_K_FLUX' if certificate['required_K_negative'] else 'INFEASIBLE_NATIVE_K_SOURCE_BOUND'),
                'feasibility_failure_kind':'' if feasible else 'BOUND_INCOMPATIBILITY; NOT_RANK_DEFICIENCY; NOT_NUMERICAL_FAILURE',
                'max_abs_complete_rhs':accepted.get('max_abs_rhs'),
                'max_abs_scaled_complete_rhs':accepted.get('max_abs_scaled_independent_rhs'),
                'max_abs_allocation_only_rhs':float(max(abs(allocated.rhs))),
                'max_abs_signed_source_equation_residual':float(max(abs(signed_residual))),
                'signed_multiplier_role':'unconstrained algebraic diagnostic; not a physical parameterization',
                'ae4_productive_fraction':abs(a.cl_cell_fmol_s)/turnover,
                'ae4_cancellation_fraction':1-abs(a.cl_cell_fmol_s)/turnover,
                **{k:v for k,v in certificate.items() if not isinstance(v,str)},
                **{'balancing_'+k+'_multiplier':float(accepted_multipliers[i]) if feasible else None for i,k in enumerate(AXES)},
                **{'signed_diagnostic_'+k+'_multiplier':float(signed[i]) for i,k in enumerate(AXES)}}
            payload={'row':row,'contract_sha256':sha256_file(OUT/'contract.json'),
                'core_state':manifest['roots'][root_id]['core_state'],
                'core_state_sha256':manifest['roots'][root_id]['core_state_sha256'],
                'complete_state':state.tolist(),'complete_state_bytes_sha256':__import__('hashlib').sha256(state.tobytes()).hexdigest(),
                'baseline_rhs':ev.rhs.tolist(),'allocation_only_rhs':allocated.rhs.tolist(),
                'baseline_residuals':residual_summary(ev),'accepted_residuals':accepted,
                'allocation_only_whole_cell_parameters':asdict(modified.parameters),
                'allocation_only_ae4_parameters':asdict(modified.ae4_parameters),
                'accepted_whole_cell_parameters':asdict(model.parameters) if feasible else None,
                'accepted_ae4_parameters':asdict(model.ae4_parameters) if feasible else None,
                'regulatory_model':asdict(model.regulatory_model),'stimulus':asdict(model.stimulus),
                'lp_diagnostics':solvers,'K_certificate':certificate,
                'full_source_matrix':matrix.tolist(),'allocation_source_change':delta.tolist(),
                'unconstrained_minimum_change_multipliers':dict(zip(AXES,signed.tolist())),
                'source_affine_family':'All signed solutions are this minimum-change vector plus null(M); the nonnegative/source-constrained family is empty for rejected endpoints.'}
            write_json(OUT/'fixed_wt_payloads'/f'{candidate_id}.json',payload)
            decisions.append(row)
    write_rows(OUT/'fixed_wt_solutions.csv',decisions)
    print(json.dumps({'conditions':len(decisions),'feasible':sum(r['fixed_wt_feasible'] for r in decisions),
        'modified_feasible':sum(r['fixed_wt_feasible'] and r['target_share'] is not None for r in decisions)}))


def freeze():
    assert not (OUT/'frozen_fixed_wt_manifest.json').exists()
    check_inputs()
    payloads=sorted((OUT/'fixed_wt_payloads').glob('*.json'))
    assert len(payloads)==30
    data=[read(p) for p in payloads]
    assert {(p['row']['root_id'],p['row']['condition']) for p in data}=={
        (r,c) for r in read(OUT/'contract.json')['root_ids'] for c in CONDITIONS}
    hashes={relative(p):sha256_file(p) for p in OUT.rglob('*') if p.is_file()}
    for p in (Path(__file__),REPO/'tests/test_task17_fixed_wt.py',ANALYSIS/'method.md'):
        hashes[relative(p)]=sha256_file(p)
    write_json(OUT/'frozen_fixed_wt_manifest.json',{
        'manifest_id':'TASK17_FIXED_WT_PRE_GENOTYPE_V1','created_utc':now(),
        'pre_freeze_head':git('rev-parse','HEAD'),'branch':BRANCH,'scientific_baseline':BASE,
        'contract_sha256':sha256_file(OUT/'contract.json'),'hashes':hashes,
        'all_decisions':[p['row'] for p in data],
        'feasible_parameterizations':[{'candidate_id':p['row']['candidate_id'],
            'payload_path':relative(path),'payload_sha256':sha256_file(path),
            'whole_cell_parameters_sha256':sha256_object(p['accepted_whole_cell_parameters']),
            'ae4_parameters_sha256':sha256_object(p['accepted_ae4_parameters']),
            'core_state_sha256':p['core_state_sha256']} for path,p in zip(payloads,data) if p['row']['fixed_wt_feasible']],
        'genotype_outcomes_used':False,'pushed_checkpoint_required_before_genotype':True,
        'complete_feasible_family':'Only the ten unchanged inherited baselines; all twenty modified endpoints certified infeasible in an outer relaxation of nuisance bounds.'})
    print('Fixed-WT checkpoint frozen; remote push required before any genotype evaluation.')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('stage',choices=('prepare','solve','freeze'))
    args=parser.parse_args()
    {'prepare':prepare,'solve':solve,'freeze':freeze}[args.stage]()


if __name__=='__main__':
    main()
