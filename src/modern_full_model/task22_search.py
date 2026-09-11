"""Bounded WT-only restoration and authoritative Task 22 production screens."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import math
import time

import numpy as np
from scipy.optimize import least_squares

from . import task22_contract as contract
from .task22_contract import (REPO, OUT, PHYSIOLOGY, CAPACITIES, CAPACITY_FOLD,
    OSMOTIC_RATIO, ANTI_CANCELLATION_MAX, read, require_stage0, write_json, now)
from .task21_calibration import WTProblem, SPEC, with_capacities, inverse
from .run_calcium_fast_screen import load_freeze
from .model import WT
from .water import compartment_osmolarity_mOsm, cell_impermeant_osmoles_fmol
from .validation import sha256_file, sha256_object


def physiology_failures(values):
    return [name for name, band in PHYSIOLOGY.items()
            if not math.isfinite(values[name]) or not band['lower'] <= values[name] <= band['upper']]


def capacity_failures(values, reference, fractions=()):
    failed = []
    for name, ref in reference.items():
        v = values[name]
        if not math.isfinite(v) or v < 0 or (ref == 0 and v != 0):
            failed.append(name)
        elif ref > 0 and not CAPACITY_FOLD[0] <= v / ref <= CAPACITY_FOLD[1]:
            failed.append(name)
    if any(not math.isfinite(v) or not 0 <= v <= 1 for v in fractions):
        failed.append('fractions')
    return failed


def osmotic_failures(osmolarities):
    bath = osmolarities['bath']
    if not math.isfinite(bath) or bath <= 0:
        return ['bath']
    return [side for side in ('cell', 'lumen') if not math.isfinite(osmolarities[side])
            or not OSMOTIC_RATIO[0] <= osmolarities[side] / bath <= OSMOTIC_RATIO[1]]


def anti_cancellation(fluxes, signed_chloride):
    pool = sum(max(v, 0.) for v in signed_chloride.values())
    finite = all(math.isfinite(v) for v in [*fluxes.values(), *signed_chloride.values()])
    ratios = {k: abs(v) / pool for k, v in fluxes.items()} if finite and pool > 0 else {}
    worst = max(ratios, key=ratios.get) if ratios else None
    return dict(positive_loading_pool_fmol_s=pool, signed_basolateral_cl_fmol_s=signed_chloride,
        absolute_flux_fmol_s={k:abs(v) for k,v in fluxes.items()}, ratios=ratios,
        worst_source=worst, maximum_ratio=ratios[worst] if worst else None,
        passes=bool(ratios and ratios[worst] <= ANTI_CANCELLATION_MAX))


def flux_ledger(ev, model):
    d = ev.diagnostics; h = d.homeostasis; a = d.ae4; m = d.membranes
    flux = {}
    for side in ('apical', 'basolateral'):
        j = getattr(m, 'pump_' + side + '_fmol_s')
        flux.update({f'pump_{side}_cycles':j, f'pump_{side}_Na':-3*j, f'pump_{side}_K':2*j})
    for name in ('na', 'k', 'cl', 'tic', 'alkalinity'):
        flux['AE4_' + name] = float(getattr(a, name + '_cell_fmol_s'))
    flux['AE4_cycles'] = a.cl_cell_fmol_s
    n, e, nh = h.nkcc1_inward_fmol_s, h.ae2_inward_fmol_s, h.nhe1_inward_fmol_s
    flux.update(NKCC1_cycles=n, NKCC1_Na=n, NKCC1_K=n, NKCC1_Cl=2*n,
        AE2_cycles=e, AE2_Cl=e, AE2_TIC=-e, AE2_alkalinity=-e,
        NHE1_cycles=nh, NHE1_Na=nh, NHE1_alkalinity=nh)
    conversion = 1e15 / model.parameters.constants.faraday_C_mol
    for name in ('k_apical', 'k_basolateral', 'cl_apical', 'para_na', 'para_k', 'para_cl', 'para_hco3'):
        valence = -1 if 'cl' in name or 'hco3' in name else 1
        flux[name] = m.currents_A[name] * conversion / valence
    flux.update({'CO2_' + k:v for k,v in d.co2_fluxes_fmol_s.items()})
    flux.update({'outflow_' + k:v for k,v in d.outflow_sources_fmol_s.items()})
    signed = dict(AE4=a.cl_cell_fmol_s, NKCC1=2*n, AE2=e)
    return dict(solute_fluxes_fmol_s=flux, water=asdict(d.water), **anti_cancellation(flux, signed))


def production_screen(problem, z):
    fast = problem.fast(z)
    model = with_capacities(problem.model, fast['capacities'])
    ev = model.evaluate(0., fast['state'], genotype=WT)
    d = ev.diagnostics; o = d.observables; ci = o.cell_concentrations_mM; li = o.lumen_concentrations_mM
    values = dict(cl_i_mM=ci['cl'], ph_i=o.cell_acid_base.ph, na_i_mM=ci['na'],
        k_i_mM=ci['k'], volume_i_pL=float(fast['state'][5]), na_l_mM=li['na'],
        k_l_mM=li['k'], cl_l_mM=li['cl'], ph_l=o.lumen_acid_base.ph,
        tic_l_mM=li['tic'], volume_l_pL=float(fast['state'][11]))
    closure = inverse.residual_summary(ev)
    failures = ['physiology:' + name for name in physiology_failures(values)]
    failures += ['osmolarity:' + side for side in osmotic_failures(o.osmolarities_mOsm)]
    fractions = [model.parameters.membranes.apical_k_fraction, model.parameters.membranes.apical_pump_fraction,
                 *fast['state'][12:]]
    failures += ['capacity:' + name for name in capacity_failures(fast['capacities'], problem.reference, fractions)]
    for name, tol in [('max_abs_scaled_independent_rhs', SPEC.root_scaled_tolerance),
                     ('max_abs_omitted_rhs', SPEC.omitted_row_raw_tolerance),
                     ('max_abs_current_A', SPEC.current_tolerance_A),
                     ('max_abs_state_charge_fmol', SPEC.charge_tolerance_fmol),
                     ('max_conservation_tolerance_ratio', 1.)]:
        if not math.isfinite(closure[name]) or closure[name] > tol:
            failures.append('closure:' + name)
    # The inherited regulatory equilibrium is fixed; every row must close.
    regulatory_max = float(np.max(np.abs(ev.rhs[12:]), initial=0.))
    if regulatory_max > SPEC.omitted_row_raw_tolerance:
        failures.append('closure:regulatory_rhs')
    a = d.ae4
    if not math.isfinite(a.diagnostics['affinity']) or a.diagnostics['affinity'] <= 0:
        failures.append('thermodynamics:nonpositive_affinity')
    if not math.isfinite(a.cl_cell_fmol_s) or a.cl_cell_fmol_s <= 0:
        failures.append('thermodynamics:nonpositive_AE4_loading')
    if a.na_cell_fmol_s * a.k_cell_fmol_s < 0:
        failures.append('mechanism:opposing_cations')
    if not np.all(np.isfinite(fast['state'])) or min(fast['state'][:12]) <= 0:
        failures.append('physical:state')
    ledger = flux_ledger(ev, model)
    if not ledger['passes']:
        failures.append('anti_cancellation')
    logs = [math.log(fast['capacities'][k] / problem.reference[k]) for k in CAPACITIES]
    return dict(complete_state=fast['state'].tolist(), whole_cell_parameters=asdict(model.parameters),
        ae4_parameters=asdict(model.ae4_parameters), observables=values,
        osmolarities_mOsm=o.osmolarities_mOsm,
        osmolarity_ratios={s:o.osmolarities_mOsm[s]/o.osmolarities_mOsm['bath'] for s in ('cell','lumen')},
        calibrated_capacities=fast['capacities'], reference_capacities=problem.reference,
        capacity_folds={k:fast['capacities'][k]/problem.reference[k] for k in CAPACITIES},
        production_diagnostics=asdict(d), complete_rhs=ev.rhs.tolist(), closure=closure,
        regulatory_max_abs_rhs=regulatory_max, flux_ledger=ledger,
        objective=[((values['cl_i_mM']-50.10)/1.50)**2+((values['ph_i']-6.91)/.07)**2,
                   max(abs(v) for v in logs), sum(v*v for v in logs)],
        resting_gate_failures=failures, resting_gate_pass=not failures,
        open_voltage_transcription_max_abs_error=float(max(abs(fast['rhs']-ev.rhs[:12]))))


def production_osmolarities(problem, fast):
    """Search calls the same osmolarity function as production, not a new sum."""
    g = problem.p.geometry; b = problem.p.bath
    result = {}
    for name in ('cell','lumen'):
        c = fast[name]
        impermeant = cell_impermeant_osmoles_fmol(other_impermeant_fmol=g.cell_impermeant_osmoles_fmol,
            finite_buffer_fmol=g.cell_buffer_total_fmol) if name == 'cell' else 0.
        result[name] = compartment_osmolarity_mOsm(na_mM=c['na'],k_mM=c['k'],cl_mM=c['cl'],
            total_carbon_mM=c['tic'], impermeant_osmoles_fmol=impermeant,volume_pL=c['volume'])
    result['bath'] = compartment_osmolarity_mOsm(na_mM=b.na_mM,k_mM=b.k_mM,cl_mM=b.cl_mM,
        total_carbon_mM=b.tic_mM,untracked_osmolyte_mM=b.untracked_osmolyte_mM)
    return result


def bounds(problem):
    b = problem.p.bath
    bath = compartment_osmolarity_mOsm(na_mM=b.na_mM,k_mM=b.k_mM,cl_mM=b.cl_mM,
        total_carbon_mM=b.tic_mM,untracked_osmolyte_mM=b.untracked_osmolyte_mM)
    # TIC_i upper and cation-ratio bounds are implied by the declared screens.
    # TIC_i > 0 has a numerical exploration floor only, never an acceptance band.
    lo = [math.log(47.10),6.77,math.log(1e-12),math.log(.3),math.log(2/200),
          math.log(80),6.8,math.log(1),math.log(.02),math.log(100/30)]
    hi = [math.log(53.10),7.05,math.log(1.2*bath),math.log(5),math.log(60/60),
          math.log(180),8.,math.log(80),math.log(.5),math.log(200/1)]
    # One-trillionth log inset avoids exp(log(endpoint)) rounding outside the hard box.
    lo += [math.log(.01)+1e-12]*15 + [-np.inf,-np.inf]
    hi += [math.log(100)-1e-12]*15 + [np.inf,np.inf]
    return np.array(lo), np.array(hi)


def inequalities(problem, fast):
    values = {f'{ion}_{s}_mM':fast[side][ion] for s,side in [('i','cell'),('l','lumen')]
              for ion in ('na','k','cl')}
    values.update(ph_i=fast['cell']['ph'],ph_l=fast['lumen']['ph'],
        volume_i_pL=fast['cell']['volume'],volume_l_pL=fast['lumen']['volume'],tic_l_mM=fast['lumen']['tic'])
    margins = []
    for name, band in PHYSIOLOGY.items():
        scale = band['upper']-band['lower'];v=values[name]
        margins.extend(((v-band['lower'])/scale,(band['upper']-v)/scale))
    osm = production_osmolarities(problem,fast)
    for s in ('cell','lumen'):
        ratio=osm[s]/osm['bath'];margins.extend(((ratio-.8)/.4,(1.2-ratio)/.4))
    margins.append(fast['affinity'])
    j,n,e,nh,pa,pb = fast['pathway_fluxes'][:6]
    pool=max(j,0)+max(2*n,0)+max(e,0)
    magnitudes=np.r_[abs(fast['pathway_fluxes']),abs(2*j),abs(2*n),abs(3*pa),abs(3*pb),
        [fast['flow']*abs(fast['lumen'][k]) for k in ('na','k','cl','tic','alkalinity')]]
    # Smooth signed rate inequality; every authoritative ledger component is also screened.
    margins.extend((100*pool-magnitudes)/(1+100*pool))
    return np.asarray(margins)


def seed_inventory():
    c = require_stage0(); manifest,_=load_freeze(); jobs=[]
    for i,root in enumerate(c['seed_root_ids']):
        jobs.append(dict(attempt_id=f'legacy_{i:02d}',root_id=root,kind='inherited_root',source=root))
    for i,path in enumerate(c['task21_witness_files']):
        witness=read(REPO/path)
        jobs.append(dict(attempt_id=f'witness_{i:02d}',root_id=witness['root_id'],kind='task21_resting_witness',
            source=path,source_sha256=sha256_file(REPO/path)))
    # Additional WT-only high-carbon, low-volume starts probe positive AE4 affinity.
    for i,root in enumerate(c['seed_root_ids']):
        jobs.append(dict(attempt_id=f'robust_{i:02d}',root_id=root,kind='WT_only_high_carbon',source=root))
    return jobs


def starting_coordinates(problem, job):
    z=problem.initial.copy()
    if job['kind']=='task21_resting_witness':
        w=read(REPO/job['source']);d=w['production_diagnostics'];o=d['observables']
        for offset,side,vi in [(0,'cell',5),(5,'lumen',11)]:
            c=o[side+'_concentrations_mM']
            z[offset:offset+5]=[math.log(c['cl']),o[side+'_acid_base']['ph'],math.log(c['tic']),
                math.log(w['complete_state'][vi]),math.log(c['na']/c['k'])]
        z[10:25]=[math.log(w['calibrated_capacities'][k]/problem.reference[k]) for k in CAPACITIES]
        z[25:]=[100*d['membranes']['v_apical_V'],100*d['membranes']['v_basolateral_V']]
    elif job['kind']=='WT_only_high_carbon':
        z[:10]=[math.log(50.1),7.04,math.log(20),math.log(1.),math.log(25/125),
            math.log(120),7.8,math.log(40),math.log(.15),math.log(140/10)]
    return z


def solve(job):
    require_stage0();manifest,_=load_freeze();p=WTProblem(manifest,job['root_id'])
    path=OUT/'attempts'/(job['attempt_id']+'.json')
    if path.exists():raise FileExistsError('Attempt already exists: '+str(path))
    lo,hi=bounds(p);original=starting_coordinates(p,job)
    start=np.clip(original,lo+1e-12,hi-1e-12)
    def residual(z):
        d=p.fast(z)
        return np.r_[d['equalities'],100*np.minimum(inequalities(p,d),0)]
    started=now();begin=time.monotonic()
    fit=least_squares(residual,start,bounds=(lo,hi),max_nfev=1200,
        ftol=1e-11,xtol=1e-11,gtol=1e-11,x_scale='jac')
    try:
        candidate=production_screen(p,fit.x);error=None
    except (ValueError,OverflowError,FloatingPointError,RuntimeError) as exc:
        candidate=None;error=f'{type(exc).__name__}: {exc}'
    result=dict(**job,started_utc=started,completed_utc=now(),elapsed_s=time.monotonic()-begin,
        stage0_checkpoint_sha=read(OUT/'stage0_checkpoint.json')['checkpoint_sha'],
        method='bounded_WT_only_least_squares_restoration',max_nfev=1200,nfev=fit.nfev,
        optimizer_success=bool(fit.success),optimizer_message=fit.message,
        original_seed_coordinates=original.tolist(),projected_seed_coordinates=start.tolist(),
        coordinates=fit.x.tolist(),restoration_norm=float(np.linalg.norm(fit.fun)),
        minimum_scaled_domain_margin=float(min(inequalities(p,p.fast(fit.x)))),
        active_coordinate_indices=np.flatnonzero((fit.x-lo<1e-7)|(hi-fit.x<1e-7)).tolist(),
        candidate=candidate,production_exception=error,genotype_evaluations=0,
        classification='RESTING_GATE_PASS' if candidate and candidate['resting_gate_pass'] else 'NUMERICAL_RESTORATION_DID_NOT_FIND_ADMISSIBLE_WT',
        global_optimality_claim=False)
    write_json(path,result)
    print(job['attempt_id'],fit.nfev,result['classification'],
        candidate['resting_gate_failures'] if candidate else error,flush=True)
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args();jobs=seed_inventory()
    write_json(OUT/'search_plan.json',dict(created_utc=now(),jobs=jobs,
        search_budget='27 independent bounded restorations, maximum 1200 function evaluations each. A proved domain-wide contradiction terminates further search; otherwise unresolved numerical failure cannot claim impossibility.',
        ranking_applies_to='Admissible solutions only; restoration merit is not biological ranking.',
        genotype_evaluations=0))
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(solve,jobs))


if __name__=='__main__':main()
