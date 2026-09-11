"""WT-only Task 21 calibration; no genotype runner or outcomes are imported.

The fast residual is an open-voltage transcription of the production balances.
Production evaluation is authoritative for all retained candidates.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
import hashlib
import json
import math
import time

import numpy as np
from scipy.optimize import least_squares, minimize
from scipy.special import expit, logit

from .acid_base import _carbon_fractions, total_alkalinity_mM
from .calibration import CORE_INDEPENDENT_ROWS, WTCalibrationSpec
from .membranes import hill_activation
from .model import ModernFullModel, WT
from .nkcc_stimulation import StimulatedNkcc1Model
from .pooled_ae4 import MODEL_ID, select_ae4_mechanism
from .run_calcium_fast_screen import load_freeze, write_json, write_rows
from .task14_blind import slug
from . import task18_fixed_wt as inverse
from .task21_contract import REPO, OUT, CAPACITIES, PHYSIOLOGY, BRANCH, git
from .validation import sha256_file, sha256_object

STAGE0_SHA = '9820a6b62b29cd9895b9020b85ae89ad7d44b3fb'
STAGE0_BLOB = '3b390f50a50fe9a6dd40fad651c19c6eb2bf332c'
SPEC = WTCalibrationSpec()


def read(path):
    return json.loads(path.read_text())


def require_stage0():
    contract = read(OUT/'precalibration_contract.json')
    receipt = read(OUT/'stage0_checkpoint.json')
    if git('branch', '--show-current') != BRANCH:
        raise AssertionError('Wrong branch')
    if (receipt['checkpoint_sha'] != STAGE0_SHA or
            receipt['remote_contract_blob_sha'] != STAGE0_BLOB or
            receipt['remote_ref_verified'] is not True):
        raise AssertionError('Stage 0 remote checkpoint is not verified')
    data = (OUT/'precalibration_contract.json').read_bytes()
    if hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest() != STAGE0_BLOB:
        raise AssertionError('Precalibration contract changed')
    for path, digest in contract['source_hashes'].items():
        if sha256_file(REPO/path) != digest:
            raise AssertionError(f'Frozen source changed: {path}')
    if git('rev-parse', 'HEAD:archive') != contract['archive_tree_sha']:
        raise AssertionError('Archive changed')
    if (OUT/'wt_manifest.json').exists():
        raise AssertionError('WT ensemble is already frozen')
    return contract


def reference_capacities(model):
    return {**inverse.reference_capacities(model),
            'AE2': model.parameters.homeostasis.ae2_capacity_fmol_s}


def with_capacities(model, capacities):
    # Use the serialized/reloaded numeric representation for authoritative gates.
    capacities={name:float(value) for name,value in capacities.items()}
    intermediate = inverse.model_with_capacities(model, capacities)
    return StimulatedNkcc1Model(ModernFullModel(
        parameters=replace(intermediate.parameters, homeostasis=replace(
            intermediate.parameters.homeostasis,
            ae2_capacity_fmol_s=capacities['AE2'])),
        stimulus=intermediate.stimulus, regulatory_model=intermediate.regulatory_model,
        ae4_parameters=intermediate.ae4_parameters, ae4_evaluator=intermediate.ae4_evaluator))


def physiology_failures(observables):
    return [name for name, band in PHYSIOLOGY.items()
            if not (math.isfinite(observables[name]) and
                    band['lower'] <= observables[name] <= band['upper'])]


class WTProblem:
    """Ten neutral state coordinates, 15 log capacities, two voltage coordinates."""
    def __init__(self, manifest, root):
        legacy, self.seed, old = inverse.baseline(manifest, root)
        self.model = select_ae4_mechanism(legacy, MODEL_ID)
        self.root = root
        self.reference = reference_capacities(self.model)
        self.ref = np.array([self.reference[k] for k in CAPACITIES])
        self.p = self.model.parameters
        self.bath_ab = old.diagnostics.observables.bath_acid_base
        self.regulation = np.array(self.seed[12:])
        self.gate = hill_activation(old.diagnostics.stimulus.calcium_uM,
            self.p.membranes.calcium_half_uM, self.p.membranes.calcium_hill)
        self.initial = self.encode(old)

    def encode(self, evaluation):
        d = evaluation.diagnostics
        o = d.observables
        core = self.model.layout.decode(self.seed)
        result = []
        for c, ab, volume in ((o.cell_concentrations_mM, o.cell_acid_base, core.cell.volume_pL),
                              (o.lumen_concentrations_mM, o.lumen_acid_base, core.lumen.volume_pL)):
            result.extend((math.log(c['cl']), ab.ph, math.log(c['tic']), math.log(volume),
                           math.log(c['na']/c['k'])))
        return np.array(result + [0.]*15 + [d.membranes.v_apical_V*100,
                                           d.membranes.v_basolateral_V*100])

    def state(self, z):
        """Coordinates enforce both compartment charge identities exactly."""
        p = self.p
        arrays = []
        concentrations = []
        for offset, buffer, bpka, fixed in (
            (0,p.geometry.cell_buffer_total_fmol,p.acid_base.cell_buffer_pka,
             p.geometry.fixed_cell_anion_equivalents_fmol),
            (5,p.geometry.lumen_buffer_total_fmol,p.acid_base.lumen_buffer_pka,0.)):
            cl, ph, tic, volume, logodds = z[offset:offset+5]
            cl,tic,volume = np.exp([cl,tic,volume])
            ta = total_alkalinity_mM(ph=ph,total_carbon_mM=tic,
                buffer_total_mM=buffer/volume,buffer_pka=bpka,parameters=p.acid_base)
            pool = cl+ta+fixed/volume
            na = pool*expit(logodds)
            k = pool-na
            co2,hco3,co3 = np.array(_carbon_fractions(ph,p.acid_base))*tic
            c = dict(na=na,k=k,cl=cl,tic=tic,alkalinity=ta,co2=co2,hco3=hco3,
                     co3=co3,ph=ph,volume=volume)
            concentrations.append(c)
            arrays.extend(np.array([na,k,cl,tic,ta,1.])*volume)
        return np.r_[arrays,self.regulation], concentrations

    def fast(self, z):
        state,(ci,li) = self.state(z)
        p = self.p; b = p.bath; mp = p.membranes; w = p.water
        cap = dict(zip(CAPACITIES,self.ref*np.exp(z[10:25])))
        va,vb = z[25:27]/100; vt = vb-va
        thermal = p.constants.thermal_voltage_V
        cf = 1e15/p.constants.faraday_C_mol
        a = math.log(b.cl_mM/ci['cl'])+math.log((ci['na']+ci['k'])/(b.na_mM+b.k_mM))+2*math.log(ci['hco3']/self.bath_ab.hco3_mM)
        j = cap['AE4']*2*math.sinh(a/2)
        partition = ci['na']/(ci['na']+ci['k']) if j >= 0 else b.na_mM/(b.na_mM+b.k_mM)
        width = p.homeostasis.thermodynamic_saturation_log_width
        n = cap['NKCC1']*math.tanh(math.log(b.na_mM*b.k_mM*b.cl_mM**2/(ci['na']*ci['k']*ci['cl']**2))/width)
        h = cap['NHE1']*math.tanh(math.log(b.na_mM/ci['na'])/width + math.log(10)*(b.ph-ci['ph'])/width)
        e = cap['AE2']*math.tanh(math.log(b.cl_mM*ci['hco3']/(ci['cl']*self.bath_ab.hco3_mM))/width)
        na_factor = ci['na']**3/(ci['na']**3+mp.nak_na_half_mM**3)
        def pump(k,capacity):
            return capacity*na_factor*k**2/(k**2+mp.nak_k_half_mM**2)
        pa,pb = pump(li['k'],cap['pump_apical']),pump(b.k_mM,cap['pump_basolateral'])
        ka = (self.gate*cap['K_apical']+mp.g_apical_background_S)*(va-thermal*math.log(li['k']/ci['k']))*cf
        kb = (self.gate*cap['K_basolateral']+mp.g_basolateral_background_S)*(vb-thermal*math.log(b.k_mM/ci['k']))*cf
        cl = -self.gate*cap['CaCC']*(va+thermal*math.log(li['cl']/ci['cl']))*cf
        para = {}
        for ion,bath_value,zion in (('Na',b.na_mM,1),('K',b.k_mM,1),('Cl',b.cl_mM,-1),('HCO3',self.bath_ab.hco3_mM,-1)):
            para[ion] = cap['paracellular_'+ion]*(vt-thermal/zion*math.log(bath_value/li[ion.lower()]))*cf/zion
        co2b = cap['CO2_basolateral']*(self.bath_ab.co2_mM-ci['co2'])
        co2a = cap['CO2_apical']*(li['co2']-ci['co2'])
        oi = sum(ci[k] for k in ('na','k','cl','tic'))+(p.geometry.cell_impermeant_osmoles_fmol+p.geometry.cell_buffer_total_fmol)/ci['volume']
        ol = sum(li[k] for k in ('na','k','cl','tic'))
        ob = b.na_mM+b.k_mM+b.cl_mM+b.tic_mM+b.untracked_osmolyte_mM
        qb = w.basolateral_hydraulic_pL_s_mOsm*(oi-ob)
        qa = w.apical_hydraulic_pL_s_mOsm*(ol-oi)
        qt = w.paracellular_hydraulic_pL_s_mOsm*(ol-ob)
        q = w.outflow_rate_s*max(li['volume']-w.lumen_dead_volume_pL,0.)
        rhs = np.array([n+h-partition*j-3*(pa+pb),n-(1-partition)*j+2*(pa+pb)-ka-kb,
            2*n+e+j-cl,-e-2*j+co2b+co2a,h-e-2*j,qb-qa,
            3*pa-para['Na']-q*li['na'],-2*pa+ka-para['K']-q*li['k'],
            cl-para['Cl']-q*li['cl'],-para['HCO3']-co2a-q*li['tic'],
            -para['HCO3']-q*li['alkalinity'],qa+qt-q])
        ipara = para['Na']+para['K']-para['Cl']-para['HCO3']
        current = np.array([ka-cl+pa-ipara,kb+pb+ipara])
        eq = np.r_[rhs[list(CORE_INDEPENDENT_ROWS)]/SPEC.independent_rhs_scales,current/.02]
        obs = dict(cl_i_mM=ci['cl'],ph_i=ci['ph'],na_i_mM=ci['na'],k_i_mM=ci['k'],volume_i_pL=ci['volume'])
        return dict(state=state,cell=ci,lumen=li,capacities=cap,rhs=rhs,equalities=eq,
                    current_fmol_s=current,affinity=a,j_ae4=j,observables=obs,flow=q,
                    pathway_fluxes=np.array([j,n,e,h,pa,pb,ka,kb,cl,co2b,co2a,
                                             para['Na'],para['K'],para['Cl'],para['HCO3']]))

    def constraints(self, z, affinity_floor=1e-6):
        d = self.fast(z)
        inequalities = []
        for name, scale in (('na_i_mM',20.),('k_i_mM',120.)):
            v=d['observables'][name]; band=PHYSIOLOGY[name]
            inequalities.extend(((v-band['lower'])/scale,(band['upper']-v)/scale))
        inequalities.append(d['affinity']-affinity_floor)
        return np.asarray(inequalities)

    def production(self, z):
        d = self.fast(z)
        model = with_capacities(self.model,d['capacities'])
        ev = model.evaluate(0.,d['state'],genotype=WT)
        summary = inverse.residual_summary(ev)
        failures = physiology_failures(d['observables'])
        for name, tol in (('max_abs_scaled_independent_rhs',SPEC.root_scaled_tolerance),
            ('max_abs_omitted_rhs',SPEC.omitted_row_raw_tolerance),
            ('max_abs_current_A',SPEC.current_tolerance_A),
            ('max_abs_state_charge_fmol',SPEC.charge_tolerance_fmol),
            ('max_conservation_tolerance_ratio',1.)):
            if not math.isfinite(summary[name]) or summary[name]>tol: failures.append(name)
        ae4=ev.diagnostics.ae4
        if ae4.diagnostics['affinity']<=0: failures.append('nonpositive_affinity')
        if ae4.cl_cell_fmol_s<=0: failures.append('nonpositive_AE4_loading')
        if ae4.na_cell_fmol_s*ae4.k_cell_fmol_s<0: failures.append('opposing_cation_sources')
        if not np.all(np.isfinite(d['state'])) or np.min(d['state'][:12])<=0: failures.append('invalid_core')
        if not all(math.isfinite(v) and v>0 for v in d['capacities'].values()): failures.append('invalid_capacity')
        return dict(complete_state=d['state'].tolist(),observables=d['observables'],
            whole_cell_parameters=asdict(model.parameters),ae4_parameters=asdict(model.ae4_parameters),
            calibrated_capacities=d['capacities'],capacity_folds={k:d['capacities'][k]/self.reference[k] for k in CAPACITIES},
            reference_capacities=self.reference,production_diagnostics=asdict(ev.diagnostics),
            complete_rhs=ev.rhs.tolist(),closure=summary,resting_gate_failures=failures,
            resting_gate_pass=not failures,
            open_voltage_transcription_max_abs_error=float(np.max(np.abs(d['rhs']-ev.rhs[:12]))),
            objective=[float(((d['cell']['cl']-50.1)/1.5)**2+((d['cell']['ph']-6.91)/.07)**2),
                       float(np.max(np.abs(z[10:25]))),float(np.sum(z[10:25]**2))])

    def bounds(self, log_limit):
        # Numerical exploration limits, not biological/capacity feasibility walls.
        # All inherited capacity boxes are discarded; successive searches expand.
        lo=[math.log(47.1),6.77,-12,math.log(.3),-16,-12,3.01,-12,-12,-16]+[-log_limit]*15+[-100,-100]
        hi=[math.log(53.1),7.05,10,math.log(5),16,12,10.99,12,12,16]+[log_limit]*15+[100,100]
        return np.array(lo),np.array(hi)


def solve_seed(root):
    require_stage0()
    manifest,_=load_freeze()
    problem=WTProblem(manifest,root)
    path=OUT/'seed_attempts'/f'{slug(root)}.json'
    if path.exists(): return read(path)
    started=inverse.now(); attempts=[]; best=None
    for number,(limit,kind) in enumerate(((8.,'legacy'),(16.,'expanded'),(32.,'high_carbon_na'),(32.,'high_carbon_k'))):
        lower,upper=problem.bounds(limit)
        start=problem.initial.copy() if best is None else np.array(best.x)
        if kind.startswith('high_carbon'):
            start[:10]=[math.log(50.1),6.91,math.log(25.),math.log(.60),
                logit(.08 if kind.endswith('na') else .23),math.log(1200.),8.5,math.log(400.),
                math.log(20.),logit(.999 if kind.endswith('na') else .001)]
            start[10:25]=0.
            start[25:]=[-4.,2.]
        start=np.clip(start,lower+1e-8,upper-1e-8)
        def restoration(z):
            d=problem.fast(z)
            fit=[(d['cell']['cl']-50.1)/1.5,(d['cell']['ph']-6.91)/.07]
            return np.r_[d['equalities'],100*np.minimum(problem.constraints(z),0.),
                         .01*np.asarray(fit),1e-5*z[10:25]]
        begin=time.monotonic()
        fit=least_squares(restoration,start,bounds=(lower,upper),max_nfev=2500,
            ftol=1e-11,xtol=1e-11,gtol=1e-11,x_scale='jac')
        if best is None or np.linalg.norm(restoration(fit.x))<np.linalg.norm(restoration(best.x)):
            best=fit
        payload=problem.production(fit.x)
        attempts.append(dict(attempt=number,kind=kind,capacity_log_exploration_limit=limit,
            search_bounds_active=np.flatnonzero((fit.x-lower<1e-5)|(upper-fit.x<1e-5)).tolist(),
            optimizer_success=bool(fit.success),optimizer_message=fit.message,nfev=fit.nfev,
            elapsed_s=time.monotonic()-begin,restoration_norm=float(np.linalg.norm(fit.fun)),
            coordinates=fit.x.tolist(),candidate=payload))
        print(root,kind,payload['objective'],payload['closure']['max_abs_scaled_independent_rhs'],
              payload['resting_gate_failures'],flush=True)
        write_json(path.with_suffix('.progress.json'),dict(root_id=root,started_utc=started,attempts=attempts))
    result=dict(root_id=root,started_utc=started,completed_utc=inverse.now(),
        stage0_checkpoint_sha=STAGE0_SHA,legacy_state_role='numerical start and capacity reference only',
        attempts=attempts,status='RESTORATION_COMPLETE_PENDING_LEXICOGRAPHIC_REFINEMENT',
        genotype_results_used=False)
    write_json(path,result)
    path.with_suffix('.progress.json').unlink(missing_ok=True)
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--root');args=parser.parse_args()
    require_stage0(); manifest,_=load_freeze()
    roots=[args.root] if args.root else sorted(manifest['roots'])
    if args.workers==1:
        for root in roots: solve_seed(root)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            list(pool.map(solve_seed,roots))


if __name__=='__main__': main()
