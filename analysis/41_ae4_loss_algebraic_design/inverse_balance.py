"""Select one shared NKCC1 capacity using algebraic equations, no trajectories.

The temporary residual substitution below is the explicit candidate equation
for inverse design. The production model is implemented only after selection.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[k]='1'
from pathlib import Path
import importlib.util
import json
import sys
import time
from dataclasses import replace, asdict

ROOT=Path(__file__).resolve().parents[2]
HERE=ROOT/'analysis/41_ae4_loss_algebraic_design'
OUT=ROOT/'results/41_ae4_loss_algebraic_design'
spec=importlib.util.spec_from_file_location('task40_common',ROOT/'analysis/40_ae4_equal_cation_routing/validation_common.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
import numpy as np
from scipy.optimize import root
from scipy.special import expit
from threadpoolctl import threadpool_info
from modern_full_model.task30_nhe1_repair import (
    state_from_coordinates,COORDINATE_LOWER,COORDINATE_UPPER,INDEPENDENT_ROWS,INDEPENDENT_RHS_SCALES)
from modern_full_model.model import WT, AE4_NULL


def read_json(path):return json.loads(path.read_text())
def write_json(path,data):path.write_text(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+'\n')


def state_coordinates(model,y,genotype):
    e=model.evaluate(600.,y,genotype=genotype)
    c,l=e.diagnostics.observables.cell_concentrations_mM,e.diagnostics.observables.lumen_concentrations_mM
    return np.array([c['na'],c['k'],c['tic'],e.diagnostics.observables.cell_acid_base.ph,y[5],
        l['na'],l['k'],l['tic'],e.diagnostics.observables.lumen_acid_base.ph,y[11]])


def main():
    start_time=time.perf_counter()
    assert not (OUT/'capacity_design.json').exists(),'Do not overwrite an attempted design.'
    old.verify_frozen_sources()
    _,_,model,_,rest=old.frozen_models_and_state()
    cases=(('wt',WT),('ae4_null',AE4_NULL))
    yold={name:np.array(read_json(ROOT/'results/40_ae4_equal_cation_routing'/f'{name}_verification.json')['last_checked_state_vector']) for name,_ in cases}
    span=COORDINATE_UPPER-COORDINATE_LOWER
    starts=[state_coordinates(model,yold[name],genotype) for name,genotype in cases]
    target_deficit=.275
    ae4=model.evaluate(600.,yold['wt']).diagnostics.ae4.cl_cell_fmol_s
    # With matched composition/return and negligible storage, D ~ A/(2C+A).
    capacity_seed=ae4*(1-target_deficit)/(2*target_deficit)
    initial=np.r_[*[np.log((x-COORDINATE_LOWER)/(COORDINATE_UPPER-x)) for x in starts],np.log(capacity_seed)]
    calls=0

    def equations(z,details=False):
        nonlocal calls
        calls+=not details
        if calls>2000:raise RuntimeError('Single inverse-design call limit exceeded')
        capacity=float(np.exp(z[-1]));rows=[];records=[]
        for j,(name,genotype) in enumerate(cases):
            coords=COORDINATE_LOWER+span*expit(z[10*j:10*(j+1)])
            y=state_from_coordinates(model,coords);y[12:]=1.
            e=model.evaluate(600.,y,genotype=genotype)
            n=e.diagnostics.homeostasis.nkcc1_inward_fmol_s
            limited=float(np.clip(n,-capacity,capacity));delta=limited-n
            rhs=np.array(e.rhs);rhs[:3]+=np.array([delta,delta,2*delta])
            rows.extend(rhs[INDEPENDENT_ROWS]/INDEPENDENT_RHS_SCALES)
            q=e.diagnostics.water.lumen_outflow_pL_s
            records.append({'case':name,'state_vector':y.tolist(),'coordinates':coords.tolist(),
                'raw_rhs':rhs.tolist(),'q_out_pL_s':q,'unrestricted_nkcc1_fmol_s':float(n),
                'limited_nkcc1_fmol_s':limited,'cap_active':abs(n)>capacity,
                'ae4_cl_fmol_s':float(e.diagnostics.ae4.cl_cell_fmol_s),
                'observables':old.parent.inherited.row_for(model,600.,y,e),
                'conservation_residuals':dict(e.diagnostics.conservation_residuals)})
        rows.append(records[1]['q_out_pL_s']/records[0]['q_out_pL_s']-(1-target_deficit))
        return (np.array(rows),records) if details else np.array(rows)

    write_json(OUT/'design_intent.json',{'target_null_sustained_deficit':target_deficit,
        'target_finite_horizon_range':[.20,.35],'design_unknowns':21,'design_equations':21,
        'capacity_seed_from_literal_A_over_2C_plus_A':capacity_seed,'solver':'root/hybr',
        'parameter_selection_uses_target':True,'trajectories_before_selection':0,
        'scientific_workers':1,'thread_pools':[{k:p[k] for k in ('internal_api','num_threads')} for p in threadpool_info()]})
    fit=root(equations,initial,method='hybr',options={'xtol':1e-11,'maxfev':1900})
    residual,states=equations(fit.x,details=True)
    capacity=float(np.exp(fit.x[-1]))
    result={'solver_success':bool(fit.success),'solver_message':str(fit.message),
        'root_calls':1,'actual_residual_calls':calls,'max_abs_scaled_design_residual':float(np.max(np.abs(residual))),
        'capacity_fmol_s':capacity,'target_sustained_deficit':target_deficit,
        'actual_sustained_deficit':1-states[1]['q_out_pL_s']/states[0]['q_out_pL_s'],
        'states':states,'same_capacity_for_all_genotypes':True,'parameter_selection_uses_target':True,
        'production_trajectories_run':0,'wall_seconds':time.perf_counter()-start_time,
        'changes_task40_rest':capacity<rest['observables']['nkcc1_cycle_inward_fmol_s'],
        'parent_rest_nkcc1_fmol_s':rest['observables']['nkcc1_cycle_inward_fmol_s']}
    write_json(OUT/'capacity_design.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
