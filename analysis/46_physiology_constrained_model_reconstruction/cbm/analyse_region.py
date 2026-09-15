"""Distributed reaction groups using the unchanged, verified native VFFVA worker."""
from __future__ import annotations
import argparse, json, time, subprocess
from pathlib import Path
import build_cbm as b
from carbonscope_flux_core.solver import prepare_highs_flux_region, _ReusableWorker, run_highs_fba, run_highs_vffva
GROUPS={'chloride':['N','A','E','J','pCl','oCl'],'sodium':['H','B','Pb','Pa','pNa','oNa'],'potassium':['Kb','Ka','pK','oK'],'carbon':['Db','Da','pHCO3','oHCO3','oCO2']}

def group(genotype,name,M=100):
    start=time.monotonic(); prepared=prepare_highs_flux_region(b.build(genotype,M,direction="maximise"),None)
    assert prepared.retention is None
    # Calls the identical endpoint implementation exercised in Task46G. No altered solver,
    # reduced matrix, index renaming, retention constraint or new generic capability.
    worker=_ReusableWorker(prepared); out={k:{} for k in GROUPS[name]}
    for sense,label in [('min','minimum'),('max','maximum')]:
        worker.begin_pass(sense)
        for k in GROUPS[name]: out[k][label]=worker.solve_endpoint(b.IDS.index(k),sense)
    return {'genotype':genotype,'group':name,'input_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=b.AREA.parents[1],text=True).strip(),'builder_sha256':b.identity(),'M':M,'fraction_of_optimum':None,'native_api':'prepare_highs_flux_region plus unchanged _ReusableWorker.begin_pass/solve_endpoint','ranges':out,'endpoints':2*len(out),'elapsed_seconds':time.monotonic()-start}

def compact_witness(genotype):
    v={k:0. for k in b.IDS}
    v.update(N=.5 if genotype=='KO' else .25,A=0 if genotype=='KO' else .5,B=0 if genotype=='KO' else .5,Pb=1/6,Kb=5/6 if genotype=='KO' else 1/3,J=1,pNa=-1,oNa=1,oCl=1)
    assert b.validate(v,genotype)['valid']
    return v

def summary(genotype):
    maximum=run_highs_fba(b.build(genotype,direction="maximise")); minimum=run_highs_fba(b.build(genotype,direction="minimise"))
    w=compact_witness(genotype)
    # Single predeclared alternative computational M, no parameter grid.
    alt=run_highs_vffva(b.build(genotype,1000,direction="maximise"),None,workers=4,chunk_size=6)
    return {'genotype':genotype,'builder_sha256':b.identity(),'maximum_truncated_J':maximum.objective_value,'minimum_J':minimum.objective_value,'witness':w,'witness_validation':b.validate(w,genotype),'solver_maximum_witness':maximum.fluxes.to_dict(),'solver_maximum_validation':b.validate(maximum.fluxes.to_dict(),genotype),'solver_minimum_witness':minimum.fluxes.to_dict(),'maximum_active_bounds':[r.reaction_id+':'+side for r in b.build(genotype,direction="maximise").reactions for side,val in [('lower',r.lower_bound),('upper',r.upper_bound)] if abs(maximum.fluxes[r.reaction_id]-val)<1e-7],'alternative_M':1000,'alternative_ranges':alt.ranges.to_dict(orient='index'),'absolute_physical_maximum':'unidentified; the uncapped homogeneous model has an output-carrying ray, not a finite capacity prediction','ray_proof':'Multiply the reported complete witness by any positive scalar after removing computational box and J reporting truncation. All homogeneous balances and signs are preserved.','finite_witness_big_M_invariance':'Every non-output absolute flux is at most1, hence the reference-output witness remains feasible for all M>=1.','FVA_endpoint_vectors_are_not_joint_states':True}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--genotype',choices=['WT','KO'],required=True);p.add_argument('--group',choices=list(GROUPS)+['summary'],required=True);p.add_argument('--output',required=True);a=p.parse_args()
    result=summary(a.genotype) if a.group=='summary' else group(a.genotype,a.group)
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'genotype':a.genotype,'group':a.group,'output':a.output}))
