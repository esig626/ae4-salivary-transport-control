"""One frozen algebraic capacity; inherited Task 40 equations otherwise."""
from inverse_balance import ROOT, HERE, OUT, old, read_json, write_json
from dataclasses import asdict
import hashlib
import json
import numpy as np
from modern_full_model.model import WT
from modern_full_model.nkcc1_capacity_limited import Nkcc1Capacity,CapacityLimitedNkcc1Model
from modern_full_model.validation import sha256_object, CONSERVATION_RESIDUAL_TOLERANCES

PREPARED_HEAD='76f3144a4e662eb12d6da26c83ebf6f2ff47a59d'
BRANCH='codex/task-40-ae4-equal-cation-routing'
MODEL_ID='TASK40_WITH_SHARED_FINITE_NKCC1_CAPACITY'
parent=old.parent


def models_and_reference():
    design=read_json(OUT/'capacity_design.json')
    assert design['solver_success'] and design['max_abs_scaled_design_residual']<1e-7
    _,rest,stim,y,record=old.frozen_models_and_state()
    p=Nkcc1Capacity(design['capacity_fmol_s'])
    return rest,stim,CapacityLimitedNkcc1Model(rest,p),CapacityLimitedNkcc1Model(stim,p),y,record


def parameter_hashes(model):
    return {**old.parameter_hashes(model),'nkcc1_finite_capacity':sha256_object(model.nkcc1_capacity)}


def diagnose(model,t,y,e):
    row,failures,ratios=old.diagnose(model,t,y,e)
    r=e.diagnostics.regulatory
    row.update(nkcc1_unrestricted_demand_fmol_s=float(r['nkcc1_unrestricted_demand_fmol_s']),
        nkcc1_finite_capacity_fmol_s=float(r['nkcc1_finite_capacity_fmol_s']),
        nkcc1_capacity_constraint_active=float(r['nkcc1_capacity_constraint_active']))
    return row,failures,ratios


def verify_frozen_sources():
    old.verify_frozen_sources()
    f=read_json(OUT/'frozen_inputs.json')
    for x in f['execution_files']:
        assert hashlib.sha256((ROOT/x['path']).read_bytes()).hexdigest()==x['sha256'],x['path']
    return f


def frozen_models_and_state():
    f=verify_frozen_sources();_,_,rest,stim,_,_=models_and_reference()
    r=read_json(OUT/'wt_rest.json')
    assert r['admissible'] and r['state_sha256']==sha256_object(r['state_vector'])==f['frozen_rest_state_sha256']
    return r,rest,stim,np.array(r['state_vector']),r
