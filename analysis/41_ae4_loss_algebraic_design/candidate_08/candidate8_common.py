"""Frozen Task 40 plus one explicit AE4-dependent CaCC recruitment law."""
from candidate8_design import ROOT,DIRECTORY,OUTPUT,old,read_json,write_json
from dataclasses import asdict
import hashlib
import json
import numpy as np
from modern_full_model.model import WT
from modern_full_model.ae4_cacc_recruitment import Ae4CaccRecruitmentParameters,Ae4DependentCaccModel
from modern_full_model.validation import sha256_object,CONSERVATION_RESIDUAL_TOLERANCES

HERE=DIRECTORY;OUT=OUTPUT
PREPARED_HEAD='76f3144a4e662eb12d6da26c83ebf6f2ff47a59d'
BRANCH='codex/task-40-ae4-equal-cation-routing'
MODEL_ID='TASK40_AE4_DEPENDENT_CACC_RECRUITMENT'
parent=old.parent


def models_and_reference():
    d=read_json(OUT/'design.json');assert d['solver_success']
    _,rest,stim,y,record=old.frozen_models_and_state()
    p=Ae4CaccRecruitmentParameters(d['null_stimulated_cacc_fraction'])
    return rest,stim,Ae4DependentCaccModel(rest,p),Ae4DependentCaccModel(stim,p),y,record


def parameter_hashes(model):
    return {**old.parameter_hashes(model),'ae4_cacc_recruitment':sha256_object(model.cacc_recruitment)}


def diagnose(model,t,y,e):
    row,failures,ratios=old.diagnose(model,t,y,e)
    r=e.diagnostics.regulatory
    row.update(cacc_recruitment_factor=float(r['ae4_dependent_cacc_recruitment_factor']),
        effective_cacc_conductance_S=float(r['effective_cacc_conductance_S']))
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
