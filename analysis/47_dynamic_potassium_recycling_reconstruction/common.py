"""Task 47: reconstruct frozen objects, never calibrate or alter production code."""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
from dataclasses import replace
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'src'))
from modern_full_model import parameters as par
from modern_full_model.model import ModernFullModel, WT, Genotype
from modern_full_model.transporters import AE4Parameters
from modern_full_model.ae4_equal_cation_routing import equal_routing_adapter
from modern_full_model.nkcc1_palk2010 import Nkcc1Kinetics
from modern_full_model.nkcc_stimulation import N1AlgebraicNkcc1, CompositeAe4NkccRegulation, StimulatedNkcc1Model
from modern_full_model.camp_pka import R1EffectiveActivation, RegulatoryGain, AE4Construct
from modern_full_model.nbc_minimal import MinimalNbcModel, MinimalNbcParameters
from modern_full_model.validation import SecretagogueProtocol, StimulusArm, sha256_object

BASE = 'acd554eb430ae25968ef56e3a9acdd816ca4e688'
BRANCH = 'analysis/task-47-dynamic-potassium-recycling-reconstruction'

def dump(path, data):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + '\n')

def csv_write(path, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def hashes(model):
    return {name: sha256_object(getattr(model, attr)) for name, attr in (
        ('whole_cell','parameters'), ('ae4','ae4_parameters'), ('nbc','nbc_parameters'),
        ('regulation','regulatory_model'), ('stimulus','stimulus'), ('nkcc1_core','nkcc1_kinetics'))}

def load_model():
    payload = json.loads((HERE/'input/active_parameters.json').read_text())
    types = dict(constants=par.PhysicalConstants, acid_base=par.AcidBaseParameters,
                 bath=par.BathParameters, geometry=par.GeometryParameters,
                 homeostasis=par.HomeostasisParameters, membranes=par.MembraneParameters,
                 water=par.WaterParameters, initial=par.InitialConditions)
    params = par.FullModelParameters(**{key: cls(**payload['parameters'][key]) for key, cls in types.items()})
    sp = dict(payload['stimulus'])
    sp['arm'] = StimulusArm(sp['arm'])
    stimulus = SecretagogueProtocol(**sp)
    rp = dict(payload['regulatory_model']['ae4_regulatory_model'])
    rp['construct'] = AE4Construct(rp['construct'])
    rp['gain'] = RegulatoryGain(**rp['gain'])
    regulation = CompositeAe4NkccRegulation(
        R1EffectiveActivation(**rp),
        N1AlgebraicNkcc1(**payload['regulatory_model']['nkcc1_regulatory_model']), stimulus)
    core = ModernFullModel(params, stimulus=stimulus, regulatory_model=regulation,
                          ae4_parameters=AE4Parameters(**payload['ae4_parameters']),
                          ae4_evaluator=equal_routing_adapter,
                          nkcc1_kinetics=Nkcc1Kinetics(**payload['nkcc1_kinetics']))
    model = MinimalNbcModel(StimulatedNkcc1Model(core), MinimalNbcParameters(**payload['nbc_parameters']))
    freeze = json.loads((ROOT/'results/40_ae4_equal_cation_routing/frozen_inputs.json').read_text())
    saved = json.loads((ROOT/'results/40_ae4_equal_cation_routing/wt_rest.json').read_text())
    assert hashes(model) == freeze['stimulus_parameter_hashes']
    y0 = np.array(saved['state_vector'])
    assert sha256_object(y0.tolist()) == freeze['frozen_rest_state_sha256']
    return model, y0

def clone_parameters(model, **membrane_changes):
    """Instantaneous derivative checks only; never a new production candidate."""
    old = model.base_model.base_model
    core = ModernFullModel(replace(old.parameters, membranes=replace(old.parameters.membranes, **membrane_changes)),
                          stimulus=old.stimulus, regulatory_model=old.regulatory_model,
                          ae4_parameters=old.ae4_parameters, ae4_evaluator=old.ae4_evaluator,
                          nkcc1_kinetics=old.nkcc1_kinetics)
    return MinimalNbcModel(StimulatedNkcc1Model(core), model.nbc_parameters)

def inherited_diagnose():
    spec = importlib.util.spec_from_file_location('task47_readonly_task40', ROOT/'analysis/40_ae4_equal_cation_routing/validation_common.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.diagnose

def source_hashes():
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted((ROOT/'src/modern_full_model').glob('*.py'))}
