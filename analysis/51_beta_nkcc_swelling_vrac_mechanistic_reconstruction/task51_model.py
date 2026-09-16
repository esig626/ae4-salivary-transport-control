"""One additive beta-NKCC activity arm plus the unchanged conserved VRAC law.

No parameter defaults license production. Existing transport, regulation states,
source routing, current and water equations are delegated unchanged.
"""
from dataclasses import dataclass, replace
from pathlib import Path
import math
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,str(ROOT/'analysis/48_joint_experimental_constraint_reconstruction'))
sys.path.insert(0,str(ROOT/'analysis/49_camp_vrac_secretory_branch_reconstruction'))
from protocol_layer import build_model, GENOTYPES, Interventions
from vrac_model import VracModel

@dataclass(frozen=True)
class BetaNkccRegulation:
    parent: object
    beta_only_multiplier: float

    def __post_init__(self):
        if not math.isfinite(self.beta_only_multiplier) or self.beta_only_multiplier<1:
            raise ValueError('beta-only NKCC multiplier must be finite and >=1')

    def __getattr__(self,name):return getattr(self.parent,name)

    def evaluate(self,time_s,state,beta_input):
        old=dict(self.parent.evaluate(time_s,state,beta_input=beta_input))
        ca=float(old['nkcc1_capacity_multiplier'])
        extra=float(beta_input)*(self.beta_only_multiplier-1.)
        return {**old,'nkcc1_capacity_multiplier':ca+extra,
            'nkcc1_ca_capacity_multiplier':ca,'nkcc1_beta_increment':extra,
            'nkcc1_beta_only_multiplier':self.beta_only_multiplier,
            'nkcc1_beta_status':'CROSS_SPECIES_SOURCE_ASSUMPTION',
            'nkcc1_combination':'ADDITIVE_ABOVE_BASELINE_NO_INTERACTION'}

class Task51Model(VracModel):
    def __init__(self,parent,*,conductance_S,reference_volume_pL,parameter_status):
        super().__init__(parent,conductance_S=conductance_S,
                         reference_volume_pL=reference_volume_pL)
        self.task51_parameter_status=parameter_status

    def evaluate(self,t,y,*,genotype=GENOTYPES['WT']):
        ev=super().evaluate(t,y,genotype=genotype)
        reg={**ev.diagnostics.regulatory,'vrac_parameter_status':self.task51_parameter_status}
        return replace(ev,diagnostics=replace(ev.diagnostics,regulatory=reg))

def build_task51(arm,*,beta_only_multiplier,conductance_S,reference_volume_pL,
                 source_bath,interventions=Interventions(),parameter_status='UNLICENSED_TEST_OR_INFERENCE'):
    parent=build_model(arm,source_bath=source_bath,interventions=interventions)
    # ProtocolModel -> MinimalNbcModel -> StimulatedNkcc1Model -> core.
    # Installing in the shared regulator makes both homeostasis evaluations
    # consume the identical total activity, including under active calcium.
    core=parent.model.base_model.base_model
    inherited_names=tuple(core.regulatory_model.state_names)
    core.regulatory_model=BetaNkccRegulation(core.regulatory_model,beta_only_multiplier)
    assert tuple(core.regulatory_model.state_names)==inherited_names
    return Task51Model(parent,conductance_S=conductance_S,
        reference_volume_pL=reference_volume_pL,parameter_status=parameter_status)

def load_rest(track,genotype_name='WT'):
    import json
    if track=='track1':
        path=ROOT/'results/40_ae4_equal_cation_routing/wt_rest.json';key='state_vector'
    elif track=='track2':
        path=ROOT/'analysis/48_joint_experimental_constraint_reconstruction/output/rests'/f'{genotype_name}.json';key='state'
    else:raise ValueError('track1 or track2 required')
    return np.array(json.loads(path.read_text())[key],dtype=float)
