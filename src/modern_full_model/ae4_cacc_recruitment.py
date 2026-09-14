"""Target-directed AE4-dependent CaCC recruitment hypothesis (Task 41).

This phenomenological dependency is not established by AE4 stoichiometry.
The conductance is changed before the existing coupled current/NBC closure.
No secretion output, ionic source vector, or current is rescaled afterward.
"""
from dataclasses import dataclass,replace
from functools import lru_cache
import math
from .model import ModernFullModel,WT
from .nbc_minimal import MinimalNbcModel,normalized_secretory_activation
from .nkcc_stimulation import StimulatedNkcc1Model


@dataclass(frozen=True)
class Ae4CaccRecruitmentParameters:
    null_stimulated_fraction: float

    def __post_init__(self):
        if not math.isfinite(self.null_stimulated_fraction) or not 0<self.null_stimulated_fraction<=1:
            raise ValueError('Residual CaCC recruitment must be in (0,1]')


def recruitment_factor(activation,ae4_expression,residual_fraction):
    if not all(math.isfinite(v) for v in (activation,ae4_expression,residual_fraction)):
        raise ValueError('Recruitment arguments must be finite')
    if not 0<=activation<=1 or ae4_expression<0 or not 0<residual_fraction<=1:
        raise ValueError('Recruitment arguments outside their domain')
    # These branches preserve exact parent nesting, including roundoff.
    if activation==0 or ae4_expression==1 or residual_fraction==1:return 1.
    return (1-activation)+activation*(residual_fraction+(1-residual_fraction)*ae4_expression)


class Ae4DependentCaccModel:
    def __init__(self,base_model,parameters):
        self.base_model=base_model
        self.cacc_recruitment=parameters

    def __getattr__(self,name):
        return getattr(self.base_model,name)

    @lru_cache(maxsize=16)
    def effective_model(self,factor):
        if factor==1:return self.base_model
        template=self.base_model;core=template.base_model.base_model
        parameters=replace(core.parameters,membranes=replace(core.parameters.membranes,
            g_cl_apical_S=core.parameters.membranes.g_cl_apical_S*factor))
        new=ModernFullModel(parameters,stimulus=core.stimulus,
            regulatory_model=core.regulatory_model,ae4_parameters=core.ae4_parameters,
            ae4_evaluator=core.ae4_evaluator,nkcc1_kinetics=core.nkcc1_kinetics)
        return MinimalNbcModel(StimulatedNkcc1Model(new),template.nbc_parameters)

    def evaluate(self,time_s,vector,*,genotype=WT):
        activation=normalized_secretory_activation(self.stimulus(time_s).calcium_uM,self.nbc_parameters)
        factor=recruitment_factor(activation,genotype.ae4_expression,self.cacc_recruitment.null_stimulated_fraction)
        model=self.effective_model(factor)
        e=model.evaluate(time_s,vector,genotype=genotype)
        r=dict(e.diagnostics.regulatory)
        r.update(ae4_dependent_cacc_activation=float(activation),
            ae4_dependent_cacc_recruitment_factor=float(factor),
            effective_cacc_conductance_S=model.parameters.membranes.g_cl_apical_S,
            ae4_cacc_hypothesis_status='TARGET_SELECTED_UNVALIDATED_RECRUITMENT_DEPENDENCY')
        return replace(e,diagnostics=replace(e.diagnostics,regulatory=r))

    def rhs(self,time_s,vector,*,genotype=WT):
        return self.evaluate(time_s,vector,genotype=genotype).rhs
