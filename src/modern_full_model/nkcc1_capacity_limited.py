"""Shared finite NKCC1 throughput hypothesis for target-directed Task 41.

The Palk law supplies the unconstrained signed demand. A shared finite cycle
capacity limits both directions and scales with NKCC1 expression. AE4
expression never enters this law. This is a phenomenological added capacity,
not an independently measured transporter parameter or an output correction.
"""
from dataclasses import dataclass, replace
import math
import numpy as np
from .model import WT


@dataclass(frozen=True)
class Nkcc1Capacity:
    cycle_capacity_fmol_s: float

    def __post_init__(self):
        if not math.isfinite(self.cycle_capacity_fmol_s) or self.cycle_capacity_fmol_s<=0:
            raise ValueError('Cycle capacity must be finite and positive')


def capacity_limited_cycles(demand,capacity,expression=1.):
    if not all(math.isfinite(x) for x in (demand,capacity,expression)) or capacity<=0 or expression<0:
        raise ValueError('Invalid capacity-limited cycle input')
    limit=capacity*expression
    return min(max(demand,-limit),limit)


class CapacityLimitedNkcc1Model:
    def __init__(self,base_model,parameters):
        self.base_model=base_model
        self.nkcc1_capacity=parameters

    def __getattr__(self,name):
        return getattr(self.base_model,name)

    def evaluate(self,time_s,vector,*,genotype=WT):
        e=self.base_model.evaluate(time_s,vector,genotype=genotype)
        h=e.diagnostics.homeostasis
        demand=h.nkcc1_inward_fmol_s
        n=capacity_limited_cycles(demand,self.nkcc1_capacity.cycle_capacity_fmol_s,genotype.nkcc1_expression)
        dn=n-demand
        h=replace(h,nkcc1_inward_fmol_s=n,
            na_cell_fmol_s=h.na_cell_fmol_s+dn,k_cell_fmol_s=h.k_cell_fmol_s+dn,
            cl_cell_fmol_s=h.cl_cell_fmol_s+2*dn)
        h=replace(h,charge_source_residual_fmol_s=math.fsum((h.na_cell_fmol_s,h.k_cell_fmol_s,-h.cl_cell_fmol_s,-h.alkalinity_cell_fmol_s)))
        rhs=np.array(e.rhs,copy=True);rhs[:3]+=[dn,dn,2*dn]
        residuals=dict(e.diagnostics.conservation_residuals)
        residuals['cell_bulk_charge_rate_fmol_s']=float(rhs[0]+rhs[1]-rhs[2]-rhs[4])
        residuals['homeostasis_charge_fmol_s']=h.charge_source_residual_fmol_s
        regulatory=dict(e.diagnostics.regulatory)
        regulatory.update(nkcc1_unrestricted_demand_fmol_s=float(demand),
            nkcc1_finite_capacity_fmol_s=self.nkcc1_capacity.cycle_capacity_fmol_s*genotype.nkcc1_expression,
            nkcc1_capacity_constraint_active=abs(demand)>self.nkcc1_capacity.cycle_capacity_fmol_s*genotype.nkcc1_expression)
        return replace(e,rhs=rhs,diagnostics=replace(e.diagnostics,homeostasis=h,
            regulatory=regulatory,conservation_residuals=residuals))

    def rhs(self,time_s,vector,*,genotype=WT):
        return self.evaluate(time_s,vector,genotype=genotype).rhs
