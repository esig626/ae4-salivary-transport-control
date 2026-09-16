"""The one predeclared VRAC law, inserted into the inherited current closure.

No conductance default is provided and no diagnostic value is licensed for
production. All inherited parameters, chemistry, water and transporter laws
are immutable. The exact Task13B gate is reused, not rediscovered.
"""
from __future__ import annotations
from dataclasses import replace
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,str(ROOT/'analysis/48_joint_experimental_constraint_reconstruction'))
from protocol_layer import build_model, GENOTYPES
from modern_full_model.membranes import hill_activation, nernst_voltage_V, current_to_fmol_s
from modern_full_model.nbc_minimal import evaluate_minimal_nbc
from modern_full_model.vbeta_diagnostic import VbetaParameters, VbetaTopology, vbeta_effective_conductance_S

FIELDS=('na','k','cl','tic','alkalinity')

class VracModel:
    def __init__(self, parent, *, conductance_S, reference_volume_pL):
        self.parent=parent
        self.vrac_parameters=VbetaParameters(
            topology=VbetaTopology.V1_RECTIFIED_SWELLING_GATE,
            maximum_conductance_S=conductance_S,
            reference_cell_volume_pL=reference_volume_pL)

    def __getattr__(self, name):
        return getattr(self.parent,name)

    def evaluate(self,t,y,*,genotype=GENOTYPES['WT']):
        base=self.parent.evaluate(t,y,genotype=genotype)
        d=base.diagnostics; old=d.membranes; p=self.parameters; mp=p.membranes
        gv= vbeta_effective_conductance_S(self.vrac_parameters,
            beta_input=d.stimulus.beta_input,cell_volume_pL=float(y[5]))
        reg=dict(d.regulatory)
        reg.update(vrac_effective_conductance_S=gv,
            vrac_rectified_swelling=max(float(y[5])/self.vrac_parameters.reference_cell_volume_pL-1,0),
            vrac_reference_volume_pL=self.vrac_parameters.reference_cell_volume_pL,
            vrac_conductance_scale_S=self.vrac_parameters.maximum_conductance_S,
            vrac_parameter_status='UNLICENSED_UNLESS_SEPARATELY_CALIBRATED')
        if gv==0.0:
            currents=dict(old.currents_A)
            currents.update(vrac_apical=0.0,cl_apical_total=currents['cl_apical'])
            reg.update(vrac_current_A=0.,vrac_chloride_export_fmol_s=0.)
            return replace(base,diagnostics=replace(d,regulatory=reg,
                membranes=replace(old,currents_A=currents)))

        obs=d.observables;ci=obs.cell_concentrations_mM
        li=obs.lumen_concentrations_mM;c=p.constants
        gate=hill_activation(d.stimulus.calcium_uM,mp.calcium_half_uM,mp.calcium_hill)
        gka=mp.g_k_total_S*gate*mp.apical_k_fraction+mp.g_apical_background_S
        gkb=mp.g_k_total_S*gate*(1-mp.apical_k_fraction)+mp.g_basolateral_background_S
        gcl=mp.g_cl_apical_S*gate
        gp_species=dict(para_na=mp.g_para_na_S,para_k=mp.g_para_k_S,
                       para_cl=mp.g_para_cl_S,para_hco3=mp.g_para_hco3_S)
        gp=sum(gp_species.values())
        ecl=nernst_voltage_V(ci['cl'],li['cl'],valence=-1,thermal_voltage_V=c.thermal_voltage_V)
        va0,vb0=old.v_apical_V,old.v_basolateral_V
        epsa,epsb=old.current_residuals_A['apical'],old.current_residuals_A['basolateral']
        denominator=gka+gcl+gv+gp
        def nbc(vb):
            return evaluate_minimal_nbc(na_i_mM=ci['na'],hco3_i_mM=ci['hco3'],
                na_o_mM=p.bath.na_mM,hco3_o_mM=obs.bath_acid_base.hco3_mM,
                v_basolateral_V=vb,thermal_voltage_V=c.thermal_voltage_V,
                faraday_C_mol=c.faraday_C_mol,
                activation_fraction=d.regulatory['minimal_nbc_activation_fraction'],
                parameters=self.nbc_parameters)
        def da(db):
            return (gp*db-gv*(va0-ecl)-epsa)/denominator
        def balance(db):
            return epsb+gkb*db+nbc(vb0+db).current_A-d.regulatory['minimal_nbc_current_A']+gp*(db-da(db))
        db=brentq(balance,-0.5-vb0,0.5-vb0,xtol=1e-14,rtol=1e-13)
        dva=da(db);dt=db-dva;va=va0+dva;vb=vb0+db
        newnbc=nbc(vb)
        iv=gv*(va-ecl)
        delta_i=dict(k_apical=gka*dva,cl_apical=gcl*dva,k_basolateral=gkb*db,
                     **{k:v*dt for k,v in gp_species.items()})
        currents=dict(old.currents_A)
        for k,v in delta_i.items():currents[k]+=v
        currents.update(vrac_apical=iv,nbc_basolateral=newnbc.current_A,
            cl_apical_total=currents['cl_apical']+iv)
        currents['apical_total']=sum(currents[k] for k in ['k_apical','cl_apical','pump_apical','vrac_apical'])
        currents['basolateral_old']=currents['k_basolateral']+currents['pump_basolateral']
        currents['basolateral_total']=currents['basolateral_old']+newnbc.current_A
        currents['paracellular_total']=sum(currents[k] for k in gp_species)
        def mol(i,z):return current_to_fmol_s(i,valence=z,faraday_C_mol=c.faraday_C_mol)
        dj={k:mol(v,-1 if k in ['cl_apical','para_cl','para_hco3'] else 1) for k,v in delta_i.items()}
        jv=mol(iv,-1)
        dB=newnbc.cycle_inward_fmol_s-d.regulatory['minimal_nbc_cycle_inward_fmol_s']
        dcell=np.array([dB,-dj['k_apical']-dj['k_basolateral'],-dj['cl_apical']-jv,2*dB,2*dB])
        dlumen=np.array([-dj['para_na'],dj['k_apical']-dj['para_k'],
            dj['cl_apical']+jv-dj['para_cl'],-dj['para_hco3'],-dj['para_hco3']])
        cs={k:old.cell_sources_fmol_s[k]+v for k,v in zip(FIELDS,dcell)}
        ls={k:old.lumen_sources_fmol_s[k]+v for k,v in zip(FIELDS,dlumen)}
        charge=lambda s:s['na']+s['k']-s['cl']-s['alkalinity']
        ia,ib,ip=[currents[k] for k in ['apical_total','basolateral_total','paracellular_total']]
        mem=replace(old,v_apical_V=va,v_basolateral_V=vb,v_transepithelial_V=vb-va,
            currents_A=currents,cell_sources_fmol_s=cs,lumen_sources_fmol_s=ls,
            current_residuals_A=dict(apical=ia-ip,basolateral=ib+ip),
            charge_residuals_fmol_s=dict(cell_source_minus_current=charge(cs)+mol(ia+ib,1),
                lumen_source_minus_current=charge(ls)-mol(ia-ip,1)))
        raw=np.array(base.rhs,copy=True);raw[:5]+=dcell;raw[6:11]+=dlumen
        carbon_external=(d.homeostasis.tic_cell_fmol_s+d.ae4.tic_cell_fmol_s+cs['tic']
            +d.co2_fluxes_fmol_s['bath_to_cell']+ls['tic']-d.outflow_sources_fmol_s['tic'])
        conservation=dict(d.conservation_residuals)
        conservation.update(cell_bulk_charge_rate_fmol_s=raw[0]+raw[1]-raw[2]-raw[4],
            lumen_bulk_charge_rate_minus_outflow_fmol_s=raw[6]+raw[7]-raw[8]-raw[10]+charge(d.outflow_sources_fmol_s),
            carbon_accounting_fmol_s=raw[3]+raw[9]-carbon_external,
            apical_current_A=ia-ip,basolateral_current_A=ib+ip,
            minimal_nbc_charge_source_fmol_s=-newnbc.cycle_inward_fmol_s,
            minimal_nbc_carbon_source_fmol_s=2*newnbc.cycle_inward_fmol_s,
            vrac_cell_plus_lumen_chloride_fmol_s=0.)
        reg.update(minimal_nbc_current_A=newnbc.current_A,
            minimal_nbc_cycle_inward_fmol_s=newnbc.cycle_inward_fmol_s,
            minimal_nbc_affinity_log=newnbc.affinity_log,
            vrac_current_A=iv,vrac_chloride_export_fmol_s=jv)
        return replace(base,rhs=raw,diagnostics=replace(d,membranes=mem,regulatory=reg,
            conservation_residuals=conservation))

    def rhs(self,t,y,*,genotype=GENOTYPES['WT']):
        return self.evaluate(t,y,genotype=genotype).rhs

    def solve_dynamics(self,time_span_s,*,initial_state,genotype=GENOTYPES['WT'],
                       method='BDF',rtol=1e-8,atol=1e-10,max_step_s=2.,t_eval=None):
        return solve_ivp(lambda t,y:self.rhs(t,y,genotype=genotype),time_span_s,
            initial_state,method=method,rtol=rtol,atol=atol,max_step=max_step_s,t_eval=t_eval)
