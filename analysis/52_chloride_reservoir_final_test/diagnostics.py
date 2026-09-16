"""Full conserved-coordinate readout and unchanged Task37 acceptance gates."""
from dataclasses import asdict
import math
import numpy as np
from task52_model import current_to_fmol_s
from task52_pre_palk.validation import CONSERVATION_RESIDUAL_TOLERANCES

LIMITS = (('na_i_mM',None,40.,False),('k_i_mM',50.,200.,True),
          ('cl_i_mM',30.,80.,True),('ph_i',6.6,7.3,True),
          ('hco3_i_mM',0.,100.,False),('tic_i_mM',0.,None,False),
          ('volume_i_pL',0.,3.,False))

def finite(value):
    if isinstance(value,dict): return all(finite(x) for x in value.values())
    if isinstance(value,(tuple,list)): return all(finite(x) for x in value)
    if isinstance(value,(float,int,np.number)): return bool(np.isfinite(value))
    return True

def flux_values(model, e):
    d=e.diagnostics; h=d.homeostasis; m=d.membranes; r=d.regulatory
    F=model.parameters.constants.faraday_C_mol
    flux=lambda k,z:current_to_fmol_s(m.currents_A[k],valence=z,faraday_C_mol=F)
    return dict(N=h.nkcc1_inward_fmol_s,E=h.ae2_inward_fmol_s,A=d.ae4.cl_cell_fmol_s,
        H=h.nhe1_inward_fmol_s,B=r['minimal_nbc_cycle_inward_fmol_s'],
        P=m.pump_apical_fmol_s+m.pump_basolateral_fmol_s,
        J=flux('cl_apical_total',-1),J_tmem16a=flux('cl_tmem16a',-1),
        J_aux=flux('cl_auxiliary',-1),K_apical=flux('k_apical',1),
        K_basolateral=flux('k_basolateral',1),J_para_Cl=flux('para_cl',-1),
        C=sum(d.co2_fluxes_fmol_s.values()),Q=d.water.lumen_outflow_pL_s)

def row_for(model,t,y,e):
    d=e.diagnostics; m=d.membranes; r=d.regulatory; f=flux_values(model,e)
    row={'time_s':float(t)}
    for prefix,c,ab in (('i',d.observables.cell_concentrations_mM,d.observables.cell_acid_base),
                         ('l',d.observables.lumen_concentrations_mM,d.observables.lumen_acid_base)):
        for key,value in c.items(): row[f'{key}_{prefix}_mM']=float(value)
        row['ph_'+prefix]=float(ab.ph)
        for key,value in asdict(ab).items(): row[f'speciation_{prefix}_{key}']=float(value)
    names=getattr(model,'single_state_names',getattr(model,'state_names',()))
    for name,value in zip(names,y): row[name]=float(value)
    for name,value in zip(names,e.rhs): row['d_'+name+'_dt']=float(value)
    row.update({k+'_fmol_s':float(v) for k,v in f.items() if k!='Q'})
    row.update(q_out_pL_s=float(f['Q']),volume_i_pL=float(y[5]),volume_l_pL=float(y[11]),
        E_Cl_V=float(r['task52_E_Cl_V']),V_a_V=float(m.v_apical_V),
        V_b_V=float(m.v_basolateral_V),V_trans_V=float(m.v_transepithelial_V),
        chloride_driving_force_V=float(r['task52_Cl_driving_force_V']),
        G_aux_effective_S=float(r['task52_g_aux_effective_S']),
        G_tmem16a_effective_S=float(r['task52_g_tmem16a_effective_S']),
        I_aux_A=float(m.currents_A['cl_auxiliary']),
        I_tmem16a_A=float(m.currents_A['cl_tmem16a']),
        I_total_Cl_A=float(m.currents_A['cl_apical_total']),
        I_nbc_A=float(r['minimal_nbc_current_A']),
        cell_charge_fmol=float(y[0]+y[1]-y[2]-y[4]-model.parameters.geometry.fixed_cell_anion_equivalents_fmol),
        cell_osmolarity_mOsm=float(d.observables.osmolarities_mOsm['cell']),
        lumen_osmolarity_mOsm=float(d.observables.osmolarities_mOsm['lumen']),
        bath_osmolarity_mOsm=float(d.observables.osmolarities_mOsm['bath']),
        chloride_concentration_slope_mM_s=float((e.rhs[2]-(y[2]/y[5])*e.rhs[5])/y[5]),
        auxiliary_component_residual_A=float(r['task52_apical_current_component_residual_A']))
    for key,value in asdict(d.water).items(): row['water_'+key]=float(value)
    for key,value in d.conservation_residuals.items(): row['residual_'+key]=float(value)
    for key,value in d.co2_fluxes_fmol_s.items(): row['co2_'+key+'_fmol_s']=float(value)
    for key,value in d.outflow_sources_fmol_s.items(): row['outflow_'+key+'_fmol_s']=float(value)
    for key,value in m.currents_A.items(): row['current_'+key+'_A']=float(value)
    for key in ('na_cell_fmol_s','k_cell_fmol_s','cl_cell_fmol_s','tic_cell_fmol_s','alkalinity_cell_fmol_s'):
        row['ae4_'+key]=float(getattr(d.ae4,key))
    for key in ('minimal_nbc_activation_fraction','minimal_nbc_affinity_log','nhe1_stimulation_multiplier',
                'task52_local_unimposed_nkcc1_cycles','task52_local_unimposed_ae2_cycles'):
        if key in r: row[key]=float(r[key])
    # Same numeric schema for both cells, with explicitly local unconstrained
    # law values retained as diagnostics only, never used to drive matched KO.
    row.setdefault('task52_local_unimposed_nkcc1_cycles',float(f['N']))
    row.setdefault('task52_local_unimposed_ae2_cycles',float(f['E']))
    expected=np.array([f['N']+f['H']+f['B']+d.ae4.na_cell_fmol_s-3*f['P'],
        f['N']+d.ae4.k_cell_fmol_s+2*f['P']-f['K_apical']-f['K_basolateral'],
        2*f['N']+f['A']+f['E']-f['J'],
        f['C']+2*f['B']-f['E']+d.ae4.tic_cell_fmol_s,
        f['H']+2*f['B']-f['E']+d.ae4.alkalinity_cell_fmol_s])
    row['independent_cell_assembly_max_abs_fmol_s']=float(np.max(np.abs(e.rhs[:5]-expected)))
    row['chloride_equation_residual_fmol_s']=float(e.rhs[2]-expected[2])
    row['tracked_cell_osmole_source_fmol_s']=float(4*f['N']+f['H']+3*f['B']-2*f['A']-f['P']-f['K_apical']-f['K_basolateral']-f['J']+f['C'])
    row['tracked_cell_osmole_residual_fmol_s']=float(sum(e.rhs[:4])-row['tracked_cell_osmole_source_fmol_s'])
    return row

def diagnose(model,t,y,e):
    row=row_for(model,t,y,e); failures=[]; ratios={}
    if not(np.all(np.isfinite(y)) and np.all(np.isfinite(e.rhs)) and finite(asdict(e.diagnostics))):
        failures.append({'gate':'finite_states_rhs_diagnostics'})
    for i,v in enumerate(y[:12]):
        if v<=0: failures.append({'gate':'positive_core','index':i,'value':float(v)})
    for key,low,high,inclusive in LIMITS:
        v=row[key]
        if (low is not None and (v<low if inclusive else v<=low)) or (high is not None and (v>high if inclusive else v>=high)):
            failures.append({'gate':key,'value':v,'lower':low,'upper':high,'inclusive':inclusive})
    residuals=e.diagnostics.conservation_residuals
    assert set(residuals)==set(CONSERVATION_RESIDUAL_TOLERANCES)|{'minimal_nbc_charge_source_fmol_s','minimal_nbc_carbon_source_fmol_s'}
    for key,tol in CONSERVATION_RESIDUAL_TOLERANCES.items():
        ratio=abs(float(residuals[key]))/tol;ratios[key]=ratio
        if not math.isfinite(ratio) or ratio>1:
            failures.append({'gate':key,'value':float(residuals[key]),'tolerance':tol})
    B=row['B_fmol_s']
    for key,value in (('minimal_nbc_charge_source_fmol_s',-B),('minimal_nbc_carbon_source_fmol_s',2*B)):
        if residuals[key]!=value: failures.append({'gate':key,'value':float(residuals[key]),'expected':value})
    for key,tol in (('auxiliary_component_residual_A',1e-20),('independent_cell_assembly_max_abs_fmol_s',1e-10),
                    ('tracked_cell_osmole_residual_fmol_s',1e-10)):
        if abs(row[key])>tol: failures.append({'gate':key,'value':row[key],'tolerance':tol})
    return row,failures,ratios
