"""Algebraic calculations for the two analysis figures."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parent
P=json.loads((ROOT/'inputs.json').read_text())
VT=P['R_J_mol_K']*P['T_K']/P['F_C_mol']
JSTAR=P['c_star_mM']*P['V_star_pL']/P['t_star_s']
ISTAR=P['F_C_mol']*1e-15*JSTAR
h=1/(1+(P['calcium_half_uM']/P['calcium_uM'])**P['calcium_hill'])
GP=sum(P['G_para_S'].values())*VT/ISTAR
GKA=P['f_K']*P['G_K_S']*h*VT/ISTAR
GKB=(1-P['f_K'])*P['G_K_S']*h*VT/ISTAR
GCL=(P['G_Cl_S']*h+P['G_aux_S'])*VT/ISTAR
GA=GKA+GCL
GB_WT=P['WT']['J_B_fmol_s']/np.tanh(P['WT']['A_B']/P['w_B'])
GB_KO=P['KO']['J_B_fmol_s']/np.tanh(P['KO']['A_B']/P['w_B'])
GAMMA_B=GB_WT/JSTAR
H=10**-P['bath']['ph']; K1=10**-6.1; K2=10**-10.3
HCOE=P['bath']['tic']*K1*H/(H*H+K1*H+K1*K2)
BATH={**P['bath'],'hco3':HCOE}
LOWER_SLOPE=GKB+GP*GA/(GA+GP)

def coefficients(case):
    s=P[case]
    e_ka=np.log(s['k_l']/s['k_i'])
    e_kb=np.log(BATH['k']/s['k_i'])
    e_cl=-np.log(s['cl_l']/s['cl_i'])
    pa=s['I_pump_a_A']/ISTAR
    pb=s['I_pump_b_A']/ISTAR
    ca=-GKA*e_ka-GCL*e_cl+pa
    cp=0.0
    for ion,z in (('na',1),('k',1),('cl',-1),('hco3',-1)):
        ep=np.log(BATH[ion]/s[ion+'_l'])/z
        cp-=P['G_para_S'][ion]*VT/ISTAR*ep
    ab0=np.log(BATH['na']*HCOE**2/(s['na_i']*s['hco3_i']**2))
    return ca,cp,e_kb,pb,ab0

def psi_a(case,psi_b):
    ca,cp,_,_,_=coefficients(case)
    return (GP*np.asarray(psi_b)-ca+cp)/(GA+GP)

def residual(case,psi_b,recruitment=1.0):
    _,cp,e_kb,pb,ab0=coefficients(case)
    ps=np.asarray(psi_b)
    jb=recruitment*GAMMA_B*np.tanh((ab0+ps)/P['w_B'])
    return GKB*(ps-e_kb)+pb+jb+GP*(ps-psi_a(case,ps))+cp

def derivative(case,psi_b,recruitment=1.0):
    ab0=coefficients(case)[-1]
    x=(ab0+np.asarray(psi_b))/P['w_B']
    return LOWER_SLOPE+recruitment*GAMMA_B/P['w_B']*(1-np.tanh(x)**2)

def root(case,recruitment=1.0):
    return brentq(lambda ps:residual(case,ps,recruitment),-20,20,xtol=1e-14)

def carrier_flux(carrier_factor,rate_factor,case='WT'):
    s=P[case]; c=P['carrier']
    lcl=np.log(BATH['cl']/s['cl_i'])
    lna=np.log(s['na_i']/BATH['na'])+2*np.log(s['hco3_i']/HCOE)
    lk=np.log(s['k_i']/BATH['k'])+2*np.log(s['hco3_i']/HCOE)
    aplus=np.asarray(rate_factor)[...,None]*np.array([c['kCl_s_inv']*np.exp(lcl/2),c['kNa_s_inv']*np.exp(lna/2),c['kK_s_inv']*np.exp(lk/2)])
    aminus=np.asarray(rate_factor)[...,None]*np.array([c['kCl_s_inv']*np.exp(-lcl/2),c['kNa_s_inv']*np.exp(-lna/2),c['kK_s_inv']*np.exp(-lk/2)])
    num=aplus[...,0]*(aplus[...,1]+aplus[...,2])-aminus[...,0]*(aminus[...,1]+aminus[...,2])
    den=(aplus+aminus).sum(axis=-1)
    return np.asarray(carrier_factor)*num/den

def checks():
    report={'VT_V':VT,'J_star_fmol_s':JSTAR,'I_star_A':ISTAR,
      'NBC_capacity_from_WT_fmol_s':GB_WT,'NBC_capacity_from_KO_fmol_s':GB_KO,
      'NBC_capacity_recovery_difference':abs(GB_WT-GB_KO),
      'slope_lower_bound':LOWER_SLOPE}
    for case in ('WT','KO'):
        r=root(case)
        report[case]={'calculated_V_b_V':r*VT,'recorded_V_b_V':P[case]['V_b_V'],
          'V_b_error_V':abs(r*VT-P[case]['V_b_V']),
          'V_a_error_V':abs(psi_a(case,r)*VT-P[case]['V_a_V']),
          'minimum_tested_slope':float(derivative(case,np.linspace(-20,20,2001)).min())}
    factors=np.geomspace(.1,10,201)
    base=carrier_flux(1.,1.)
    report['maximum_relative_carrier_symmetry_error']=float(np.max(np.abs(carrier_flux(factors,1/factors)/base-1)))
    return report
