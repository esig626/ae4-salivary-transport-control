"""Task 44: exact coordinate changes around the unchanged Task 40 model."""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import csv, hashlib, importlib.util, json, sys
from dataclasses import asdict, replace
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from modern_full_model.model import ModernFullModel, ConstantStimulus, Genotype, WT
from modern_full_model.nbc_minimal import MinimalNbcModel
from modern_full_model.nkcc_stimulation import StimulatedNkcc1Model
from modern_full_model.ae4_catalan2025 import with_mechanism_class, CLASSES
from modern_full_model.ae4_cacc_recruitment import Ae4DependentCaccModel, Ae4CaccRecruitmentParameters
from modern_full_model.validation import sha256_object
from modern_full_model.membranes import current_to_fmol_s
KEEP=np.array([0,1,2,3,5,6,7,8,9,11,12])
BASE_SHA='c4d4f207f702d8eee09ab94d2d90ae90552dd641'

def write_json(path,obj):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n')

def write_csv(path,rows):
    if not rows:return
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,keys,lineterminator="\n");w.writeheader();w.writerows(rows)

def reference():
    spec=importlib.util.spec_from_file_location('task40_readonly',ROOT/'analysis/40_ae4_equal_cation_routing/validation_common.py')
    task=importlib.util.module_from_spec(spec);spec.loader.exec_module(task)
    saved,_,_,rest,stim,_=task.models_and_reference()
    frozen=json.loads((ROOT/'results/40_ae4_equal_cation_routing/wt_rest.json').read_text())
    manifest=json.loads((ROOT/'results/40_ae4_equal_cation_routing/frozen_inputs.json').read_text())
    assert sha256_object(frozen['state_vector'])==manifest['frozen_rest_state_sha256']
    assert task.parameter_hashes(stim)==manifest['stimulus_parameter_hashes']
    audit=[]
    for item in json.loads((ROOT/'analysis/40_ae4_equal_cation_routing/parent_source_manifest.json').read_text())['files']:
        p=ROOT/item['path'];data=p.read_bytes()
        actual=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        if p.suffix=='.py' or item['path'].startswith('src/'):
            assert actual==item['git_blob_sha1'],item['path']
        audit.append(dict(path=item['path'],matches_task40=actual==item['git_blob_sha1'],actual_blob=actual))
    return rest,stim,np.array(frozen['state_vector']),task,audit

def clone(template,*,constant=True,parameters=None,ae4_parameters=None,nbc_parameters=None,nkcc1_kinetics=None,regulation=None,stimulus=None):
    old=template.base_model.base_model
    protocol=stimulus or (ConstantStimulus(.25,1.) if constant else old.stimulus)
    reg=replace(regulation or old.regulatory_model,stimulus=protocol)
    if parameters is not None and parameters.geometry.cell_buffer_total_fmol != old.parameters.geometry.cell_buffer_total_fmol:
        # Adjust only the unused constructor seed to the same fixed charge.
        from modern_full_model.acid_base import total_alkalinity_mM
        ii=parameters.initial; gg=parameters.geometry
        ta=total_alkalinity_mM(ph=ii.cell_ph,total_carbon_mM=ii.cell_tic_mM,buffer_total_mM=gg.cell_buffer_total_fmol/ii.cell_volume_pL,buffer_pka=parameters.acid_base.cell_buffer_pka,parameters=parameters.acid_base)
        cl=ii.cell_na_mM+ii.cell_k_mM-ta-gg.fixed_cell_anion_equivalents_fmol/ii.cell_volume_pL
        parameters=replace(parameters,initial=replace(ii,cell_cl_mM=cl))
    core=ModernFullModel(parameters or old.parameters,stimulus=protocol,regulatory_model=reg,
                        ae4_parameters=ae4_parameters or old.ae4_parameters,
                        ae4_evaluator=old.ae4_evaluator,nkcc1_kinetics=nkcc1_kinetics or old.nkcc1_kinetics)
    return MinimalNbcModel(StimulatedNkcc1Model(core),nbc_parameters or template.nbc_parameters)

def cacc(template,rho):
    if not 0<=rho<1:raise ValueError('rho must be in [0,1)')
    return Ae4DependentCaccModel(template,Ae4CaccRecruitmentParameters(1.-rho))

class Scaling:
    def __init__(self,model,rest):
        self.C0=model.parameters.bath.na_mM
        self.V0=float(rest[5]);self.L0=float(rest[11])
        self.n0=self.C0*self.V0;self.nl0=self.C0*self.L0
        self.J0=float(model.evaluate(1.,rest).diagnostics.homeostasis.nkcc1_inward_fmol_s)
        self.t0=self.n0/self.J0
        self.D=np.array([self.n0]*5+[self.V0]+[self.nl0]*5+[self.L0]+[1.])
        self.fixed=model.parameters.geometry.fixed_cell_anion_equivalents_fmol
    def encode(self,y):return np.asarray(y)/self.D
    def decode(self,x):return np.asarray(x)*self.D
    def independent(self,y):return self.encode(y)[KEEP]
    def expand(self,z):
        x=np.zeros(13);x[KEEP]=z
        x[4]=x[0]+x[1]-x[2]-self.fixed/self.n0
        x[10]=x[6]+x[7]-x[8]
        return self.decode(x)
    def rhs(self,model,tau,x,g):return model.rhs(tau*self.t0,self.decode(x),genotype=g)*self.t0/self.D
    def reduced_rhs(self,model,z,g):return model.rhs(1.,self.expand(z),genotype=g)[KEEP]*self.t0/self.D[KEEP]
    def payload(self):
        return dict(C0_mM=self.C0,V0_pL=self.V0,L0_pL=self.L0,n0_fmol=self.n0,nl0_fmol=self.nl0,
                    J0_fmol_s=self.J0,t0_s=self.t0,epsilon_lumen=self.L0/self.V0,
                    state_scales=self.D.tolist(),independent_rows=KEEP.tolist())

def observe(model,t,y,g=WT):
    e=model.evaluate(float(t),np.asarray(y),genotype=g);d=e.diagnostics
    c=d.observables.cell_concentrations_mM;l=d.observables.lumen_concentrations_mM
    m=d.membranes;h=d.homeostasis;r=d.regulatory;f=model.parameters.constants.faraday_C_mol
    vals=dict(time_s=float(t),q_pL_s=d.water.lumen_outflow_pL_s,
        na_i_mM=c['na'],k_i_mM=c['k'],cl_i_mM=c['cl'],T_i_mM=c['tic'],A_i_mM=y[4]/y[5],
        ph_i=d.observables.cell_acid_base.ph,hco3_i_mM=c['hco3'],V_i_pL=y[5],V_l_pL=y[11],
        cl_l_mM=l['cl'],ph_l=d.observables.lumen_acid_base.ph,
        osm_lumen_mOsm=d.observables.osmolarities_mOsm['lumen'],
        N=h.nkcc1_inward_fmol_s,H=h.nhe1_inward_fmol_s,E=h.ae2_inward_fmol_s,
        B=r['minimal_nbc_cycle_inward_fmol_s'],J4=d.ae4.cl_cell_fmol_s,
        P=m.pump_apical_fmol_s+m.pump_basolateral_fmol_s,
        CaCC=current_to_fmol_s(m.currents_A['cl_apical'],valence=-1,faraday_C_mol=f),
        paraCl=current_to_fmol_s(m.currents_A['para_cl'],valence=-1,faraday_C_mol=f),
        Va_V=m.v_apical_V,Vb_V=m.v_basolateral_V,
        charge_max_fmol=max(abs(v) for v in d.state_charge_fmol.values()),
        current_residual_A=max(abs(v) for v in m.current_residuals_A.values()),
        max_amount_rhs=max(abs(e.rhs[i]) for i in [0,1,2,3,4,6,7,8,9,10]),
        max_volume_rhs=max(abs(e.rhs[i]) for i in [5,11]))
    return {k:float(v) for k,v in vals.items()}

def physiology(row,y):
    fail=[]
    limits={'na_i_mM':(0.,40.),'k_i_mM':(50.,200.),'cl_i_mM':(30.,80.),'ph_i':(6.6,7.3),
            'hco3_i_mM':(0.,100.),'V_i_pL':(0.,3.)}
    for name,(lo,hi) in limits.items():
        if not lo<=row[name]<=hi:fail.append(name)
    if min(y[:12])<=0:fail.append('positive_core')
    return fail

def integrate(model,y0,g=WT,*,end=600.,rtol=1e-8,max_step=2.,method='Radau',dense=True,scaled=None):
    if scaled is None:
        initial=np.r_[y0,0.]
        def rhs(t,y):
            e=model.evaluate(t,y[:-1],genotype=g)
            return np.r_[e.rhs,e.diagnostics.water.lumen_outflow_pL_s]
        sol=solve_ivp(rhs,(0.,end),initial,method=method,rtol=rtol,atol=1e-11,max_step=max_step,dense_output=dense)
    else:
        initial=np.r_[scaled.encode(y0),0.]
        def rhs(tau,x):
            e=model.evaluate(tau*scaled.t0,scaled.decode(x[:-1]),genotype=g)
            return np.r_[e.rhs*scaled.t0/scaled.D,e.diagnostics.water.lumen_outflow_pL_s*scaled.t0/scaled.V0]
        sol=solve_ivp(rhs,(0.,end/scaled.t0),initial,method=method,rtol=rtol,
                      atol=np.r_[1e-11/scaled.D,1e-11/scaled.V0],max_step=max_step/scaled.t0,dense_output=dense)
    if not sol.success:raise RuntimeError(sol.message)
    return sol

def jacobian(fun,x,step=1e-5):
    cols=[]
    for k in range(len(x)):
        delta=np.zeros(len(x));delta[k]=step*max(abs(x[k]),.1)
        # The retained regulatory coordinate is physically in [0,1]. The
        # production evaluator clips out-of-domain trials; centring at r=1
        # would halve the column. Use a second-order inward derivative.
        if len(x)==len(KEEP) and k==len(x)-1 and x[k]+delta[k]>1.:
            cols.append((3*fun(x)-4*fun(x-delta)+fun(x-2*delta))/(2*delta[k]))
        elif len(x)==len(KEEP) and k==len(x)-1 and x[k]-delta[k]<0.:
            cols.append((-3*fun(x)+4*fun(x+delta)-fun(x+2*delta))/(2*delta[k]))
        else:
            cols.append((fun(x+delta)-fun(x-delta))/(2*delta[k]))
    return np.stack(cols,axis=-1)

def stationary(model,scale,seed,g=WT):
    zseed=scale.independent(seed);zseed[-1]=1.
    # Ten independent amount/volume unknowns, with exact regulatory equilibrium.
    def fun(logx):
        z=np.r_[np.exp(logx),1.]
        return scale.reduced_rhs(model,z,g)[:10]
    fit=root(fun,np.log(zseed[:10]),method='hybr',options={'xtol':1e-10,'maxfev':1500})
    z=np.r_[np.exp(fit.x),1.];y=scale.expand(z)
    residual=scale.reduced_rhs(model,z,g)
    good=np.max(np.abs(residual))<1e-7
    if not good:raise RuntimeError(f'Root residual {np.max(np.abs(residual)):.3g}: {fit.message}')
    j1=jacobian(lambda zz:scale.reduced_rhs(model,zz,g),z)
    j2=jacobian(lambda zz:scale.reduced_rhs(model,zz,g),z,5e-6)
    lam,vec=np.linalg.eig(j2/scale.t0)
    row=observe(model,1,y,g);fail=physiology(row,y)
    singular=np.linalg.svd(j2,compute_uv=False)
    result=dict(observables=row,state=y.tolist(),z=z.tolist(),residual_max=float(np.max(np.abs(residual))),
                root_success=bool(fit.success),root_message=str(fit.message),nfev=int(fit.nfev),
                physiological=not fail,failed_gates=fail,
                eigenvalues_per_s=[[float(v.real),float(v.imag)] for v in lam],
                max_real_eigenvalue_per_s=float(max(lam.real)),
                slowest_decay_s=float(-1./max(lam.real)) if max(lam.real)<0 else None,
                jacobian_relative_step_error=float(np.linalg.norm(j1-j2)/np.linalg.norm(j2)),
                jacobian_dimensionless=j2.tolist(),condition_number=float(singular[0]/singular[-1]),
                relative_smallest_singular_value=float(singular[-1]/singular[0]))
    return result
