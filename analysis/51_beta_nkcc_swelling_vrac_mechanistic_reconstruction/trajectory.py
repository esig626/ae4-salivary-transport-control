"""Authoritative orchestrator trajectory runner with immutable case records."""
import csv
import json
import time
import traceback
from dataclasses import asdict,replace
import numpy as np
import scipy
from scipy.integrate import solve_ivp,cumulative_trapezoid
from threadpoolctl import threadpool_info
from checkpoint import HERE,ROOT,OUT,dump,sha,integrity
from task51_model import build_task51,load_rest,GENOTYPES,Interventions
from modern_full_model.membranes import current_to_fmol_s
from modern_full_model.validation import CONSERVATION_RESIDUAL_TOLERANCES

def run_case(name,*,track,genotype_name,arm,beta_multiplier,conductance_S,
             tmem16a_off=False,purpose='INDEPENDENT_IDENTIFICATION',parameter_status='UNLICENSED_INFERENCE_TRIAL'):
    integrity()
    dest=OUT/'trajectories';dest.mkdir(exist_ok=True)
    path=dest/f'{name}.json'
    if path.exists():raise RuntimeError('Case already recorded; reuse, do not overwrite')
    start=time.perf_counter()
    y0=load_rest(track,'WT' if genotype_name=='AE4_5PCT' else genotype_name)
    genotype=replace(GENOTYPES['WT'],name='AE4_5PCT',ae4_expression=.05) if genotype_name=='AE4_5PCT' else GENOTYPES[genotype_name]
    solver=dict(method='Radau' if track=='track1' else 'BDF',rtol=1e-7 if track=='track1' else 1e-8,
                atol=[1e-10]*13,max_step=2.)
    solver['atol'][5]=solver['atol'][11]=1e-12
    record=dict(name=name,purpose=purpose,track=track,genotype=asdict(genotype),arm=arm,
        beta_multiplier=beta_multiplier,conductance_S=conductance_S,reference_volume_pL=float(y0[5]),
        initial_state=y0.tolist(),tmem16a_off=tmem16a_off,solver=solver,status='STARTED',
        parameter_status=parameter_status,model_sha256=sha(HERE/'task51_model.py'),
        runner_sha256=sha(HERE/'trajectory.py'),software=dict(numpy=np.__version__,scipy=scipy.__version__),
        blas_threads=[p['num_threads'] for p in threadpool_info()])
    assert all(n==1 for n in record['blas_threads'])
    dump(path,record)
    count=0
    try:
        m=build_task51(arm,beta_only_multiplier=beta_multiplier,conductance_S=conductance_S,
            reference_volume_pL=float(y0[5]),source_bath=track=='track2',
            interventions=Interventions(t16ainh_a01=tmem16a_off),parameter_status=parameter_status)
        record.update(parent_payload_hash=m.task48_payload_hash,stimulus=asdict(m.stimulus),bath=asdict(m.parameters.bath))
        original=m.evaluate
        def counted(*a,**kw):
            nonlocal count
            count+=1
            return original(*a,**kw)
        m.evaluate=counted
        sol=solve_ivp(lambda t,y:m.rhs(t,y,genotype=genotype),(1e-6,600.),y0,
            dense_output=True,**solver)
        record.update(solver_success=bool(sol.success),solver_message=sol.message,
            nfev=sol.nfev,njev=sol.njev,nlu=sol.nlu,accepted_steps=len(sol.t)-1)
        if not sol.success:raise RuntimeError(sol.message)
        times=np.r_[0.,1e-6,np.arange(1.,601.)]
        states=np.column_stack([y0,sol.sol(times[1:])])
        all_times=np.unique(np.r_[times,sol.t])
        maxima={k:0. for k in CONSERVATION_RESIDUAL_TOLERANCES}
        physical={k:[float('inf'),float('-inf')] for k in ['na_i_mM','k_i_mM','cl_i_mM','ph_i','volume_i_pL']}
        rows={};current_charge_max=0.;identity_max=0.
        F=m.parameters.constants.faraday_C_mol
        def mol(i,z=1):return current_to_fmol_s(i,valence=z,faraday_C_mol=F)
        for t in all_times:
            y=y0 if t==0 else sol.sol(t)
            if not np.all(np.isfinite(y)) or np.min(y[:12])<=0:raise RuntimeError('Nonpositive/nonfinite state')
            e=m.evaluate(float(t),y,genotype=genotype);d=e.diagnostics;o=d.observables;r=d.regulatory;mem=d.membranes
            ci=o.cell_concentrations_mM;c=mem.currents_A;h=d.homeostasis
            P=mem.pump_apical_fmol_s+mem.pump_basolateral_fmol_s
            row=dict(time_s=float(t),nkcc_multiplier=r['nkcc1_capacity_multiplier'],
                nkcc_ca_multiplier=r['nkcc1_ca_capacity_multiplier'],nkcc_beta_increment=r['nkcc1_beta_increment'],
                nkcc_cycles_fmol_s=h.nkcc1_inward_fmol_s,nkcc_cl_fmol_s=2*h.nkcc1_inward_fmol_s,
                ae4_cl_fmol_s=d.ae4.cl_cell_fmol_s,ae2_cl_fmol_s=h.ae2_inward_fmol_s,
                nhe_fmol_s=h.nhe1_inward_fmol_s,nbc_cycles_fmol_s=r['minimal_nbc_cycle_inward_fmol_s'],
                vrac_conductance_S=r['vrac_effective_conductance_S'],vrac_current_A=c['vrac_apical'],
                vrac_cl_fmol_s=r['vrac_chloride_export_fmol_s'],tmem16a_current_A=c['cl_apical'],
                tmem16a_cl_fmol_s=mol(c['cl_apical'],-1),total_cl_export_fmol_s=mol(c['cl_apical_total'],-1),
                na_i_mM=ci['na'],k_i_mM=ci['k'],cl_i_mM=ci['cl'],tic_i_mM=ci['tic'],hco3_i_mM=ci['hco3'],
                ph_i=o.cell_acid_base.ph,volume_i_pL=float(y[5]),relative_swelling=float(y[5]/y0[5]-1),
                v_apical_V=mem.v_apical_V,v_basolateral_V=mem.v_basolateral_V,v_trans_V=mem.v_transepithelial_V,
                pump_cycles_fmol_s=P,k_apical_fmol_s=mol(c['k_apical']),k_basolateral_fmol_s=mol(c['k_basolateral']),
                para_cl_return_fmol_s=mol(c['para_cl'],-1),lumen_osm_mOsm=o.osmolarities_mOsm['lumen'],
                volume_l_pL=float(y[11]),flow_pL_s=d.water.lumen_outflow_pL_s,
                chloride_derivative_mM_s=float((e.rhs[2]-ci['cl']*e.rhs[5])/y[5]))
            for k in physical:physical[k]=[min(physical[k][0],row[k]),max(physical[k][1],row[k])]
            for k,tol in CONSERVATION_RESIDUAL_TOLERANCES.items():
                v=abs(float(d.conservation_residuals[k]));maxima[k]=max(maxima[k],v)
                if v>tol:raise RuntimeError(f'Conservation failed {k}: {v} > {tol}')
                row[k]=float(d.conservation_residuals[k])
            current_charge_max=max(current_charge_max,*map(abs,mem.charge_residuals_fmol_s.values()))
            identity_max=max(identity_max,abs(row['total_cl_export_fmol_s']-(6*P-h.nhe1_inward_fmol_s+2*e.rhs[0]-e.rhs[4]-e.rhs[2])))
            if t in times:rows[float(t)]=row
        ordered=[rows[float(t)] for t in times]
        values=lambda k:np.array([r[k] for r in ordered])
        cumulative=cumulative_trapezoid(values('flow_pL_s'),times,initial=0)
        for row,q in zip(ordered,cumulative):row['cumulative_pL']=float(q)
        late=times>=300;window=times>=60
        ranges=dict(na_i_mM=(0,40),k_i_mM=(50,200),cl_i_mM=(30,80),ph_i=(6.6,7.3),volume_i_pL=(0,3))
        failures={k:physical[k] for k,(lo,hi) in ranges.items() if physical[k][0]<lo or physical[k][1]>hi}
        integrals={k:float(np.trapezoid(np.maximum(values(k)[window],0),times[window]))
            for k in ['nkcc_cl_fmol_s','ae4_cl_fmol_s','ae2_cl_fmol_s']}
        integrals.update({k:float(np.trapezoid(values(k)[window],times[window])) for k in
            ['vrac_cl_fmol_s','tmem16a_cl_fmol_s','total_cl_export_fmol_s','pump_cycles_fmol_s','nhe_fmol_s','nbc_cycles_fmol_s']})
        keep=[i for i,t in enumerate(times) if t in [0,1e-6,1,2,5] or t%10==0]
        csvpath=dest/f'{name}.csv'
        with csvpath.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(ordered[0]));writer.writeheader();writer.writerows(ordered[i] for i in keep)
        statepath=dest/f'{name}_states.csv'
        with statepath.open('w',newline='') as f:
            writer=csv.writer(f);writer.writerow(['time_s',*m.state_names]);writer.writerows([float(times[i]),*states[:,i]] for i in keep)
        record.update(status='COMPLETE',numerically_valid=True,physiologically_admissible=not failures,
            physiology_ranges=physical,physiology_departures=failures,checked_samples=len(all_times),
            conservation_maxima=maxima,conservation_tolerances=dict(CONSERVATION_RESIDUAL_TOLERANCES),
            membrane_charge_max_fmol_s=current_charge_max,total_export_identity_max_fmol_s=identity_max,
            swelling_mean_300_600=float(np.trapezoid(values('relative_swelling')[late],times[late])/300),
            swelling_endpoint=float(values('relative_swelling')[-1]),swelling_max=float(max(values('relative_swelling'))),
            volume_positive_samples=int(np.count_nonzero(values('relative_swelling')>0)),
            vrac_positive_samples=int(np.count_nonzero(values('vrac_conductance_S')>0)),
            secretion_cumulative_pL=float(cumulative[-1]),mean_flow_60_600_pL_s=float(np.trapezoid(values('flow_pL_s')[window],times[window])/540),
            initial_chloride_derivative_mM_s=float(values('chloride_derivative_mM_s')[1]),
            chloride_change_mM=float(values('cl_i_mM')[-1]-values('cl_i_mM')[0]),
            conditional_cl_slope_60_600_mM_s=float(np.polyfit(times[window]-330,values('cl_i_mM')[window],1)[0]),
            identified_spq_uptake=None,integrals_60_600=integrals,
            observations_60s=[r for r in ordered if r['time_s']%60==0],initial_right_limit=ordered[1],
            endpoint_state=states[:,-1].tolist(),csv_sha256=sha(csvpath),states_sha256=sha(statepath))
    except Exception:
        record.update(status='FAILED',numerically_valid=False,error=traceback.format_exc())
    record.update(core_evaluations=count+1,elapsed_seconds=time.perf_counter()-start)
    dump(path,record)
    print(json.dumps({k:record.get(k) for k in ['name','status','swelling_mean_300_600','physiologically_admissible','core_evaluations']}),flush=True)
    return record
