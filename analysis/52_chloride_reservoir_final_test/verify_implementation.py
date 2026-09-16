"""One recorded fixed-state Task52 verification; no trajectory or root solve."""
from pathlib import Path
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import json
import hashlib
import traceback
import time
TASK=Path(__file__).resolve().parent

def write(path,value): path.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')

def main():
    path=TASK/'output/verification_attempt_01.json'
    assert not path.exists(), 'Do not silently repeat an existing verification attempt'
    record={'status':'RUNNING','production_trajectories':0,'stationary_solves':0,'fits':0,'tests':[]}
    write(path,record);start=time.perf_counter()
    try:
        import numpy as np
        import scipy
        from threadpoolctl import threadpool_info
        import task52_model as m
        from project_onset import project_onset
        from diagnostics import diagnose,flux_values
        from dataclasses import asdict
        pools=threadpool_info();assert all(x['num_threads']==1 for x in pools)
        record['software']={'numpy':np.__version__,'scipy':scipy.__version__}
        record['thread_pools']=[{k:x[k] for k in ('internal_api','num_threads','version') if k in x} for x in pools]
        def passed(name,**evidence): record['tests'].append(dict(name=name,pass_=True,**evidence))
        audit=m.dependency_audit();assert audit['pass'];record['dependency_audit']=audit
        snapshot=json.loads((TASK/'output/frozen_task37_snapshot.json').read_text())
        for name,r in snapshot.items():
            assert hashlib.sha256((TASK.parents[1]/name).read_bytes()).hexdigest()==r['sha256']
        passed('exact_snapshot_and_no_forbidden_import',files=len(snapshot))
        original=m.accepted_state();parent=m.make_parent()
        state_records={}
        for mode in ('central','alternate'):
            state_records[mode]={}
            for genotype,Cl,pH in (('WT',50.10,6.91),('KO',36.50,6.89)):
                projection=project_onset(parent.parameters,original,Cl,pH,mode)
                y=np.array(projection['state']);local=m.SharedAuxiliaryModel(parent,0.)
                e=local.evaluate(0.,y,genotype=m.WT if genotype=='WT' else m.AE4_NULL)
                row,failures,ratios=diagnose(local,0.,y,e)
                assert not failures, (mode,genotype,failures)
                assert abs(row['cell_charge_fmol'])<=1e-9
                assert abs(row['cell_osmolarity_mOsm']-row['bath_osmolarity_mOsm'])<=1e-9
                assert abs(row['ph_i']-pH)<=1e-10
                assert abs(row['cl_i_mM']-Cl)<=1e-12
                assert projection['unchanged_lumen_and_regulatory_coordinates']
                if mode=='central': assert y[0]==original[0] and y[5]==original[5]
                else:
                    assert abs(y[0]/y[5]-original[0]/original[5])<1e-12
                    assert abs(y[3]/y[5]-original[3]/original[5])<1e-12
                projection['recovered_ph_i']=row['ph_i']
                projection['normal_speciation_and_physical_gates']='PASS'
                projection['initial_diagnostics']=row
                projection['max_conservation_ratio']=max(ratios.values())
                state_records[mode][genotype]=projection
                passed('projector_'+mode+'_'+genotype,cl_i_mM=row['cl_i_mM'],ph_i=row['ph_i'],
                       k_i_mM=row['k_i_mM'],tic_i_mM=row['tic_i_mM'],volume_i_pL=row['volume_i_pL'])
        write(TASK/'output/projected_onsets.json',state_records)
        ypair=np.r_[state_records['central']['WT']['state'],state_records['central']['KO']['state']]
        # Zero auxiliary: WT inherits the exact parent vector field.
        noaux=m.PairedModel(G_aux_S=0.)
        for t,a in ((0.,0.),(1e-6,0.),(600.,1.)):
            y=original.copy();y[12]=a
            ref=parent.evaluate(t,y,genotype=m.WT)
            wt,ko=noaux.evaluate(t,np.r_[y,y])
            assert np.array_equal(ref.rhs,wt.rhs)
            assert ref.diagnostics.water==wt.diagnostics.water
            assert ref.diagnostics.ae4==wt.diagnostics.ae4
            assert all(getattr(ko.diagnostics.ae4,k)==0. for k in ('na_cell_fmol_s','k_cell_fmol_s','cl_cell_fmol_s','tic_cell_fmol_s','alkalinity_cell_fmol_s'))
        passed('exact_WT_parent_nesting_and_KO_AE4_zero',comparisons=3)
        for scale in m.ALLOWED_SUPPLY:
            pair=m.PairedModel(G_aux_S=2.32e-9,ko_supply_multiplier=scale)
            wt,ko=pair.evaluate(1e-6,ypair)
            unforced=pair.local.evaluate(1e-6,ypair[13:],genotype=m.AE4_NULL)
            hw,hk,hl=wt.diagnostics.homeostasis,ko.diagnostics.homeostasis,unforced.diagnostics.homeostasis
            assert hk.nkcc1_inward_fmol_s==scale*hw.nkcc1_inward_fmol_s
            assert hk.ae2_inward_fmol_s==scale*hw.ae2_inward_fmol_s
            dN=hk.nkcc1_inward_fmol_s-hl.nkcc1_inward_fmol_s
            dE=hk.ae2_inward_fmol_s-hl.ae2_inward_fmol_s
            expected=np.r_[dN,dN,2*dN+dE,-dE,-dE,np.zeros(8)]
            assert np.allclose(ko.rhs-unforced.rhs,expected,rtol=0,atol=1e-14)
            assert ko.diagnostics.membranes==unforced.diagnostics.membranes
            assert ko.diagnostics.water==unforced.diagnostics.water
            assert hk.nhe1_inward_fmol_s==hl.nhe1_inward_fmol_s
            assert ko.diagnostics.regulatory['minimal_nbc_cycle_inward_fmol_s']==unforced.diagnostics.regulatory['minimal_nbc_cycle_inward_fmol_s']
            for label,y,e in (('WT',ypair[:13],wt),('KO',ypair[13:],ko)):
                _,failures,_=diagnose(pair,1e-6,y,e);assert not failures,(scale,label,failures)
            passed('paired_supply_exact_stoichiometry_and_local_laws_'+str(scale),N_WT=hw.nkcc1_inward_fmol_s,
                   N_KO=hk.nkcc1_inward_fmol_s,E_WT=hw.ae2_inward_fmol_s,E_KO=hk.ae2_inward_fmol_s)
        for kwargs in ({'ko_supply_multiplier':1.02},{'ko_supply_multiplier':1.2},{'G_aux_S':1e-9}):
            try: m.PairedModel(**kwargs)
            except ValueError: pass
            else: raise AssertionError('Unlicensed parameter accepted')
        passed('reject_unlicensed_supply_and_auxiliary_values')
        for G in m.ALLOWED_G:
            local=m.SharedAuxiliaryModel(parent,G)
            ew=local.evaluate(1e-6,ypair[:13],genotype=m.WT)
            ek=local.evaluate(1e-6,ypair[:13],genotype=m.AE4_NULL)
            assert ew.diagnostics.membranes==ek.diagnostics.membranes
            assert ew.diagnostics.regulatory['task52_g_aux_effective_S']==G
            for e in (ew,ek):
                r=e.diagnostics.regulatory; currents=e.diagnostics.membranes.currents_A
                assert abs(currents['cl_apical_total']-currents['cl_tmem16a']-currents['cl_auxiliary'])<1e-20
                assert currents['cl_auxiliary']==G*(e.diagnostics.membranes.v_apical_V-r['task52_E_Cl_V'])
                assert flux_values(local,e)['J_aux']*r['task52_Cl_driving_force_V']>=0
            off=local.evaluate(0.,ypair[:13]);ref=parent.evaluate(0.,ypair[:13])
            assert np.array_equal(off.rhs,ref.rhs) and off.diagnostics.membranes.currents_A['cl_auxiliary']==0.
            passed('shared_auxiliary_law_split_and_beta_zero_'+str(G))
        cch=m.PairedModel(protocol='CCH_ONLY',G_aux_S=2.32e-9)
        ew,ek=cch.evaluate(1e-6,ypair)
        cch0=m.PairedModel(protocol='CCH_ONLY',G_aux_S=0.)
        zw,zk=cch0.evaluate(1e-6,ypair)
        assert np.array_equal(ew.rhs,zw.rhs) and np.array_equal(ek.rhs,zk.rhs)
        passed('CCh_beta_zero_exact_paired_nesting')
        ipr=m.PairedModel(protocol='IPR_ONLY',G_aux_S=2.32e-9)
        ew,ek=ipr.evaluate(1e-6,ypair)
        for y,e in ((ypair[:13],ew),(ypair[13:],ek)):
            assert e.diagnostics.regulatory['minimal_nbc_activation_fraction']==0.
            assert e.diagnostics.membranes.currents_A['cl_auxiliary']!=0.
            _,failures,_=diagnose(ipr,1e-6,y,e);assert not failures,failures
        passed('IPR_auxiliary_current_in_zero_NBC_closure')
        # Two synthetic fixed states exercise both signs; no physiological
        # claim or production trajectory is made for these fixtures.
        signs=[]
        local=m.SharedAuxiliaryModel(parent,2.32e-9)
        for lumen_cl in (1.,1000.):
            y=ypair[:13].copy();y[8]=lumen_cl*y[11]
            e=local.evaluate(1e-6,y);d=e.diagnostics;f=flux_values(local,e)
            signs.append(float(np.sign(f['J_aux'])))
            total=f['J'];para=f['J_para_Cl'];mem=d.membranes
            assert abs(mem.cell_sources_fmol_s['cl']+total)<1e-12
            assert abs(mem.lumen_sources_fmol_s['cl']-total+para)<1e-12
            assert abs(mem.cell_sources_fmol_s['cl']+mem.lumen_sources_fmol_s['cl']+para)<1e-12
        assert set(signs)=={-1.,1.}
        passed('both_current_signs_and_cell_lumen_Cl_cancellation',signs=signs)
        record['dependency_audit']=m.dependency_audit();assert record['dependency_audit']['pass']
        record['counts']=dict(m.COUNTS)
        record['count_scope']='Top-level local/paired calls plus explicitly counted historical constructor core calls; nested wrapper reevaluations are not labelled as independent trajectories.'
        record['status']='PASS'
    except Exception as exc:
        record['status']='FAIL';record['error']=type(exc).__name__+': '+str(exc);record['traceback']=traceback.format_exc()
        if 'm' in locals():record['counts']=dict(m.COUNTS)
    record['elapsed_s']=time.perf_counter()-start
    record['source_sha256']={str(q.relative_to(TASK)):hashlib.sha256(q.read_bytes()).hexdigest()
                              for q in TASK.glob('*.py')}
    write(path,record)
    print(json.dumps({'status':record['status'],'tests_passed':len(record['tests']),'counts':record.get('counts'),
                      'error':record.get('error'),'elapsed_s':record['elapsed_s']}))
    if record['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
