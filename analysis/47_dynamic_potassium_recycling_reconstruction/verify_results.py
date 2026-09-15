"""Independent read only checks of archived states, scalar laws and balances.

Does not import or evaluate the production model. No integration or fitting.
"""
from pathlib import Path
import csv
import gzip
import hashlib
import json
import math
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def read(case,name):
    with gzip.open(HERE/'output'/case/(name+'.csv.gz'),'rt') as f:
        return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]

def main():
    freeze=json.loads((HERE/'PREDICTION_MANIFEST.json').read_text())
    for path,sha in freeze['sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha,path
    payload=json.loads((HERE/'input/active_parameters.json').read_text());p=payload['parameters']
    mp=p['membranes'];bath=p['bath'];pc=p['constants'];f=pc['faraday_C_mol']*1e-15
    vt=pc['gas_constant_J_mol_K']*pc['temperature_K']/pc['faraday_C_mol']
    compression=json.loads((HERE/'output/compression_manifest.json').read_text())
    for path,record in compression.items():
        blob=(HERE/path).read_bytes()
        assert hashlib.sha256(blob).hexdigest()==record['compressed_sha256']
        raw=gzip.decompress(blob)
        assert hashlib.sha256(raw).hexdigest()==record['uncompressed_sha256']
        assert len(raw)==record['uncompressed_bytes']
    report={}
    for case in ['wt','ae4_null','ae4_5pct']:
        rr,ss=read(case,'trajectory'),read(case,'states')
        assert len(rr)==len(ss)==602
        assert np.array_equal([r['time_s'] for r in rr],np.r_[0.,1e-6,np.arange(1.,601.)])
        maxerr={}; min_affinity=math.inf;min_reserve=math.inf
        for r,s in zip(rr,ss):
            assert r['time_s']==s['time_s']
            assert all(math.isfinite(v) for v in list(r.values())+list(s.values()))
            assert all(v>0 for key,v in s.items() if key!='time_s' and ('fmol' in key or 'pL' in key))
            ni,ki,cli=(s[ion+'_i_fmol']/s['volume_i_pL'] for ion in ['na','k','cl'])
            kl=s['k_l_fmol']/s['volume_l_pL']
            ca=payload['stimulus']['resting_calcium_uM' if r['time_s']==0 else 'stimulated_calcium_uM']
            gate=ca**mp['calcium_hill']/(ca**mp['calcium_hill']+mp['calcium_half_uM']**mp['calcium_hill'])
            ga=mp['g_k_total_S']*mp['apical_k_fraction']*gate
            gb=mp['g_k_total_S']*(1-mp['apical_k_fraction'])*gate
            ka=ga*(r['v_apical_V']-vt*math.log(kl/ki))/f
            kb=gb*(r['v_basolateral_V']-vt*math.log(bath['k_mM']/ki))/f
            ng=ni**3/(ni**3+mp['nak_na_half_mM']**3)
            cea=mp['nak_capacity_fmol_s']*mp['apical_pump_fraction']*kl**2/(kl**2+mp['nak_k_half_mM']**2)
            ceb=mp['nak_capacity_fmol_s']*(1-mp['apical_pump_fraction'])*bath['k_mM']**2/(bath['k_mM']**2+mp['nak_k_half_mM']**2)
            pump=(cea+ceb)*ng
            xx=ni*ki*cli**2
            nn=payload['nkcc1_kinetics']['alpha_eff_fmol_s']*r['nkcc1_activity_multiplier']*(157.5-2.0096e-5*xx)/(1.0306+1.3852e-6*xx)
            jj=r['cacc_cl_export_fmol_s'];aa=r['ae4_cl_inward_fmol_s'];hh=r['nhe1_inward_fmol_s'];bb=r['nbc_cycle_inward_fmol_s']
            errors=dict(na_concentration=ni-r['na_i_mM'],k_concentration=ki-r['k_i_mM'],cl_concentration=cli-r['cl_i_mM'],
                k_apical=ka-r['k_apical_efflux_fmol_s'],k_basal=kb-r['k_basolateral_efflux_fmol_s'],
                pump=pump-r['pump_total_cycles_fmol_s'],nkcc=nn-r['nkcc1_cycle_inward_fmol_s'],
                k_storage=nn+2*pump-aa/2-ka-kb-r['k_storage_fmol_s'],
                na_storage=nn+hh+bb-aa/2-3*pump-r['na_storage_fmol_s'],
                electrical_sum=ka+kb+pump+bb-jj,
                dynamic_cbm=ka+kb-(jj/3+nn+hh/3-aa/2)-( -2*r['na_storage_fmol_s']+r['ta_storage_fmol_s']+r['cl_storage_fmol_s'])/3+r['k_storage_fmol_s'])
            for key,val in errors.items():maxerr[key]=max(maxerr.get(key,0),abs(val))
            assert max(abs(v) for v in errors.values())<1e-10
            assert 0<pump<cea+ceb<=mp['nak_capacity_fmol_s']
            assert 0<ni<40 and 50<=ki<=200 and 30<=cli<=80 and 6.6<=r['ph_i']<=7.3
            assert 0<r['volume_i_pL']<3 and 0<r['hco3_i_mM']<100 and r['tic_i_mM']>0
            ideal=math.log(bath['na_mM']*bath['k_mM']*bath['cl_mM']**2/xx)
            assert nn*ideal>=0
            min_affinity=min(min_affinity,ideal);min_reserve=min(min_reserve,cea+ceb-pump)
        integral=float(np.trapezoid([r['q_out_pL_s'] for r in rr],[r['time_s'] for r in rr]))
        assert abs(integral-rr[-1]['cumulative_outflow_pL'])<1e-12
        summary=json.loads((HERE/'output'/case/'summary.json').read_text())
        assert summary['parameter_hashes']==freeze['parameter_hashes']
        assert summary['initial_state_sha256']==freeze['initial_state_sha256']
        report[case]=dict(saved_states_checked=len(rr),max_errors=maxerr,
                         min_ideal_nkcc_affinity=min_affinity,min_pump_reserve_fmol_s=min_reserve,
                         cumulative_quadrature_error_pL=abs(integral-rr[-1]['cumulative_outflow_pL']))
    derivative_errors=[]
    for case in ['wt','ae4_null']:
        with (HERE/'output'/case/'local_derivatives.csv').open() as handle:
            derivative_errors.extend(float(r['max_scaled_derivative_error']) for r in csv.DictReader(handle))
    assert len(derivative_errors)==24 and max(derivative_errors)<2e-5
    result=dict(status='PASS',frozen_files_unchanged=len(freeze['sha256']),lossless_archives_checked=len(compression),
                independent_saved_states_checked=1806,cases=report,
                local_derivative_checks=24,max_scaled_local_derivative_error=max(derivative_errors),
                production_integrations_in_this_verifier=0,production_model_imported=False)
    (HERE/'output/final_verification.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ('status','frozen_files_unchanged','independent_saved_states_checked','local_derivative_checks','max_scaled_local_derivative_error')}))

if __name__=='__main__':main()
