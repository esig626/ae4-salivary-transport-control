"""Read-only analysis of completed cohorts; never calls an ODE/root solver."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import csv,json,hashlib,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
HERE=ROOT/'analysis/41_ae4_loss_algebraic_design'
OUT=ROOT/'results/41_ae4_loss_algebraic_design'
FINAL=OUT/'candidate_08'
CASES=('wt','ae4_5pct','ae4_null')
LABELS=('WT','AE4 5%','AE4 null')


def read(p):return json.loads(p.read_text())
def write(p,x):p.write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n')
def csvread(p):
    with p.open() as f:return list(csv.DictReader(f))
def csvwrite(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)


def main():
    v={c:read(FINAL/f'{c}_verification.json') for c in CASES}
    assert all(x['status']=='PASS' and x['complete_600_s'] for x in v.values())
    g={c:{k:np.array([float(r[k]) for r in csvread(FINAL/f'{c}_dense_timeseries.csv')])
        for k in csvread(FINAL/f'{c}_dense_timeseries.csv')[0]} for c in CASES}
    times=g['wt']['time_s'];assert all(np.array_equal(g[c]['time_s'],times) for c in CASES)
    f={c:{r['quantity']:float(r['value']) for r in csvread(FINAL/f'{c}_integrated_fluxes.csv')} for c in CASES}
    summary=[];curve=[];fluxes=[];endpoints=[];conservation=[]
    for c in CASES:
        s=v[c]['secretion'];wc=v['wt']['secretion'];Q=s['cumulative_0_600_pL']
        D=100*(1-Q/wc['cumulative_0_600_pL']);R=Q/wc['cumulative_0_600_pL']
        dNa=g[c]['na_i_mM']-g['wt']['na_i_mM']
        dCl=g[c]['cl_i_mM']-g['wt']['cl_i_mM']
        positive=times>0
        running_D=100*(1-g[c]['cumulative_outflow_0_t_pL'][positive]/g['wt']['cumulative_outflow_0_t_pL'][positive])
        peak=int(np.argmax(running_D))
        summary.append({'case':c,'ae4_expression':v[c]['ae4_expression'],
            'cacc_stimulated_fraction':g[c]['cacc_recruitment_factor'][-1],
            'cumulative_0_600_pL':Q,'ratio_to_WT':R,'deficit_percent':D,
            'absolute_deficit_pL':wc['cumulative_0_600_pL']-Q,
            'endpoint_flow_pL_s':s['endpoint_pL_s'],'endpoint_deficit_percent':100*(1-s['endpoint_pL_s']/wc['endpoint_pL_s']),
            'mean_60_600_pL_s':s['mean_60_600_pL_s'],
            'mean_flow_relative_to_own_WT_rest':s['mean_60_600_pL_s']/s['q_rest_pL_s'],
            'NKCC1_compensation_percent':100*(f[c]['positive_nkcc1_cl_loading']/f['wt']['positive_nkcc1_cl_loading']-1),
            'AE2_positive_loading_share_percent':100*v[c]['chloride_partition']['positive_shares']['ae2'],
            'Na_i_difference_at_600_mM':float(dNa[-1]),
            'Na_i_difference_max_on_saved_grid_mM':float(max(dNa)),
            'Na_i_higher_than_WT_on_saved_grid':bool(np.any(dNa>1e-9)),
            'Cl_i_difference_at_600_mM':float(dCl[-1]),
            'peak_cumulative_deficit_percent':float(running_D[peak]),
            'peak_cumulative_deficit_time_s':float(times[positive][peak]),
            'within_20_35_percent':bool(20<=D<=35) if c!='wt' else None})
        for idx,t in enumerate(times):
            if t%60:continue
            W=g['wt']['cumulative_outflow_0_t_pL'][idx];C=g[c]['cumulative_outflow_0_t_pL'][idx]
            curve.append({'case':c,'time_s':float(t),'q_out_pL_s':float(g[c]['q_out_pL_s'][idx]),
                'cumulative_pL':float(C),'WT_cumulative_pL':float(W),'deficit_pL':float(W-C),
                'deficit_percent':float(100*(1-C/W)) if W>0 else None})
        for k,val in f[c].items():
            fluxes.append({'quantity':k,'case':c,'integral_60_600':val,
                'difference_from_WT':val-f['wt'][k],
                'percent_change_from_WT':100*(val/f['wt'][k]-1) if f['wt'][k]!=0 else None})
        for field in ('na_i_mM','k_i_mM','cl_i_mM','tic_i_mM','hco3_i_mM','ph_i','volume_i_pL',
            'v_apical_V','v_basolateral_V','v_transepithelial_V'):
            endpoints.append({'quantity':field,'case':c,'initial':float(g[c][field][0]),
                'final_600':float(g[c][field][-1]),'minimum_all_monitored':v[c]['minima'][field],
                'maximum_all_monitored':v[c]['maxima'][field]})
        for field,x in v[c]['max_abs_conservation_residuals'].items():
            conservation.append({'case':c,'quantity':field,'max_abs_residual':x,
                'tolerance':v[c]['conservation_tolerances'][field],
                'max_ratio':v[c]['max_conservation_ratios'][field]})
    assert all(s['within_20_35_percent'] for s in summary[1:])
    for field in ('cumulative_0_600_pL','endpoint_pL_s','mean_60_600_pL_s'):
        assert v['wt']['secretion'][field]==read(ROOT/'results/40_ae4_equal_cation_routing/wt_verification.json')['secretion'][field]
    assert all(v[c]['initial_state_sha256']==v['wt']['initial_state_sha256'] for c in CASES)
    b=read(FINAL/'design.json')['null_stimulated_cacc_fraction']
    cell_balance=[]
    for c in CASES:
        grid=g[c];ncl=grid['cl_i_mM']*grid['volume_i_pL']
        source=grid['nkcc1_cl_inward_fmol_s']+grid['ae4_cl_inward_fmol_s']+grid['ae2_cl_inward_signed_fmol_s']-grid['cacc_cl_export_fmol_s']
        quad=float(np.trapezoid(source,times));delta=float(ncl[-1]-ncl[0])
        cell_balance.append({'case':c,'initial_Cl_i_amount_fmol':float(ncl[0]),
            'final_Cl_i_amount_fmol':float(ncl[-1]),'delta_Cl_i_amount_fmol':delta,
            'integrated_net_Cl_i_source_fmol':quad,'one_second_quadrature_difference_fmol':quad-delta,
            'Na_i_amount_initial_fmol':float(grid['na_i_mM'][0]*grid['volume_i_pL'][0]),
            'Na_i_amount_final_fmol':float(grid['na_i_mM'][-1]*grid['volume_i_pL'][-1])})
    result={'status':'TASK41_TARGET_REACHED','selected_candidate':'08','target_used_for_design':True,
        'independently_validated_mechanism':False,'target_metric':'0-600 s cumulative secretion versus same-model WT',
        'target_range_percent':[20,35],'selected_null_cacc_recruitment_fraction':b,
        'mechanism':'Explicit AE4-dependent stimulated CaCC recruitment; Task 40 otherwise',
        'frozen_rest_sha256':v['wt']['initial_state_sha256'],'WT_bitwise_nests_Task40':True,
        'cases':summary,'cell_chloride_balance_quadrature':cell_balance,
        'endpoint':{c:{k:float(val[-1]) for k,val in g[c].items()} for c in CASES},
        'long_run_validity_established':False,
        'final_trajectories':3,'final_numerical_retries':0,'new_final_WT_rest_solves':0}
    write(OUT/'final_summary.json',result)
    csvwrite(OUT/'final_secretion_comparison.csv',summary)
    csvwrite(OUT/'final_cumulative_deficit.csv',curve)
    csvwrite(OUT/'final_flux_comparison.csv',fluxes)
    csvwrite(OUT/'final_intracellular_comparison.csv',endpoints)
    csvwrite(OUT/'final_conservation_summary.csv',conservation)
    csvwrite(OUT/'final_cell_chloride_balance.csv',cell_balance)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
        'axes.spines.right':False,'svg.fonttype':'none','svg.hashsalt':'task41-final'})
    fig,axs=plt.subplots(3,2,figsize=(11,10),layout='constrained')
    colors=('#233b56','#cf7b26','#16877f')
    for c,label,color in zip(CASES,LABELS,colors):
        x=g[c];axs[0,0].plot(times,1e3*x['q_out_pL_s'],label=label,color=color,lw=2)
        if c!='wt':
            use=times>=1;deficit=100*(1-x['cumulative_outflow_0_t_pL'][use]/g['wt']['cumulative_outflow_0_t_pL'][use])
            axs[0,1].plot(times[use],deficit,label=label,color=color,lw=2)
            axs[0,1].annotate(f'{deficit[-1]:.2f}%',(600,deficit[-1]),xytext=(-5,5),textcoords='offset points',ha='right',color=color)
        for ax,key in ((axs[1,0],'na_i_mM'),(axs[1,1],'cl_i_mM'),(axs[2,0],'ph_i'),(axs[2,1],'nkcc1_cl_inward_fmol_s')):
            ax.plot(times,x[key],color=color,lw=2,label=label)
    axs[0,1].axhspan(20,35,color='#16877f',alpha=.1,label='Target range')
    axs[1,1].axhline(30,color='#777777',ls=':',lw=1);axs[1,1].set_ylim(28,63)
    axs[2,0].axhline(7.3,color='#777777',ls=':',lw=1);axs[2,0].set_ylim(6.88,7.32)
    titles=['Fluid outflow','Cumulative secretion deficit','Intracellular sodium',
        'Intracellular chloride','Intracellular pH','NKCC1 chloride loading']
    yl=['Outflow (fL/s)','Deficit versus WT (%)','Na_i (mM)','Cl_i (mM)','pH_i','Cl influx (fmol/s)']
    for ax,title,ylabel in zip(axs.ravel(),titles,yl):
        ax.set(title=title,xlabel='Time (s)',ylabel=ylabel,xlim=(0,600));ax.grid(alpha=.15)
    axs[0,0].legend(frameon=False,ncol=3,fontsize=9)
    axs[0,1].set_ylim(0,5*np.ceil(max(s['peak_cumulative_deficit_percent'] for s in summary)/5))
    fig.suptitle('AE4-dependent CaCC recruitment: a target-selected model',fontsize=15)
    fig.savefig(OUT/'final_trajectories.svg',metadata={'Date':None})
    fig.savefig(OUT/'final_trajectories.png',dpi=160)
    plt.close(fig)
    print(json.dumps({'cases':summary,'endpoint':result['endpoint'],'cell_balance':cell_balance},indent=2))

if __name__=='__main__':main()
