"""Read completed Task 39 text results only; no model import or integration."""
from pathlib import Path
import csv
import json

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/39_palk_nkcc1_full_validation'
CASES=('wt','ae4_5pct','ae4_null')


def read_json(p): return json.loads(p.read_text())
def write_json(p,v): p.write_text(json.dumps(v,indent=2,sort_keys=True,allow_nan=False)+'\n')
def series(case):
    with (OUT/(case+'_timeseries.csv')).open() as f:
        return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
def integrated(case):
    with (OUT/(case+'_integrated_fluxes.csv')).open() as f:
        return {r['quantity']:(float(r['value']),r['unit']) for r in csv.DictReader(f)}
def write_csv(name,rows):
    with (OUT/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n')
        w.writeheader();w.writerows(rows)


def main():
    budget=read_json(OUT/'budget.json')
    assert budget['status']=='PALK_NKCC1_FULL_VALIDATION_COMPLETE'
    assert budget['no_further_scientific_execution']
    vs={c:read_json(OUT/(c+'_verification.json')) for c in CASES}
    ts={c:series(c) for c in CASES};ints={c:integrated(c) for c in CASES}
    assert all(v['status']=='PASS' and v['complete_600_s'] for v in vs.values())
    assert len({v['initial_state_sha256'] for v in vs.values()})==1
    assert all(v['parameter_hashes']==vs['wt']['parameter_hashes'] for v in vs.values())
    old=read_json(ROOT/'results/38_ae4_perturbation_validation/comparison_summary.json')
    flow=[];intracellular=[];flux=[];comparison={};compensation=[]
    for j,wr in enumerate(ts['wt']):
        r={'time_s':wr['time_s'],'wt_q_out_pL_s':wr['q_out_pL_s'],
           'wt_cumulative_pL':wr['cumulative_outflow_0_t_pL']}
        ir={'time_s':wr['time_s']}
        for field in ('na_i_mM','k_i_mM','cl_i_mM','ph_i','hco3_i_mM','tic_i_mM','volume_i_pL','nkcc1_x_i_mM4'):
            ir['wt_'+field]=wr[field]
        for case in CASES[1:]:
            cr=ts[case][j];assert cr['time_s']==wr['time_s']
            r[case+'_q_out_pL_s']=cr['q_out_pL_s']
            r[case+'_flow_ratio']=cr['q_out_pL_s']/wr['q_out_pL_s']
            r[case+'_cumulative_pL']=cr['cumulative_outflow_0_t_pL']
            r[case+'_cumulative_minus_wt_pL']=cr['cumulative_outflow_0_t_pL']-wr['cumulative_outflow_0_t_pL']
            r[case+'_cumulative_ratio']=(cr['cumulative_outflow_0_t_pL']/wr['cumulative_outflow_0_t_pL'] if j else None)
            for field in ('na_i_mM','k_i_mM','cl_i_mM','ph_i','hco3_i_mM','tic_i_mM','volume_i_pL','nkcc1_x_i_mM4'):
                ir[case+'_'+field]=cr[field]
                ir[case+'_minus_wt_'+field]=cr[field]-wr[field]
        flow.append(r);intracellular.append(ir)
    for quantity,(wt,unit) in ints['wt'].items():
        r={'quantity':quantity,'unit':unit,'window_start_s':60,'window_end_s':600,'wt':wt}
        for case in CASES[1:]:
            value=ints[case][quantity][0]
            r[case]=value;r[case+'_minus_wt']=value-wt
            r[case+'_ratio_to_wt']=value/wt if wt!=0 else None
        flux.append(r)
    for case in CASES[1:]:
        s=vs[case]['secretion'];w=vs['wt']['secretion']
        ratio=s['cumulative_0_600_pL']/w['cumulative_0_600_pL']
        mean_ratio=s['mean_60_600_pL_s']/w['mean_60_600_pL_s']
        end_ratio=s['endpoint_pL_s']/w['endpoint_pL_s']
        nkcc=ints[case]['positive_nkcc1_cl_loading'][0]/ints['wt']['positive_nkcc1_cl_loading'][0]
        old_nkcc=old['perturbations'][case]['nkcc1_positive_loading_change_percent']
        percentages=[(r['time_s'],100*(1-r[case+'_cumulative_ratio'])) for r in flow[1:]]
        flow_deficits=[(r['time_s'],100*(1-r[case+'_flow_ratio'])) for r in flow]
        cl_differences=[(r['time_s'],r[case+'_minus_wt_cl_i_mM']) for r in intracellular]
        comparison[case]={
            'cumulative_ratio':ratio,'cumulative_deficit_percent':100*(1-ratio),
            'mean_flow_ratio_60_600':mean_ratio,'mean_flow_deficit_percent':100*(1-mean_ratio),
            'endpoint_flow_ratio':end_ratio,'endpoint_flow_deficit_percent':100*(1-end_ratio),
            'integrated_positive_nkcc1_loading_ratio':nkcc,
            'nkcc1_compensation_percent':100*(nkcc-1),
            'task38_nkcc1_compensation_percent':old_nkcc,
            'compensation_change_percentage_points':100*(nkcc-1)-old_nkcc,
            'Palk_reduces_relative_compensation':nkcc-1<old_nkcc/100,
            'endpoint_differences':{k:ts[case][-1][k]-ts['wt'][-1][k]
                for k in ('na_i_mM','k_i_mM','cl_i_mM','ph_i','hco3_i_mM','tic_i_mM','volume_i_pL','nkcc1_x_i_mM4')},
            'endpoint_X_ratio':ts[case][-1]['nkcc1_x_i_mM4']/ts['wt'][-1]['nkcc1_x_i_mM4'],
            'flux_ratios':{k:ints[case][k][0]/v[0] if v[0] else None for k,v in ints['wt'].items()},
            'largest_tabulated_cumulative_deficit_time_and_percent':max(percentages,key=lambda x:x[1]),
            'largest_tabulated_flow_deficit_time_and_percent':max(flow_deficits,key=lambda x:x[1]),
            'most_negative_tabulated_Cl_difference_time_and_mM':min(cl_differences,key=lambda x:x[1]),
            'time_grid_caveat':'Extrema and crossings in these comparisons use the published 60-second readouts; quadratures retain the one-second grid.',
        }
        compensation.append({'case':case,'task38_compensation_percent':old_nkcc,
            'task39_compensation_percent':100*(nkcc-1),'change_percentage_points':100*(nkcc-1)-old_nkcc,
            'task38_cumulative_deficit_percent':100*(1-old['perturbations'][case]['secretion']['ratio_to_frozen_wt']),
            'task39_cumulative_deficit_percent':100*(1-ratio),
            'compensation_definition':'60-600 s positive NKCC1 chloride loading relative to the WT with the same core law'})
    report={'status':budget['status'],'scientific_hypothesis':'NOT_SUPPORTED_IN_FROZEN_TASK38_ARCHITECTURE',
        'perturbations':comparison,'wt':vs['wt']['secretion'],
        'chloride_partitions':{c:vs[c]['chloride_partition'] for c in CASES},
        'held_out_validation':{'experimental_null_deficit_percent':35,'original_2018_model_flow_deficit_percent':24,
            'predicted_null_cumulative_deficit_percent':comparison['ae4_null']['cumulative_deficit_percent'],
            'predicted_null_endpoint_flow_deficit_percent':comparison['ae4_null']['endpoint_flow_deficit_percent'],
            'used_for_fitting':False,'full_2018_model_reproduction_claimed':False,
            'comparison_scope':'The 24% value is a historical model flow-rate deficit; protocols and whole-gland/whole-cell observables differ. Both present cumulative and endpoint results are reported.'},
        'budget':budget}
    write_csv('flow_and_cumulative_comparison.csv',flow)
    write_csv('intracellular_comparison.csv',intracellular)
    write_csv('integrated_flux_comparison.csv',flux)
    write_csv('nkcc1_compensation_comparison.csv',compensation)
    write_json(OUT/'comparison_summary.json',report)
    print(json.dumps({'status':report['status'],'perturbations':comparison},indent=2))


if __name__=='__main__': main()
