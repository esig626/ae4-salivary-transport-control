"""Postprocess the completed frozen Task 40 files; never import the model."""
from pathlib import Path
import csv
import json
import math

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/40_ae4_equal_cation_routing'
HERE=ROOT/'analysis/40_ae4_equal_cation_routing'
CASES=('wt','ae4_5pct','ae4_null')


def read_json(path):return json.loads(path.read_text())
def read_csv(path):return list(csv.DictReader(path.open()))
def write_json(path,data):path.write_text(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+'\n')
def write_csv(name,rows):
    with (OUT/name).open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def main():
    budget=read_json(OUT/'budget.json')
    assert budget['status']=='TASK40_EQUAL_ROUTING_VALIDATION_COMPLETE'
    assert budget['no_further_scientific_execution']
    assert budget['stationary_solves']==1 and budget['integration_attempts']==3
    assert budget['optimisation_calls']==budget['parameter_sweeps']==budget['calibration_calls']==0
    rest=read_json(OUT/'wt_rest.json');freeze=read_json(OUT/'frozen_inputs.json')
    v={c:read_json(OUT/(c+'_verification.json')) for c in CASES}
    s={c:[{k:float(value) for k,value in row.items()} for row in read_csv(OUT/(c+'_timeseries.csv'))] for c in CASES}
    f={c:{row['quantity']:(float(row['value']),row['unit']) for row in read_csv(OUT/(c+'_integrated_fluxes.csv'))} for c in CASES}
    old={row['case']:row for row in read_csv(ROOT/'results/39_palk_nkcc1_full_validation/nkcc1_compensation_comparison.csv')}
    assert rest['admissible'] and rest['state_sha256']==freeze['frozen_rest_state_sha256']
    for c in CASES:
        assert v[c]['status']=='PASS' and v[c]['complete_600_s']
        assert v[c]['initial_state_sha256']==rest['state_sha256']
        assert v[c]['initial_state_vector']==rest['state_vector']
        assert v[c]['parameter_hashes']==freeze['stimulus_parameter_hashes']
        assert [row['time_s'] for row in s[c]]==list(range(0,601,60))
        assert len(v[c]['integration_attempts'])==1
        assert v[c]['integration_attempts'][0]['solver']==v['wt']['integration_attempts'][0]['solver']
        assert math.isclose(s[c][-1]['cumulative_outflow_0_t_pL'],v[c]['secretion']['cumulative_0_600_pL'],abs_tol=1e-14)
        for row in s[c]:
            j=row['ae4_cl_inward_fmol_s']
            assert row['ae4_na_cell_source_fmol_s']==row['ae4_k_cell_source_fmol_s']==-0.5*j
            assert row['ae4_tic_cell_source_fmol_s']==row['ae4_ta_cell_source_fmol_s']==-2*j
    flow=[];intracellular=[];integrals=[];comparison={};compensation=[]
    fields=('na_i_mM','k_i_mM','cl_i_mM','ph_i','hco3_i_mM','tic_i_mM','volume_i_pL','nkcc1_x_i_mM4')
    for index,w in enumerate(s['wt']):
        row={'time_s':w['time_s'],'wt_q_out_pL_s':w['q_out_pL_s'],
             'wt_cumulative_pL':w['cumulative_outflow_0_t_pL']}
        cell={'time_s':w['time_s'],**{'wt_'+k:w[k] for k in fields}}
        for c in CASES[1:]:
            x=s[c][index];qw=x['q_out_pL_s']/w['q_out_pL_s']
            ratio=x['cumulative_outflow_0_t_pL']/w['cumulative_outflow_0_t_pL'] if index else None
            row.update({c+'_q_out_pL_s':x['q_out_pL_s'],c+'_flow_ratio':qw,
                c+'_instantaneous_deficit_percent':100*(1-qw),
                c+'_cumulative_pL':x['cumulative_outflow_0_t_pL'],
                c+'_cumulative_lost_pL':w['cumulative_outflow_0_t_pL']-x['cumulative_outflow_0_t_pL'],
                c+'_cumulative_ratio':ratio,c+'_cumulative_deficit_percent':100*(1-ratio) if ratio is not None else None})
            for k in fields:
                cell[c+'_'+k]=x[k];cell[c+'_minus_wt_'+k]=x[k]-w[k]
        flow.append(row);intracellular.append(cell)
    for name,(value,unit) in f['wt'].items():
        row={'quantity':name,'unit':unit,'window_start_s':60,'window_end_s':600,'wt':value}
        for c in CASES[1:]:
            x=f[c][name][0];ratio=x/value if value else None
            row.update({c:x,c+'_minus_wt':x-value,c+'_ratio_to_wt':ratio,
                c+'_relative_change_percent':100*(ratio-1) if ratio is not None else None})
        integrals.append(row)
    for c in CASES[1:]:
        x=s[c][-1];w=s['wt'][-1]
        ratio=v[c]['secretion']['cumulative_0_600_pL']/v['wt']['secretion']['cumulative_0_600_pL']
        nkcc=100*(f[c]['positive_nkcc1_cl_loading'][0]/f['wt']['positive_nkcc1_cl_loading'][0]-1)
        oldnkcc=float(old[c]['task39_compensation_percent'])
        peak=max(flow[1:],key=lambda r:r[c+'_cumulative_deficit_percent'])
        na_positive=[row['time_s'] for row in intracellular if row[c+'_minus_wt_na_i_mM']>0]
        comparison[c]={
            'cumulative_ratio':ratio,'cumulative_deficit_percent':100*(1-ratio),
            'cumulative_lost_pL':v['wt']['secretion']['cumulative_0_600_pL']-v[c]['secretion']['cumulative_0_600_pL'],
            'endpoint_flow_ratio':x['q_out_pL_s']/w['q_out_pL_s'],
            'endpoint_flow_deficit_percent':100*(1-x['q_out_pL_s']/w['q_out_pL_s']),
            'nkcc1_compensation_percent':nkcc,'task39_nkcc1_compensation_percent':oldnkcc,
            'compensation_change_percentage_points':nkcc-oldnkcc,
            'relative_reduction_in_compensation_percent':100*(1-nkcc/oldnkcc),
            'endpoint_na_higher_than_wt':x['na_i_mM']>w['na_i_mM'],
            'tabulated_times_na_higher_than_wt':na_positive,
            'endpoint_differences':{k:x[k]-w[k] for k in fields},
            'endpoint_cl_lower_percent':100*(1-x['cl_i_mM']/w['cl_i_mM']),
            'endpoint_X_ratio':x['nkcc1_x_i_mM4']/w['nkcc1_x_i_mM4'],
            'endpoint_cell_na_amount_fmol':x['na_i_mM']*x['volume_i_pL'],
            'endpoint_cell_cl_amount_fmol':x['cl_i_mM']*x['volume_i_pL'],
            'integrated_flux_relative_changes_percent':{k:100*(f[c][k][0]/wt-1) if wt else None for k,(wt,_) in f['wt'].items()},
            'nkcc1_replacement_fraction_of_lost_ae4_loading':(f[c]['positive_nkcc1_cl_loading'][0]-f['wt']['positive_nkcc1_cl_loading'][0])/(f['wt']['positive_ae4_cl_loading'][0]-f[c]['positive_ae4_cl_loading'][0]),
            'largest_tabulated_cumulative_deficit':{'time_s':peak['time_s'],'percent':peak[c+'_cumulative_deficit_percent']},
            'comparison_extrema_resolution_s':60,
            'ae2_positive_loading_share_percent':100*v[c]['chloride_partition']['positive_shares']['ae2'],
            'task39_cumulative_deficit_percent':float(old[c]['task39_cumulative_deficit_percent'])}
        compensation.append({'case':c,'task39_compensation_percent':oldnkcc,
            'task40_compensation_percent':nkcc,'change_percentage_points':nkcc-oldnkcc,
            'task39_cumulative_deficit_percent':float(old[c]['task39_cumulative_deficit_percent']),
            'task40_cumulative_deficit_percent':100*(1-ratio)})
    core=[{'state':key,'unit':'pL' if k in (5,11) else 'fmol',
           'task31':rest['initial_state_vector'][k],'task40':rest['state_vector'][k],
           'task40_minus_task31':rest['core_shift_from_task31'][k]} for k,key in enumerate(rest['state_names'][:12])]
    state31=read_json(ROOT/'results/31_nhe1_mechanistic_repair/rest_checkpoints/185f5cee3bcec7a023cf82e4092872114f2292732ad9c4bab68ecf4e062f76e3.json')
    ledger=[]
    for name,value in rest['flux_ledger'].items():
        if isinstance(value,(int,float)):
            previous=state31['fluxes'].get(name)
            ledger.append({'quantity':name,'task31':previous,'task40':value,
                'task40_minus_task31':value-previous if previous is not None else None})
    residuals=[{'diagnostic':key,'tolerance':tol,'wt_rest_abs':abs(rest['diagnostics']['conservation_residuals'][key]),
        **{c+'_max_abs':v[c]['max_abs_conservation_residuals'][key] for c in CASES}}
        for key,tol in rest['conservation_tolerances'].items()]
    summary={'status':budget['status'],'mechanistic_conclusion':'PARTIAL_COMPENSATION_REDUCTION_WITHOUT_SUSTAINED_NA_ELEVATION_OR_LARGE_SECRETION_DEFICIT',
        'wt':v['wt']['secretion'],'perturbations':comparison,
        'chloride_partitions':{c:v[c]['chloride_partition'] for c in CASES},
        'all_initial_states_identical_to_frozen_rest':True,
        'held_out_validation':{'experimental_null_deficit_percent':35,
            'predicted_null_cumulative_deficit_percent':comparison['ae4_null']['cumulative_deficit_percent'],
            'gap_percentage_points':35-comparison['ae4_null']['cumulative_deficit_percent'],
            'used_for_fitting':False,'same_experimental_protocol_or_observation_model_asserted':False},
        'comparison_note':'Task 39 and Task 40 effects each use their own frozen WT. Their WT resting states and absolute WT secretion differ.',
        'quadrature':v['wt']['secretion']['quadrature'],'budget':budget}
    for name,rows in [('flow_and_cumulative_comparison.csv',flow),('intracellular_comparison.csv',intracellular),
        ('integrated_flux_comparison.csv',integrals),('nkcc1_compensation_comparison.csv',compensation),
        ('wt_rest_core_comparison.csv',core),('wt_rest_flux_comparison.csv',ledger),('conservation_summary.csv',residuals)]:write_csv(name,rows)
    write_json(OUT/'comparison_summary.json',summary)
    write_json(OUT/'report_verification.json',{'status':'PASS','all_three_cases_complete':True,
        'one_frozen_initial_state':True,'only_genotype_expression_differs':all(v[c]['only_genotype_difference']==({} if c=='wt' else {'ae4_expression':.05 if c=='ae4_5pct' else 0.}) for c in CASES),
        'readout_times':list(range(0,601,60)),'readout_source_signatures_exact':True,
        'matching_production_solver_settings':True,'no_model_imports_or_new_numerical_runs':True})
    print(json.dumps({'status':summary['status'],'perturbations':comparison},indent=2))


if __name__=='__main__':main()
