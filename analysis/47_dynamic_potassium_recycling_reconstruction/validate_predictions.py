"""Postfreeze validation of fixed predictions, with no model execution."""
from common import *

def main():
    receipt=json.loads((HERE/'REMOTE_PREDICTION_RECEIPT.json').read_text())
    assert receipt['remote_verified']
    freeze=json.loads((HERE/'PREDICTION_MANIFEST.json').read_text())
    for path,sha in freeze['sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha,path
    cases={case:json.loads((HERE/'output'/case/'summary.json').read_text()) for case in ['wt','ae4_null','ae4_5pct']}
    wt=cases['wt'];qwt=wt['endpoint']['cumulative_outflow_pL']
    rows=[]
    for case,s in cases.items():
        assert s['parameter_hashes']==freeze['parameter_hashes']
        assert s['initial_state_sha256']==freeze['initial_state_sha256']
        assert s['all_inherited_gates_passed']
        row=dict(case=case,ae4_expression=s['expression'],
            cumulative_outflow_pL=s['endpoint']['cumulative_outflow_pL'],
            cumulative_deficit_percent=100*(1-s['endpoint']['cumulative_outflow_pL']/qwt),
            early_0_60_deficit_percent=100*(1-s['integrals']['0_60']['q_out_pL_s']/wt['integrals']['0_60']['q_out_pL_s']),
            late_60_600_deficit_percent=100*(1-s['integrals']['60_600']['q_out_pL_s']/wt['integrals']['60_600']['q_out_pL_s']),
            endpoint_flow_deficit_percent=100*(1-s['endpoint']['q_out_pL_s']/wt['endpoint']['q_out_pL_s']),
            nkcc_compensation_percent=100*(s['integrals']['60_600']['k_nkcc_influx_fmol_s']/wt['integrals']['60_600']['k_nkcc_influx_fmol_s']-1),
            all_inherited_gates_passed=True)
        rows.append(row)
    null=next(r for r in rows if r['case']=='ae4_null')
    assert abs(null['cumulative_deficit_percent']-freeze['predictions']['null_cumulative_deficit_percent'])<1e-12
    # The experimental magnitude enters only here, after verified publication.
    data=dict(prediction_commit=receipt['prediction_commit'],
        unchanged_frozen_files=len(freeze['sha256']),retuned_after_freeze=False,
        benchmark=dict(central_null_deficit_percent=35.0,
            source='User supplied central benchmark; Peña-Münzenmayer et al. 2015, DOI 10.1074/jbc.M114.612895',
            comparison='Descriptive whole gland benchmark versus primary acinar model; no new fitted observation model',
            five_percent_condition='Computational expression intervention, not a separately established experimental magnitude'),
        null_direction_recovered=null['cumulative_deficit_percent']>0,
        null_magnitude_gap_percentage_points=35.0-null['cumulative_deficit_percent'],
        wt_acceptable_under_inherited_gates=all(wt['wt_activation'].values()),
        nkcc_compensation_reduced_by_new_mechanism=False,
        correction_made=False,results=rows,
        biological_scope='The binding ceiling is absent from the frozen model; biological capacities and missing kinetic dependencies remain unidentifiable from current evidence.')
    dump(HERE/'output/final_validation.json',data)
    csv_write(HERE/'output/final_validation.csv',rows)
    print(json.dumps(dict(prediction_commit=receipt['prediction_commit'],null_deficit_percent=null['cumulative_deficit_percent'],five_percent_deficit_percent=rows[2]['cumulative_deficit_percent'],all_frozen_files_unchanged=True)))

if __name__=='__main__':main()
