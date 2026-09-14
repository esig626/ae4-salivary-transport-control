"""One closed-form finite-horizon correction using two existing null runs."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from inverse_balance import ROOT,old,read_json,write_json
import time,json

DIRECTORY=Path(__file__).resolve().parent
OUTPUT=ROOT/'results/41_ae4_loss_algebraic_design/candidate_08'


def main():
    start=time.perf_counter();assert not (OUTPUT/'design.json').exists()
    p=ROOT/'results/41_ae4_loss_algebraic_design/candidate_07'
    previous=read_json(p/'design.json');b=previous['null_stimulated_cacc_fraction']
    Qb=read_json(p/'ae4_null_verification.json')['secretion']['cumulative_0_600_pL']
    Q1=read_json(ROOT/'results/40_ae4_equal_cation_routing/ae4_null_verification.json')['secretion']['cumulative_0_600_pL']
    Qwt=read_json(p/'wt_verification.json')['secretion']['cumulative_0_600_pL']
    k=b*(Q1-Qb)/(Qb-b*Q1);Qinf=Q1*(1+k);Qt=.7*Qwt
    correction=k*Qt/(Qinf-Qt)
    assert 0<b<correction<1 and k>0
    result={'solver_success':True,'root_calls':0,'actual_residual_calls':0,
        'method':'Closed-form two-point saturating-conductance surrogate; not an exact ODE identity',
        'inputs':{'previous_b':b,'previous_null_Q_pL':Qb,'task40_null_Q_pL':Q1,'wt_Q_pL':Qwt},
        'surrogate':{'k':k,'Qinf_pL':Qinf,'target_Q_pL':Qt,'target_null_cumulative_deficit':.30},
        'null_stimulated_cacc_fraction':correction,'parameter_selection_uses_trajectories':True,
        'additional_trajectories_during_selection':0,'wall_seconds':time.perf_counter()-start,
        'prior_stationary_null_ph_limitation':'Prior b=0.08140757 has limiting pH=7.3191; no long-run validation claimed at corrected b'}
    write_json(OUTPUT/'design.json',result);print(json.dumps(result,indent=2))

if __name__=='__main__':main()
