"""Reports from saved records only; no model import or phenotype comparison."""
from pathlib import Path
import argparse
import csv
import json

ROOT = Path(__file__).resolve().parents[2]
HERE = ROOT/'analysis/42_catalan_2025_ae4_mechanism_classes'
OUT = ROOT/'results/42_catalan_2025_ae4_mechanism_classes'
PANEL = ('C0', 'C1', 'C2', 'C3a', 'C3b', 'C4a', 'C4b')
CASES = ('wt', 'ae4_5pct', 'ae4_null')


def read(path): return json.loads(path.read_text())
def number(v): return 'unavailable' if v is None else f'{v:.10g}'
def table(headers, rows):
    return '\n'.join(['| '+' | '.join(headers)+' |', '| '+' | '.join(['---']*len(headers))+' |',
                      *['| '+' | '.join(map(str,row))+' |' for row in rows]])


def write_csv(path, rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def audit():
    p=read(OUT/'preflight.json')
    rows=[[c['class_id'],*map(number,c['source_coefficients_na_k_cl_tic_ta']),
           number(c['affinity_over_rt']),c['forward_direction'],str(c['scalar_direction_agrees_at_reference'])]
          for c in p['classes']]
    text='''# Task 42 source, charge and thermodynamic audit

This execution uses exactly C0, C1, C2, C3a, C3b, C4a and C4b. All seven focused tests pass before any stationary solve or production integration. The frozen reference is the accepted Task 40 WT rest. Every inherited scientific input is checked against the remote preparation tree.

The source layer replaces only the AE4 contribution object. The complete model already consumes separate TIC and TA properties, so carbonate contributes one carbon and two alkalinity equivalents without duplicating or editing the whole cell equations. C0 gives a bitwise identical complete RHS to Task 40 in the tested resting and stimulated states and all three expressions. Tests also check every other class against the literal declared source difference, invariant scalar J4, unchanged non AE4 diagnostics and independent cellular balances.

The table gives cellular source coefficients per positive inward chloride cycle and forward affinity divided by RT. Concentrations approximate activities.

'''+table(['Class','Na','K','Cl','TIC','TA','A/RT','Forward direction','J4 agrees'],rows)+'''

For each row, charge is exactly `nu_Na + nu_K - nu_Cl - nu_TA = 0`, using rational arithmetic. For every transported quantity, `S(-J4) = -S(J4)`. Source coefficients are fixed; there is no adjustable cation fraction.

## Whole cell amount balances

Let N be inward NKCC1 cycles, H inward NHE1 cycles, B inward NBC cycles, E inward AE2 chloride cycles and P the sum of both Na/K pump cycle rates. Let K_a and K_b be outward apical and basolateral potassium fluxes, C_a outward CaCC chloride flux and D the sum of inward CO2 fluxes. All rates below use fmol/s; TA counts equivalents. For each class, insert its five coefficients from the table:

```
dNa_i/dt  = N + H + B - 3P + nu_Na J4
dK_i/dt   = N + 2P - K_a - K_b + nu_K J4
dCl_i/dt  = 2N + E - C_a + nu_Cl J4
dTIC_i/dt = D + 2B - E + nu_TIC J4
dTA_i/dt  = H + 2B - E + nu_TA J4
```

Volume evolves under the unchanged bath to cell and cell to lumen water fluxes. Concentration changes follow `dc/dt = (dn/dt - c dV/dt)/V`. Paracellular transport enters the inherited lumen balances and current closure, not these cellular transport sources. Multiplying the amount balances by charge coefficients gives `d(Na_i+K_i-Cl_i-TA_i)/dt = -B-P-K_a-K_b+C_a`; AE4, NKCC1, NHE1 and AE2 cancel exactly. The unchanged membrane current solve closes the remaining electrogenic terms.

## Thermodynamic calculation

For signed cellular species coefficients nu_s, the forward free energy divided by RT is `sum_s nu_s log(c_i,s/c_o,s) + (V_b/V_T) sum_s nu_s z_s`. Forward affinity is its negative. Each individual chemical and electrical term is saved in preflight.json. The electrical contribution cancels because every class is electroneutral, including carbonate valence minus two. The nonzero membrane potential is included before this cancellation.

Carbonate is taken from the inherited acid base speciation in the cell and bath, with `CO3/HCO3 = 10^(pH-pKa2)`. The inherited pKa2, temperature, reference concentrations and potential are recorded in preflight.json. Carbonate is never replaced by bicarbonate in the audit or cellular sources.

Positive A supports the proposed inward chloride direction, negative A opposes it, and values within 1e-12 of zero are labelled indeterminate. This direction diagnostic does not remove any algebraically valid class from the prescribed panel. Keeping J4(state) fixed is an intentional source comparison assumption. A mismatch between its direction and a class affinity is recorded, not repaired by changing rates. These models do not claim a thermodynamically derived kinetic law for every new class.

## Execution and publication boundaries

The seven source classes, parameter values, solver tolerances, physiology gates and stimulus are fixed before execution. Each class receives one WT root solve from the accepted Task 40 state. Every class passing rest receives WT, 5% and null trajectories from that identical class rest, with the inherited 600 s protocol. No root search, parameter fitting, sweep, added class or genotype rest is permitted. Numerical retry allowances remain those in the supplied task; no retry follows a physiology failure.

The source audit is published first, then all seven rests, then each class immediately after its three trajectories. The complete prediction checkpoint is published last and execution stops before phenotype reveal. Source and execution hashes are frozen in preflight.json and each class manifest. Publication history records confirmed remote commits.

Previous Task 42 source scripts were reused after review; previous numerical outputs were not copied into this execution. Changes before the audit comprise publication stage controls, a verifier that does not assume previous outcomes, compact reporting and accurate class completion status when a trajectory fails. The scientific source layer and all inherited mechanisms remain unchanged.

The required instructions and Task 40 report contain a phenotype reference. No experimental phenotype source or Task 41 content is inspected for this execution. No target value is used by source code, gates, execution or reports. This checkpoint makes no claim about agreement with the held out phenotype.
'''
    (HERE/'source_and_balance_audit.md').write_text(text)
    print('Source and balance audit written.')


def class_metrics(name):
    directory=OUT/name;rest=read(directory/'wt_rest.json')
    affinity=next(c for c in read(OUT/'preflight.json')['classes'] if c['class_id']==name)
    row={'class_id':name,'wt_rest_status':rest['status'],'affinity_over_rt':affinity['affinity_over_rt'],
         'forward_direction':affinity['forward_direction']}
    records={};integrals={};failure=[]
    for case in CASES:
        f=directory/(case+'_verification.json')
        v=read(f) if f.exists() else None; records[case]=v
        passed=v is not None and v['status']=='PASS'
        row[case+'_status']=v['status'] if v else 'NOT_RUN_REST_FAILED'
        row[case+'_cumulative_0_600_pL']=v['secretion']['cumulative_0_600_pL'] if passed else None
        row[case+'_endpoint_pL_s']=v['secretion']['endpoint_pL_s'] if passed else None
        ip=directory/(case+'_all_window_integrals.csv')
        with ip.open() if ip.exists() else open('/dev/null') as f:
            integrals[case]={(r['quantity'],int(r['window_start_s']),int(r['window_end_s'])):float(r['value']) for r in csv.DictReader(f)}
        if v and not passed:
            first=v['first_failure']
            gates=', '.join(g['gate'] for g in first['failed_gates']) if first else ('numerical failure' if v['solver_exception'] else 'WT activation or solver gate')
            at=number(first['time_s'] if first else v['last_checked_time_s'])
            failure.append(f'{case}: {gates} at {at} s')
    for case in CASES[1:]:
        valid=records[case] is not None and records['wt'] is not None and records[case]['status']==records['wt']['status']=='PASS'
        for label,lo,hi in (('early',0,60),('late',60,600),('total',0,600)):
            key=('water_outflow',lo,hi)
            row[case+'_'+label+'_deficit_pct']=100*(1-integrals[case][key]/integrals['wt'][key]) if valid else None
        key=('positive_nkcc1_cl_loading',60,600)
        row[case+'_nkcc1_compensation_pct']=100*(integrals[case][key]/integrals['wt'][key]-1) if valid and integrals['wt'][key]!=0 else None
    row['failure_reason']='; '.join(failure) if rest['admissible'] else '; '.join(g['gate'] for g in rest['failed_gates'])
    return row,rest,records


def summarise_class(name):
    row,rest,records=class_metrics(name); directory=OUT/name
    write_csv(directory/'class_summary.csv',[row])
    cols=['case','status','Q(0,600) pL','q(600) pL/s','early deficit %','late deficit %','total deficit %','NKCC1 change %']
    data=[]
    for case in CASES:
        data.append([case,row[case+'_status'],number(row[case+'_cumulative_0_600_pL']),number(row[case+'_endpoint_pL_s']),
                     *[number(row.get(case+'_'+suffix)) for suffix in ('early_deficit_pct','late_deficit_pct','total_deficit_pct','nkcc1_compensation_pct')]])
    states=['na_i_mM','k_i_mM','cl_i_mM','ph_i','volume_i_pL','q_out_pL_s']
    rest_rows=[[key,number(rest.get('observables',{}).get(key))] for key in states]
    gate_rows=[[k,str(v)] for k,v in rest.get('gates',{}).items()]
    residual_rows=[]
    for case,v in records.items():
        if v:
            for key,value in v['max_abs_conservation_residuals'].items():
                residual_rows.append([case,key,number(value),number(v['conservation_tolerances'][key])])
            for key in ('source_signature_residual_fmol_s','ae4_charge_residual_fmol_s','cell_balance_residual_fmol_s'):
                residual_rows.append([case,key,number(max(abs(v['minima'][key]),abs(v['maxima'][key]))),'0' if key=='source_signature_residual_fmol_s' else '1e-10'])
    text=f'# Task 42 {name}: frozen predictions\n\nWT rest: {rest["status"]}. Forward affinity A/RT: {number(row["affinity_over_rt"])} ({row["forward_direction"]}).\n\n'
    text+=table(cols,data)+'\n\n'
    text+='Early means 0 to 60 s, late means 60 to 600 s. Deficits and NKCC1 compensation use this class\'s WT. A failed WT trajectory leaves these ratios unavailable. All successful genotype outputs remain saved. No extrapolation of an incomplete trajectory is used.\n\n'
    text+='Failure record: '+(row['failure_reason'] or 'None')+'.\n\n'
    text+='## WT rest\n\n'+table(['Observable','Value'],rest_rows)+'\n\n'+table(['Gate','Pass'],gate_rows)+'\n\n'
    text+='## Conserved sources and currents\n\n'+table(['Case','Residual','Maximum absolute value','Tolerance'],residual_rows)+'\n\n'
    text+='Full WT state vectors, root diagnostics and exact failed gates are in wt_rest.json. Each saved timeseries contains intracellular Na, K, Cl, pH, volume, secretion and transport fluxes. Files ending in all_window_integrals.csv contain AE4, NKCC1, AE2, NBC, NHE1, pump, CaCC and paracellular integrals over all three declared windows. Verification files retain endpoints, monitoring details, software, source hashes and failures.\n\nThese are predictions under the inherited Task 40 cycle law. Phenotype interpretation remains locked.\n'
    (directory/'class_report.md').write_text(text)
    print(f'{name}: saved prediction summary written.')


def final():
    complete=read(OUT/'panel_complete.json');verified=read(OUT/'saved_prediction_verification.json')
    assert verified['status']=='SAVED_PREDICTIONS_VERIFIED'
    rows=[class_metrics(n)[0] for n in PANEL]
    write_csv(OUT/'mechanism_comparison.csv',rows)
    def pair(r,suffix): return ' / '.join(number(r[c+'_'+suffix]) for c in CASES[1:])
    data=[[r['class_id'],r['wt_rest_status'],number(r['affinity_over_rt']),
           ' / '.join(number(r[c+'_cumulative_0_600_pL']) for c in CASES),
           pair(r,'early_deficit_pct'),pair(r,'late_deficit_pct'),pair(r,'total_deficit_pct'),
           pair(r,'nkcc1_compensation_pct'),r['failure_reason'] or 'None'] for r in rows]
    text='# Task 42 complete prediction checkpoint\n\nTASK42_PRE_REVEAL_PREDICTIONS_COMPLETE\n\n'
    text+=f'All seven predeclared mechanism classes are recorded. Stationary solves: {complete["stationary_solves"]}. Production integrations: {complete["integration_attempts"]}. Numerical retries: {complete["numerical_retries"]}. Parameter fitting, sweeps and genotype resting solves: zero. One worker and one BLAS thread. All class records were published before this complete checkpoint.\n\n'
    text+=table(['Class','WT rest','A/RT','Q WT / 5% / null (pL)','Early deficit 5% / null (%)','Late deficit 5% / null (%)','Total deficit 5% / null (%)','NKCC1 change 5% / null (%)','Failure'],data)+'\n\n'
    text+='All ratios use an admissible completed WT trajectory from the same class. Unavailable means a required reference or trajectory failed; an incomplete WT is not assigned a 600 s secretion prediction. Each passing rest received the three prescribed genotype attempts. Physiological gate failures are recorded without rescue, and successful outputs from that class are retained. Early and late mean 0 to 60 s and 60 to 600 s. NKCC1 compensation is the change in integrated positive chloride loading over 60 to 600 s.\n\n'
    text+='The source audit records charge identities, the explicit carbonate calculation and all forward affinities. Every class keeps the inherited Task 40 scalar J4(state) even where its direction disagrees with a class affinity. This intentionally isolates transported source coefficients and is not a derivation of class specific kinetics.\n\n'
    text+='The source layer, seven focused tests, inherited source manifest, execution scripts, all WT rest records, frozen genotype inputs, compact trajectories, window integrals, source and current residuals, failure records and publication history accompany this checkpoint. Per class reports are saved beside the numerical records. The saved prediction verifier performs no model integration and confirms that the frozen execution hashes and scientific parameter hashes are unchanged.\n\n'
    text+='No held out phenotype comparison, Task 41 comparison, class selection or phenotype conclusion is made. Execution stops here before reveal.\n'
    (HERE/'prediction_checkpoint.md').write_text(text)
    print('Complete prediction table and checkpoint report written; no phenotype comparison.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=('audit','class','final'))
    parser.add_argument('class_id',nargs='?',choices=PANEL)
    args=parser.parse_args()
    if args.stage=='audit': audit()
    elif args.stage=='class':
        assert args.class_id is not None
        summarise_class(args.class_id)
    else: final()
