"""Freeze the complete empty WT evidence only when its exclusion is proved."""
from __future__ import annotations

from collections import Counter
import json

from .task22_contract import *
from .run_calcium_fast_screen import write_rows
from .validation import sha256_file


def freeze():
    contract=require_stage0();plan=read(OUT/'search_plan.json')
    proof=read(OUT/'structural_certificate.json')
    if not proof['all_ten_root_backgrounds_covered']:
        raise AssertionError('This empty-ensemble writer requires the domain-wide certificate')
    jobs=plan['jobs'];expected={j['attempt_id'] for j in jobs}
    paths=sorted((OUT/'attempts').glob('*.json'))
    if {p.stem for p in paths}!=expected:
        raise AssertionError('WT search incomplete')
    attempts=[read(p) for p in paths]
    if {a['root_id'] for a in attempts if a['kind']=='inherited_root'}!=set(contract['seed_root_ids']):
        raise AssertionError('Missing inherited root')
    if {a['source'] for a in attempts if a['kind']=='task21_resting_witness'}!=set(contract['task21_witness_files']):
        raise AssertionError('Missing Task 21 witness')
    if any(a['candidate'] and a['candidate']['resting_gate_pass'] for a in attempts):
        raise AssertionError('Resting candidate contradicts certificate; resolve before any freeze')
    if any(a['genotype_evaluations']!=0 for a in attempts):
        raise AssertionError('Phenotype firewall breached')
    for name in ('validation_focused.txt','validation_broad_wt.txt','validation_execution.json'):
        if not (OUT/name).is_file():raise AssertionError('Missing validation evidence: '+name)
    execution=read(OUT/'validation_execution.json')
    if any(run['exit_code']!=0 for run in execution['runs']):
        raise AssertionError('Required validation failed')
    decisions=[];calibrations=[];physiology=[];osmotic=[];capacities=[];fluxes=[];anti=[]
    for a in attempts:
        p=a['candidate'];identifier=a['attempt_id']
        base=dict(attempt_id=identifier,root_id=a['root_id'],seed_kind=a['kind'])
        observed=[]
        failures=p['resting_gate_failures'] if p else ['production_exception']
        if any(v.startswith('physiology:') or v.startswith('osmolarity:') for v in failures):
            observed.append('PHYSIOLOGICAL_BOUND_CONFLICT_AT_ENDPOINT')
        if any(v.startswith('capacity:') for v in failures):observed.append('CAPACITY_BOUND_CONFLICT_AT_ENDPOINT')
        if any(v.startswith('thermodynamics:') for v in failures):observed.append('THERMODYNAMIC_CONFLICT_AT_ENDPOINT')
        if any(v.startswith('closure:') for v in failures):observed.append('NUMERICAL_NONCONVERGENCE_OF_RESTING_BALANCES')
        decisions.append(dict(**base,source=a['source'],nfev=a['nfev'],optimizer_success=a['optimizer_success'],
            optimizer_message=a['optimizer_message'],elapsed_s=a['elapsed_s'],restoration_norm=a['restoration_norm'],
            resting_admissible=False,observed_failure_classes=json.dumps(observed),gate_failures=json.dumps(failures),
            domain_wide_failure_class='PROVEN_STRUCTURAL_CONTRADICTION',
            active_coordinate_indices=json.dumps(a['active_coordinate_indices']),genotype_evaluations=0))
        if p is None:continue
        calibrations.append(dict(**base,resting_gate_pass=False,admissible_WT=False,
            measured_WT_objective=p['objective'][0],maximum_abs_log_fold=p['objective'][1],
            sum_squared_log_fold=p['objective'][2],**p['closure'],
            affinity=p['production_diagnostics']['ae4']['diagnostics']['affinity'],
            AE4_Cl_fmol_s=p['production_diagnostics']['ae4']['cl_cell_fmol_s']))
        physiology.append(dict(**base,**p['observables'],resting_gate_pass=False,
            physiology_gate_failures=json.dumps([f for f in failures if f.startswith('physiology:')])) )
        osmotic.append(dict(**base,**{f'{k}_mOsm':v for k,v in p['osmolarities_mOsm'].items()},
            **{f'{k}_bath_ratio':v for k,v in p['osmolarity_ratios'].items()},
            osmotic_gate_pass=not any(f.startswith('osmolarity:') for f in failures)))
        for name in CAPACITIES:
            capacities.append(dict(**base,capacity=name,reference=p['reference_capacities'][name],
                candidate=p['calibrated_capacities'][name],fold=p['capacity_folds'][name],
                lower_fold=.01,upper_fold=100.,domain_pass=.01<=p['capacity_folds'][name]<=100.))
        ledger=p['flux_ledger']
        for name,v in ledger['solute_fluxes_fmol_s'].items():
            fluxes.append(dict(**base,source=name,signed_fmol_s=v,absolute_fmol_s=abs(v),
                positive_loading_pool_fmol_s=ledger['positive_loading_pool_fmol_s'],
                absolute_over_positive_pool=ledger['ratios'].get(name),
                ratio_gate_pass=bool(ledger['ratios'] and ledger['ratios'][name]<=100.)))
        anti.append(dict(**base,positive_loading_pool_fmol_s=ledger['positive_loading_pool_fmol_s'],
            maximum_ratio=ledger['maximum_ratio'],worst_source=ledger['worst_source'],
            anti_cancellation_gate_pass=ledger['passes'],
            **{'Cl_'+k+'_fmol_s':v for k,v in ledger['signed_basolateral_cl_fmol_s'].items()},
            **ledger['water']))
    for filename,rows in [('wt_search_decisions.csv',decisions),('wt_calibration_table.csv',calibrations),
        ('physiology_table.csv',physiology),('osmolarity_diagnostics.csv',osmotic),
        ('capacity_folds.csv',capacities),('resting_flux_diagnostics.csv',fluxes),
        ('anti_cancellation_diagnostics.csv',anti)]:write_rows(OUT/filename,rows)
    write_json(OUT/'wt_solution_payloads.json',dict(created_utc=now(),admissible_solution_ids=[],
        solutions=[],status='EMPTY_PROVED_DOMAIN_CONTRADICTION',
        failed_endpoint_payloads=[str(p.relative_to(REPO))+'#candidate' for p in paths],
        lexicographic_ranking='Empty admissible set; no endpoint selected or ranked as WT.'))
    write_json(OUT/'wt_dynamic_results.json',[])
    write_rows(OUT/'wt_dynamic_results.csv',[dict(root_id=root,calcium_uM=ca,duration_s=600.,
        status='NOT_RUN_NO_ADMISSIBLE_RESTING_WT',dynamic_gate_pass=False,integrations_run=0)
        for root in contract['seed_root_ids'] for ca in CALCIUM])
    summary=dict(created_utc=now(),classification='PROVEN_STRUCTURAL_CONTRADICTION',
        subtype='PHYSIOLOGICAL_THERMODYNAMIC_CONFLICT_FROM_CARBON_ALKALINITY_BALANCES',
        proof_scope=proof['root_certificates'][0]['proof_scope'],search_complete=True,
        inherited_root_starts=10,task21_witness_starts=7,additional_WT_only_starts=10,
        attempts=len(attempts),resting_admissible_count=0,admissible_WT_count=0,
        task22_WT_integrations=0,genotype_evaluations=0,historical_genotype_regressions=0,
        genotype_refits=0,experimental_loss_magnitude_compared=False,
        endpoint_failure_counts=dict(Counter(f for a in attempts for f in a['candidate']['resting_gate_failures'])),
        minimum_scaled_independent_rhs=min(a['candidate']['closure']['max_abs_scaled_independent_rhs'] for a in attempts),
        no_capacity_domain_violations=all(r['domain_pass'] for r in capacities),
        minimum_capacity_fold=min(r['fold'] for r in capacities),maximum_capacity_fold=max(r['fold'] for r in capacities),
        immutable_archive_tree_sha=contract['archive_tree_sha'],
        inherited_source_hashes_unchanged=True,
        disposition='STOP before genotype evaluation. No claim that pooled AE4 is impossible outside this domain/chassis.')
    write_json(OUT/'wt_validation_summary.json',summary)
    artifacts=[p for p in OUT.rglob('*') if p.is_file() and p.name not in ('wt_manifest.json','wt_checkpoint_receipt.json')]
    artifacts+=list(ANALYSIS.glob('*.md'))
    artifacts+=list((REPO/'src/modern_full_model').glob('task22*.py'))
    artifacts+=[REPO/'tests/test_task22_physiological_domain.py']
    write_json(OUT/'wt_manifest.json',dict(created_utc=now(),branch=BRANCH,
        stage0_checkpoint_sha=read(OUT/'stage0_checkpoint.json')['checkpoint_sha'],
        classification=summary['classification'],search_complete=True,
        seed_root_ids=contract['seed_root_ids'],task21_witness_files=contract['task21_witness_files'],
        attempt_ids=sorted(expected),retained_solution_ids=[],rest_valid_solution_ids=[],
        wt_dynamic_conditions_run=0,genotype_evaluations=0,historical_genotype_regressions=0,
        structural_certificate='results/22_physiological_wt_domain/structural_certificate.json',
        file_sha256={str(p.relative_to(REPO)):sha256_file(p) for p in sorted(artifacts)},
        source_hashes=contract['source_hashes'],archive_tree_sha=contract['archive_tree_sha'],
        no_genotype_evaluation_authorized=True))
    print('Complete WT evidence frozen: 27 attempts, 0 admissible WT; push and remotely verify the manifest.')


if __name__=='__main__':freeze()
