"""Rerun the frozen WT/5% comparison with only AE4 cation routing changed.

Uses the same saved WT and genotype initial cores as Task 14B, without a
resting solve or recalibration. Initial drift is recorded explicitly because
an inherited root need not remain stationary after the routing intervention.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
import csv
import hashlib
import json
from pathlib import Path
import statistics
import time

import numpy as np
import scipy

from .ae4_routing_only import MODEL_ID, with_routing_only
from .genotype_evaluation import (genotype_with_expression, simulate_genotype,
                                 _root_state_with_basal_regulation)
from .model import WT
from .run_calcium_fast_screen import write_json, write_rows, MANIFEST_SHA256
from .task14_blind import build_model, CALCIUM, routing
from .validation import PRODUCTION_RADAU, PRODUCTION_BDF, sha256_file, sha256_object

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / 'results/ae4_routing_only'
FROZEN = REPO / 'results/14B_five_percent_ae4_check/frozen_inputs.json'
MANIFEST = REPO / 'results/13B_modern_full_model/native_dynamic_contract_manifest.json'
SOURCE_COMMIT = '3f41579ef81f1e73cc37dc1ce4ab3380f98c4659'


def inputs():
    if sha256_file(MANIFEST) != MANIFEST_SHA256:
        raise AssertionError('The original ten root manifest changed')
    manifest = json.loads(MANIFEST.read_text())
    frozen = json.loads(FROZEN.read_text())
    cores = {r['root_id']: r['core_state'] for r in frozen['selected_five_percent_states']}
    if set(cores) != set(manifest['roots']) or len(cores) != 10:
        raise AssertionError('The full ten root panel is required')
    for root, payload in manifest['roots'].items():
        for name in ('native_root_object', 'core_state', 'whole_cell_parameters', 'ae4_parameters'):
            if sha256_object(payload[name]) != payload[name + '_sha256']:
                raise AssertionError((root, name))
    for row in frozen['selected_five_percent_states']:
        if sha256_object(row['core_state']) != row['core_state_sha256']:
            raise AssertionError('The inherited 5% initial state changed')
    return manifest, frozen, cores


def output_stem(root, calcium, architecture, expression, solver):
    return f'{root}_ca{calcium:.2f}_{architecture}_e{expression:.2f}_{solver}'


def run_one(request):
    root, calcium, architecture, expression, solver_label = request
    start = time.perf_counter()
    manifest, frozen, cores = inputs()
    original = build_model(manifest, root, calcium)
    model = with_routing_only(original) if architecture == 'routing_only' else original
    core = manifest['roots'][root]['core_state'] if expression == 1. else cores[root]
    genotype = WT if expression == 1. else genotype_with_expression('AE4', expression)
    solver = PRODUCTION_RADAU if solver_label == 'production_radau' else PRODUCTION_BDF
    y0 = _root_state_with_basal_regulation(model, core)
    initial = model.evaluate(0., y0, genotype=genotype)
    trace = simulate_genotype(model, core, transporter='AE4', genotype=genotype,
                              solver=solver, time_s=frozen['time_grid_s'])
    exact_grid = np.array_equal(trace.time_s, frozen['time_grid_s'])
    valid = bool(trace.numerical_gate_pass and exact_grid)
    n = len(trace.time_s)
    ae4 = np.zeros((n, 4))
    reference = np.zeros((n, 3))
    max_net_cl_change = max_anion_change = max_other_rhs_change = 0.
    opposite_count = reference_opposite_count = 0
    negative_cycles = 0
    for i, (t, y) in enumerate(zip(trace.time_s, trace.states.T)):
        new_eval = model.evaluate(float(t), y, genotype=genotype)
        old_eval = original.evaluate(float(t), y, genotype=genotype)
        a, b = new_eval.diagnostics.ae4, old_eval.diagnostics.ae4
        ae4[i] = (a.na_cell_fmol_s, a.k_cell_fmol_s, a.cl_cell_fmol_s, a.hco3_cell_fmol_s)
        reference[i] = (b.na_cell_fmol_s, b.k_cell_fmol_s, b.cl_cell_fmol_s)
        max_net_cl_change = max(max_net_cl_change, abs(a.cl_cell_fmol_s - b.cl_cell_fmol_s))
        max_anion_change = max(max_anion_change, abs(a.hco3_cell_fmol_s - b.hco3_cell_fmol_s))
        max_other_rhs_change = max(max_other_rhs_change, float(np.max(np.abs(new_eval.rhs[2:] - old_eval.rhs[2:]))))
        opposite_count += int(a.na_cell_fmol_s * a.k_cell_fmol_s < 0.)
        reference_opposite_count += int(b.na_cell_fmol_s * b.k_cell_fmol_s < 0.)
        negative_cycles += int(a.cl_cell_fmol_s < 0.)
    if architecture == 'routing_only' and (max_net_cl_change != 0. or max_anion_change != 0.
                                          or max_other_rhs_change != 0. or opposite_count):
        raise AssertionError('Routing intervention changed an excluded component')
    stem = output_stem(*request)
    path = OUT / 'trajectories' / (stem + '.npz')
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {k:v for k,v in asdict(trace).items() if isinstance(v, np.ndarray)}
    arrays.update(ae4_sources_na_k_cl_hco3_fmol_s=ae4,
                  reference_sources_na_k_cl_fmol_s=reference)
    np.savez_compressed(path, **arrays)
    charge = ae4[:,0] + ae4[:,1] - ae4[:,2] - ae4[:,3]
    row = dict(root_id=root, calcium_uM=calcium, architecture=architecture,
               ae4_expression=expression, routing_family=routing(manifest, root),
               solver=solver_label, valid=valid, solver_success=bool(trace.success),
               solver_message=trace.message, exact_time_grid=bool(exact_grid),
               total_0_600_pL=float(trace.cumulative_flow_pL[-1]) if valid else None,
               endpoint_flow_pL_s=float(trace.flow_pL_s[-1]) if valid else None,
               positive_core=bool(trace.positive_core), nonnegative_flow=bool(trace.all_flow_nonnegative),
               conservation_ratio=float(trace.max_dimensionless_conservation_ratio),
               initial_max_amount_rhs_fmol_s=float(np.max(np.abs(initial.rhs[[0,1,2,3,4,6,7,8,9,10]]))),
               initial_max_volume_rhs_pL_s=float(np.max(np.abs(initial.rhs[[5,11]]))),
               max_pointwise_chloride_law_change_fmol_s=max_net_cl_change,
               max_pointwise_anion_law_change_fmol_s=max_anion_change,
               max_pointwise_other_rhs_change=max_other_rhs_change,
               opposing_cation_samples=opposite_count, reference_opposing_cation_samples=reference_opposite_count,
               negative_net_cycle_samples=negative_cycles,
               max_ae4_charge_residual_fmol_s=float(np.max(np.abs(charge))),
               sample_count=n, initial_core_sha256=sha256_object(core),
               whole_cell_parameters_sha256=sha256_object(model.parameters),
               ae4_parameters_sha256=sha256_object(model.ae4_parameters),
               trajectory_path=str(path.relative_to(REPO)), wall_seconds=time.perf_counter()-start)
    for ion in ('na', 'k', 'cl', 'ph'):
        values = getattr(trace, 'cell_' + ion + ('_mM' if ion != 'ph' else ''))
        row['initial_cell_' + ion] = float(values[0])
        row['endpoint_cell_' + ion] = float(values[-1])
    for j, ion in enumerate(('na', 'k', 'cl', 'hco3')):
        row['initial_ae4_' + ion + '_source_fmol_s'] = float(ae4[0,j])
        row['endpoint_ae4_' + ion + '_source_fmol_s'] = float(ae4[-1,j])
    write_json(OUT / 'trajectories' / (stem + '.json'), row)
    return row


def summarise(rows):
    rows = [r for r in rows if r['solver'] == 'production_radau']
    index = {(r['root_id'], r['calcium_uM'], r['architecture'], r['ae4_expression']):r for r in rows}
    manifest, _, _ = inputs()
    old_saved = {(r['root_id'],float(r['calcium_uM'])):r for r in csv.DictReader(
        (REPO / 'results/14B_five_percent_ae4_check/five_percent_ae4_flow.csv').open())}
    pairs = []
    for root in sorted(manifest['roots']):
        for ca in CALCIUM:
            keys = [(root,ca,a,e) for a in ('original','routing_only') for e in (1.,.05)]
            if not all(k in index for k in keys):
                continue
            pair = dict(root_id=root, calcium_uM=ca, routing_family=routing(manifest,root))
            for arch in ('original','routing_only'):
                wt, low = index[root,ca,arch,1.], index[root,ca,arch,.05]
                valid = wt['valid'] and low['valid'] and wt['total_0_600_pL'] > 0.
                ratio = low['total_0_600_pL'] / wt['total_0_600_pL'] if valid else None
                pair.update({arch+'_valid':valid, arch+'_wt_pL':wt['total_0_600_pL'],
                             arch+'_5pct_pL':low['total_0_600_pL'], arch+'_ratio':ratio,
                             arch+'_change_percent':100*(ratio-1.) if valid else None})
            saved = float(old_saved[root,ca]['ae4_5pct_to_wt_total_ratio'])
            pair['saved_original_ratio'] = saved
            pair['original_reproduction_absolute_ratio_error'] = abs(pair['original_ratio']-saved) if pair['original_valid'] else None
            pair['changes_increase_to_decrease'] = (pair['original_valid'] and pair['routing_only_valid']
                and pair['original_ratio'] > 1. and pair['routing_only_ratio'] < 1.)
            pairs.append(pair)
    write_rows(OUT/'comparison.csv', pairs)
    summary = dict(pair_count=len(pairs), valid_pairs=sum(p['original_valid'] and p['routing_only_valid'] for p in pairs),
                   changes_increase_to_decrease=sum(p['changes_increase_to_decrease'] for p in pairs),
                   per_calcium=[])
    for ca in CALCIUM:
        group = [p for p in pairs if p['calcium_uM'] == ca]
        item = dict(calcium_uM=ca, pairs=len(group))
        for arch in ('original','routing_only'):
            values = [p[arch+'_change_percent'] for p in group if p[arch+'_valid']]
            if values:
                item.update({arch+'_change_percent_min':min(values), arch+'_change_percent_max':max(values),
                             arch+'_change_percent_median':statistics.median(values)})
        summary['per_calcium'].append(item)
    summary['max_original_ratio_reproduction_error'] = max((p['original_reproduction_absolute_ratio_error']
        for p in pairs if p['original_valid']), default=None)
    summary['max_conservation_ratio'] = max((r['conservation_ratio'] for r in rows), default=None)
    summary['initial_state_convention'] = 'Saved Task 14B WT and 5% cores unchanged; no new resting solve'
    write_json(OUT/'summary.json', summary)
    return summary


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--pilot', action='store_true')
    parser.add_argument('--confirm', action='store_true')
    args=parser.parse_args()
    manifest, frozen, _=inputs()
    requests=[(r,c,a,e,'production_radau') for r in sorted(manifest['roots']) for c in CALCIUM
              for a in ('original','routing_only') for e in (1.,.05)]
    if args.pilot:
        requests=requests[:4]
    if args.confirm:
        # Crosscheck the extrema of the modified total secretion ratio.
        pairs=list(csv.DictReader((OUT/'comparison.csv').open()))
        chosen={p['root_id']+':'+p['calcium_uM']:p for p in (
            min(pairs,key=lambda p:float(p['routing_only_ratio'])),
            max(pairs,key=lambda p:float(p['routing_only_ratio'])))}
        requests=[(p['root_id'],float(p['calcium_uM']),'routing_only',e,'production_bdf')
                  for p in chosen.values() for e in (1.,.05)]
    files=[MANIFEST,FROZEN,REPO/'results/14_scale_free_genotype_holdout/prediction_contract.json']
    files += sorted((REPO/'src').rglob('*.py'))
    hashes={str(p.relative_to(REPO)):sha256_file(p) for p in files}
    contract=dict(source_commit=SOURCE_COMMIT,model_id=MODEL_ID,root_ids=sorted(manifest['roots']),
                  calcium_uM=CALCIUM,ae4_expression=[1.,.05],initial_states='Inherited WT and 5% cores',
                  no_new_resting_solves=True, no_parameter_fits=True, no_capacity_changes=True,
                  cation_partition='Na/(Na+K) on the cation donor side, equal weighting',
                  net_chloride_law='Exactly inherited evaluate_ae4_qss source_cl_i_fmol_s',
                  production_solver=asdict(PRODUCTION_RADAU), time_grid_s=frozen['time_grid_s'],
                  source_sha256=hashes,numpy_version=np.__version__,scipy_version=scipy.__version__)
    contract_path=OUT/'contract.json'
    if contract_path.exists():
        if json.loads(contract_path.read_text()) != json.loads(json.dumps(contract)):
            raise AssertionError('Frozen run contract changed')
    else:
        write_json(contract_path,contract)
    rows=[]; pending=[]
    for req in requests:
        p=OUT/'trajectories'/(output_stem(*req)+'.json')
        if p.exists():
            rows.append(json.loads(p.read_text()))
        else:
            pending.append(req)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        jobs={pool.submit(run_one,req):req for req in pending}
        for job in as_completed(jobs):
            row=job.result()
            rows.append(row)
            write_rows(OUT/('confirmation_trajectories.csv' if args.confirm else 'trajectory_summary.csv'),rows)
            print(json.dumps({k:row[k] for k in ('root_id','calcium_uM','architecture','ae4_expression',
                                               'valid','total_0_600_pL','wall_seconds')}),flush=True)
    if {str(p.relative_to(REPO)):sha256_file(p) for p in files} != hashes:
        raise AssertionError('Model or frozen inputs changed during integration')
    if not args.confirm:
        print(json.dumps(summarise(rows),indent=2),flush=True)


if __name__=='__main__':
    main()
