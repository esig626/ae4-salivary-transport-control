#!/usr/bin/env python3
"""Build manuscript time series from recorded files; never evaluate a model.

Each plot has one time axis. Saved cumulative integrals are used where available.
EPS, PDF and PNG versions, plotting CSVs and an input hash manifest are emitted.
"""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'manuscript/figures'
DATA = ROOT / 'manuscript/data/time_series'
for directory in (OUT, DATA):
    directory.mkdir(parents=True, exist_ok=True)
INPUTS: dict[str, str] = {}
AUDIT: dict = {'method': 'File-only arithmetic on recorded data; no model evaluation', 'central': {}}
FIGURES = []


def register(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f'Recorded source missing: {path}')
    INPUTS[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path


def read_csv(relative: str) -> dict[str, np.ndarray]:
    with register(ROOT / relative).open(newline='') as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    return {key: np.asarray([float(row[key]) for row in rows]) for key in reader.fieldnames or ()}


def store(name: str, columns: dict[str, np.ndarray]) -> None:
    lengths = {len(np.atleast_1d(v)) for v in columns.values()}
    if len(lengths) != 1:
        raise ValueError(f'Inconsistent plotting columns in {name}')
    with (DATA / (name.replace('/', '_') + '.csv')).open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        writer.writerows(zip(*columns.values()))


def axes(ylabel: str):
    fig, ax = plt.subplots(figsize=(7.0, 3.9))
    ax.set_xlabel('Time after stimulation (s)', fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_xlim(0, 600)
    ax.tick_params(labelsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    return fig, ax


def save(fig, ax, name: str, legend: bool = True, **legend_kw):
    if legend:
        ax.legend(frameon=False, fontsize=legend_kw.pop('fontsize', 9), **legend_kw)
    fig.tight_layout(pad=1.1)
    for extension in ('pdf', 'eps', 'png'):
        kw = {'dpi': 190} if extension == 'png' else {}
        fig.savefig(OUT / f'{name}.{extension}', **kw)
    plt.close(fig)
    FIGURES.append(name)


# 1. Historical chloride loading, signed fluxes rather than endpoint bars.
w = read_csv('results/37_wt_nbc_validation/wt_timeseries.csv')
name = 'fig01_chloride_loading'
fig, ax = axes('Signed chloride flux (fmol/s)')
for field, label, style in [
    ('nkcc1_cl_inward_fmol_s', 'NKCC1', '-'),
    ('ae4_cl_inward_fmol_s', 'AE4', '--'),
    ('ae2_cl_inward_signed_fmol_s', 'AE2', ':'),
]:
    ax.plot(w['time_s'], w[field], style, marker='o', markersize=3, linewidth=1.6, label=label)
store(name, {k: w[k] for k in ('time_s', 'nkcc1_cl_inward_fmol_s', 'ae4_cl_inward_fmol_s', 'ae2_cl_inward_signed_fmol_s')})
save(fig, ax, name)

# 2. Later model compensation, not a claimed reproduction of the isolated assay.
fig, ax = axes('NKCC1 change after AE4 loss (%)')
compensation = {}
for folder, label, style in [
    ('39_palk_nkcc1_full_validation', 'Palk/Benjamin, donor routing', '-'),
    ('40_ae4_equal_cation_routing', 'Palk/Benjamin, equal routing', '--'),
    ('41_ae4_loss_algebraic_design/candidate_08', 'Calibrated export coupling', ':'),
]:
    wt = read_csv(f'results/{folder}/wt_timeseries.csv')
    ko = read_csv(f'results/{folder}/ae4_null_timeseries.csv')
    if not np.array_equal(wt['time_s'], ko['time_s']):
        raise ValueError('Historical genotype grids differ')
    change = 100 * (ko['nkcc1_cl_inward_fmol_s'] / wt['nkcc1_cl_inward_fmol_s'] - 1)
    ax.plot(wt['time_s'], change, style, marker='o', markersize=3, linewidth=1.6, label=label)
    compensation[folder] = change
    store(f'fig02_{folder}', {'time_s': wt['time_s'], 'NKCC1_change_percent': change})
save(fig, ax, 'fig02_nkcc_compensation')

# 3. Recorded class trajectories end where the run ended; no extension of failures.
fig, ax = axes('Wild type intracellular pH')
styles = ['-', '--', ':', '-.', '-', '--', ':']
for cls, style in zip(('C0','C1','C2','C3a','C3b','C4a','C4b'), styles):
    values = read_csv(f'results/42_catalan_2025_ae4_mechanism_classes/{cls}/wt_timeseries.csv')
    ax.plot(values['time_s'], values['ph_i'], style, linewidth=1.5, marker='o', markersize=2.5, label=cls)
    store(f'fig03_{cls}', {'time_s': values['time_s'], 'WT_pH': values['ph_i']})
fig.set_size_inches(7.0, 4.3)
save(fig, ax, 'fig03_class_ph', ncol=4, loc='upper center', bbox_to_anchor=(0.5, -0.19))

# Central paired results, including recorded diagnostic masks and saved quadrature.
central = {}
base = ROOT / 'analysis/52_chloride_reservoir_final_test/output/reporter_recovery_52D/cases'
for case in ('case_01','case_02','case_03'):
    npz_path = register(base / case / 'Radau_trajectory.npz')
    summary = json.loads(register(base / case / 'Radau_summary.json').read_text())
    actual_hash = INPUTS[str(npz_path.relative_to(ROOT))]
    if actual_hash != summary['trajectory_sha256']:
        raise ValueError(f'Trajectory hash mismatch: {case}')
    with np.load(npz_path, allow_pickle=False) as archive:
        arrays = {key: archive[key].copy() for key in archive.files}
    names = arrays['columns'].tolist()
    time = arrays['time_s']
    def column(field, genotype, arrays=arrays, names=names):
        index = names.index(field)
        present = arrays['observations_present'][:, genotype, index]
        values = arrays['observations'][:, genotype, index]
        if not np.all(present & np.isfinite(values)):
            raise ValueError(f'Unavailable plotted diagnostic: {field}')
        return values
    if not summary['completed_600s'] or not summary['reservoir_budget']['pass']:
        raise ValueError(f'Central case is not valid: {case}')
    # Integer-second headlines include the separately stored epsilon onset.
    grid = (time == 0) | (time == 1e-6) | (time == np.round(time))
    if not np.all(np.isin(np.arange(601), time[grid])):
        raise ValueError('Incomplete integer-second grid')
    trapz = getattr(np, 'trapezoid', None) or np.trapz
    calc = [float(trapz(column('q_out_pL_s', g)[grid], time[grid])) for g in (0,1)]
    saved = summary['headline_0_600']['fluid_pL']
    if not np.allclose(calc, [saved['WT'],saved['KO']], rtol=0, atol=2e-12):
        raise ValueError('Recorded fluid integral mismatch')
    delta = arrays['paired_states'][:,2] - arrays['paired_states'][:,15]
    if not np.allclose(delta, arrays['delta_chloride_fmol'], rtol=0, atol=1e-12):
        raise ValueError('Chloride amount difference mismatch')
    qtime = arrays['quadrature_time_s']
    if not np.array_equal(qtime, arrays['raw_accepted_time_s']):
        raise ValueError('Quadrature/accepted-state grids differ')
    delta_q = arrays['raw_accepted_paired_states'][:,2] - arrays['raw_accepted_paired_states'][:,15]
    flux_names = arrays['quadrature_flux_names'].tolist()
    cumulative = arrays['quadrature8_cumulative']
    onset = summary['reservoir_budget']['onset_right_limit_integrals']
    onset_difference = np.array([onset['WT'][str(k)]-onset['KO'][str(k)] for k in flux_names])
    diff = cumulative[:,0,:] - cumulative[:,1,:] + onset_difference
    idx = {key:flux_names.index(key) for key in ('N','E','A','J')}
    rhs = delta[0] + 2*diff[:,idx['N']] + diff[:,idx['E']] + diff[:,idx['A']] - delta_q
    error = diff[:,idx['J']] - rhs
    # Add the independently recorded epsilon onset term to match the reported budget.
    if np.max(np.abs(error)) > 1e-5:
        raise ValueError('Reservoir identity fails on saved quadrature grid')
    if not np.isclose(error[-1], summary['reservoir_budget']['identity_error_fmol'], rtol=0, atol=1e-11):
        raise ValueError('Saved endpoint budget does not match the report')
    audit = {
        'fluid_WT_pL': calc[0], 'fluid_KO_pL': calc[1],
        'fluid_deficit_percent': 100*(1-calc[1]/calc[0]),
        'max_saved_grid_budget_error_fmol': float(np.max(np.abs(error))),
        'endpoint_budget_error_fmol': float(error[-1]),
        'saved_endpoint_budget_error_fmol': summary['reservoir_budget']['identity_error_fmol'],
        'NPZ_SHA256_verified': True,
    }
    AUDIT['central'][case] = audit
    central[case] = (arrays, summary, column, qtime, delta_q, diff, rhs, idx)
    output_columns = {'time_s': time}
    for field in ('q_out_pL_s','cl_i_mM','ph_i','volume_i_pL','chloride_driving_force_V','N_fmol_s','E_fmol_s','A_fmol_s','J_fmol_s','J_aux_fmol_s'):
        for g,label in enumerate(('WT','KO')):
            output_columns[f'{label}_{field}'] = column(field,g)
    store(f'central_{case}', output_columns)

# 4. Cumulative deficits at the same ten declared observation times.
fig, ax = axes('Cumulative secretion deficit (%)')
for folder, label, style in [
    ('40_ae4_equal_cation_routing', 'Equal routing, compensating supply', '--'),
    ('41_ae4_loss_algebraic_design/candidate_08', 'Target calibrated coupling', ':'),
]:
    wt=read_csv(f'results/{folder}/wt_timeseries.csv'); ko=read_csv(f'results/{folder}/ae4_null_timeseries.csv')
    keep=wt['time_s']>0
    deficit=100*(1-ko['cumulative_outflow_0_t_pL'][keep]/wt['cumulative_outflow_0_t_pL'][keep])
    ax.plot(wt['time_s'][keep],deficit,style,marker='s',markersize=3,linewidth=1.5,label=label)
    store(f'fig04_{folder}',{'time_s':wt['time_s'][keep],'cumulative_deficit_percent':deficit})
for case, g, style, mark in [('case_01',0.0,'-','o'),('case_02',2.32,'--','^'),('case_03',4.49,'-.','v')]:
    summary=central[case][1]; t=np.arange(60,601,60,dtype=float)
    values=np.array([100*summary['cumulative_observations'][str(int(x))]['fluid_pL']['KO_deficit_fraction'] for x in t])
    ax.plot(t,values,style,marker=mark,markersize=3.5,linewidth=1.3,label=f'Reservoir, auxiliary {g:g} nS')
    store(f'fig04_{case}',{'time_s':t,'cumulative_deficit_percent':values})
ax.set_ylim(0,42)
fig.set_size_inches(7.0, 4.7)
save(fig,ax,'fig04_cumulative_deficits',loc='upper center',bbox_to_anchor=(0.5,-0.19),ncol=2,fontsize=8)

# 5, 6 and 8. Central WT/KO traces: absolute output, substrate and driving force.
a, summary, col, qt, dq, diff, rhs, index = central['case_02']
for number,field,ylabel,scale in [
    ('05','q_out_pL_s','Fluid outflow (fL/s)',1000),
    ('06','cl_i_mM','Intracellular chloride (mM)',1),
    ('08','chloride_driving_force_V','Outward chloride driving force (mV)',1000),
]:
    fig,ax=axes(ylabel)
    for g,label,style in [(0,'Wild type','-'),(1,'AE4 knockout','--')]:
        ax.plot(a['time_s'],scale*col(field,g),style,linewidth=1.8,label=label)
    save(fig,ax,{'05':'fig05_fluid_flow','06':'fig06_chloride','08':'fig08_driving_force'}[number])

# 7. Exact budget on the original Gauss/accepted-state grid.
fig,ax=axes('WT minus KO chloride amount (fmol)')
ax.plot(qt,diff[:,index['J']],linewidth=2.0,label='Extra apical chloride export')
ax.plot(qt,rhs,'--',linewidth=1.6,label='Initial + AE4 supply − remaining')
ax.plot(qt,diff[:,index['A']],':',linewidth=1.7,label='Cumulative AE4 supply difference')
ax.plot(qt,dq,'-.',linewidth=1.6,label='Remaining intracellular difference')
store('fig07_reservoir_budget',{'time_s':qt,'apical_export_difference_fmol':diff[:,index['J']],'mass_balance_reconstruction_fmol':rhs,'AE4_supply_difference_fmol':diff[:,index['A']],'remaining_reservoir_difference_fmol':dq,'identity_error_fmol':diff[:,index['J']]-rhs})
save(fig,ax,'fig07_reservoir_budget',loc='upper left')

manifest={
    'source_baseline':'764e32a648ee3b69265364788df5a59f76cb2843',
    'figures': sorted(FIGURES), 'input_sha256':dict(sorted(INPUTS.items())),
    'all_axes':'Time after stimulation in seconds',
    'historical_sampling':'Stored CSV grids: 60-second comparative samples and the finer saved class grids; joining lines are visual guides, not newly solved trajectories.',
    'central_sampling':'Saved accepted-step and integer-second diagnostics; independently saved Gauss quadrature for the budget.',
    'reported_sensitivities':'Not plotted: only summary values are present in the manuscript decision record.',
}
(DATA/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
(DATA/'arithmetic_verification.json').write_text(json.dumps(AUDIT,indent=2)+'\n')
print(json.dumps(AUDIT,indent=2))
print(f'Wrote {len(FIGURES)} time series figures in EPS, PDF and PNG formats.')
