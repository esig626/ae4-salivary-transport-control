"""Plot principal numerical results from archived trajectories, without simulation."""
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=('case_01', 'case_02', 'case_03'), default='case_02')
    parser.add_argument('--data', type=Path, default=ROOT / 'data/central')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if any(output == p or p in output.parents for p in
           (ROOT / 'data', ROOT / 'model', ROOT / 'parameters', ROOT / 'provenance')):
        parser.error('Write figures outside the archived input directories.')
    directory = args.data / args.case
    summary = json.loads((directory / 'Radau_summary.json').read_text())
    if not summary['completed_600s']:
        raise ValueError('Only a complete archived principal case can be plotted.')
    with np.load(directory / 'Radau_trajectory.npz', allow_pickle=False) as archive:
        a = {key: archive[key] for key in archive.files}
    time = a['time_s']
    states = a['paired_states']
    columns = a['columns'].tolist()
    if states.shape != (len(time), 26) or not np.all(np.diff(time) > 0):
        raise ValueError('Unexpected state dimensions or time ordering.')
    present = a.get('observations_present', np.isfinite(a['observations']))

    def values(name):
        k = columns.index(name)
        v = a['observations'][:, :, k]
        if not present[:, :, k].all() or not np.isfinite(v).all():
            raise ValueError('Missing required diagnostic: ' + name)
        return v

    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'serif', 'font.size': 12,
                         'axes.labelsize': 12, 'xtick.labelsize': 10,
                         'ytick.labelsize': 10, 'legend.fontsize': 9})

    def finish(fig, ax, filename, ylabel):
        ax.set(xlabel='Time after stimulation (s)', ylabel=ylabel, xlim=(0, 600))
        ax.legend(frameon=False, loc='best')
        fig.tight_layout()
        for extension in ('eps', 'pdf', 'png'):
            fig.savefig(output / (filename + '.' + extension), dpi=300)
        plt.close(fig)

    for name, column, scale, ylabel in (
        ('Fig_1', 'q_out_pL_s', 1000, 'Fluid flow (fL/s)'),
        ('Fig_2', 'cl_i_mM', 1, 'Intracellular chloride (mM)'),
        ('Fig_3', 'chloride_driving_force_V', 1000, 'Chloride driving force (mV)'),
    ):
        fig, ax = plt.subplots(figsize=(7, 3.9))
        v = values(column) * scale
        ax.plot(time, v[:, 0], label='Control', linewidth=1.5)
        ax.plot(time, v[:, 1], label='Ae4 knockout', linewidth=1.5, linestyle='--')
        finish(fig, ax, name, ylabel)
        np.savetxt(output / (name + '_data.csv'), np.column_stack((time, v)),
                   delimiter=',', header='time_s,control,knockout', comments='')

    # Use saved eight-node integrals, not reconstructed or interpolated curves.
    qt = a['quadrature_time_s']
    integral = a['quadrature8_cumulative']
    names = a['quadrature_flux_names'].tolist()
    idx = np.searchsorted(time, qt)
    if np.any(idx >= len(time)) or not np.allclose(time[idx], qt, rtol=0, atol=1e-12):
        raise ValueError('Quadrature times must coincide with recorded states.')
    onset = summary['reservoir_budget']['onset_right_limit_integrals']
    initial = states[0, 2] - states[0, 15]
    remaining = states[idx, 2] - states[idx, 15]
    uptake = integral[:, 0, names.index('A')] + onset['WT']['A']
    export = (integral[:, 0, names.index('J')] - integral[:, 1, names.index('J')]
              + onset['WT']['J'] - onset['KO']['J'])
    rhs = initial + uptake - remaining
    plot_time = np.r_[0., qt]
    series = ((np.r_[0., export], 'Additional control export'),
              (np.r_[0., uptake], 'Continued Ae4 uptake'),
              (np.full(len(plot_time), initial), 'Initial intracellular difference'),
              (np.r_[initial, remaining], 'Remaining intracellular difference'),
              (np.r_[0., rhs], 'Initial difference + uptake - remaining'))
    fig, ax = plt.subplots(figsize=(7, 3.9))
    for k, (v, label) in enumerate(series):
        ax.plot(plot_time, v, label=label, linewidth=1.5,
                linestyle='--' if k == len(series) - 1 else '-')
    finish(fig, ax, 'Fig_4', 'Chloride amount (fmol)')
    np.savetxt(output / 'Fig_4_data.csv', np.column_stack((plot_time, *[v for v, _ in series])),
               delimiter=',', header='time_s,additional_export,Ae4_uptake,initial_difference,remaining_difference,equation_rhs', comments='')
    print('Plotted Fig_1 to Fig_4 from stored data. No ODE was integrated.')
    print('Fig_5 is not generated because the full continuation dataset is not supplied.')


if __name__ == '__main__':
    main()
