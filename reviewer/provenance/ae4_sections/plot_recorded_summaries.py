"""Regenerate the two summary plots without evaluating the dynamical model."""
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent
FIGURES = ROOT/'figures'
FIGURES.mkdir(exist_ok=True)


def save(fig, name):
    fig.tight_layout()
    for extension in ('pdf', 'png', 'svg'):
        fig.savefig(FIGURES/f'{name}.{extension}', dpi=240)
    plt.close(fig)


def main():
    central = pd.read_csv(ROOT/'data/fig04_case_02.csv')
    fig, ax = plt.subplots(figsize=(7.0, 3.9))
    ax.plot(central['time_s'], central['cumulative_deficit_percent'],
            marker='o', markersize=4, linewidth=1.5, label='Model')
    ax.errorbar([600], [35.0], yerr=[4.7], fmt='s', markersize=5,
                capsize=4, label='Experiment, mean and standard error')
    ax.set(xlabel='Time after stimulation (s)',
           ylabel='Cumulative secretion deficit (%)', xlim=(0, 630), ylim=(0, 43))
    ax.legend(frameon=False, loc='lower right', fontsize=9)
    save(fig, 'cumulative_central')

    results = pd.read_csv(ROOT/'data/central_results.csv')
    conductance = results['G_aux_nS'].to_numpy()
    control = results['WT_fluid_pL'].to_numpy()
    knockout = results['KO_fluid_pL'].to_numpy()
    ratio = knockout/control
    fig, ax = plt.subplots(figsize=(7.0, 3.9))
    for values, label, marker, linestyle in (
        (control/control[0], 'Control secretion', 'o', '-'),
        (knockout/knockout[0], 'Knockout secretion', 's', '--'),
        (ratio/ratio[0], 'Knockout to control ratio', '^', ':'),
    ):
        ax.plot(conductance, values, label=label, marker=marker,
                linestyle=linestyle, linewidth=1.5, markersize=5)
    ax.set(xlabel='Shared auxiliary chloride conductance (nS)',
           ylabel='Value relative to zero auxiliary conductance',
           xlim=(-0.15, 4.65), ylim=(0.998, 1.014))
    ax.set_xticks(conductance)
    ax.ticklabel_format(axis='y', style='plain', useOffset=False)
    ax.legend(frameon=False, fontsize=9)
    save(fig, 'conductance_normalised')


if __name__ == '__main__':
    main()
