# Chloride availability limits Slc4a9-dependent salivary secretion

This repository contains the model, numerical simulations, parameter records and plotting material accompanying the manuscript **Chloride availability limits Slc4a9-dependent salivary secretion**.

The current manuscript model is a conservation-explicit dynamical model of primary salivary secretion that compares a control acinar cell with an Ae4 knockout under matched non-Ae4 chloride supply.

## Quick start

Clone the repository and create a Python environment.

```bash
git clone https://github.com/esig626/ae4-salivary-transport-control.git
cd ae4-salivary-transport-control

python -m venv .venv
source .venv/bin/activate
python -m pip install -r supplementary_material/requirements.txt
```

On Windows, activate the environment with:

```text
.venv\Scripts\activate
```

Run the principal simulation:

```bash
python supplementary_material/model/run_model.py
```

The default run is the combined CCh + IPR control versus Ae4 knockout comparison used for the main manuscript result. Output is written to `ae4_output/`.

## The manuscript model

The model used in the article is:

**[`supplementary_material/model/ae4_salivary_model.py`](supplementary_material/model/ae4_salivary_model.py)**

The public model class is `PairedModel`. It contains the paired 26-variable control and knockout system and exposes the model right-hand side through `rhs(t, y)`.

The directly executable entry point is:

**[`supplementary_material/model/run_model.py`](supplementary_material/model/run_model.py)**

The underlying transport, membrane, acid-base and water submodels are retained in:

`supplementary_material/model/frozen_task37/src/modern_full_model/`

These include the reversible NKCC1 and Ae2 transport laws, NHE1, NBCe1-B, the Ae4 carrier, membrane currents, carbonate chemistry and osmotic water transport.

## Run the simulations reported in the article

All commands below are run from the repository root. Use a different `--output` directory for each run because existing output directories are deliberately not overwritten.

### Principal control and Ae4 knockout comparison

Combined CCh + IPR stimulation, central auxiliary conductance and matched NKCC1/Ae2 supply:

```bash
python supplementary_material/model/run_model.py \
  --simulation principal \
  --output output_principal
```

This is the simulation underlying the principal secretion, intracellular chloride, chloride driving-force and chloride-accounting results.

### Auxiliary apical chloride conductance sensitivity

Run the three shared auxiliary conductances used in the manuscript, 0, 2.32 and 4.49 nS:

```bash
python supplementary_material/model/run_model.py \
  --simulation auxiliary \
  --output output_auxiliary
```

### CCh-only simulation

```bash
python supplementary_material/model/run_model.py \
  --simulation cch-only \
  --output output_cch
```

This reproduces the stimulus-specific comparison in which most of the secretion deficit persists without the beta-adrenergic input.

### IPR-only simulation

```bash
python supplementary_material/model/run_model.py \
  --simulation ipr-only \
  --output output_ipr
```

The knockout trajectory reaches the prescribed upper intracellular pH acceptance limit before 600 s, as reported in the manuscript. The script therefore records the accepted trajectory up to the stopping point rather than treating the missing 600 s endpoint as a valid result.

### Compensation calculation

Evaluate any knockout NKCC1/Ae2 supply multiplier (sigma) between 1 and 2.01:

```bash
python supplementary_material/model/run_model.py \
  --simulation compensation \
  --sigma 1.49 \
  --output output_sigma_149
```

Evaluate the manuscript compensation root directly:

```bash
python supplementary_material/model/run_model.py \
  --simulation compensation-root \
  --output output_sigma_root
```

Run the 41-point nonlinear compensation scan over (1.00 \leq \sigma \leq 2.00) in steps of 0.025:

```bash
python supplementary_material/model/run_model.py \
  --simulation compensation-scan \
  --output output_compensation
```

The scan also writes `compensation_scan.csv` with the control secretion, knockout secretion and cumulative secretion deficit at each value of `sigma`.

## Plot the numerical figures

The stored principal trajectories can be plotted without rerunning the ODE system:

```bash
python supplementary_material/scripts/plot_figures.py \
  --output article_figures
```

The script exports EPS, PDF and PNG figures together with the numerical values used for plotting.

## Model-analysis figures

The two mathematical figures used in the **Model analysis** section have their own reproducible source folder:

[`supplementary_material/model_analysis/`](supplementary_material/model_analysis/README.md)

Run both figure scripts with:

```bash
python supplementary_material/model_analysis/plot_model_analysis.py
```

The membrane-potential uniqueness figure uses the archived **principal central simulation** in `supplementary_material/data/central/case_02/` at `t = 600 s`. That source simulation can be rerun with:

```bash
python supplementary_material/model/run_model.py --simulation principal --output output_principal
```

The Ae4 parameter-equivalence figure is analytical and does not use a dynamical trajectory. Its README links directly to the Ae4 carrier implementation and records the exact inputs used for the published plot.

## Inspect the parameters

A readable parameter document is provided at:

[`supplementary_material/PARAMETERS.md`](supplementary_material/PARAMETERS.md)

The full-precision machine-readable records are:

- `supplementary_material/parameters/parameter_and_case_freeze.json`
- `supplementary_material/parameters/projected_onsets.json`

To export the constructed model parameters:

```bash
python supplementary_material/scripts/export_parameters.py \
  --output ae4_parameters.json
```

The tables in the manuscript are rounded for presentation. The model uses the full-precision values retained in these records.

## Simulation output

Each completed simulation writes a directory containing:

- `Radau_summary.json` with the numerical summary, endpoint values and checks
- `Radau_trajectory.npz` with the paired state trajectory and named diagnostics
- `run_summary.json` at the top level with a compact summary of the requested simulations

The principal archived trajectories used in the manuscript are also retained under:

`supplementary_material/data/central/`

## Use the model from Python

The model can also be imported directly:

```python
import sys
sys.path.insert(0, "supplementary_material/model")

from ae4_salivary_model import PairedModel

model = PairedModel(
    protocol="CCH_IPR",
    G_aux_S=2.32e-9,
    ko_supply_multiplier=1.0,
)
```

`model.rhs(t, y)` evaluates the derivatives of the paired control and Ae4 knockout state. `model.evaluate(t, y)` additionally returns the transport and diagnostic quantities used in the analysis.

## Supplementary material

The full supplementary directory is here:

[`supplementary_material/`](supplementary_material/README.md)

It contains the model, parameter records, numerical data, plotting utilities and the material supplied in `ae4_sections.zip`.

## Repository structure

| Directory | Purpose |
| --- | --- |
| `supplementary_material/` | Model, data, parameters and scripts accompanying the manuscript |
| `manuscript/` | Manuscript source and figure files |
| `analysis/` | Numbered scientific analyses used during model development |
| `archive/` | Earlier model and manuscript material retained for provenance |
| `docs/` | Project ledgers and supporting documentation |

The historical analysis directories are retained for provenance. For the current manuscript, start with `supplementary_material/model/ae4_salivary_model.py` and `supplementary_material/model/run_model.py`.
