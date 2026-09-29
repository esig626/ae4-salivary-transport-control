# Supplementary material

This directory accompanies **Intracellular chloride availability limits Slc4a9-dependent salivary secretion**.

## Manuscript model

**The model used for the principal control and Ae4 knockout simulations is
[`model/ae4_salivary_model.py`](model/ae4_salivary_model.py).**

The public entry point is `PairedModel`. It evaluates the 26-variable paired control and knockout system. The transport, membrane, acid-base and water submodels used by this entry point are retained under [`model/frozen_task37/src/modern_full_model/`](model/frozen_task37/src/modern_full_model/).

Use this model entry point for the manuscript simulations rather than selecting one of the historical model variants elsewhere in the repository.

## Contents

| Location | Contents |
| --- | --- |
| [`model/ae4_salivary_model.py`](model/ae4_salivary_model.py) | **Manuscript model and principal public entry point** |
| `model/frozen_task37/` | Transport, membrane, acid-base and water model dependencies |
| `parameters/parameter_and_case_freeze.json` | Full-precision model inputs, cases and numerical settings |
| `parameters/projected_onsets.json` | Control and Ae4 knockout onset states |
| [`PARAMETERS.md`](PARAMETERS.md) | Parameter document with definitions, units and source locations |
| `data/central/` | Stored numerical trajectories and summaries for the principal simulations |
| `data/summary_52D.csv` | Compact principal simulation summary |
| `scripts/reproduce_principal.py` | Python entry point for rerunning the principal simulations |
| `scripts/plot_figures.py` | Python plotting script for the manuscript numerical figures |
| `scripts/export_parameters.py` | Export the constructed model parameters |
| `scripts/verify_files.py` | Check the supplied scientific files against the manifest |
| `provenance/ae4_sections/` | Data and Python plotting material supplied in `ae4_sections.zip` |
| [`SOURCE_MANIFEST.json`](SOURCE_MANIFEST.json) | Source revision and file identifiers |

## Running the principal simulations

From the repository root, using Python 3.12,

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -r supplementary_material/requirements.txt
python supplementary_material/scripts/verify_files.py
python supplementary_material/scripts/reproduce_principal.py --case case_02 --output /tmp/ae4_reproduction
```

Use `--case all` to run the three principal combined-stimulation cases.

## Plotting the numerical figures

```bash
python supplementary_material/scripts/plot_figures.py --output /tmp/ae4_figures
```

The plotting script reads the stored numerical trajectories and exports the plotted values alongside EPS, PDF and PNG versions.

## Parameters and data

The appendix tables contain rounded values for presentation. The numerical implementation uses the full-precision values retained in the model and parameter records. The control and knockout onset states, prescribed inputs, solver settings and case definitions are in `parameters/`.

The material supplied in `ae4_sections.zip` is retained under `provenance/ae4_sections/`, including its data files and Python summary plotting script.
