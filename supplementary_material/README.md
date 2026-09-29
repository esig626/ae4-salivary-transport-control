# Supplementary material

This directory accompanies **Chloride availability limits Slc4a9-dependent salivary secretion**.

## Start here

The manuscript model is **[`model/ae4_salivary_model.py`](model/ae4_salivary_model.py)** and the executable simulation script is **[`model/run_model.py`](model/run_model.py)**.

From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r supplementary_material/requirements.txt
python supplementary_material/model/run_model.py
```

The default command runs the principal combined CCh + IPR control versus Ae4 knockout comparison.

## Article simulations

```bash
# Principal comparison
python supplementary_material/model/run_model.py --simulation principal --output output_principal

# Auxiliary conductance sensitivity
python supplementary_material/model/run_model.py --simulation auxiliary --output output_auxiliary

# CCh alone
python supplementary_material/model/run_model.py --simulation cch-only --output output_cch

# IPR alone
python supplementary_material/model/run_model.py --simulation ipr-only --output output_ipr

# Compensation at sigma = 1.49
python supplementary_material/model/run_model.py --simulation compensation --sigma 1.49 --output output_sigma_149

# Manuscript compensation root
python supplementary_material/model/run_model.py --simulation compensation-root --output output_sigma_root

# 41-point compensation scan
python supplementary_material/model/run_model.py --simulation compensation-scan --output output_compensation
```

The IPR-only knockout reaches the prescribed intracellular pH acceptance limit before 600 s, so that run records the accepted trajectory up to the stopping point.

## Contents

| Location | Contents |
| --- | --- |
| [`model/ae4_salivary_model.py`](model/ae4_salivary_model.py) | **Manuscript model and public model class** |
| [`model/run_model.py`](model/run_model.py) | **Directly executable article simulations** |
| `model/frozen_task37/` | Transport, membrane, acid-base and water model dependencies |
| `parameters/parameter_and_case_freeze.json` | Full-precision model inputs and numerical settings |
| `parameters/projected_onsets.json` | Control and Ae4 knockout onset states |
| [`PARAMETERS.md`](PARAMETERS.md) | Parameter definitions, units and source locations |
| `data/central/` | Stored principal numerical trajectories and summaries |
| `data/summary_52D.csv` | Compact principal simulation summary |
| `scripts/plot_figures.py` | Plot numerical manuscript figures from stored trajectories |
| `scripts/export_parameters.py` | Export constructed model parameters |
| `scripts/verify_files.py` | Check supplied scientific files against the manifest |
| `provenance/ae4_sections/` | Data and Python plotting material supplied in `ae4_sections.zip` |
| [`SOURCE_MANIFEST.json`](SOURCE_MANIFEST.json) | Source revision and file identifiers |

## Plot figures

```bash
python supplementary_material/scripts/plot_figures.py --output article_figures
```

## Export parameters

```bash
python supplementary_material/scripts/export_parameters.py --output ae4_parameters.json
```

The manuscript tables contain rounded values for readability. The simulations use the full-precision parameter and onset records retained in `parameters/`.
