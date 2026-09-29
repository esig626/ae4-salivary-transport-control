# Ae4 salivary secretion model

This folder contains the model used for the simulations reported in **Chloride availability limits Slc4a9-dependent salivary secretion**.

## Files to use

**Model implementation:** [`ae4_salivary_model.py`](ae4_salivary_model.py)

**Executable simulations:** [`run_model.py`](run_model.py)

The public model class is `PairedModel`. Its `rhs(t, y)` method evaluates the derivatives of the 26-variable paired control and Ae4 knockout system.

## Quick run

From the repository root:

```bash
python supplementary_material/model/run_model.py
```

This runs the principal combined CCh + IPR simulation.

Other article simulations are available through `--simulation`:

```bash
python supplementary_material/model/run_model.py --simulation auxiliary --output output_auxiliary
python supplementary_material/model/run_model.py --simulation cch-only --output output_cch
python supplementary_material/model/run_model.py --simulation ipr-only --output output_ipr
python supplementary_material/model/run_model.py --simulation compensation-root --output output_sigma_root
python supplementary_material/model/run_model.py --simulation compensation-scan --output output_compensation
```

Use `python supplementary_material/model/run_model.py --help` for the complete command-line interface.

The underlying transport, membrane, acid-base and water equations are provided by the preserved modules under `frozen_task37/src/modern_full_model/`. Full-precision onset states and numerical settings are in `../parameters/`, and stored principal simulation data are in `../data/central/`.
