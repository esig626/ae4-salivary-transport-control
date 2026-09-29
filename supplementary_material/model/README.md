# Ae4 salivary secretion model

This folder contains the model used for the simulations reported in the manuscript.

## Run it

From the repository root,

```bash
python supplementary_material/model/run_model.py --output ae4_output
```

The default is `case_02`, the principal manuscript simulation. Use

```bash
python supplementary_material/model/run_model.py --case all --output ae4_output
```

to run all three principal auxiliary-conductance cases.

## Model implementation

**[`ae4_salivary_model.py`](ae4_salivary_model.py) is the manuscript model.**

Its public entry point is `PairedModel`. The `rhs(t, y)` method returns the derivatives of the 26-variable paired control and Ae4 knockout system.

The underlying transport, membrane, acid-base and water equations are provided by the preserved modules under `frozen_task37/src/modern_full_model/`.

Full-precision onset states and numerical settings are in `../parameters/`. Stored principal simulation data are in `../data/central/`.
