# Ae4 salivary secretion model

**Use [`ae4_salivary_model.py`](ae4_salivary_model.py). This is the model entry point for the principal simulations reported in the manuscript.**

The `PairedModel` class evaluates the paired control and Ae4 knockout system. Its `rhs(t, y)` method returns the derivatives of the 26 dynamical variables. The underlying transport, membrane, acid-base and water equations are provided by the preserved modules under `frozen_task37/src/modern_full_model/`.

The principal simulation data are stored in `../data/central/`, while full-precision onset states and numerical settings are stored in `../parameters/`.
