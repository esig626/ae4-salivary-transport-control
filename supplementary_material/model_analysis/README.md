# Model-analysis figures

This folder contains the Python code and source mapping for the two mathematical figures used in the **Model analysis** section of *Chloride availability limits Slc4a9-dependent salivary secretion*.

The plotting entry point is:

`plot_model_analysis.py`

Run it from the repository root with:

```bash
python supplementary_material/model_analysis/plot_model_analysis.py
```

It regenerates both EPS files under `supplementary_material/model_analysis/figures/` and writes the plotted numerical values to CSV files in this directory.

## Membrane-potential uniqueness figure

**Generated figure:** `figures/Ae4_voltage_uniqueness.eps`

**Code:** `plot_model_analysis.py` and `analysis_math.py`

**Simulation source:** the principal central simulation

`../data/central/case_02/Radau_summary.json`

with the corresponding stored trajectory

`../data/central/case_02/Radau_trajectory.npz`.

This is the same principal comparison used in the manuscript: combined CCh + IPR stimulation, `G_aux = 2.32 nS`, central projected onset states, and matched knockout NKCC1/Ae2 supply (`sigma = 1`). The figure uses the control and knockout ionic compositions at `t = 600 s` from that simulation.

The figure itself is **not a second ODE simulation**. After loading the archived `t = 600 s` compositions, the script holds those compositions fixed and varies the dimensionless basolateral membrane potential `psi_b` in the scalar algebraic membrane-potential equation. It plots `dH/dpsi_b` and marks the algebraic solution corresponding to each archived model state.

To rerun the source simulation from the manuscript model:

```bash
python supplementary_material/model/run_model.py \
  --simulation principal \
  --output output_principal
```

The archived manuscript record remains in `supplementary_material/data/central/case_02/`; the rerun is written separately and does not overwrite it.

## Ae4 parameter-equivalence figure

**Generated figure:** `figures/Ae4_parameter_equivalence.eps`

**Code:** `plot_model_analysis.py` and `analysis_math.py`

This figure does **not** come from a dynamical simulation. It is an analytical plot of the exact carrier scaling symmetry derived in the manuscript. With

```text
Lambda_4 = N_4 k_0 / J_*
```

and the relative sodium and potassium rate ratios fixed, constant `Lambda_4` gives

```text
k_0 / k_0_ref = (Lambda_4 / Lambda_4_ref) / (N_4 / N_4_ref).
```

Every point on a given curve therefore has the same effective Ae4 carrier combination. In particular, `(0.5, 2)`, `(1, 1)`, and `(2, 0.5)` lie on the reference curve and give identical Ae4 cycle rates at every ionic state when the other model quantities are unchanged.

The underlying Ae4 carrier implementation is in:

`../model/frozen_task37/src/modern_full_model/ae4_routing_only.py`

and the publication-facing model entry point is:

`../model/ae4_salivary_model.py`.

## Exact inputs and provenance

`inputs.json` records the numerical values used to generate the published figure files, including the source simulation path and source commit. These values were extracted from the archived principal endpoint and the model parameter record. `numerical_checks.json` verifies that the algebraic membrane-potential roots recover the archived model voltages and that the Ae4 parameter scaling symmetry is satisfied numerically.

No parameter fitting, image tracing, or new trajectory integration is performed by `plot_model_analysis.py`.
