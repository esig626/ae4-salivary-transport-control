# Ae4 manuscript review package

This folder accompanies **Intracellular chloride availability limits Slc4a9-dependent salivary secretion**. It separates the principal model and its complete archived simulations from older summaries and figure assets.

## Start with the model

The principal model entry is [`model/task52_model.py`](model/task52_model.py), specifically `PairedModel`. Its cell equations, chemical equilibria, transport laws and membrane potential calculations are in [`model/frozen_task37/src/modern_full_model/`](model/frozen_task37/src/modern_full_model/). These are unchanged copies of the scientific files at repository revision `c44160e1849828cde887a5e24618335c88f97af7`.

Use this entry rather than selecting a different model from the repository's older analyses. The control uses locally evaluated NKCC1 and Ae2 rates. The knockout receives the same complete cycle rates in the principal comparison, has no Ae4 transport, and starts from its own projected onset state. Both cases have the same apical conductances. No Ae4 dependent apical conductance multiplier is used.

## Contents

| Location | Contents |
| --- | --- |
| `model/` | Principal paired model, unchanged dependencies, original numerical runner and its recorded output repair |
| `parameters/parameter_and_case_freeze.json` | Original full precision onset states, prescribed protocols, case definitions, numerical settings and input hashes |
| `parameters/projected_onsets.json` | Original onset construction records |
| `PARAMETERS.md` | Parameter definitions, units and source locations |
| `data/central/case_01` to `case_03` | Complete archived trajectories and summaries for the three shared auxiliary conductances |
| `data/summary_52D.csv` | Principal secretion and chloride summaries |
| `scripts/` | Portable principal simulation launcher, parameter exporter, data integrity check and plots from archived trajectories |
| `provenance/ae4_sections/` | Selected data and the original Python plotting script extracted unchanged from the uploaded `ae4_sections.zip` |
| `SOURCE_MANIFEST.json` | Source revision, original archive digest and Git blob hashes for copied scientific files |
| `REPRODUCIBILITY_STATUS.md` | Exact coverage and unresolved deposition requirements |

The supplied ZIP was an earlier manuscript package, not a complete model release. Its PDF figures use an older numbering scheme. Those PDFs, the old manuscript ledger and its bibliography are not duplicated here. In particular, the ZIP's old Fig_5 is not the current compensation figure, and its source class figure is not part of the present principal analysis. The original ZIP remains a separate supplied archive, identified by its SHA256 in the manifest.

## Read the archived results

`case_01`, `case_02` and `case_03` use shared auxiliary conductances of 0, 2.32 and 4.49 nS, respectively. All three use combined CCh and IPR stimulation, the central onset construction and knockout supply multiplier 1. Their cumulative fluid deficits are approximately 20.04%, 20.08% and 20.10%. `case_02` is the principal central case.

Each `Radau_trajectory.npz` contains the stored time grid, 26 state coordinates, named diagnostics, diagnostic presence masks and cumulative flux quadrature. Missing optional diagnostics are not zero and must be read with their presence mask. The neighbouring `Radau_summary.json` documents the method, accepted endpoint, physiological checks, integral definitions and numerical residuals.

## Reproduce and inspect

From the repository root, using Python 3.12, create an isolated environment and install the listed dependencies.

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -r reviewer/requirements.txt
python reviewer/scripts/verify_files.py
python reviewer/scripts/plot_figures.py --output /tmp/ae4_figures
python reviewer/scripts/export_parameters.py --output /tmp/ae4_parameters.json
python reviewer/scripts/reproduce_principal.py --case case_02 --output /tmp/ae4_reproduction
```

On Windows, activate the environment with `.venv\Scripts\activate` and choose suitable output paths. Select `--case all` to run all three principal cases. Existing output directories are rejected by the simulation launcher. The parameter exporter constructs the original model but does not integrate a trajectory or fit a parameter.

The portable launcher reuses the archived `run_attempt` routine and its recorded diagnostic serialization repair. It replaces only the original launcher checks tied to a historical branch, commit and one-time output directory. Equations, parameter values, initial states, solver settings, acceptance checks and quadrature are not retuned. Fresh outputs are written separately from archived data. The new adapters were syntax checked during packaging; a fresh full numerical reproduction was not performed during this publication step.

`plot_figures.py` generates the numerical panels corresponding to current Fig_1 through Fig_4 and exports their plotted values. It uses stored trajectories and stored flux quadrature, not image tracing or generated curves. It does not promise identical typography to every earlier figure export. It deliberately does not generate the current Fig_5 because the full original compensation dataset and implementation have not been located in the supplied archive or inspected repository snapshot.

## Scope of the release

The complete principal trajectories are available. The newer 41 point compensation continuation, its tangent and root verification files, and complete CCh only and IPR only output records are not supplied in this package. Older reported summaries are preserved and labelled as such. Do not treat them as independently reproduced trajectories. See [`REPRODUCIBILITY_STATUS.md`](REPRODUCIBILITY_STATUS.md) before making a claim that every result in the current manuscript is reproducible from this release.
