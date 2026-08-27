# Research plan

## Current objective

Develop a compact mathematical paper on identifiability and discrimination of transporter mechanisms in a reduced pump leak model of epithelial secretion.

The earlier sensitivity centred plan is superseded by `docs/QUICK_PAPER_PLAN.md`.

## Core question

What can be inferred about AE2 and AE4 transporter activities, or about competing transporter mechanisms, from steady state physiological measurements?

The intended mathematical object is a parameter to observation map obtained from an independently reconstructed steady state reduction of the published 2018 model.

## Priority stages

### Stage A. Independent reconstruction

Reconstruct the published baseline model and derive the smallest useful steady state system without copying the historical unpublished reduction.

### Stage B. Identifiability

Study the map from transporter activity parameters to measurable outputs. Establish local rank conditions, observational equivalence sets, and minimal measurement panels.

### Stage C. Observation geometry

Use computation to illustrate the analytical results, including equivalence curves and separation gained by adding observables. Numerical geometry supports the mathematics but is not the headline contribution.

### Stage D. Optional discrimination

Only if justified by the observation geometry, compare explicit transporter mechanism classes and quantify discrimination under a small explicit noise model.

### Stage E. Viability audit and paper

Assess whether the strongest analytical result is enough for a short mathematical biology paper. Draft only after the result has been checked and frozen.

## Methods that are not the project

- conventional sensitivity analysis
- parameter sweeps as conclusions
- AUTO or MATCONT bifurcation analysis as the main paper
- a large new composite testing theory unless forced by the problem

## Stop conditions

Reconsider rather than force a paper if

- the published model cannot be independently reconstructed
- the identifiability result reduces to a trivial parameter count with no model specific structure
- all biologically plausible measurement panels are trivially full rank and there is no meaningful observational equivalence
- the mechanism classes overlap everywhere in realistic observation space
- the only new result is numerical conditioning or sensitivity

See `prompts/10_identifiability_discrimination_quickpaper.md` for the executable Codex task.
