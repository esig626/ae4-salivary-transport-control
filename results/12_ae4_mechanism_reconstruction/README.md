# Phase 12 machine-readable evidence

All simulated secretion values are dimensionless within-protocol ratios;
absolute physical flow and time are not certified on this comparison chassis.

| File | Role |
| --- | --- |
| `phenotype_targets.csv` | Primary experimental measurements, uncertainty, protocol, and direct-versus-interpreted status |
| `model_comparison.csv`, `model_comparison.json` | Frozen C1--C8 calibration, WT closure, knockout predictions, phenotype gates, conservation, and thermodynamic status |
| `round2_envelope.csv` | Minimal structural interpolation grid; a generous diagnostic domain, not a source-admissible family |
| `resting_equilibria.csv`, `resting_equilibria.json` | KO-blind moving-WT-root search and resting chloride/pH/volume gates |
| `genotype_resting_equilibria.csv`, `genotype_resting_equilibria.json` | Candidate-independent AE4-null plus candidate/capacity-specific AE2-null resting-state falsification tests |
| `thermodynamic_rest_profile.csv`, `thermodynamic_rest_profile.json` | Candidate affinity, equilibrium, and flux-direction checks at resting roots |
| `rest_capacity_envelope.csv`, `rest_capacity_envelope.json` | Finite `1e-8`--`1e4` capacity continuation and refined connected pass bands |
| `pka_activation_envelope.csv`, `pka_activation_envelope.json` | Independent beta/cAMP/PKA-fold diagnostics across frozen capacity points |
| `trajectories.csv` | Sparse, solver-labelled BDF/Radau landmarks (`t=0,50,100/right-limit,150,200`, plus both sides of other input changes); dense internal grids are used for integration |
| `pathway_localization.csv` | Exact one-component residual and two-direction closure geometry |
| `calibration_manifest.json` | Chassis fingerprint, calibration inputs, held-out targets, conventions, and solver protocol |
| `solver_crosscheck.json` | Independent BDF/Radau and conservation comparisons |
| `summary.json` | Controlling result summary and quantified bounds by evidentiary domain |
| `artifact_sha256.json` | Final source/artifact hashes and scoped A/B reproducibility proof |

The moving-rest capacity continuation and the +/-10% volume band are explicit
post-hoc, KO-blind sensitivity tests. They are neither experimental uncertainty
nor a global proof about disconnected or remote roots. Regenerate every model
output with the command in the Phase-12 analysis README.
