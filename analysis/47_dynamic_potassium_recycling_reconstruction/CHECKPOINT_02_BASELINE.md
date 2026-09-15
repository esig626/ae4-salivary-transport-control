# Checkpoint 02: frozen baseline dynamics

Previous published milestone: `0df7a90a27a263757dbffc3917d00482749617b5`.

WT and AE4 null integrations are complete. The detailed current, concentration, amount, storage, voltage, pump and NKCC diagnostics are saved with lossless CSV compression. Baseline interpretation and the complete equations are separate from the numerical outputs. All inherited gates and independent amount/current/storage identities passed. Production source and scientific parameters are unchanged.

The WT summary serialization incident and exact recovery limits are preserved in `output/EXECUTION_NOTES.md`. No trajectory was rerun. There have been two integrations, no resting solves, no optimisation and no parameter sweep. The local derivative checks use saved states only.

Finite recycling is insufficient as a binding mechanism in this frozen model. No independently specified numerical repair has been identified, so there is no reconstruction candidate to tune or select. Biological extrapolation is limited by unmeasured effective conductance and pump kinetics. The next step is to publish the immutable prediction manifest, then run only the deferred 5% case and the phenotype comparison.
