# Executable protocol and observation layer

`protocol_layer.py` reconstructs the inherited model from the immutable Task 47 input without running Task 47. Every public evaluation, RHS and integration call applies only the explicitly requested acute inhibitor targets. Inhibitors do not change the genotype resting state. The cache identity contains parameter vector, genotype, bath, stimulus and inhibitors. NHE aliases map to combined agonists with or without EIPA.

The source B+ transported ions are represented using the qualified spectator accounting frozen in PREDECLARED_EXECUTION.md. All numerical outputs must retain this qualification. Intracellular equations and inherited transporter laws remain unchanged. The legacy bath builder exists only for an exact nesting regression check; it is not the experimental bath.

The catalogue in output/protocol_coverage.json contains every requested genotype and assay family. Blocked protocols raise UnrepresentableProtocol. This is not a claim that all source experiments are executable. Missing carbonic anhydrase kinetics, extracellular NKCC dependence, sodium zero limits, optical calibration and source regression windows remain scientific blockers. No arbitrary concentration floor or timing window is introduced.

Observation functions preserve concentration dilution, inverse SPQ normalisation, explicit recovery windows, integrated gland water and bicarbonate, and uncertainty in both ratio groups. Source regression windows and dye coefficients must be supplied explicitly; otherwise those observations stay unavailable.

Nine focused regression tests pass, including unchanged inherited RHS, targeted drug action through the public RHS, actual gene deletions, charge accounting, source bicarbonate preservation, independent agonist arms, blocker handling, cache separation, and measurement units. Independent protocol review found two API defects (mask bypass and cache collisions); both were corrected and checked before this checkpoint. No trajectories, resting roots or parameter inference were run for 48B.

Run: `python analysis/48_joint_experimental_constraint_reconstruction/test_protocol_layer.py`.
