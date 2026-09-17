# Active time series figure plan

All eight figures have time in seconds on the horizontal axis. EPS, PDF and PNG are emitted together. No model solves, new fits, synthetic interpolation models or invented source traces are allowed.

1. `fig01_chloride_loading`: signed historical WT NKCC1, AE4 and AE2 chloride fluxes. Saved 60 s samples; the original integrated loading fractions remain text.
2. `fig02_nkcc_compensation`: pointwise knockout/WT NKCC1 change in historical Palk donor, equal-routing and accepted calibrated candidate_08 models. Distinguish pointwise from integrated compensation.
3. `fig03_class_ph`: saved WT pH traces for all seven molecularly motivated source classes, stopped at their actual failure/endpoints.
4. `fig04_cumulative_deficits`: saved cumulative fluid deficits at ten observation times for equal routing, calibrated coupling and three central reservoir cases. No experimental final-time band across all times.
5. `fig05_fluid_flow`: central 2.32 nS WT/KO outflow, displayed in fL/s.
6. `fig06_chloride`: central WT/KO intracellular chloride concentrations. Onset values are imposed, later values predicted.
7. `fig07_reservoir_budget`: independent extra apical export, mass-balance reconstruction, cumulative AE4 input and remaining intracellular difference. Use original Gauss integrals plus the recorded epsilon contribution.
8. `fig08_driving_force`: central WT/KO outward chloride driving force, E_Cl minus V_a, in mV.

Scalar original-model outcomes and summary-only continuation sensitivities remain tables. Every underlying source hash and all plotted values are recorded in `data/time_series/`.
