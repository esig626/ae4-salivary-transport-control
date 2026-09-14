# Task 42 C2: frozen predictions

WT rest: PASS. Forward affinity A/RT: 6.871880179 (supported).

| case | status | Q(0,600) pL | q(600) pL/s | early deficit % | late deficit % | total deficit % | NKCC1 change % |
| --- | --- | --- | --- | --- | --- | --- | --- |
| wt | TASK42_SOURCE_CLASS_WT_FAILED | unavailable | unavailable | unavailable | unavailable | unavailable | unavailable |
| ae4_5pct | PASS | 0.9354967649 | 0.001519632568 | unavailable | unavailable | unavailable | unavailable |
| ae4_null | PASS | 0.9437588582 | 0.001549339821 | unavailable | unavailable | unavailable | unavailable |

Early means 0 to 60 s, late means 60 to 600 s. Deficits and NKCC1 compensation use this class's WT. A failed WT trajectory leaves these ratios unavailable. All successful genotype outputs remain saved. No extrapolation of an incomplete trajectory is used.

Failure record: wt: ph_i at 190 s.

## WT rest

| Observable | Value |
| --- | --- |
| na_i_mM | 12.2969325 |
| k_i_mM | 115.2110929 |
| cl_i_mM | 59.23722178 |
| ph_i | 6.833672113 |
| volume_i_pL | 1.424732964 |
| q_out_pL_s | 0.001066652925 |

| Gate | Pass |
| --- | --- |
| bulk_charge | True |
| current | True |
| independent_jacobian_full_rank | True |
| independent_rhs | True |
| no_coordinate_boundary_hits | True |
| omitted_charge_rows | True |
| production_physiology_and_conservation | True |
| regulatory_rhs | True |
| rest_protocol_time_invariance | True |
| root_solver_success | True |

## Conserved sources and currents

| Case | Residual | Maximum absolute value | Tolerance |
| --- | --- | --- | --- |
| wt | ae4_charge_fmol_s | 0 | 1e-10 |
| wt | apical_current_A | 2.358757215e-25 | 1e-20 |
| wt | basolateral_current_A | 2.42596564e-23 | 1e-20 |
| wt | buffer_site_accounting_fmol_s | 0 | 1e-10 |
| wt | carbon_accounting_fmol_s | 5.551115123e-17 | 1e-10 |
| wt | cell_bulk_charge_rate_fmol_s | 2.51271226e-13 | 1e-10 |
| wt | cell_speciation_alkalinity_mM | 6.021849686e-12 | 1e-09 |
| wt | homeostasis_charge_fmol_s | 5.551115123e-17 | 1e-10 |
| wt | lumen_bulk_charge_rate_minus_outflow_fmol_s | 2.438153845e-15 | 1e-10 |
| wt | lumen_speciation_alkalinity_mM | 6.110667528e-13 | 1e-09 |
| wt | water_volume_accounting_pL_s | 2.168404345e-19 | 1e-12 |
| wt | source_signature_residual_fmol_s | 0 | 0 |
| wt | ae4_charge_residual_fmol_s | 0 | 1e-10 |
| wt | cell_balance_residual_fmol_s | 8.326672685e-17 | 1e-10 |
| ae4_5pct | ae4_charge_fmol_s | 0 | 1e-10 |
| ae4_5pct | apical_current_A | 2.77880987e-25 | 1e-20 |
| ae4_5pct | basolateral_current_A | 2.493174065e-23 | 1e-20 |
| ae4_5pct | buffer_site_accounting_fmol_s | 0 | 1e-10 |
| ae4_5pct | carbon_accounting_fmol_s | 5.551115123e-17 | 1e-10 |
| ae4_5pct | cell_bulk_charge_rate_fmol_s | 2.606734273e-13 | 1e-10 |
| ae4_5pct | cell_speciation_alkalinity_mM | 6.270539643e-12 | 1e-09 |
| ae4_5pct | homeostasis_charge_fmol_s | 6.678685383e-17 | 1e-10 |
| ae4_5pct | lumen_bulk_charge_rate_minus_outflow_fmol_s | 2.837140245e-15 | 1e-10 |
| ae4_5pct | lumen_speciation_alkalinity_mM | 6.026290578e-13 | 1e-09 |
| ae4_5pct | water_volume_accounting_pL_s | 4.33680869e-19 | 1e-12 |
| ae4_5pct | source_signature_residual_fmol_s | 0 | 0 |
| ae4_5pct | ae4_charge_residual_fmol_s | 0 | 1e-10 |
| ae4_5pct | cell_balance_residual_fmol_s | 6.938893904e-17 | 1e-10 |
| ae4_null | ae4_charge_fmol_s | 0 | 1e-10 |
| ae4_null | apical_current_A | 2.6495629e-25 | 1e-20 |
| ae4_null | basolateral_current_A | 2.478633781e-23 | 1e-20 |
| ae4_null | buffer_site_accounting_fmol_s | 0 | 1e-10 |
| ae4_null | carbon_accounting_fmol_s | 5.551115123e-17 | 1e-10 |
| ae4_null | cell_bulk_charge_rate_fmol_s | 2.561076351e-13 | 1e-10 |
| ae4_null | cell_speciation_alkalinity_mM | 5.950795412e-12 | 1e-09 |
| ae4_null | homeostasis_charge_fmol_s | 6.830473687e-17 | 1e-10 |
| ae4_null | lumen_bulk_charge_rate_minus_outflow_fmol_s | 2.829333989e-15 | 1e-10 |
| ae4_null | lumen_speciation_alkalinity_mM | 6.05293593e-13 | 1e-09 |
| ae4_null | water_volume_accounting_pL_s | 4.33680869e-19 | 1e-12 |
| ae4_null | source_signature_residual_fmol_s | 0 | 0 |
| ae4_null | ae4_charge_residual_fmol_s | 0 | 1e-10 |
| ae4_null | cell_balance_residual_fmol_s | 5.551115123e-17 | 1e-10 |

Full WT state vectors, root diagnostics and exact failed gates are in wt_rest.json. Each saved timeseries contains intracellular Na, K, Cl, pH, volume, secretion and transport fluxes. Files ending in all_window_integrals.csv contain AE4, NKCC1, AE2, NBC, NHE1, pump, CaCC and paracellular integrals over all three declared windows. Verification files retain endpoints, monitoring details, software, source hashes and failures.

These are predictions under the inherited Task 40 cycle law. Phenotype interpretation remains locked.
