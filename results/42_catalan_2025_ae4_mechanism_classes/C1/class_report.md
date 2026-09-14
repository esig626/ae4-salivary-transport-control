# Task 42 C1: frozen predictions

WT rest: PASS. Forward affinity A/RT: 4.845371436 (supported).

| case | status | Q(0,600) pL | q(600) pL/s | early deficit % | late deficit % | total deficit % | NKCC1 change % |
| --- | --- | --- | --- | --- | --- | --- | --- |
| wt | PASS | 0.8562568246 | 0.001238723645 | unavailable | unavailable | unavailable | unavailable |
| ae4_5pct | PASS | 0.9429438823 | 0.001534055855 | 0.2456131549 | -11.58997669 | -10.12395524 | 87.82711688 |
| ae4_null | PASS | 0.948424455 | 0.001553579362 | 0.2613580431 | -12.32275463 | -10.76401702 | 95.48207056 |

Early means 0 to 60 s, late means 60 to 600 s. Deficits and NKCC1 compensation use this class's WT. A failed WT trajectory leaves these ratios unavailable. All successful genotype outputs remain saved. No extrapolation of an incomplete trajectory is used.

Failure record: None.

## WT rest

| Observable | Value |
| --- | --- |
| na_i_mM | 12.14141637 |
| k_i_mM | 116.1167131 |
| cl_i_mM | 59.46498148 |
| ph_i | 6.915079955 |
| volume_i_pL | 1.437672416 |
| q_out_pL_s | 0.001070463877 |

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
| wt | apical_current_A | 2.6495629e-25 | 1e-20 |
| wt | basolateral_current_A | 2.418857057e-23 | 1e-20 |
| wt | buffer_site_accounting_fmol_s | 0 | 1e-10 |
| wt | carbon_accounting_fmol_s | 2.775557562e-17 | 1e-10 |
| wt | cell_bulk_charge_rate_fmol_s | 2.524369602e-13 | 1e-10 |
| wt | cell_speciation_alkalinity_mM | 6.227907079e-12 | 1e-09 |
| wt | homeostasis_charge_fmol_s | 5.551115123e-17 | 1e-10 |
| wt | lumen_bulk_charge_rate_minus_outflow_fmol_s | 2.745199901e-15 | 1e-10 |
| wt | lumen_speciation_alkalinity_mM | 6.488143356e-13 | 1e-09 |
| wt | water_volume_accounting_pL_s | 2.168404345e-19 | 1e-12 |
| wt | source_signature_residual_fmol_s | 0 | 0 |
| wt | ae4_charge_residual_fmol_s | 0 | 1e-10 |
| wt | cell_balance_residual_fmol_s | 5.551115123e-17 | 1e-10 |
| ae4_5pct | ae4_charge_fmol_s | 0 | 1e-10 |
| ae4_5pct | apical_current_A | 2.77880987e-25 | 1e-20 |
| ae4_5pct | basolateral_current_A | 2.492850948e-23 | 1e-20 |
| ae4_5pct | buffer_site_accounting_fmol_s | 0 | 1e-10 |
| ae4_5pct | carbon_accounting_fmol_s | 5.551115123e-17 | 1e-10 |
| ae4_5pct | cell_bulk_charge_rate_fmol_s | 2.586125758e-13 | 1e-10 |
| ae4_5pct | cell_speciation_alkalinity_mM | 6.160405519e-12 | 1e-09 |
| ae4_5pct | homeostasis_charge_fmol_s | 6.808789643e-17 | 1e-10 |
| ae4_5pct | lumen_bulk_charge_rate_minus_outflow_fmol_s | 2.846681224e-15 | 1e-10 |
| ae4_5pct | lumen_speciation_alkalinity_mM | 5.324629626e-13 | 1e-09 |
| ae4_5pct | water_volume_accounting_pL_s | 4.33680869e-19 | 1e-12 |
| ae4_5pct | source_signature_residual_fmol_s | 0 | 0 |
| ae4_5pct | ae4_charge_residual_fmol_s | 0 | 1e-10 |
| ae4_5pct | cell_balance_residual_fmol_s | 6.938893904e-17 | 1e-10 |
| ae4_null | ae4_charge_fmol_s | 0 | 1e-10 |
| ae4_null | apical_current_A | 2.972680326e-25 | 1e-20 |
| ae4_null | basolateral_current_A | 2.4938203e-23 | 1e-20 |
| ae4_null | buffer_site_accounting_fmol_s | 0 | 1e-10 |
| ae4_null | carbon_accounting_fmol_s | 5.551115123e-17 | 1e-10 |
| ae4_null | cell_bulk_charge_rate_fmol_s | 2.574954139e-13 | 1e-10 |
| ae4_null | cell_speciation_alkalinity_mM | 5.99698069e-12 | 1e-09 |
| ae4_null | homeostasis_charge_fmol_s | 6.743737513e-17 | 1e-10 |
| ae4_null | lumen_bulk_charge_rate_minus_outflow_fmol_s | 3.051378594e-15 | 1e-10 |
| ae4_null | lumen_speciation_alkalinity_mM | 5.253575353e-13 | 1e-09 |
| ae4_null | water_volume_accounting_pL_s | 4.33680869e-19 | 1e-12 |
| ae4_null | source_signature_residual_fmol_s | 0 | 0 |
| ae4_null | ae4_charge_residual_fmol_s | 0 | 1e-10 |
| ae4_null | cell_balance_residual_fmol_s | 5.551115123e-17 | 1e-10 |

Full WT state vectors, root diagnostics and exact failed gates are in wt_rest.json. Each saved timeseries contains intracellular Na, K, Cl, pH, volume, secretion and transport fluxes. Files ending in all_window_integrals.csv contain AE4, NKCC1, AE2, NBC, NHE1, pump, CaCC and paracellular integrals over all three declared windows. Verification files retain endpoints, monitoring details, software, source hashes and failures.

These are predictions under the inherited Task 40 cycle law. Phenotype interpretation remains locked.
