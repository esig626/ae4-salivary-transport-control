# Task 42 C3b: frozen predictions

WT rest: PASS. Forward affinity A/RT: -4.949908274 (opposed).

| case | status | Q(0,600) pL | q(600) pL/s | early deficit % | late deficit % | total deficit % | NKCC1 change % |
| --- | --- | --- | --- | --- | --- | --- | --- |
| wt | PASS | 1.035524581 | 0.001714057909 | unavailable | unavailable | unavailable | unavailable |
| ae4_5pct | PASS | 0.9644678367 | 0.001576204583 | 1.566863442 | 7.488332246 | 6.861908006 | 10.22835141 |
| ae4_null | PASS | 0.956128264 | 0.0015579507 | 1.677376096 | 8.375881551 | 7.667255646 | 11.30431823 |

Early means 0 to 60 s, late means 60 to 600 s. Deficits and NKCC1 compensation use this class's WT. A failed WT trajectory leaves these ratios unavailable. All successful genotype outputs remain saved. No extrapolation of an incomplete trajectory is used.

Failure record: None.

## WT rest

| Observable | Value |
| --- | --- |
| na_i_mM | 11.3476966 |
| k_i_mM | 117.2274766 |
| cl_i_mM | 60.89846289 |
| ph_i | 6.912172389 |
| volume_i_pL | 1.459366845 |
| q_out_pL_s | 0.001080066708 |

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
| wt | apical_current_A | 2.132575017e-25 | 1e-20 |
| wt | basolateral_current_A | 2.435012928e-23 | 1e-20 |
| wt | buffer_site_accounting_fmol_s | 0 | 1e-10 |
| wt | carbon_accounting_fmol_s | 5.551115123e-17 | 1e-10 |
| wt | cell_bulk_charge_rate_fmol_s | 2.51271226e-13 | 1e-10 |
| wt | cell_speciation_alkalinity_mM | 5.751843446e-12 | 1e-09 |
| wt | homeostasis_charge_fmol_s | 6.591949209e-17 | 1e-10 |
| wt | lumen_bulk_charge_rate_minus_outflow_fmol_s | 2.282028733e-15 | 1e-10 |
| wt | lumen_speciation_alkalinity_mM | 5.83977311e-13 | 1e-09 |
| wt | water_volume_accounting_pL_s | 4.33680869e-19 | 1e-12 |
| wt | source_signature_residual_fmol_s | 0 | 0 |
| wt | ae4_charge_residual_fmol_s | 0 | 1e-10 |
| wt | cell_balance_residual_fmol_s | 8.326672685e-17 | 1e-10 |
| ae4_5pct | ae4_charge_fmol_s | 0 | 1e-10 |
| ae4_5pct | apical_current_A | 2.584939414e-25 | 1e-20 |
| ae4_5pct | basolateral_current_A | 2.476371959e-23 | 1e-20 |
| ae4_5pct | buffer_site_accounting_fmol_s | 0 | 1e-10 |
| ae4_5pct | carbon_accounting_fmol_s | 5.551115123e-17 | 1e-10 |
| ae4_5pct | cell_bulk_charge_rate_fmol_s | 2.576792946e-13 | 1e-10 |
| ae4_5pct | cell_speciation_alkalinity_mM | 5.911715562e-12 | 1e-09 |
| ae4_5pct | homeostasis_charge_fmol_s | 6.938893904e-17 | 1e-10 |
| ae4_5pct | lumen_bulk_charge_rate_minus_outflow_fmol_s | 2.705301261e-15 | 1e-10 |
| ae4_5pct | lumen_speciation_alkalinity_mM | 5.884182031e-13 | 1e-09 |
| ae4_5pct | water_volume_accounting_pL_s | 4.33680869e-19 | 1e-12 |
| ae4_5pct | source_signature_residual_fmol_s | 0 | 0 |
| ae4_5pct | ae4_charge_residual_fmol_s | 0 | 1e-10 |
| ae4_5pct | cell_balance_residual_fmol_s | 8.326672685e-17 | 1e-10 |
| ae4_null | ae4_charge_fmol_s | 0 | 1e-10 |
| ae4_null | apical_current_A | 2.746498128e-25 | 1e-20 |
| ae4_null | basolateral_current_A | 2.484773012e-23 | 1e-20 |
| ae4_null | buffer_site_accounting_fmol_s | 0 | 1e-10 |
| ae4_null | carbon_accounting_fmol_s | 5.551115123e-17 | 1e-10 |
| ae4_null | cell_bulk_charge_rate_fmol_s | 2.56177024e-13 | 1e-10 |
| ae4_null | cell_speciation_alkalinity_mM | 6.128431096e-12 | 1e-09 |
| ae4_null | homeostasis_charge_fmol_s | 6.591949209e-17 | 1e-10 |
| ae4_null | lumen_bulk_charge_rate_minus_outflow_fmol_s | 2.886579864e-15 | 1e-10 |
| ae4_null | lumen_speciation_alkalinity_mM | 6.039613254e-13 | 1e-09 |
| ae4_null | water_volume_accounting_pL_s | 4.33680869e-19 | 1e-12 |
| ae4_null | source_signature_residual_fmol_s | 0 | 0 |
| ae4_null | ae4_charge_residual_fmol_s | 0 | 1e-10 |
| ae4_null | cell_balance_residual_fmol_s | 5.551115123e-17 | 1e-10 |

Full WT state vectors, root diagnostics and exact failed gates are in wt_rest.json. Each saved timeseries contains intracellular Na, K, Cl, pH, volume, secretion and transport fluxes. Files ending in all_window_integrals.csv contain AE4, NKCC1, AE2, NBC, NHE1, pump, CaCC and paracellular integrals over all three declared windows. Verification files retain endpoints, monitoring details, software, source hashes and failures.

These are predictions under the inherited Task 40 cycle law. Phenotype interpretation remains locked.
