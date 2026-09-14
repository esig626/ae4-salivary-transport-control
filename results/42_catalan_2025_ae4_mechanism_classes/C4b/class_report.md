# Task 42 C4b: frozen predictions

WT rest: PASS. Forward affinity A/RT: -4.501077859 (opposed).

| case | status | Q(0,600) pL | q(600) pL/s | early deficit % | late deficit % | total deficit % | NKCC1 change % |
| --- | --- | --- | --- | --- | --- | --- | --- |
| wt | TASK42_SOURCE_CLASS_WT_FAILED | unavailable | unavailable | unavailable | unavailable | unavailable | unavailable |
| ae4_5pct | PASS | 0.962415362 | 0.001574130214 | unavailable | unavailable | unavailable | unavailable |
| ae4_null | PASS | 0.9541312327 | 0.001555640569 | unavailable | unavailable | unavailable | unavailable |

Early means 0 to 60 s, late means 60 to 600 s. Deficits and NKCC1 compensation use this class's WT. A failed WT trajectory leaves these ratios unavailable. All successful genotype outputs remain saved. No extrapolation of an incomplete trajectory is used.

Failure record: wt: ph_i at 317 s.

## WT rest

| Observable | Value |
| --- | --- |
| na_i_mM | 11.40177147 |
| k_i_mM | 116.6443434 |
| cl_i_mM | 60.99517069 |
| ph_i | 6.840843855 |
| volume_i_pL | 1.456772585 |
| q_out_pL_s | 0.001079940004 |

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
| wt | basolateral_current_A | 2.415625883e-23 | 1e-20 |
| wt | buffer_site_accounting_fmol_s | 0 | 1e-10 |
| wt | carbon_accounting_fmol_s | 5.551115123e-17 | 1e-10 |
| wt | cell_bulk_charge_rate_fmol_s | 2.511602037e-13 | 1e-10 |
| wt | cell_speciation_alkalinity_mM | 5.988098906e-12 | 1e-09 |
| wt | homeostasis_charge_fmol_s | 6.938893904e-17 | 1e-10 |
| wt | lumen_bulk_charge_rate_minus_outflow_fmol_s | 2.220012368e-15 | 1e-10 |
| wt | lumen_speciation_alkalinity_mM | 8.935074902e-13 | 1e-09 |
| wt | water_volume_accounting_pL_s | 4.33680869e-19 | 1e-12 |
| wt | source_signature_residual_fmol_s | 0 | 0 |
| wt | ae4_charge_residual_fmol_s | 0 | 1e-10 |
| wt | cell_balance_residual_fmol_s | 8.326672685e-17 | 1e-10 |
| ae4_5pct | ae4_charge_fmol_s | 0 | 1e-10 |
| ae4_5pct | apical_current_A | 2.584939414e-25 | 1e-20 |
| ae4_5pct | basolateral_current_A | 2.494466535e-23 | 1e-20 |
| ae4_5pct | buffer_site_accounting_fmol_s | 0 | 1e-10 |
| ae4_5pct | carbon_accounting_fmol_s | 5.551115123e-17 | 1e-10 |
| ae4_5pct | cell_bulk_charge_rate_fmol_s | 2.58466859e-13 | 1e-10 |
| ae4_5pct | cell_speciation_alkalinity_mM | 6.206590797e-12 | 1e-09 |
| ae4_5pct | homeostasis_charge_fmol_s | 6.808789643e-17 | 1e-10 |
| ae4_5pct | lumen_bulk_charge_rate_minus_outflow_fmol_s | 2.66973943e-15 | 1e-10 |
| ae4_5pct | lumen_speciation_alkalinity_mM | 6.403766406e-13 | 1e-09 |
| ae4_5pct | water_volume_accounting_pL_s | 2.168404345e-19 | 1e-12 |
| ae4_5pct | source_signature_residual_fmol_s | 0 | 0 |
| ae4_5pct | ae4_charge_residual_fmol_s | 0 | 1e-10 |
| ae4_5pct | cell_balance_residual_fmol_s | 8.326672685e-17 | 1e-10 |
| ae4_null | ae4_charge_fmol_s | 0 | 1e-10 |
| ae4_null | apical_current_A | 2.617251157e-25 | 1e-20 |
| ae4_null | basolateral_current_A | 2.491881595e-23 | 1e-20 |
| ae4_null | buffer_site_accounting_fmol_s | 0 | 1e-10 |
| ae4_null | carbon_accounting_fmol_s | 2.775557562e-17 | 1e-10 |
| ae4_null | cell_bulk_charge_rate_fmol_s | 2.583905312e-13 | 1e-10 |
| ae4_null | cell_speciation_alkalinity_mM | 6.142641951e-12 | 1e-09 |
| ae4_null | homeostasis_charge_fmol_s | 6.722053469e-17 | 1e-10 |
| ae4_null | lumen_bulk_charge_rate_minus_outflow_fmol_s | 2.693158196e-15 | 1e-10 |
| ae4_null | lumen_speciation_alkalinity_mM | 6.394884622e-13 | 1e-09 |
| ae4_null | water_volume_accounting_pL_s | 4.33680869e-19 | 1e-12 |
| ae4_null | source_signature_residual_fmol_s | 0 | 0 |
| ae4_null | ae4_charge_residual_fmol_s | 0 | 1e-10 |
| ae4_null | cell_balance_residual_fmol_s | 6.938893904e-17 | 1e-10 |

Full WT state vectors, root diagnostics and exact failed gates are in wt_rest.json. Each saved timeseries contains intracellular Na, K, Cl, pH, volume, secretion and transport fluxes. Files ending in all_window_integrals.csv contain AE4, NKCC1, AE2, NBC, NHE1, pump, CaCC and paracellular integrals over all three declared windows. Verification files retain endpoints, monitoring details, software, source hashes and failures.

These are predictions under the inherited Task 40 cycle law. Phenotype interpretation remains locked.
