import math

from src.modern_full_model.nbc_minimal import (
    REFERENCE_AFFINITY_LOG,
    REFERENCE_REQUIRED_NBC_FLUX_FMOL_S,
    MinimalNbcParameters,
    evaluate_minimal_nbc,
)


def test_reference_task31_state_gives_derived_required_flux():
    flux = evaluate_minimal_nbc(
        na_i_mM=11.636125748680639,
        hco3_i_mM=4.909374551508538,
        na_o_mM=145.0,
        hco3_o_mM=23.77832554186389,
    )
    assert math.isclose(flux.affinity_log, REFERENCE_AFFINITY_LOG, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(
        flux.inward_fmol_s,
        REFERENCE_REQUIRED_NBC_FLUX_FMOL_S,
        rel_tol=1e-12,
        abs_tol=1e-15,
    )


def test_nbc_is_exactly_electroneutral_and_chloride_neutral():
    flux = evaluate_minimal_nbc(
        na_i_mM=12.0,
        hco3_i_mM=5.0,
        na_o_mM=145.0,
        hco3_o_mM=24.0,
    )
    assert flux.cl_cell_fmol_s == 0.0
    assert flux.transported_charge_equivalents_fmol_s == 0.0
    assert flux.charge_source_residual_fmol_s == 0.0
    assert flux.na_cell_fmol_s == flux.tic_cell_fmol_s
    assert flux.na_cell_fmol_s == flux.alkalinity_cell_fmol_s


def test_nbc_reverses_when_intracellular_activity_product_exceeds_bath():
    flux = evaluate_minimal_nbc(
        na_i_mM=200.0,
        hco3_i_mM=30.0,
        na_o_mM=100.0,
        hco3_o_mM=10.0,
        parameters=MinimalNbcParameters(capacity_fmol_s=0.2, log_width=2.0),
    )
    assert flux.affinity_log < 0.0
    assert flux.inward_fmol_s < 0.0


def test_zero_activity_scale_nests_exactly_to_no_nbc_source():
    flux = evaluate_minimal_nbc(
        na_i_mM=12.0,
        hco3_i_mM=5.0,
        na_o_mM=145.0,
        hco3_o_mM=24.0,
        activity_scale=0.0,
    )
    assert flux.inward_fmol_s == 0.0
    assert flux.na_cell_fmol_s == 0.0
    assert flux.tic_cell_fmol_s == 0.0
    assert flux.alkalinity_cell_fmol_s == 0.0
