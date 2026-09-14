import math

from src.modern_full_model.membranes import (
    ElectricalEnvironment,
    evaluate_membrane_closure,
)
from src.modern_full_model.nbc_minimal import (
    DERIVED_AE4_CL_FMOL_S,
    DERIVED_REFERENCE_NBC_CAPACITY_FMOL_S,
    DERIVED_REQUIRED_NBC_CYCLE_FMOL_S,
    DERIVED_STIMULATED_NHE1_FMOL_S,
    FULLY_RECRUITED_CALCIUM_UM,
    RESTING_CALCIUM_UM,
    SOURCE_FIXED_NHE1_STIMULATED_MULTIPLIER,
    TARGET_NKCC1_POSITIVE_CL_SHARE,
    TASK31_R09_NKCC1_CL_FMOL_S,
    MinimalNbcParameters,
    evaluate_membrane_closure_with_nbc,
    evaluate_minimal_nbc,
    nbc_affinity_log,
    nhe1_stimulation_multiplier,
    normalized_secretory_activation,
)
from src.modern_full_model.task31_nhe1_mechanistic import (
    load_background,
    task31_parameters,
)


TASK31_CARRIER_FMOL = 2.339370005697548e-05


def _task31_parameters():
    background = load_background("R09")
    return task31_parameters(
        background,
        carrier_amount_fmol=TASK31_CARRIER_FMOL,
    )


def _task31_reference_environment(calcium_uM):
    return ElectricalEnvironment(
        na_i_mM=11.636125748680639,
        k_i_mM=116.76320932871538,
        cl_i_mM=60.30496692587399,
        hco3_i_mM=4.909374551508538,
        na_l_mM=124.68503924961041,
        k_l_mM=25.67998080632129,
        cl_l_mM=145.91441438704814,
        hco3_l_mM=4.447012642264572,
        na_e_mM=145.0,
        k_e_mM=5.0,
        cl_e_mM=126.16159280366088,
        hco3_e_mM=23.77832554186389,
        calcium_uM=calcium_uM,
    )


def test_stoichiometry_is_one_na_two_bicarbonate_and_electrogenic():
    p = _task31_parameters()
    flux = evaluate_minimal_nbc(
        na_i_mM=12.0,
        hco3_i_mM=5.0,
        na_o_mM=145.0,
        hco3_o_mM=24.0,
        v_basolateral_V=-0.070,
        thermal_voltage_V=p.constants.thermal_voltage_V,
        faraday_C_mol=p.constants.faraday_C_mol,
        activation_fraction=1.0,
    )
    assert flux.cycle_inward_fmol_s > 0.0
    assert flux.cl_cell_fmol_s == 0.0
    assert flux.na_cell_fmol_s == flux.cycle_inward_fmol_s
    assert flux.tic_cell_fmol_s == 2.0 * flux.cycle_inward_fmol_s
    assert flux.alkalinity_cell_fmol_s == 2.0 * flux.cycle_inward_fmol_s
    assert flux.transported_charge_equivalents_fmol_s == -flux.cycle_inward_fmol_s
    assert flux.charge_source_residual_fmol_s == -flux.cycle_inward_fmol_s
    assert flux.current_A > 0.0
    expected_current = (
        p.constants.faraday_C_mol * flux.cycle_inward_fmol_s * 1.0e-15
    )
    assert math.isclose(flux.current_A, expected_current, rel_tol=1e-14, abs_tol=0.0)


def test_voltage_enters_affinity_with_correct_reversal_sign():
    p = _task31_parameters()
    chemical = math.log((145.0 * 24.0**2) / (12.0 * 5.0**2))
    v_reversal = -p.constants.thermal_voltage_V * chemical
    affinity = nbc_affinity_log(
        na_i_mM=12.0,
        hco3_i_mM=5.0,
        na_o_mM=145.0,
        hco3_o_mM=24.0,
        v_basolateral_V=v_reversal,
        thermal_voltage_V=p.constants.thermal_voltage_V,
    )
    assert abs(affinity) < 1e-14

    forward = evaluate_minimal_nbc(
        na_i_mM=12.0,
        hco3_i_mM=5.0,
        na_o_mM=145.0,
        hco3_o_mM=24.0,
        v_basolateral_V=v_reversal + 0.020,
        thermal_voltage_V=p.constants.thermal_voltage_V,
        faraday_C_mol=p.constants.faraday_C_mol,
        activation_fraction=1.0,
    )
    reverse = evaluate_minimal_nbc(
        na_i_mM=12.0,
        hco3_i_mM=5.0,
        na_o_mM=145.0,
        hco3_o_mM=24.0,
        v_basolateral_V=v_reversal - 0.020,
        thermal_voltage_V=p.constants.thermal_voltage_V,
        faraday_C_mol=p.constants.faraday_C_mol,
        activation_fraction=1.0,
    )
    assert forward.cycle_inward_fmol_s > 0.0
    assert reverse.cycle_inward_fmol_s < 0.0


def test_rest_is_exact_zero_nbc_and_unit_nhe1_multiplier():
    p = _task31_parameters()
    activation = normalized_secretory_activation(RESTING_CALCIUM_UM)
    assert activation == 0.0
    assert nhe1_stimulation_multiplier(activation) == 1.0

    flux = evaluate_minimal_nbc(
        na_i_mM=11.636125748680639,
        hco3_i_mM=4.909374551508538,
        na_o_mM=145.0,
        hco3_o_mM=23.77832554186389,
        v_basolateral_V=-0.06710109365819838,
        thermal_voltage_V=p.constants.thermal_voltage_V,
        faraday_C_mol=p.constants.faraday_C_mol,
        activation_fraction=activation,
    )
    assert flux.cycle_inward_fmol_s == 0.0
    assert flux.current_A == 0.0
    assert flux.na_cell_fmol_s == 0.0
    assert flux.tic_cell_fmol_s == 0.0
    assert flux.alkalinity_cell_fmol_s == 0.0


def test_central_stimulus_is_full_nbc_and_source_fixed_2p3_nhe1():
    activation = normalized_secretory_activation(FULLY_RECRUITED_CALCIUM_UM)
    assert activation == 1.0
    assert nhe1_stimulation_multiplier(activation) == SOURCE_FIXED_NHE1_STIMULATED_MULTIPLIER


def test_reference_resting_current_closure_is_exactly_nested():
    p = _task31_parameters()
    environment = _task31_reference_environment(RESTING_CALCIUM_UM)
    canonical = evaluate_membrane_closure(environment, p)
    augmented, nbc = evaluate_membrane_closure_with_nbc(
        environment,
        p,
        activation_fraction=0.0,
    )
    assert augmented == canonical
    assert nbc.cycle_inward_fmol_s == 0.0
    assert nbc.current_A == 0.0


def test_stimulated_reference_current_closure_hits_derived_cycle_scale():
    p = _task31_parameters()
    environment = _task31_reference_environment(FULLY_RECRUITED_CALCIUM_UM)
    membranes, nbc = evaluate_membrane_closure_with_nbc(
        environment,
        p,
        activation_fraction=1.0,
    )
    assert math.isclose(
        MinimalNbcParameters().capacity_fmol_s,
        DERIVED_REFERENCE_NBC_CAPACITY_FMOL_S,
        rel_tol=0.0,
        abs_tol=0.0,
    )
    assert math.isclose(
        nbc.cycle_inward_fmol_s,
        DERIVED_REQUIRED_NBC_CYCLE_FMOL_S,
        rel_tol=2e-12,
        abs_tol=1e-14,
    )
    assert abs(membranes.current_residuals_A["apical"]) < 1e-18
    assert abs(membranes.current_residuals_A["basolateral"]) < 1e-18
    assert abs(membranes.charge_residuals_fmol_s["cell_source_minus_current"]) < 1e-10
    assert abs(membranes.charge_residuals_fmol_s["lumen_source_minus_current"]) < 1e-10


def test_70_30_algebra_and_ae4_expression_predictions_close_exactly():
    c = TASK31_R09_NKCC1_CL_FMOL_S
    a = DERIVED_AE4_CL_FMOL_S
    h = DERIVED_STIMULATED_NHE1_FMOL_S
    b = DERIVED_REQUIRED_NBC_CYCLE_FMOL_S

    assert math.isclose(c / (c + a), TARGET_NKCC1_POSITIVE_CL_SHARE, rel_tol=0.0, abs_tol=1e-15)
    assert math.isclose(h + 2.0 * b - 2.0 * a, 0.0, rel_tol=0.0, abs_tol=1e-15)

    ratio_5pct = (c + 0.05 * a) / (c + a)
    ratio_null = c / (c + a)
    assert math.isclose(ratio_5pct, 0.715, rel_tol=0.0, abs_tol=1e-15)
    assert math.isclose(ratio_null, 0.70, rel_tol=0.0, abs_tol=1e-15)
