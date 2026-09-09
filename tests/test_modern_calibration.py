"""Focused WT-only tests for the Task-13B calibration/root layer."""

from __future__ import annotations

from dataclasses import replace
import unittest

import numpy as np

from modern_full_model import calibration as calibration_module
from modern_full_model.camp_pka import R1EffectiveActivation
from modern_full_model.calibration import (
    BoundedValue,
    ModelVariant,
    WTCalibrationParameters,
    WTCalibrationSpec,
    build_wt_model,
    decode_rest_coordinates,
    encode_electroneutral_coordinates,
    fit_wt_source_compatible_variant,
    solve_wt_resting_branches,
)


REFERENCE_COORDINATES = np.asarray(
    (
        14.756773495778988,
        111.10374637122908,
        6.194122572604099,
        6.91,
        1.3,
        146.13568410951905,
        4.087361120509596,
        12.891629639536427,
        7.276453648930861,
        0.10164106693766255,
    )
)
REFERENCE_PARAMETERS = WTCalibrationParameters(
    ae4_carrier_amount_fmol_at_unit_rate_gauge=0.15012642650706573,
    nhe1_capacity_fmol_s=0.007440980590414442,
    common_membrane_conductance_scale=10.0,
    cell_other_impermeant_osmoles_fmol=124.71681310732363,
)
REFERENCE_VARIANT = ModelVariant(
    "G2_BALANCED_APICAL_K_BIASED_G10_TEST",
    mixed_bath_na_attempt_fraction=0.10,
    apical_pump_fraction=0.10,
    apical_k_fraction=0.10,
)


def _narrow_reference_spec() -> WTCalibrationSpec:
    base = WTCalibrationSpec()
    coordinate_spans = (2.0, 8.0, 1.0, 0.05, 0.08, 8.0, 1.5, 2.0, 0.15, 0.025)
    coordinate_bounds = tuple(
        BoundedValue(
            item.name,
            float(target - span),
            float(target),
            float(target + span),
            item.units,
            item.tier,
            item.provenance,
            "test-local narrow multistart domain around the audited reference",
        )
        for item, target, span in zip(
            base.coordinate_bounds, REFERENCE_COORDINATES, coordinate_spans
        )
    )
    parameter_values = (
        REFERENCE_PARAMETERS.ae4_carrier_amount_fmol_at_unit_rate_gauge,
        REFERENCE_PARAMETERS.nhe1_capacity_fmol_s,
        REFERENCE_PARAMETERS.common_membrane_conductance_scale,
        REFERENCE_PARAMETERS.cell_other_impermeant_osmoles_fmol,
    )
    parameter_factors = (1.25, 1.35, 2.0, 1.20)
    parameter_bounds = tuple(
        BoundedValue(
            item.name,
            float(target / factor),
            float(target),
            float(target * factor),
            item.units,
            item.tier,
            item.provenance,
            "test-local narrow multistart domain around the audited reference",
        )
        for item, target, factor in zip(
            base.parameter_bounds, parameter_values, parameter_factors
        )
    )
    return replace(
        base,
        coordinate_bounds=coordinate_bounds,
        parameter_bounds=parameter_bounds,
    )


def test_electroneutral_parameterization_enforces_both_charge_manifolds() -> None:
    model = build_wt_model(REFERENCE_VARIANT, REFERENCE_PARAMETERS)
    state = encode_electroneutral_coordinates(model, REFERENCE_COORDINATES)
    evaluation = model.evaluate(0.0, state)

    assert max(abs(value) for value in evaluation.diagnostics.state_charge_fmol.values()) < 1e-12
    np.testing.assert_allclose(
        decode_rest_coordinates(model, state), REFERENCE_COORDINATES, rtol=0.0, atol=2e-13
    )
    raw = evaluation.rhs
    assert abs(raw[0] + raw[1] - raw[2] - raw[4]) < 1e-15
    assert abs(raw[6] + raw[7] - raw[8] - raw[10]) < 1e-15


def test_minimal_wt_calibration_recovers_interior_full_rank_root() -> None:
    report = fit_wt_source_compatible_variant(
        REFERENCE_VARIANT,
        spec=_narrow_reference_spec(),
        fixed_conductance_scale=10.0,
        fixed_other_impermeant_start_fmol=(
            REFERENCE_PARAMETERS.cell_other_impermeant_osmoles_fmol
        ),
        start_count=3,
        seed=13034,
        max_nfev=3000,
    )

    assert len(report.attempts) == 3
    assert len(report.branches) == 1
    branch = report.branches[0]
    assert branch.passes_numerical_gate
    assert branch.passes_wt_gate
    assert branch.independent_jacobian_rank == 13
    assert branch.independent_jacobian_nullity == 0
    assert branch.boundary_hits == ()
    assert abs(branch.observables["cl_i_mM"] - 50.10) < 1e-9
    assert abs(branch.observables["ph_i"] - 6.91) < 1e-10
    assert abs(branch.observables["volume_i_pL"] - 1.30) < 1e-10


def test_fixed_parameter_geometry_returns_all_detected_root_clusters() -> None:
    spec = _narrow_reference_spec()
    model = build_wt_model(REFERENCE_VARIANT, REFERENCE_PARAMETERS)
    report = solve_wt_resting_branches(
        model,
        variant=REFERENCE_VARIANT,
        calibration=REFERENCE_PARAMETERS,
        spec=spec,
        start_count=3,
        seed=13091,
        max_nfev=2000,
        reference_coordinates=REFERENCE_COORDINATES,
    )

    assert len(report.attempts) == 4  # reference plus three named design starts
    assert len(report.branches) == 1
    branch = report.branches[0]
    assert branch.passes_wt_gate
    assert branch.independent_jacobian_rank == 10
    assert branch.independent_jacobian_nullity == 0
    assert set(branch.member_start_ids).issubset(
        {attempt.start_id for attempt in report.attempts}
    )


class _OneStateRegulator:
    state_names = ("test_regulator",)

    def initial_state(self, *, beta_input: float) -> tuple[float]:
        return (0.2,)

    def evaluate(
        self, time_s: float, state: tuple[float, ...], beta_input: float
    ) -> dict[str, float]:
        return {"capacity_multiplier": 1.0, "test_regulator": float(state[0])}

    def rhs(
        self, time_s: float, state: tuple[float, ...], beta_input: float
    ) -> tuple[float]:
        return (-float(state[0]),)


def test_dynamic_regulatory_root_requires_explicit_s_inverse_residual_scale() -> None:
    dynamic_variant = replace(
        REFERENCE_VARIANT,
        variant_id="DYNAMIC_SCALE_TEST",
        regulatory_family="TEST_DYNAMIC",
    )
    model = build_wt_model(
        dynamic_variant,
        REFERENCE_PARAMETERS,
        regulatory_model=_OneStateRegulator(),
    )
    coordinates = np.r_[REFERENCE_COORDINATES, 0.2]

    try:
        calibration_module._independent_scaled_rhs(
            model, coordinates, WTCalibrationSpec()
        )
    except ValueError as error:
        assert "explicit numerical s^-1 scale" in str(error)
    else:
        raise AssertionError("missing regulatory residual scale was accepted")

    scaled = calibration_module._independent_scaled_rhs(
        model,
        coordinates,
        replace(WTCalibrationSpec(), regulatory_rhs_scales_s_inv=(0.5,)),
    )
    assert scaled.shape == (11,)
    assert abs(scaled[-1] + 0.4) < 1e-14


def test_r1_basal_zero_is_reported_as_source_supported_boundary() -> None:
    regulator = R1EffectiveActivation()
    dynamic_variant = replace(
        REFERENCE_VARIANT,
        variant_id="G3_R1_REST_TEST",
        regulatory_family="R1_EFFECTIVE_BETA_PKA_AE4",
    )
    model = build_wt_model(
        dynamic_variant,
        REFERENCE_PARAMETERS,
        regulatory_model=regulator,
    )
    spec = replace(
        _narrow_reference_spec(),
        regulatory_rhs_scales_s_inv=(1.0 / regulator.tau_activation_s,),
    )
    report = solve_wt_resting_branches(
        model,
        variant=dynamic_variant,
        calibration=REFERENCE_PARAMETERS,
        spec=spec,
        start_count=3,
        seed=13131,
        max_nfev=3000,
        reference_coordinates=np.r_[REFERENCE_COORDINATES, 1.0e-8],
    )

    assert len(report.branches) == 1
    branch = report.branches[0]
    name = "ae4_effective_activation_fraction"
    assert branch.passes_wt_gate
    assert branch.independent_jacobian_rank == 11
    assert name in branch.boundary_hits
    assert branch.source_supported_boundary_hits == (name,)
    assert branch.disallowed_boundary_hits == ()


class TestModernCalibration(unittest.TestCase):
    """Expose the focused checks to the repository's unittest runner."""

    def test_electroneutral_parameterization(self) -> None:
        test_electroneutral_parameterization_enforces_both_charge_manifolds()

    def test_minimal_wt_calibration(self) -> None:
        test_minimal_wt_calibration_recovers_interior_full_rank_root()

    def test_fixed_parameter_root_geometry(self) -> None:
        test_fixed_parameter_geometry_returns_all_detected_root_clusters()

    def test_regulatory_rhs_scaling(self) -> None:
        test_dynamic_regulatory_root_requires_explicit_s_inverse_residual_scale()

    def test_r1_basal_boundary_rule(self) -> None:
        test_r1_basal_zero_is_reported_as_source_supported_boundary()
