"""Fixed seven-state comparison chassis for AE4 mechanism reconstruction.

The module is an independent Python implementation of the *documented*
Task-11 seven-state historical branch.  It does not import or execute anything
under :mod:`archive`, and it deliberately separates all AE4-specific kinetics
and stoichiometry behind one callable boundary.

Concentrations are in mM, cell height is in micrometres, voltage is in mV,
and fluxes retain the internally consistent but physically undocumented
historical chassis unit.  Consequently, absolute flow and physical time must
not be inferred from this implementation; ratios under a common protocol are
the defensible comparison quantities.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import hashlib
import json
import math
from typing import Any, Callable, Mapping, Protocol, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp


STATE_ORDER = (
    "na_l",
    "k_l",
    "height",
    "na_i",
    "k_i",
    "cl_i",
    "hco3_i",
)

BASELINE_STATE = np.array((118.7, 5.6, 28.7, 25.0, 120.0, 50.0, 10.0))


@dataclass(frozen=True)
class AE4Environment:
    """Concentrations presented to an AE4 law, all in mM.

    This small structural interface is intentionally independent of a
    particular mechanism implementation.  Candidate modules may accept this
    object directly or copy its eight named attributes into their own context.
    """

    na_i: float
    k_i: float
    cl_i: float
    hco3_i: float
    na_e: float
    k_e: float
    cl_e: float
    hco3_e: float
    # Dimensionless beta-adrenergic/cAMP/PKA stimulus occupancy.  It is kept
    # separate from calcium because native AE4 is activated by the former but
    # not by muscarinic/Ca signalling.  The historical/default protocol is 0.
    pka_activation: float = 0.0


@dataclass(frozen=True)
class AE4Balance:
    """AE4 source terms for intracellular amount balances.

    Positive ``cycle_flux`` denotes the candidate's declared forward cycle.
    The four species fields are signed *intracellular amount-rate* terms and
    are inserted directly into the corresponding chassis balances.  For the
    1 Cl : 1 cation : 2 HCO3 historical direction, for example, they are
    ``(+Cl, -cation, -2 HCO3)``.
    """

    cycle_flux: float
    na_i: float
    k_i: float
    cl_i: float
    hco3_i: float
    transported_charge: float = 0.0
    affinity: float = math.nan
    valid: bool = True
    diagnostics: Mapping[str, Any] = field(default_factory=dict)

    @property
    def charge_rate(self) -> float:
        """Net inward charge-equivalent rate implied by the source vector."""

        return self.na_i + self.k_i - self.cl_i - self.hco3_i


class AE4Law(Protocol):
    """Structural protocol accepted by :class:`FixedChassis`."""

    def __call__(
        self, environment: AE4Environment, activity_scale: float
    ) -> AE4Balance:
        ...


@dataclass(frozen=True)
class HistoricalSodiumOnlyAE4:
    """Documented Task-11 Na-only negative control.

    This law exists solely to reproduce the internally closed historical
    comparison baseline.  It is historical implementation evidence, not a
    primary-experimental or published-mechanism claim.
    """

    capacity: float = 18.9502607221547
    forward_rate: float = 0.0192
    reverse_rate: float = 0.00159305881398359

    def __call__(
        self, environment: AE4Environment, activity_scale: float
    ) -> AE4Balance:
        forward = (
            self.forward_rate
            * environment.cl_e
            * environment.hco3_i**2
            * environment.na_i
        )
        reverse = (
            self.reverse_rate
            * environment.cl_i
            * environment.hco3_e**2
            * environment.na_e
        )
        flux = float(activity_scale * self.capacity * (forward - reverse))
        affinity = math.log(forward / reverse) if forward > 0 and reverse > 0 else math.nan
        return AE4Balance(
            cycle_flux=flux,
            na_i=-flux,
            k_i=0.0,
            cl_i=flux,
            hco3_i=-2.0 * flux,
            transported_charge=0.0,
            affinity=affinity,
            diagnostics={
                "mechanism_id": "C1_historical_na_only",
                "forward_term": forward,
                "reverse_term": reverse,
                "provenance": "historical implementation evidence",
            },
        )


@dataclass(frozen=True)
class ZeroAE4:
    """Zero source used to expose the AE4 demand of the fixed chassis."""

    def __call__(
        self, environment: AE4Environment, activity_scale: float
    ) -> AE4Balance:
        del environment, activity_scale
        return AE4Balance(
            cycle_flux=0.0,
            na_i=0.0,
            k_i=0.0,
            cl_i=0.0,
            hco3_i=0.0,
            transported_charge=0.0,
            diagnostics={"mechanism_id": "zero_ae4"},
        )


@dataclass(frozen=True)
class ChassisParameters:
    """Fixed non-AE4 constants of the Task-11 seven-state branch.

    Names are descriptive replacements for the abbreviated historical names.
    Values are frozen here so simulations never read a mutable historical
    parameter file.  The source/status of each family is documented in
    ``analysis/12_ae4_mechanism_reconstruction/chassis.md``.
    """

    # Algebraic acid/base subsystem.
    buffer_scale: float = 0.199652558943136
    co2_transport_scale: float = 1970.0
    buffer_forward: float = 0.132
    buffer_reverse: float = 312.0
    impermeant_charge_amount: float = 2439.50353087137
    co2_elimination_denominator: float = 3939.9145458622197
    co2_external_lumen_sum: float = 13.2

    # Na/K pump and physical/electrical constants.
    nak_capacity: float = 0.00138020583534274
    nak_rate: float = 1_305_000.0
    k_e: float = 5.3
    nak_alpha: float = 0.641
    mm_to_molar: float = 0.001
    thermal_voltage_mv: float = 26.7137302360965
    faraday: float = 96485.3365

    # Historical channel gates and fitted conductances.
    cacc_half_activation: float = 100.26
    kcc_half_activation: float = 100.3
    cacc_hill: float = 1.49
    kcc_hill: float = 1.7
    tight_na_conductance: float = 357.863312079523
    tight_k_conductance: float = 25.893178360928
    cacc_conductance: float = 20_468_119.8347071
    kcc_conductance: float = 87_495_178.2340436

    # Fixed non-AE4 transporters.
    nkcc_capacity: float = 0.0063812
    nkcc_forward: float = 157.55
    nkcc_reverse: float = 2.0096e-5
    nkcc_denominator_forward: float = 1.0306
    nkcc_denominator_reverse: float = 1.3852e-6
    nhe1_capacity: float = 0.674737308094733
    ae2_capacity: float = 0.399084872489615

    # Water subsystem and geometry.
    water_apical: float = 0.000579346121973713
    water_basolateral: float = 5.45083716910752e-5
    water_paracellular: float = 7.31969562708724e-6
    lumen_impermeant_osmolyte: float = 48.8001357236744
    external_osmotic_sum: float = 288.1
    area_to_litre_per_height: float = 4.52961672473868e-14

    # Fixed baths and exchanger constants.
    cl_e: float = 102.6
    na_e: float = 140.2
    hco3_e_transport: float = 21.0
    ae2_cl_half: float = 5.6
    ae2_hco3_half: float = 0.0001
    nhe1_na_half: float = 15.0
    nhe1_h_half: float = 0.00045
    h_e: float = 3.8904514499428e-5

    # The historical seven-state time scaling is retained only for like-for-like
    # trajectories.  It has no established physical unit.
    rhs_scale: float = 10_000.0

    def fingerprint(self) -> str:
        """Stable SHA-256 fingerprint of every fixed non-AE4 number."""

        payload = json.dumps(asdict(self), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class Scenario:
    """Genotype/activity switches; no chassis parameter override is allowed."""

    name: str = "wt"
    ae4_activity_scale: float = 1.0
    ae2_activity_scale: float = 1.0

    def __post_init__(self) -> None:
        if self.ae4_activity_scale < 0 or self.ae2_activity_scale < 0:
            raise ValueError("activity scales must be nonnegative")


WT = Scenario("wt", 1.0, 1.0)
AE4_KNOCKOUT = Scenario("ae4_knockout", 0.0, 1.0)
AE2_KNOCKOUT = Scenario("ae2_knockout", 1.0, 0.0)


def historical_calcium_step(t: float) -> float:
    """Historical 0.05 -> 0.55 step; numeric unit is micromolar."""

    return 0.05 if t <= 100.0 else 0.55


def no_pka_activation(_t: float) -> float:
    """Default PKA protocol for exact historical-lineage reproduction."""

    return 0.0


@dataclass(frozen=True)
class SimulationResult:
    """Trajectory and reproducibility metadata returned by ``simulate``."""

    time: NDArray[np.float64]
    state: NDArray[np.float64]
    diagnostics: tuple[Mapping[str, Any], ...]
    scenario: Scenario
    parameter_fingerprint: str
    method: str
    rtol: float
    atol: float

    @property
    def endpoint(self) -> Mapping[str, Any]:
        return self.diagnostics[-1]


class FixedChassis:
    """Fixed non-AE4 salivary chassis with a swappable AE4 source law."""

    def __init__(
        self,
        ae4_law: Any | None = None,
        *,
        parameters: ChassisParameters | None = None,
        calcium_protocol: Callable[[float], float] = historical_calcium_step,
        pka_protocol: Callable[[float], float] = no_pka_activation,
        pka_discontinuities: Sequence[float] = (),
        enforce_electroneutral_ae4: bool = True,
        charge_tolerance: float = 1e-10,
    ) -> None:
        self.parameters = parameters or ChassisParameters()
        self.ae4_law = ae4_law or HistoricalSodiumOnlyAE4()
        self.calcium_protocol = calcium_protocol
        self.pka_protocol = pka_protocol
        self.pka_discontinuities = tuple(
            sorted({float(value) for value in pka_discontinuities})
        )
        self.enforce_electroneutral_ae4 = enforce_electroneutral_ae4
        self.charge_tolerance = charge_tolerance
        self._parameter_fingerprint = self.parameters.fingerprint()

    @property
    def parameter_fingerprint(self) -> str:
        return self._parameter_fingerprint

    @staticmethod
    def initial_state() -> NDArray[np.float64]:
        return BASELINE_STATE.copy()

    def _ae4_environment(
        self, y: NDArray[np.float64], *, pka_activation: float = 0.0
    ) -> AE4Environment:
        p = self.parameters
        return AE4Environment(
            na_i=float(y[3]),
            k_i=float(y[4]),
            cl_i=float(y[5]),
            hco3_i=float(y[6]),
            na_e=p.na_e,
            k_e=p.k_e,
            cl_e=p.cl_e,
            hco3_e=p.hco3_e_transport,
            pka_activation=float(pka_activation),
        )

    def _ae4_balance(
        self, y: NDArray[np.float64], scale: float, *, pka_activation: float = 0.0
    ) -> Any:
        environment = self._ae4_environment(y, pka_activation=pka_activation)
        law = self.ae4_law
        if hasattr(law, "evaluate"):
            contribution = law.evaluate(environment, activity_scale=scale)
        else:
            contribution = law(environment, scale)

        required = ("cycle_flux", "na_i", "k_i", "cl_i", "hco3_i")
        missing = [name for name in required if not hasattr(contribution, name)]
        if missing:
            raise TypeError(f"AE4 contribution lacks fields: {', '.join(missing)}")
        if not bool(getattr(contribution, "valid", True)):
            raise ValueError("AE4 candidate marked its contribution invalid")

        charge_rate = (
            float(contribution.na_i)
            + float(contribution.k_i)
            - float(contribution.cl_i)
            - float(contribution.hco3_i)
        )
        declared_charge = float(getattr(contribution, "transported_charge", 0.0))
        scale_for_tol = max(
            1.0,
            abs(float(contribution.na_i)),
            abs(float(contribution.k_i)),
            abs(float(contribution.cl_i)),
            abs(float(contribution.hco3_i)),
        )
        if self.enforce_electroneutral_ae4 and (
            abs(charge_rate) > self.charge_tolerance * scale_for_tol
            or abs(declared_charge) > self.charge_tolerance
        ):
            raise ValueError(
                "AE4 source is not electroneutral: "
                f"source charge={charge_rate:.6g}, declared={declared_charge:.6g}"
            )
        return contribution

    @staticmethod
    def _state(
        y: ArrayLike, *, require_all_positive: bool = True
    ) -> NDArray[np.float64]:
        state = np.asarray(y, dtype=float)
        if state.shape != (7,):
            raise ValueError(f"expected seven-state vector, got {state.shape}")
        if not np.all(np.isfinite(state)):
            raise ValueError("all seven chassis states must be finite")
        positive_indices = np.arange(7) if require_all_positive else np.arange(6)
        if np.any(state[positive_indices] <= 0):
            requirement = "all seven" if require_all_positive else "geometric and Nernst"
            raise ValueError(f"{requirement} chassis states must be positive")
        return state

    def _evaluate(
        self,
        y: ArrayLike,
        *,
        calcium: float,
        pka_activation: float = 0.0,
        scenario: Scenario = WT,
        include_diagnostics: bool,
    ) -> tuple[NDArray[np.float64], Mapping[str, Any]]:
        """Core evaluator; diagnostics can be skipped inside solver calls."""

        state = self._state(y, require_all_positive=include_diagnostics)
        if not math.isfinite(calcium) or calcium <= 0:
            raise ValueError("calcium must be finite and positive")
        if (
            not math.isfinite(pka_activation)
            or pka_activation < 0.0
            or pka_activation > 1.0
        ):
            raise ValueError("pka_activation must be finite and in [0, 1]")
        p = self.parameters
        na_l, k_l, height, na_i, k_i, cl_i, hco3_i = map(float, state)

        # Exact algebraic electroneutrality and historical CO2 closure.
        h_i = hco3_i + cl_i + p.impermeant_charge_amount / height - na_i - k_i
        co2_i = (
            p.co2_transport_scale * p.co2_external_lumen_sum
            - p.buffer_scale * p.buffer_reverse * hco3_i * h_i
        ) / p.co2_elimination_denominator
        cl_l = na_l + k_l

        # Na/K pump using the documented M-valued kinetic convention.
        na_m = p.mm_to_molar * na_i
        ke_m = p.mm_to_molar * p.k_e
        nak_turnover = p.nak_rate * ke_m**2 * na_m**3 / (
            ke_m**2 + p.nak_alpha * na_m**3
        )
        j_nak = p.nak_capacity * nak_turnover

        # Historical channel gates and Nernst potentials.
        p_cl = 1.0 / (1.0 + (p.cacc_half_activation / calcium) ** p.cacc_hill)
        p_k = 1.0 / (1.0 + (p.kcc_half_activation / calcium) ** p.kcc_hill)
        g_cl = p.cacc_conductance * p_cl
        g_k = p.kcc_conductance * p_k
        v_cl_mag = p.thermal_voltage_mv * math.log(cl_l / cl_i)
        v_k = p.thermal_voltage_mv * math.log(p.k_e / k_i)
        v_t_na = p.thermal_voltage_mv * math.log(na_l / p.na_e)
        v_t_k = p.thermal_voltage_mv * math.log(k_l / p.k_e)

        # Independently solve the two linear QSS current equations rather than
        # reproducing the historical closed-form voltage expressions.
        g_t = p.tight_na_conductance + p.tight_k_conductance
        matrix = np.array(((g_cl + g_t, -g_t), (-g_t, g_k + g_t)), dtype=float)
        rhs = np.array(
            (
                -g_cl * v_cl_mag
                + p.tight_na_conductance * v_t_na
                + p.tight_k_conductance * v_t_k,
                -p.faraday * j_nak
                + g_k * v_k
                - p.tight_na_conductance * v_t_na
                - p.tight_k_conductance * v_t_k,
            ),
            dtype=float,
        )
        try:
            v_a, v_b = np.linalg.solve(matrix, rhs)
        except np.linalg.LinAlgError as exc:
            raise ValueError("singular membrane-current closure") from exc
        v_t = float(v_a - v_b)
        j_cl_signed = g_cl * (v_a + v_cl_mag) / p.faraday
        j_k = g_k * (v_b - v_k) / p.faraday
        j_t_na = p.tight_na_conductance * (v_t - v_t_na) / p.faraday
        j_t_k = p.tight_k_conductance * (v_t - v_t_k) / p.faraday

        # Water and total lumen flow in the chassis-native scale.
        q_a = p.water_apical * (
            2.0 * (na_l + k_l - na_i - k_i - h_i)
            - co2_i
            + p.lumen_impermeant_osmolyte
        )
        q_b = p.water_basolateral * (
            2.0 * (na_i + k_i + h_i) + co2_i - p.external_osmotic_sum
        )
        q_t = p.water_paracellular * (
            2.0 * (na_l + k_l)
            + p.lumen_impermeant_osmolyte
            - p.external_osmotic_sum
        )
        q_total = q_a + q_t
        volume_rate = q_b - q_a

        # Fixed non-AE4 transporter laws.
        fourth_order = na_i * k_i * cl_i**2
        j_nkcc = p.nkcc_capacity * (
            (p.nkcc_forward - p.nkcc_reverse * fourth_order)
            / (
                p.nkcc_denominator_forward
                + p.nkcc_denominator_reverse * fourth_order
            )
        )
        j_ae2 = scenario.ae2_activity_scale * p.ae2_capacity * (
            (p.cl_e / (p.cl_e + p.ae2_cl_half))
            * (hco3_i / (hco3_i + p.ae2_hco3_half))
            - (cl_i / (cl_i + p.ae2_cl_half))
            * (p.hco3_e_transport / (p.hco3_e_transport + p.ae2_hco3_half))
        )
        j_nhe1 = p.nhe1_capacity * (
            (p.na_e / (p.na_e + p.nhe1_na_half))
            * (h_i / (p.nhe1_h_half + h_i))
            - (na_i / (na_i + p.nhe1_na_half))
            * (p.h_e / (p.nhe1_h_half + p.h_e))
        )
        j_buffer = p.buffer_scale * (
            p.buffer_forward * co2_i - p.buffer_reverse * hco3_i * h_i
        )
        ae4 = self._ae4_balance(
            state,
            scenario.ae4_activity_scale,
            pka_activation=pka_activation,
        )

        # Amount-rate numerators before concentration dilution.
        amount_na = j_nkcc - 3.0 * j_nak + j_nhe1 + float(ae4.na_i)
        amount_k = j_nkcc + 2.0 * j_nak - j_k + float(ae4.k_i)
        amount_cl = 2.0 * j_nkcc + j_ae2 + j_cl_signed + float(ae4.cl_i)
        amount_hco3 = j_buffer - j_ae2 + float(ae4.hco3_i)

        raw_rhs = np.array(
            (
                j_t_na - q_total * na_l,
                j_t_k - q_total * k_l,
                volume_rate,
                (amount_na - volume_rate * na_i) / height,
                (amount_k - volume_rate * k_i) / height,
                (amount_cl - volume_rate * cl_i) / height,
                (amount_hco3 - volume_rate * hco3_i) / height,
            ),
            dtype=float,
        )

        if not include_diagnostics:
            return raw_rhs, {}

        # Exact-reduction checks.  The inferred H amount balance follows by
        # differentiating electroneutrality; it must equal buffer - NHE1.
        inferred_h_amount = -amount_na - amount_k + amount_cl + amount_hco3
        ae4_charge = (
            float(ae4.na_i)
            + float(ae4.k_i)
            - float(ae4.cl_i)
            - float(ae4.hco3_i)
        )
        apical_current = j_cl_signed + j_t_na + j_t_k
        basolateral_current = j_nak + j_k - j_t_na - j_t_k
        lumen_cl_redundancy = (
            raw_rhs[0] + raw_rhs[1] - (-j_cl_signed - q_total * cl_l)
        )
        checks = {
            "cell_electroneutrality_mM": na_i
            + k_i
            + h_i
            - cl_i
            - hco3_i
            - p.impermeant_charge_amount / height,
            "lumen_electroneutrality_mM": na_l + k_l - cl_l,
            "ae4_charge_rate": ae4_charge,
            "apical_current_residual": apical_current,
            "basolateral_current_residual": basolateral_current,
            "proton_balance_identity_residual": inferred_h_amount
            - (j_buffer - j_nhe1),
            "lumen_chloride_redundancy_residual": lumen_cl_redundancy,
            "water_balance_definition_residual": raw_rhs[2] - (q_b - q_a),
        }

        diagnostics: dict[str, Any] = {
            "calcium": float(calcium),
            "pka_activation": float(pka_activation),
            "h_i": h_i,
            "pH_i": math.log10(1000.0 / h_i) if h_i > 0 else math.nan,
            "co2_i": co2_i,
            "cl_l": cl_l,
            "cell_volume_pL": p.area_to_litre_per_height * height * 1e12,
            "v_a": float(v_a),
            "v_b": float(v_b),
            "v_t": v_t,
            "q_a": q_a,
            "q_b": q_b,
            "q_t": q_t,
            "q_total": q_total,
            "volume_rate": volume_rate,
            "j_nak": j_nak,
            "j_nkcc1": j_nkcc,
            "j_ae2": j_ae2,
            "j_nhe1": j_nhe1,
            "j_buffer": j_buffer,
            "j_cl_signed": j_cl_signed,
            "j_k": j_k,
            "j_t_na": j_t_na,
            "j_t_k": j_t_k,
            "ae4_cycle_flux": float(ae4.cycle_flux),
            "ae4_na_source": float(ae4.na_i),
            "ae4_k_source": float(ae4.k_i),
            "ae4_cl_source": float(ae4.cl_i),
            "ae4_hco3_source": float(ae4.hco3_i),
            "ae4_affinity": (
                math.nan
                if getattr(ae4, "affinity", None) is None
                else float(ae4.affinity)
            ),
            "ae4_valid": bool(getattr(ae4, "valid", True)),
            "ae4_diagnostics": dict(getattr(ae4, "diagnostics", {})),
            "amount_rate_na": amount_na,
            "amount_rate_k": amount_k,
            "amount_rate_cl": amount_cl,
            "amount_rate_hco3": amount_hco3,
            "checks": checks,
            "scenario": scenario.name,
            "parameter_fingerprint": self.parameter_fingerprint,
            "raw_rhs": dict(zip(STATE_ORDER, map(float, raw_rhs), strict=True)),
            "units": {
                "concentration": "mM",
                "height": "micrometre",
                "voltage": "mV",
                "cell_volume": "pL",
                "flux": "historical chassis-native amount/area/time unit",
                "flow": "historical chassis-native height/time unit",
                "time": "historical scaled coordinate (not physical time)",
            },
        }
        return raw_rhs, diagnostics

    def evaluate(
        self,
        y: ArrayLike,
        *,
        calcium: float,
        pka_activation: float = 0.0,
        scenario: Scenario = WT,
    ) -> tuple[NDArray[np.float64], Mapping[str, Any]]:
        """Return the unscaled RHS and a complete diagnostic ledger."""

        return self._evaluate(
            y,
            calcium=calcium,
            pka_activation=pka_activation,
            scenario=scenario,
            include_diagnostics=True,
        )

    def ae4_source_demand(
        self,
        y: ArrayLike | None = None,
        *,
        calcium: float = 0.05,
        ae2_activity_scale: float = 1.0,
    ) -> Mapping[str, Any]:
        """Return the AE4 vector required to close four cell balances.

        The calculation sets AE4 identically to zero while retaining every
        non-AE4 number.  It therefore exposes a WT-closure target without using
        any knockout phenotype.  Values are amount-rate sources, before height
        division and dilution.
        """

        state = self.initial_state() if y is None else self._state(y)
        zero = FixedChassis(
            ZeroAE4(),
            parameters=self.parameters,
            calcium_protocol=self.calcium_protocol,
            enforce_electroneutral_ae4=self.enforce_electroneutral_ae4,
            charge_tolerance=self.charge_tolerance,
        )
        scenario = Scenario(
            "ae4_demand_probe",
            ae4_activity_scale=0.0,
            ae2_activity_scale=ae2_activity_scale,
        )
        raw, diagnostics = zero.evaluate(state, calcium=calcium, scenario=scenario)
        required = {
            "na_i": -float(diagnostics["amount_rate_na"]),
            "k_i": -float(diagnostics["amount_rate_k"]),
            "cl_i": -float(diagnostics["amount_rate_cl"]),
            "hco3_i": -float(diagnostics["amount_rate_hco3"]),
        }
        charge = required["na_i"] + required["k_i"] - required["cl_i"] - required["hco3_i"]
        return {
            "required_ae4_source": required,
            "required_charge_rate": charge,
            "raw_rhs_without_ae4": dict(zip(STATE_ORDER, map(float, raw), strict=True)),
            "non_ae4_amount_rates": {
                "na_i": float(diagnostics["amount_rate_na"]),
                "k_i": float(diagnostics["amount_rate_k"]),
                "cl_i": float(diagnostics["amount_rate_cl"]),
                "hco3_i": float(diagnostics["amount_rate_hco3"]),
            },
            "calcium": calcium,
            "state": dict(zip(STATE_ORDER, map(float, state), strict=True)),
            "parameter_fingerprint": self.parameter_fingerprint,
            "calibration_information": "WT baseline closure only; no KO datum",
        }

    def baseline_partition_mismatch(self) -> Mapping[str, Any]:
        """Mismatch caused by inserting the published Na/K fraction at WT.

        Cl and HCO3 signatures are left identical to the required historical
        cycle.  This isolates the cation-allocation mismatch without changing
        a chassis parameter or using a knockout observation.
        """

        demand = self.ae4_source_demand()
        required = demand["required_ae4_source"]
        flux = float(required["cl_i"])
        alpha_na = BASELINE_STATE[3] / (BASELINE_STATE[3] + BASELINE_STATE[4])
        alpha_k = 1.0 - alpha_na
        published_source = {
            "na_i": -alpha_na * flux,
            "k_i": -alpha_k * flux,
            "cl_i": flux,
            "hco3_i": -2.0 * flux,
        }
        residual = {
            species: float(published_source[species] - required[species])
            for species in published_source
        }
        return {
            "cycle_flux": flux,
            "na_fraction": alpha_na,
            "k_fraction": alpha_k,
            "required_source": required,
            "published_partition_source": published_source,
            "amount_balance_residual": residual,
            "expected_na_k_magnitude": alpha_k * flux,
            "parameter_fingerprint": self.parameter_fingerprint,
        }

    def raw_rhs(
        self, t: float, y: ArrayLike, scenario: Scenario = WT
    ) -> NDArray[np.float64]:
        raw, _ = self._evaluate(
            y,
            calcium=float(self.calcium_protocol(float(t))),
            pka_activation=float(self.pka_protocol(float(t))),
            scenario=scenario,
            include_diagnostics=False,
        )
        return raw

    def rhs(
        self, t: float, y: ArrayLike, scenario: Scenario = WT
    ) -> NDArray[np.float64]:
        """Scaled RHS used for like-for-like historical trajectories."""

        return self.parameters.rhs_scale * self.raw_rhs(t, y, scenario)

    def simulate(
        self,
        scenario: Scenario = WT,
        *,
        t_span: tuple[float, float] = (0.0, 200.0),
        t_eval: Sequence[float] | None = None,
        method: str = "BDF",
        rtol: float = 1e-6,
        atol: float = 1e-6,
        y0: ArrayLike | None = None,
    ) -> SimulationResult:
        """Integrate with an explicit split at the historical calcium step."""

        start, end = map(float, t_span)
        if not start < end:
            raise ValueError("t_span must be increasing")
        initial = self._state(self.initial_state() if y0 is None else y0)
        if t_eval is None:
            requested = np.arange(math.ceil(start), math.floor(end) + 1, dtype=float)
            if requested.size == 0 or requested[0] != start:
                requested = np.insert(requested, 0, start)
            if requested[-1] != end:
                requested = np.append(requested, end)
        else:
            requested = np.asarray(t_eval, dtype=float)
        if (
            requested.ndim != 1
            or requested.size == 0
            or np.any(np.diff(requested) <= 0)
            or requested[0] < start
            or requested[-1] > end
        ):
            raise ValueError("t_eval must be strictly increasing inside t_span")

        boundaries = [start]
        discontinuities = {100.0, *self.pka_discontinuities}
        boundaries.extend(
            value for value in sorted(discontinuities) if start < value < end
        )
        boundaries.append(end)
        pieces_t: list[NDArray[np.float64]] = []
        pieces_y: list[NDArray[np.float64]] = []
        current = initial
        for index, (left, right) in enumerate(zip(boundaries[:-1], boundaries[1:])):
            mask = (requested >= left) & (requested <= right)
            segment_eval = requested[mask]
            integration_eval = segment_eval
            if integration_eval.size == 0 or integration_eval[-1] != right:
                integration_eval = np.append(integration_eval, right)

            def segment_rhs(_t: float, state: NDArray[np.float64]) -> NDArray[np.float64]:
                raw, _ = self._evaluate(
                    state,
                    calcium=float(self.calcium_protocol(float(_t))),
                    pka_activation=float(self.pka_protocol(float(_t))),
                    scenario=scenario,
                    include_diagnostics=False,
                )
                return self.parameters.rhs_scale * raw

            solution = solve_ivp(
                segment_rhs,
                (left, right),
                current,
                method=method,
                t_eval=integration_eval,
                rtol=rtol,
                atol=atol,
            )
            if not solution.success:
                raise RuntimeError(f"chassis integration failed: {solution.message}")
            current = solution.y[:, -1]
            keep = np.isin(solution.t, segment_eval)
            segment_t = solution.t[keep]
            segment_y = solution.y[:, keep].T
            if index and segment_t.size and pieces_t[-1][-1] == segment_t[0]:
                segment_t = segment_t[1:]
                segment_y = segment_y[1:]
            if segment_t.size:
                pieces_t.append(segment_t)
                pieces_y.append(segment_y)

        time = np.concatenate(pieces_t)
        states = np.concatenate(pieces_y)
        diagnostic_rows: list[Mapping[str, Any]] = []
        for t, state in zip(time, states, strict=True):
            calcium = float(self.calcium_protocol(float(t)))
            _, diagnostic = self.evaluate(
                state,
                calcium=calcium,
                pka_activation=float(self.pka_protocol(float(t))),
                scenario=scenario,
            )
            diagnostic_rows.append(diagnostic)
        return SimulationResult(
            time=time,
            state=states,
            diagnostics=tuple(diagnostic_rows),
            scenario=scenario,
            parameter_fingerprint=self.parameter_fingerprint,
            method=method,
            rtol=rtol,
            atol=atol,
        )


def dimension_ledger() -> Mapping[str, str]:
    """Machine-readable statement of the defensible dimensional boundary."""

    return {
        "state_concentrations": "mM",
        "cell_height": "micrometre",
        "calcium_numeric": "micromolar",
        "potential": "mV",
        "cell_volume_observable": "pL",
        "all_balance_fluxes": "common historical chassis-native amount/area/time unit",
        "water_and_total_flow": "historical chassis-native height/time unit",
        "integration_time": "historical scaled coordinate; physical conversion unresolved",
        "absolute_flow": "not certified",
        "valid_primary_comparison": "dimensionless WT/KO ratios under one fixed chassis",
    }
