"""Task 42 source classes with the unchanged Task 40 scalar AE4 cycle.

The inherited AE4Contribution aliases TIC and TA to bicarbonate. This module
replaces only that contribution object, so carbonate has its own carbon and
alkalinity signatures. It does not supply class specific kinetics.
"""
from dataclasses import dataclass
from fractions import Fraction
import math
from types import MappingProxyType

from .ae4_equal_cation_routing import equal_routing_adapter
from .model import ModernFullModel
from .nbc_minimal import MinimalNbcModel
from .nkcc_stimulation import StimulatedNkcc1Model


@dataclass(frozen=True)
class MechanismClass:
    name: str
    na: float
    k: float
    cl: float
    hco3: float
    co3: float

    @property
    def conserved_coefficients(self):
        return (self.na, self.k, self.cl, self.hco3 + self.co3,
                self.hco3 + 2 * self.co3)

    @property
    def exact_charge(self):
        na, k, cl, _, ta = map(Fraction, self.conserved_coefficients)
        return na + k - cl - ta

    @property
    def species_coefficients(self):
        return dict(zip(('na', 'k', 'cl', 'hco3', 'co3'),
                        (self.na, self.k, self.cl, self.hco3, self.co3)))


CLASSES = MappingProxyType({c.name: c for c in (
    MechanismClass('C0', -.5, -.5, 1., -2., 0.),
    MechanismClass('C1', 1., -1., 1., -1., 0.),
    MechanismClass('C2', 1., -2., 1., 0., -1.),
    MechanismClass('C3a', 0., -1., 1., -2., 0.),
    MechanismClass('C3b', -1., 0., 1., -2., 0.),
    MechanismClass('C4a', 0., -1., 1., 0., -1.),
    MechanismClass('C4b', -1., 0., 1., 0., -1.),
)})


@dataclass(frozen=True)
class ClassContribution:
    mechanism: MechanismClass
    j4_fmol_s: float
    diagnostics: dict

    @property
    def na_cell_fmol_s(self): return self.mechanism.na * self.j4_fmol_s
    @property
    def k_cell_fmol_s(self): return self.mechanism.k * self.j4_fmol_s
    @property
    def cl_cell_fmol_s(self): return self.mechanism.cl * self.j4_fmol_s
    @property
    def hco3_cell_fmol_s(self): return self.mechanism.hco3 * self.j4_fmol_s
    @property
    def co3_cell_fmol_s(self): return self.mechanism.co3 * self.j4_fmol_s
    @property
    def tic_cell_fmol_s(self): return self.mechanism.conserved_coefficients[3] * self.j4_fmol_s
    @property
    def alkalinity_cell_fmol_s(self): return self.mechanism.conserved_coefficients[4] * self.j4_fmol_s
    @property
    def transported_charge_fmol_s(self): return float(self.mechanism.exact_charge) * self.j4_fmol_s
    @property
    def charge_source_residual_fmol_s(self):
        return math.fsum((self.na_cell_fmol_s, self.k_cell_fmol_s,
                          -self.cl_cell_fmol_s, -self.alkalinity_cell_fmol_s))


class CatalanSourceModel(ModernFullModel):
    def __init__(self, *args, mechanism_class, **kwargs):
        self.mechanism = CLASSES[mechanism_class]
        if kwargs.get('ae4_evaluator') is not equal_routing_adapter:
            raise ValueError('Task 42 requires the unchanged Task 40 scalar cycle')
        super().__init__(*args, **kwargs)

    def _evaluate_ae4(self, observables, *, regulatory_multiplier, expression_scale):
        inherited = super()._evaluate_ae4(observables,
            regulatory_multiplier=regulatory_multiplier, expression_scale=expression_scale)
        return ClassContribution(self.mechanism, inherited.cl_cell_fmol_s, {
            'model_id': 'CATALAN_2025_SOURCE_' + self.mechanism.name,
            'scalar_cycle_law': 'unchanged Task 40 QSS cycle',
            'regulation_gain': regulatory_multiplier,
            'expression_scale': expression_scale,
            'source_coefficients_na_k_cl_tic_ta': self.mechanism.conserved_coefficients,
            'class_specific_kinetics': False,
        })


def with_mechanism_class(template, class_id):
    """Replace the AE4 contribution without copying any whole cell equations."""
    old = template.base_model.base_model
    core = CatalanSourceModel(old.parameters, mechanism_class=class_id,
        stimulus=old.stimulus, regulatory_model=old.regulatory_model,
        ae4_parameters=old.ae4_parameters, ae4_evaluator=old.ae4_evaluator,
        nkcc1_kinetics=old.nkcc1_kinetics)
    return MinimalNbcModel(StimulatedNkcc1Model(core), template.nbc_parameters)


def forward_affinity(mechanism, concentrations_i, concentrations_o, *,
                     voltage_V, thermal_voltage_V):
    """Return A/(RT) = -sum(nu log(ci/co)) - V/Vt sum(nu z).

nu is a signed cellular species source per positive inward chloride cycle.
Concentrations approximate activities. Positive A supports the proposed
forward direction; this diagnostic does not replace the inherited J4 law.
"""
    valences = {'na': 1, 'k': 1, 'cl': -1, 'hco3': -1, 'co3': -2}
    terms = {}
    for ion, nu in mechanism.species_coefficients.items():
        ci, co = concentrations_i[ion], concentrations_o[ion]
        if not (math.isfinite(ci) and math.isfinite(co) and ci > 0 and co > 0):
            raise ValueError('Affinity requires finite positive species concentrations')
        terms[ion] = {'coefficient': nu,
                     'chemical_delta_g_over_rt': nu * math.log(ci / co),
                     'electrical_delta_g_over_rt': nu * valences[ion] * voltage_V / thermal_voltage_V}
    chemical = math.fsum(v['chemical_delta_g_over_rt'] for v in terms.values())
    electrical = math.fsum(v['electrical_delta_g_over_rt'] for v in terms.values())
    affinity = -chemical - electrical
    return {'species_terms': terms, 'chemical_delta_g_over_rt': chemical,
            'electrical_delta_g_over_rt': electrical, 'affinity_over_rt': affinity,
            'forward_direction': ('supported' if affinity > 1e-12 else
                                  'opposed' if affinity < -1e-12 else 'indeterminate')}
