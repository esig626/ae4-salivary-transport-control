"""Fixed beta conditioned effective coupling: TARGET-CALIBRATED CONSTRUCTION.

The parameter is inherited from Task 41's displayed phenotype calibration.
It is not a molecular AE4 to TMEM16A interaction or an independent estimate.
Only apical chloride conductance is modified, before inherited current closure.
"""
from dataclasses import replace
import math

from .model import WT
from .nbc_minimal import normalized_secretory_activation
from .task41_selected import Task41SelectedModel


FIXED_LAMBDA = 0.89488127156712


def effective_coupling_factor(activation, beta, ae4_expression):
    """Return the fixed law; protocol and expression are existing coordinates."""
    if not all(math.isfinite(x) for x in (activation, beta, ae4_expression)):
        raise ValueError("Coupling inputs must be finite")
    if not (0 <= activation <= 1 and 0 <= beta <= 1 and ae4_expression >= 0):
        raise ValueError("Coupling inputs outside their parent domains")
    if activation == 0 or beta == 0 or ae4_expression == 1:
        return 1.0
    return 1.0 - FIXED_LAMBDA * activation * beta * (1.0 - ae4_expression)


class Task50EffectiveCouplingModel(Task41SelectedModel):
    """Reuse the frozen conductance closure, RHS and integration dispatch.

    Supply the unchanged Task 40 MinimalNbcModel template. Initial states and
    protocols remain the caller's explicit inherited choices. No new parameter
    object, recruitment coefficient or dynamic state is created.
    """

    def __init__(self, task40_model):
        # The inherited effective_model/RHS/integration methods need only this
        # template. Task 41's separate recruitment parameter is intentionally
        # absent: evaluate below uses exactly the single fixed Task 50 lambda.
        self.base_model = task40_model

    def evaluate(self, time_s, vector, *, genotype=WT):
        stimulus = self.stimulus(time_s)
        activation = normalized_secretory_activation(
            stimulus.calcium_uM, self.nbc_parameters
        )
        factor = effective_coupling_factor(
            activation, stimulus.beta_input, genotype.ae4_expression
        )
        model = self.effective_model(factor)
        evaluation = model.evaluate(time_s, vector, genotype=genotype)
        regulatory = dict(evaluation.diagnostics.regulatory)
        regulatory.update(
            task50_calcium_activation=float(activation),
            task50_beta_input=float(stimulus.beta_input),
            task50_effective_coupling_factor=float(factor),
            task50_fixed_lambda=FIXED_LAMBDA,
            task50_effective_cl_conductance_S=model.parameters.membranes.g_cl_apical_S,
            task50_classification="TARGET-CALIBRATED CONSTRUCTION",
        )
        return replace(
            evaluation,
            diagnostics=replace(evaluation.diagnostics, regulatory=regulatory),
        )
