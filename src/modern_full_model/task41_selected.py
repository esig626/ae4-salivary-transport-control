"""Public driver for the frozen Task 41 recruitment hypothesis.

Pass the Task 40 MinimalNbcModel template and the saved shared WT rest to
solve_dynamics. The numerical RHS is the exact engine used in candidate 08.
This explicit driver prevents delegation to the parent model's integrator,
which would omit the added recruitment dependency.
"""
import numpy as np
from scipy.integrate import solve_ivp
from .model import WT
from .ae4_cacc_recruitment import Ae4DependentCaccModel,Ae4CaccRecruitmentParameters

SELECTED_NULL_CACC_RECRUITMENT=0.10511872843288446


class Task41SelectedModel(Ae4DependentCaccModel):
    def __init__(self,task40_model):
        super().__init__(task40_model,Ae4CaccRecruitmentParameters(SELECTED_NULL_CACC_RECRUITMENT))

    def solve_dynamics(self,time_span_s,*,initial_state,genotype=WT,method='Radau',
        rtol=1e-7,atol=1e-10,t_eval=None,max_step_s=2.):
        """Integrate this model's RHS; require an explicit frozen initial state.

        The recorded candidate-08 runner additionally applies physiological
        and conservation monitoring to accepted endpoints and dense samples.
        This convenience method does not substitute for that gated runner.
        """
        if method not in ('Radau','BDF'):raise ValueError('Use Radau or BDF')
        y0=self.layout.validate(initial_state)
        return solve_ivp(lambda t,y:self.rhs(t,y,genotype=genotype),time_span_s,y0,
            method=method,rtol=rtol,atol=atol,t_eval=None if t_eval is None else np.asarray(t_eval),max_step=max_step_s)
