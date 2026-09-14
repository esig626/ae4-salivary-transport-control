"""The public integration helper must call the amended, not parent, RHS."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis/41_ae4_loss_algebraic_design/candidate_08'))
from candidate8_common import models_and_reference
from modern_full_model.task41_selected import Task41SelectedModel
from modern_full_model.model import AE4_NULL
import numpy as np
import unittest
from unittest.mock import patch


class DriverTest(unittest.TestCase):
    def test_integrator_receives_selected_rhs(self):
        _,base,_,engine,y,_=models_and_reference()
        m=Task41SelectedModel(base)
        with patch('modern_full_model.task41_selected.solve_ivp') as solve:
            m.solve_dynamics((1e-6,600.),initial_state=y,genotype=AE4_NULL)
            fun=solve.call_args.args[0]
            np.testing.assert_array_equal(fun(1e-6,y),engine.rhs(1e-6,y,genotype=AE4_NULL))
            self.assertFalse(np.array_equal(fun(1e-6,y),base.rhs(1e-6,y,genotype=AE4_NULL)))
        # No ODE integration is performed by this test.


if __name__=='__main__':unittest.main(verbosity=2)
