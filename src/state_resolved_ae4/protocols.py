"""Declared stimulation protocols for fixed-chassis Stage A."""

from __future__ import annotations

import numpy as np


CALCIUM_STEP_TIME = 100.0


def historical_calcium_step(t: float) -> float:
    """Task-11/12 comparison Ca step; its physical time unit is uncertified."""

    return 0.05 if t <= CALCIUM_STEP_TIME else 0.55


def immediate_pka_step(t: float) -> float:
    """Fast-QSS PKA input sensitivity aligned to the Ca discontinuity.

    No source measured a PKA time constant, so this is not relabeled as the
    experimental 2--3 minute delay.  It is the minimal no-free-delay input.
    """

    return 0.0 if t <= CALCIUM_STEP_TIME else 1.0


def no_pka(_t: float) -> float:
    return 0.0


def stage_a_time_grid() -> np.ndarray:
    pre = np.arange(0.0, CALCIUM_STEP_TIME + 0.5, 2.0)
    post = np.arange(CALCIUM_STEP_TIME + 0.25, 200.0 + 0.125, 0.25)
    return np.concatenate((pre, np.array((CALCIUM_STEP_TIME + 1e-7,)), post))
