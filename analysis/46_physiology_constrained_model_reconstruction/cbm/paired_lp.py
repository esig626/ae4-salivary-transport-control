"""Task 46 paired WT/KO LP built above the frozen local flux core."""
from __future__ import annotations
import json
import numpy as np
import build_cbm as b
from carbonscope_flux_core.solver import compile_flux_lp, _highspy


def _csr_rows(matrix):
    starts = [0]
    indices = []
    values = []
    for row in matrix:
        nz = np.flatnonzero(row)
        indices.extend(int(i) for i in nz)
        values.extend(float(row[i]) for i in nz)
        starts.append(len(indices))
    return starts, indices, values


def solve_pair(*, wt_J=1.0, nkcc_equality=False, M=100.0):
    wt = compile_flux_lp(b.build('WT', M, direction='maximise'))
    ko = compile_flux_lp(b.build('KO', M, direction='maximise'))
    assert wt.reaction_ids == ko.reaction_ids == tuple(b.IDS)
    n = len(b.IDS)
    matrix = np.zeros((len(wt.balance_matrix) + len(ko.balance_matrix), 2 * n))
    matrix[:len(wt.balance_matrix), :n] = np.asarray(wt.balance_matrix)
    matrix[len(wt.balance_matrix):, n:] = np.asarray(ko.balance_matrix)
    rows = [row.copy() for row in matrix]
    lower = [0.0] * len(rows)
    upper = [0.0] * len(rows)

    row = np.zeros(2 * n)
    row[b.IDS.index('J')] = 1.0
    rows.append(row)
    lower.append(float(wt_J))
    upper.append(float(wt_J))

    if nkcc_equality:
        row = np.zeros(2 * n)
        row[b.IDS.index('N')] = 1.0
        row[n + b.IDS.index('N')] = -1.0
        rows.append(row)
        lower.append(0.0)
        upper.append(0.0)

    costs = np.zeros(2 * n)
    costs[n + b.IDS.index('J')] = 1.0
    lo = list(wt.lower_bounds) + list(ko.lower_bounds)
    hi = list(wt.upper_bounds) + list(ko.upper_bounds)
    starts, indices, values = _csr_rows(np.asarray(rows))

    highspy = _highspy()
    solver = highspy.Highs()
    solver.setOptionValue('output_flag', False)
    solver.setOptionValue('threads', 1)
    solver.setOptionValue('solver', 'simplex')
    solver.setOptionValue('primal_feasibility_tolerance', 1e-9)
    solver.setOptionValue('dual_feasibility_tolerance', 1e-9)
    solver.addCols(2 * n, list(costs), lo, hi, 0, [0] * (2 * n + 1), [], [])
    solver.addRows(len(rows), lower, upper, len(indices), starts, indices, values)
    solver.setMaximize()
    solver.run()
    if solver.getModelStatus() != highspy.HighsModelStatus.kOptimal:
        raise RuntimeError(solver.modelStatusToString(solver.getModelStatus()))
    x = np.asarray(solver.getSolution().col_value, float)
    activity = np.asarray(rows) @ x
    violation = max(
        float(np.max(np.maximum(np.asarray(lower) - activity, 0.0), initial=0.0)),
        float(np.max(np.maximum(activity - np.asarray(upper), 0.0), initial=0.0)),
    )
    return {
        'WT_J': float(x[b.IDS.index('J')]),
        'KO_J': float(x[n + b.IDS.index('J')]),
        'WT_NKCC': float(x[b.IDS.index('N')]),
        'KO_NKCC': float(x[n + b.IDS.index('N')]),
        'nkcc_equality': bool(nkcc_equality),
        'M': float(M),
        'row_violation': violation,
        'WT_validation': b.validate(dict(zip(b.IDS, x[:n])), 'WT', M),
        'KO_validation': b.validate(dict(zip(b.IDS, x[n:])), 'KO', M),
    }


def assay_match_diagnostic():
    # In the inhibited uptake assay, bicarbonate-dependent AE4/AE2 and CaCC
    # are absent from the initial chloride-uptake bookkeeping. The frozen CBM
    # has zero intracellular chloride storage, so its Cl_i row becomes 2*N=0.
    return {
        'assay_conditions_encoded': {'A': 0.0, 'E': 0.0, 'J': 0.0},
        'frozen_Cl_i_identity': '2*N + A + E - J = 0',
        'reduced_identity': '2*N = 0',
        'positive_initial_NKCC_uptake_representable': False,
        'reason': 'initial uptake requires transient intracellular chloride storage, absent from the sustained CBM',
    }


if __name__ == '__main__':
    out = {
        'defensible_pair': solve_pair(wt_J=1.0, nkcc_equality=False),
        'assay_match_diagnostic': assay_match_diagnostic(),
        'exact_normal_stimulation_NKCC_equality_run': False,
        'reason': 'preserved evidence supplies no justified numerical equality/tolerance for normal stimulation',
    }
    print(json.dumps(out, indent=2))
