"""Independent exact audit of the frozen H01 stoichiometric matrix."""
from fractions import Fraction
from pathlib import Path
import csv
import hashlib
import json
import subprocess

ROOT = Path(__file__).resolve().parents[3]
CBM = ROOT / 'analysis/46_physiology_constrained_model_reconstruction/cbm'
ORIGINAL = 'cc8920b6f64f0bb279363efcdcb8ce4f0077394c66487a21704d497c9b9ef76f'
EXPECTED = '561b6f17dd710dfa49ca0dfb611fcae36d1bb30840cd9cb668f2c70c94fb3552'
BASE = '5f791e3907fa6769364c5fc65ccb31881b32c006'

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def reduce_exact(matrix):
    a = [list(row) for row in matrix]
    pivots = []
    active = 0
    for col in range(len(a[0])):
        candidate = next((i for i in range(active, len(a)) if a[i][col]), None)
        if candidate is None:
            continue
        a[active], a[candidate] = a[candidate], a[active]
        scale = a[active][col]
        a[active] = [value / scale for value in a[active]]
        for i in range(len(a)):
            if i != active and a[i][col]:
                multiple = a[i][col]
                a[i] = [x - multiple*y for x, y in zip(a[i], a[active])]
        pivots.append(col)
        active += 1
        if active == len(a):
            break
    return a, pivots

def null_basis(matrix):
    reduced, pivots = reduce_exact(matrix)
    free = [i for i in range(len(matrix[0])) if i not in pivots]
    result = []
    for col in free:
        vector = [Fraction(0) for _ in matrix[0]]
        vector[col] = Fraction(1)
        for row, pivot in enumerate(pivots):
            vector[pivot] = -reduced[row][col]
        result.append(vector)
    return result, pivots, free

def zero_product(matrix, vectors):
    return all(sum(x*y for x, y in zip(row, vector)) == 0
               for row in matrix for vector in vectors)

def sparse(vector, names):
    return {name: str(value) for name, value in zip(names, vector) if value}

builder_hash = digest(CBM / 'build_cbm.py')
assert builder_hash == EXPECTED
head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
assert head == BASE, (head, BASE)
with (CBM / 'stoichiometric_matrix.csv').open() as stream:
    rows = list(csv.reader(stream))
reactions = rows[0][1:]
names = [row[0] for row in rows[1:]]
matrix = [[Fraction(value) for value in row[1:]] for row in rows[1:]]
basis, column_pivots, free_columns = null_basis(matrix)
transpose = [list(col) for col in zip(*matrix)]
left_basis, row_pivots, _ = null_basis(transpose)
assert zero_product(matrix, basis)
assert zero_product(transpose, left_basis)
assert len(reduce_exact(basis)[1]) == len(basis)
rank = len(column_pivots)
independent_names = [names[i] for i in row_pivots]
without_neutral = [row for name, row in zip(names, matrix) if name != 'neutral_outflow']
without_neutral_rank = len(reduce_exact(without_neutral)[1])
identity = json.loads((CBM / 'network_identity.json').read_text())
checks = {
    'published_H00_head_matches': head == BASE,
    'frozen_builder_matches': builder_hash == EXPECTED,
    'full_rank_matches': rank == identity['full_rank'],
    'nullity_matches': len(basis) == identity['nullity'],
    'independent_rows_match': independent_names == identity['independent_rows'],
    'all_right_basis_residuals_exactly_zero': zero_product(matrix, basis),
    'right_basis_linearly_independent': len(reduce_exact(basis)[1]) == len(basis),
    'all_left_basis_residuals_exactly_zero': zero_product(transpose, left_basis),
    'neutral_outflow_increases_rank_by_one': rank - without_neutral_rank == 1,
    'builder_still_matches_after_audit': digest(CBM / 'build_cbm.py') == EXPECTED,
}
assert all(checks.values()), checks
report = {
    'audit': 'H01 independent rational rank and nullspace audit',
    'method': 'CSV parsed as fractions.Fraction; standalone Gauss-Jordan elimination; no builder import or scientific LP',
    'published_base': head,
    'builder_sha256': builder_hash,
    'original_builder_sha256_observed_before_audit': ORIGINAL,
    'metadata_change': 'Orchestrator corrected only the AE4 bound provenance label after the original hash was independently observed. Exact matrix supplied to this audit is unchanged according to orchestrator; this report records its own matrix digest.',
    'matrix_sha256': digest(CBM / 'stoichiometric_matrix.csv'),
    'auditor_script_sha256': digest(Path(__file__)),
    'matrix_shape': [len(matrix), len(reactions)],
    'full_rank': rank,
    'right_nullity': len(basis),
    'left_nullity': len(left_basis),
    'independent_rows': independent_names,
    'pivot_reactions': [reactions[i] for i in column_pivots],
    'free_reactions': [reactions[i] for i in free_columns],
    'rank_without_neutral_outflow': without_neutral_rank,
    'right_nullspace_basis': [sparse(vector, reactions) for vector in basis],
    'left_nullspace_basis': [sparse(vector, names) for vector in left_basis],
    'checks': checks,
    'findings': [
        'Both charge rows are exact linear combinations of compartment Na, K, Cl and TA balance rows.',
        'Neutral outflow is independent of the first ten component balance rows and must be retained.',
        'The real stoichiometric solution space has dimension ten before inequality constraints. This dimension alone makes no claim about feasible flux directions or thermodynamics.',
    ],
    'pass': True,
}
import sys
destination = Path(sys.argv[1]) if len(sys.argv)>1 else Path('/tmp/task46h_rank_review.json')
destination.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({'pass': report['pass'], 'shape': report['matrix_shape'], 'rank': rank, 'nullity': len(basis), 'report': str(destination)}))
