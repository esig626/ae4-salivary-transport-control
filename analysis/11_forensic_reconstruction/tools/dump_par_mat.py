#!/usr/bin/env python3
"""Reproduce and verify the immutable historical Par.mat inventory.

This is forensic tooling, not a scientific model implementation.  It reads the
archived files without modifying them and writes only when --output is given.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any

import scipy
from scipy.io import loadmat, whosmat


ROOT = Path(__file__).resolve().parents[3]
MAT_PATH = ROOT / "archive/legacy-2017/Ae4_Dynamics_Project/Par.mat"
SCRIPT_PATH = ROOT / "archive/legacy-2017/Ae4_Dynamics_Project/Par.m"
EXPECTED_PATH = ROOT / "results/11_forensic_reconstruction/par_mat_inventory.json"


def _eval_par_expression(node: ast.AST, values: dict[str, float]) -> float:
    """Evaluate the numeric expression subset used by archived Par.m."""

    if isinstance(node, ast.Expression):
        return _eval_par_expression(node.body, values)
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
        return float(node.value)
    if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name):
        if node.value.id != "par" or node.attr not in values:
            raise ValueError(f"unknown parameter reference in Par.m: {ast.unparse(node)}")
        return values[node.attr]
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
        value = _eval_par_expression(node.operand, values)
        return value if isinstance(node.op, ast.UAdd) else -value
    if isinstance(node, ast.BinOp):
        left = _eval_par_expression(node.left, values)
        right = _eval_par_expression(node.right, values)
        if isinstance(node.op, ast.Add):
            return left + right
        if isinstance(node.op, ast.Sub):
            return left - right
        if isinstance(node.op, ast.Mult):
            return left * right
        if isinstance(node.op, ast.Div):
            return left / right
        if isinstance(node.op, ast.Pow):
            return left**right
    raise ValueError(f"unsupported expression in Par.m: {ast.dump(node)}")


def read_par_script() -> tuple[dict[str, float], int]:
    """Evaluate archived scalar assignments, preserving MATLAB script order."""

    assignment = re.compile(r"^\s*par\.([A-Za-z][A-Za-z0-9_]*)\s*=\s*(.*?)\s*;\s*$")
    values: dict[str, float] = {}
    statement_count = 0
    for line_number, raw_line in enumerate(SCRIPT_PATH.read_text().splitlines(), 1):
        source = raw_line.split("%", 1)[0]
        match = assignment.match(source)
        if not match:
            continue
        name, expression = match.groups()
        parsed = ast.parse(expression, mode="eval")
        values[name] = _eval_par_expression(parsed, values)
        statement_count += 1
    if statement_count != 46 or len(values) != 45:
        raise ValueError(
            f"unexpected Par.m assignment inventory: {statement_count} statements, "
            f"{len(values)} distinct fields"
        )
    return values, statement_count


def build_inventory() -> dict[str, Any]:
    """Return the deterministic machine-readable inventory."""

    raw = MAT_PATH.read_bytes()
    loaded = loadmat(
        MAT_PATH,
        struct_as_record=False,
        squeeze_me=False,
        mat_dtype=True,
    )
    stored = whosmat(MAT_PATH)
    if stored != [("par", (1, 1), "struct")]:
        raise ValueError(f"unexpected top-level MAT inventory: {stored!r}")

    par_array = loaded["par"]
    par = par_array[0, 0]
    field_names = list(par._fieldnames)
    script_values, assignment_count = read_par_script()
    mat_values: dict[str, float] = {}
    fields: dict[str, dict[str, Any]] = {}

    for name in field_names:
        value = getattr(par, name)
        if value.shape != (1, 1) or str(value.dtype) != "float64":
            raise ValueError(
                f"unexpected Par.mat field representation for {name}: "
                f"shape={value.shape}, dtype={value.dtype}"
            )
        scalar = float(value[0, 0])
        mat_values[name] = scalar
        fields[name] = {
            "shape": [1, 1],
            "matlab_class": "double",
            "numpy_dtype": "float64",
            "value": scalar,
        }

    missing = sorted(set(script_values) - set(mat_values))
    extra = sorted(set(mat_values) - set(script_values))
    mismatches = [
        {
            "field": name,
            "par_m_value": script_values[name],
            "par_mat_value": mat_values[name],
        }
        for name in script_values.keys() & mat_values.keys()
        if script_values[name] != mat_values[name]
    ]

    header_bytes = loaded["__header__"]
    header = header_bytes.decode("ascii") if isinstance(header_bytes, bytes) else str(header_bytes)
    header_match = re.fullmatch(
        r"MATLAB 5\.0 MAT-file, Platform: ([^,]+), Created on: (.+)", header
    )
    if not header_match:
        raise ValueError(f"unexpected MAT header: {header!r}")
    platform, created = header_match.groups()

    return {
        "schema_version": 1,
        "source_path": str(MAT_PATH.relative_to(ROOT)),
        "sha256": hashlib.sha256(raw).hexdigest(),
        "mat_file": {
            "format": "MATLAB 5.0 MAT-file",
            "platform": platform,
            "created": created,
            "version": loaded["__version__"],
            "global_variables": loaded["__globals__"],
            "stored_variables": [
                {
                    "name": "par",
                    "shape": list(par_array.shape),
                    "matlab_class": "struct",
                    "field_count": len(field_names),
                }
            ],
        },
        "reader": {
            "library": "scipy.io.loadmat",
            "scipy_version": scipy.__version__,
            "options": {
                "struct_as_record": False,
                "squeeze_me": False,
                "mat_dtype": True,
            },
        },
        "fields": fields,
        "comparison_to_par_m": {
            "source_path": str(SCRIPT_PATH.relative_to(ROOT)),
            "assignment_statement_count": assignment_count,
            "field_count": len(script_values),
            "missing_from_mat": missing,
            "extra_in_mat": extra,
            "numeric_mismatches": mismatches,
            "status": (
                "exact numerical match after evaluating 46 assignment statements "
                "(45 distinct fields) in script order"
                if not missing and not extra and not mismatches
                else "mismatch"
            ),
        },
        "notes": [
            "All 45 struct fields are MATLAB double 1x1 scalars; integer-looking values remain double when mat_dtype=True.",
            "The only stored top-level workspace variable is par. __header__, __version__, and __globals__ are scipy metadata, not additional stored scientific variables.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="compare the live archive-derived inventory with the checked-in JSON",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="write canonical JSON to this path instead of standard output",
    )
    args = parser.parse_args()

    inventory = build_inventory()
    rendered = json.dumps(inventory, indent=2, allow_nan=False) + "\n"

    if args.check:
        expected = json.loads(EXPECTED_PATH.read_text())
        if inventory != expected:
            print("Par.mat inventory differs from checked-in JSON", file=sys.stderr)
            return 1
        print(
            "Par.mat inventory verified: 1 struct, 45 double scalar fields, "
            "exact Par.m value match"
        )

    if args.output:
        args.output.write_text(rendered)
    elif not args.check:
        sys.stdout.write(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
