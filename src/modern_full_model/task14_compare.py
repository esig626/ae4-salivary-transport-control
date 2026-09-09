"""Comparison utilities. Inputs are saved predictions and explicit target records.

This module never imports the model, runs a solver, or writes a blind artifact.
The target records are provided only after the remote blind checkpoint exists.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from pathlib import Path

FINAL_NAMES = frozenset({"reveal_manifest.json", "ae4_holdout_comparison.csv",
    "ae2_validation_comparison.csv", "final_classification.json", "final_artifact_hashes.csv"})


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path):
    with Path(path).open("r", newline="") as handle:
        return list(csv.DictReader(handle))


def number(value):
    if value is None or value == "":
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("Comparison inputs must be finite")
    return result


def compare_prediction(row, target):
    """Use the target's stated resolution and uncertainty, without fitting."""
    out = dict(row)
    prediction = number(row.get("D_total"))
    observed = number(target.get("D_total"))
    interval = target.get("D_total_interval")
    out.update(target_D_total=observed, target_interval_json=json.dumps(interval),
        target_evidence_status=target.get("evidence_status"),
        comparison_status="UNAVAILABLE" if prediction is None or observed is None else "EVALUATED",
        signed_error_D_total=None, absolute_error_D_total=None, inside_reported_interval=None,
        distance_to_reported_interval=None, time_course_RMSE=None, time_points_compared=0,
        early_ratio_error=None, sustained_ratio_error=None, delta_emergence_direction_agrees=None)
    if prediction is not None and observed is not None:
        out["signed_error_D_total"] = prediction - observed
        out["absolute_error_D_total"] = abs(prediction - observed)
        if interval is not None:
            low, high = map(float, interval)
            if low > high:
                raise ValueError("Reversed experimental interval")
            out["inside_reported_interval"] = low <= prediction <= high
            out["distance_to_reported_interval"] = max(low - prediction, prediction - high, 0.0)
    residuals = []
    for point in target.get("relative_time_points", []):
        field = point["prediction_field"]
        value = number(row.get(field))
        if value is None:
            continue
        residual = value - float(point["ratio"])
        residuals.append(residual)
        out["residual_" + field] = residual
        if point.get("interval") is not None:
            low, high = map(float, point["interval"])
            out["inside_interval_" + field] = low <= value <= high
            out["interval_distance_" + field] = max(low - value, value - high, 0.0)
    if residuals:
        out["time_course_RMSE"] = math.sqrt(sum(e * e for e in residuals) / len(residuals))
        out["time_points_compared"] = len(residuals)
    for predicted, observed_key, output_key in (("R_early", "R_early", "early_ratio_error"),
            ("R_late", "R_late", "sustained_ratio_error")):
        if number(row.get(predicted)) is not None and number(target.get(observed_key)) is not None:
            out[output_key] = number(row[predicted]) - number(target[observed_key])
    if number(row.get("delta_emergence")) is not None and number(target.get("delta_emergence")) is not None:
        sign = lambda x: (x > 0) - (x < 0)
        out["delta_emergence_direction_agrees"] = sign(number(row["delta_emergence"])) == sign(number(target["delta_emergence"]))
    return out


def evaluate_rows(rows, target):
    """Keep every predeclared row, including unavailable predictions."""
    return [compare_prediction(row, target) for row in rows]


def render_csv(rows):
    columns = sorted({key for row in rows for key in row})
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode()


def write_final(results_directory, name, data):
    """The writer has an exact allowlist and rejects symlinks and traversal."""
    if name not in FINAL_NAMES:
        raise ValueError("Writing a blind prediction or unlisted file is prohibited")
    directory = Path(results_directory).resolve()
    path = directory / name
    if path.is_symlink() or path.resolve().parent != directory:
        raise ValueError("Comparison output must be a direct nonsymlink file")
    path.write_bytes(data)


def verify_blind_ledger(repository, ledger):
    repository = Path(repository).resolve()
    for row in read_csv(ledger):
        path = (repository / row["path"]).resolve()
        if not path.is_relative_to(repository):
            raise ValueError("Blind path is outside the repository")
        if file_hash(path) != row["sha256"]:
            raise ValueError("Blind artifact changed: " + row["path"])


def comparison_summary(rows):
    """Descriptive summaries only, with no best case selection."""
    groups = {}
    for row in rows:
        key = (str(row["calcium_uM"]), row["routing_family"])
        groups.setdefault(key, []).append(row)
    result = []
    for (calcium, family), group in sorted(groups.items()):
        record = {"calcium_uM": float(calcium), "routing_family": family, "row_count": len(group)}
        for metric in ("R_total", "D_total", "absolute_error_D_total", "time_course_RMSE"):
            values = sorted(number(r.get(metric)) for r in group if number(r.get(metric)) is not None)
            record[metric + "_count"] = len(values)
            record[metric + "_min"] = min(values) if values else None
            record[metric + "_max"] = max(values) if values else None
            n = len(values)
            record[metric + "_median"] = (values[(n-1)//2] + values[n//2]) / 2 if n else None
        result.append(record)
    return result
