"""Generate the checked machine outputs and structural-geometry figure."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ae4_identifiability import (
    BaselineState,
    activity_signature_matrix,
    complete_baseline_audit,
    flux_aggregate_matrix,
    flux_aggregates,
)


def write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def generate_geometry_csv(path: Path, grid_size: int = 41) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "ae2_normalized_cycle_rate",
                "ae4_normalized_cycle_rate",
                "A_shared_chloride_cycles",
                "B_bicarbonate_cycles",
                "ae4_fraction_when_A_positive",
            ]
        )
        for j2 in np.linspace(0.0, 1.0, grid_size):
            for j4 in np.linspace(0.0, 1.0, grid_size):
                a, b = flux_aggregates(float(j2), float(j4))
                fraction = j4 / a if a > 0 else ""
                writer.writerow([j2, j4, a, b, fraction])


def generate_rank_csv(path: Path) -> None:
    state = BaselineState()
    rows = [
        ("A shared chloride-cycle aggregate", 1, "exact scalar projection"),
        ("B bicarbonate-cycle aggregate", 1, "exact scalar projection"),
        ("(A,B) aggregate pair", int(np.linalg.matrix_rank(flux_aggregate_matrix())), "exact"),
        (
            "full intracellular state forcing",
            int(np.linalg.matrix_rank(activity_signature_matrix(state.na_i, state.k_i))),
            "exact where both turnover factors are nonzero",
        ),
        ("Q alone", "at most 1", "dimension bound only; physiological Jacobian not reconstructed"),
        ("Q plus one intracellular ion", "not available", "not certified: full published Jacobian unavailable"),
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["panel", "certified_rank_or_bound", "scope"])
        writer.writerows(rows)


def generate_figure(path: Path) -> None:
    values = np.linspace(0.0, 1.0, 101)
    j2, j4 = np.meshgrid(values, values)
    a = j2 + j4
    b = j2 + 2.0 * j4

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4), constrained_layout=True)
    contour = axes[0].contour(
        j2, j4, a, levels=np.arange(0.25, 2.0, 0.25), cmap="viridis"
    )
    axes[0].clabel(contour, inline=True, fontsize=8)
    axes[0].set(
        xlabel="normalized AE2 cycle rate",
        ylabel="normalized AE4 cycle rate",
        title="One scalar aggregate A: equivalence lines",
        xlim=(0, 1),
        ylim=(0, 1),
    )
    axes[0].set_aspect("equal")

    for index, fixed in enumerate(np.linspace(0.0, 1.0, 9)):
        a_line = fixed + values
        b_line = fixed + 2.0 * values
        axes[1].plot(
            a_line,
            b_line,
            color="#0072B2",
            alpha=0.55,
            linewidth=0.9,
            label="fixed AE2 rate" if index == 0 else None,
        )
        a_line = values + fixed
        b_line = values + 2.0 * fixed
        axes[1].plot(
            a_line,
            b_line,
            color="#D55E00",
            alpha=0.55,
            linewidth=0.9,
            label="fixed AE4 rate" if index == 0 else None,
        )
    axes[1].plot([0, 2], [0, 2], color="black", linewidth=1.4, label="pure AE2: B=A")
    axes[1].plot([0, 1.5], [0, 3], color="black", linestyle="--", linewidth=1.4, label="pure AE4: B=2A")
    axes[1].set(
        xlabel="A = J2 + J4",
        ylabel="B = J2 + 2 J4",
        title="Two aggregates: invertible wedge geometry",
        xlim=(0, 2),
        ylim=(0, 3),
    )
    axes[1].legend(frameon=False, fontsize=8, loc="upper left")
    fig.suptitle(
        "Certified stoichiometric geometry (flux coordinates; not a Q/state simulation)",
        fontsize=11,
    )
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    result_dir = root / "results" / "10_identifiability_discrimination"
    figure_dir = root / "figures"
    result_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    audit_path = result_dir / "baseline_printed_equation_audit.json"
    geometry_path = result_dir / "stoichiometric_geometry.csv"
    rank_path = result_dir / "panel_rank_summary.csv"
    figure_path = figure_dir / "identifiability_stoichiometric_geometry.png"

    write_json(audit_path, complete_baseline_audit())
    generate_geometry_csv(geometry_path)
    generate_rank_csv(rank_path)
    generate_figure(figure_path)

    checksums = {}
    for path in (audit_path, geometry_path, rank_path, figure_path):
        checksums[str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
    write_json(
        result_dir / "run_metadata.json",
        {
            "command": "python src/run_identifiability_audit.py",
            "date": "2026-08-27",
            "deterministic": True,
            "scope": "structural identities and printed-equation diagnostics only",
            "activity_domain": "normalized knockout-to-WT cycle rates [0,1]^2",
            "checksums_sha256": checksums,
        },
    )


if __name__ == "__main__":
    main()
