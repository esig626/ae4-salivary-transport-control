"""Generate the machine-readable diagnostics and figures for the audit."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .reduced_model import LinearizedPumpLeakModel, printed_activity_factor_check


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "10_identifiability_discrimination"
FIGURES = ROOT / "figures" / "10_identifiability_discrimination"


def json_ready(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, default=json_ready) + "\n", encoding="utf-8")


def write_panel_metrics(model: LinearizedPumpLeakModel) -> list[dict]:
    rows = model.all_panel_metrics()
    path = RESULTS / "panel_metrics.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def write_parameter_map(model: LinearizedPumpLeakModel, n: int = 41) -> list[dict]:
    rows = []
    for g2 in np.linspace(0.0, 2.0, n):
        for g4 in np.linspace(0.0, 2.0, n):
            theta = np.array([g2, g4])
            values = model.observation_values(theta)
            x = model.steady_state(theta)
            rows.append(
                {
                    "G2": g2,
                    "G4": g4,
                    "min_state_ratio": float(np.min(x)),
                    **values,
                }
            )
    path = RESULTS / "parameter_observations.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def write_equivalence(model: LinearizedPumpLeakModel) -> np.ndarray:
    points = model.q_equivalence_segment()
    q_reference = model.observation_values((1.0, 1.0))["Q_star"]
    with (RESULTS / "q_equivalence_set.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["G2", "G4", "Q_star"])
        for g2, g4 in points:
            writer.writerow([g2, g4, model.observation_values((g2, g4))["Q_star"]])
    assert max(abs(model.observation_values(p)["Q_star"] - q_reference) for p in points) < 1e-12
    return points


def plot_geometry(model: LinearizedPumpLeakModel, parameter_rows: list[dict], equivalence: np.ndarray) -> None:
    g2 = np.array([row["G2"] for row in parameter_rows])
    g4 = np.array([row["G4"] for row in parameter_rows])
    q = np.array([row["Q_star"] for row in parameter_rows])
    na = np.array([row["Na_i"] for row in parameter_rows])

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), constrained_layout=True)
    contour = axes[0].tricontourf(g2, g4, q, levels=18, cmap="viridis")
    axes[0].plot(equivalence[:, 0], equivalence[:, 1], color="white", lw=2.0, label="Q*=baseline")
    axes[0].plot([1], [1], "o", color="red", ms=5, label="anchor")
    axes[0].set(xlabel="G2 activity multiplier", ylabel="G4 activity multiplier", title="Scalar Q* equivalence")
    axes[0].legend(frameon=False, fontsize=8)
    fig.colorbar(contour, ax=axes[0], label="Q* (diagnostic proxy)")

    scatter = axes[1].scatter(q, na, c=g4, s=8, cmap="plasma", alpha=0.75)
    axes[1].set(xlabel="Q* (diagnostic proxy)", ylabel="Na_i / Na_i,0", title="Two-output parameter image")
    fig.colorbar(scatter, ax=axes[1], label="G4 activity multiplier")
    fig.suptitle("Stoichiometric diagnostic only — not a reproduced physiological simulation", fontsize=10)
    fig.savefig(FIGURES / "observation_geometry.png", dpi=180)
    plt.close(fig)


def plot_panels(panel_rows: list[dict]) -> None:
    rows = [row for row in panel_rows if row["n_observations"] == 2]
    labels = [row["panel"].replace("_star", "*") for row in rows]
    values = [row["sigma_min"] if row["rank"] == 2 else 0.0 for row in rows]
    colors = ["#2b8cbe" if row["rank"] == 2 else "#d7301f" for row in rows]
    fig, ax = plt.subplots(figsize=(9.5, 5.3), constrained_layout=True)
    y = np.arange(len(rows))
    ax.barh(y, values, color=colors)
    ax.set_yticks(y, labels=labels, fontsize=7)
    ax.invert_yaxis()
    ax.set(xlabel="smallest singular value (declared diagnostic scaling)", title="Two-observable local rank diagnostic")
    ax.axvline(0, color="black", lw=0.5)
    fig.savefig(FIGURES / "panel_rank_diagnostics.png", dpi=180)
    plt.close(fig)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    model = LinearizedPumpLeakModel()
    theta0 = np.ones(2)
    x0 = model.steady_state(theta0)
    implicit = model.state_sensitivity
    finite_difference = model.finite_difference_sensitivity()
    panel_rows = write_panel_metrics(model)
    parameter_rows = write_parameter_map(model)
    equivalence = write_equivalence(model)

    write_json(
        RESULTS / "baseline_check.json",
        {
            "check_kind": "published-baseline anchoring check; not full-model reproduction",
            "full_published_model_reproduced": False,
            "published_resting_concentrations_mM": dict(zip(("Na_i", "K_i", "Cl_i", "HCO3_i"), model.baseline_concentrations_mM)),
            "activity_anchor": {"G2": 1.0, "G4": 1.0},
            "solved_state_ratios": x0,
            "max_abs_reduced_residual": float(np.max(np.abs(model.residual(x0, theta0)))),
            "printed_activity_factor_check": printed_activity_factor_check(),
            "activity_factor_warning": "Printed AE2 and AE4 products have incompatible units and are not compared or inserted into one balance.",
            "unresolved_full_model_inputs": [
                "AE4 density-to-whole-cell flux scaling and inconsistent units",
                "current-to-molar-flux and membrane area/density conversions",
                "numerical calcium forcing used for the stimulated steady state",
                "cell-volume water-balance sign/orientation conflict",
                "net CO2 flux sign",
                "luminal acid-base closure and DAE solver details",
            ],
        },
    )
    write_json(
        RESULTS / "derivative_verification.json",
        {
            "implicit_state_sensitivity": implicit,
            "central_finite_difference_state_sensitivity": finite_difference,
            "finite_difference_step": 1e-6,
            "max_abs_error": float(np.max(np.abs(implicit - finite_difference))),
            "consistent_row_sign_flip_max_abs_error": model.row_sign_invariance_error((1, 1, 1, -1)),
            "restoring_matrix_eigenvalues": np.linalg.eigvalsh(model.restoring_matrix),
            "stoichiometric_matrix": model.stoichiometric_matrix,
        },
    )
    write_json(
        RESULTS / "domain_summary.json",
        {
            "parameter_domain": {"G2": [0.0, 2.0], "G4": [0.0, 2.0]},
            "grid_points": len(parameter_rows),
            "minimum_state_ratio_over_grid": min(row["min_state_ratio"] for row in parameter_rows),
            "maximum_Q_star_over_grid": max(row["Q_star"] for row in parameter_rows),
            "minimum_Q_star_over_grid": min(row["Q_star"] for row in parameter_rows),
            "q_equivalence_points": len(equivalence),
            "scope": "affine stoichiometric diagnostic; global statements apply only to this diagnostic map",
        },
    )
    plot_geometry(model, parameter_rows, equivalence)
    plot_panels(panel_rows)


if __name__ == "__main__":
    main()
