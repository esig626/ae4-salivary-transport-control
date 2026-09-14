#!/usr/bin/env python3
"""Build the manuscript figures from the frozen result values.

Run from the repository root:

    python manuscript/scripts/build_figures.py

The script writes PDF figures for LaTeX and PNG previews. It does not rerun the
model and therefore cannot change any scientific result.
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "manuscript" / "data" / "figure_source_values.csv"
OUT = ROOT / "manuscript" / "figures"
OUT.mkdir(parents=True, exist_ok=True)


def rows() -> List[dict[str, str]]:
    with DATA.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def save(fig: plt.Figure, stem: str) -> None:
    fig.tight_layout()
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def figure_1() -> None:
    fig, ax = plt.subplots(figsize=(12, 6.3))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 7)
    ax.axis("off")

    def box(x: float, y: float, w: float, h: float, text: str, *, lw: float = 1.2) -> None:
        patch = FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.08", linewidth=lw, fill=False
        )
        ax.add_patch(patch)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=9)

    def arrow(x1: float, y1: float, x2: float, y2: float, text: str = "") -> None:
        patch = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="->", mutation_scale=12)
        ax.add_patch(patch)
        if text:
            ax.text((x1 + x2) / 2, (y1 + y2) / 2 + 0.14, text, ha="center", fontsize=8)

    ax.text(3, 6.6, "2018 model", ha="center", fontsize=13, fontweight="bold")
    ax.text(9, 6.6, "Reconstructed model", ha="center", fontsize=13, fontweight="bold")

    box(0.3, 2.0, 1.4, 2.4, "Bath")
    box(2.1, 1.3, 2.1, 3.8, "Cell\nNa, K, Cl,\nHCO$_3$, H, CO$_2$")
    box(4.7, 2.0, 1.4, 2.4, "Lumen")
    arrow(1.7, 3.9, 2.1, 3.9, "NKCC1")
    arrow(2.1, 3.3, 1.7, 3.3, "AE4")
    arrow(1.7, 2.8, 2.1, 2.8, "AE2/NHE1")
    arrow(4.2, 3.9, 4.7, 3.9, "CaCC")
    arrow(2.8, 1.3, 2.8, 0.55, "water")
    ax.text(3.15, 5.45, "CO$_2$ hydration and an\neffective intracellular buffer", ha="center", fontsize=8)
    ax.text(3.05, 0.18, "No NBC; pH represented through H/HCO$_3$", ha="center", fontsize=8)

    box(6.25, 2.0, 1.25, 2.4, "Bath")
    box(8.0, 1.05, 2.0, 4.3, "Cell\namounts: Na, K, Cl,\nTIC, TA, volume")
    box(10.55, 2.0, 1.2, 2.4, "Lumen")
    arrow(7.5, 4.2, 8.0, 4.2, "NKCC1")
    arrow(8.0, 3.55, 7.5, 3.55, "AE4")
    arrow(7.5, 2.9, 8.0, 2.9, "NBC/NHE1/AE2")
    arrow(10.0, 4.15, 10.55, 4.15, "CaCC")
    arrow(8.65, 1.05, 8.65, 0.38, "water")
    box(8.18, 5.65, 1.65, 0.55, "β/cAMP/PKA")
    arrow(9.0, 5.65, 9.0, 5.35, "AE4 gain")
    ax.text(9.0, 0.1, "Finite carbonate chemistry; pH recovered from TIC and TA", ha="center", fontsize=8)
    ax.text(9.0, 6.35, "Source-fixed NKCC1 and NHE1; electrogenic NBC-like input", ha="center", fontsize=8)
    ax.text(6.0, 6.85, "A decade of model reconstruction", ha="center", fontsize=14, fontweight="bold")
    save(fig, "fig01_model_evolution")


def figure_2(data: List[dict[str, str]]) -> None:
    selected = [r for r in data if r["figure"] == "2"]
    names = [r["x"] for r in selected]
    values = [float(r["value"]) for r in selected]
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    bars = ax.bar(names, values)
    ax.set_ylabel("Positive chloride loading, 60-600 s (%)")
    ax.set_ylim(0, 100)
    ax.set_title("Wild-type stimulated chloride loading after NBC reconstruction")
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 2, f"{value:.1f}%", ha="center")
    ax.text(
        1.0,
        48,
        "NBC imports two bicarbonate\nequivalents per inward cycle,\nallowing productive AE4 loading\nwithout loss of pH control.",
        ha="center",
        va="center",
        fontsize=9,
    )
    save(fig, "fig02_wt_chloride_partition")


def figure_3(data: List[dict[str, str]]) -> None:
    selected = [r for r in data if r["figure"] == "3"]
    grouped: Dict[str, list[tuple[float, float]]] = {}
    for row in selected:
        grouped.setdefault(row["series"], []).append((float(row["x"]), float(row["value"])))
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    for label, pairs in grouped.items():
        pairs.sort()
        x = [pair[0] for pair in pairs]
        y = [pair[1] for pair in pairs]
        line_style = "--" if "Task 40" in label else "-"
        marker = "o" if "5%" in label else "s"
        ax.plot(x, y, linestyle=line_style, marker=marker, markersize=3.5, label=label)
    ax.axhspan(30.3, 39.7, alpha=0.12, label="Experimental null range: 35 ± 4.7%")
    ax.set_xlabel("Time after stimulation (s)")
    ax.set_ylabel("Cumulative secretion deficit relative to WT (%)")
    ax.set_xlim(50, 610)
    ax.set_ylim(0, 45)
    ax.set_title("Equal routing gives a small late deficit; target-selected CaCC coupling gives an early deficit")
    ax.legend(fontsize=8, ncol=2)
    save(fig, "fig03_secretion_deficit_time_course")


def figure_4(data: List[dict[str, str]]) -> None:
    selected = [r for r in data if r["figure"] == "4"]
    tasks = ["Task 39", "Task 40", "Task 41"]
    five = [float(next(r["value"] for r in selected if r["series"] == t and r["x"] == "AE4 5%")) for t in tasks]
    null = [float(next(r["value"] for r in selected if r["series"] == t and r["x"] == "AE4 null")) for t in tasks]
    x = np.arange(len(tasks))
    width = 0.36
    fig, ax = plt.subplots(figsize=(7.4, 4.9))
    bars1 = ax.bar(x - width / 2, five, width, label="AE4 5%")
    bars2 = ax.bar(x + width / 2, null, width, label="AE4 null")
    ax.set_xticks(x, tasks)
    ax.set_ylabel("NKCC1 compensation, 60-600 s (%)")
    ax.set_title("The constructive CaCC amendment suppresses compensatory NKCC1 loading")
    ax.legend()
    for bars in (bars1, bars2):
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.8, f"{bar.get_height():.1f}", ha="center", fontsize=8)
    ax.set_ylim(0, 38)
    save(fig, "fig04_nkcc1_compensation")


def figure_5(data: List[dict[str, str]]) -> None:
    selected = [r for r in data if r["figure"] == "5"]
    classes = ["C0", "C1", "C2", "C3a", "C3b", "C4a", "C4b"]
    values = []
    notes = []
    for cls in classes:
        row = next(r for r in selected if r["x"] == cls)
        values.append(np.nan if row["value"] == "" else float(row["value"]))
        notes.append(row["notes"])
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    positions = np.arange(len(classes))
    finite = np.isfinite(values)
    ax.bar(positions[finite], np.asarray(values)[finite])
    ax.axhspan(30.3, 39.7, alpha=0.12, label="Experimental null range: 35 ± 4.7%")
    ax.axhline(0, linewidth=0.8)
    ax.set_xticks(positions, classes)
    ax.set_ylabel("AE4-null cumulative secretion deficit at 600 s (%)")
    ax.set_title("Catalán 2025 source classes do not recover the held-out phenotype")
    for idx, (value, note) in enumerate(zip(values, notes)):
        if np.isfinite(value):
            ax.text(idx, value + (0.9 if value >= 0 else -1.8), f"{value:.1f}%", ha="center", fontsize=8)
        else:
            fail_time = note.split("at ")[-1]
            ax.plot(idx, 1.0, marker="x", markersize=9)
            ax.text(idx, 2.5, f"WT pH fail\n{fail_time}", ha="center", fontsize=8)
    ax.text(1, -15.5, "source supported", ha="center", fontsize=8)
    ax.text(3, -15.5, "source supported", ha="center", fontsize=8)
    ax.text(4, -15.5, "affinity opposed", ha="center", fontsize=8)
    ax.set_ylim(-18, 43)
    ax.legend(loc="upper right", fontsize=8)
    save(fig, "fig05_catalan_classes")


def figure_6(data: List[dict[str, str]]) -> None:
    selected = [r for r in data if r["figure"] == "6"]
    names = [r["x"] for r in selected]
    values = [float(r["value"]) for r in selected]
    errors = [0, 0, 0, 0, 4.7]
    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    bars = ax.bar(names, values, yerr=errors, capsize=5)
    ax.set_ylabel("AE4-null secretion deficit (%)")
    ax.set_ylim(0, 45)
    ax.set_title("What survives from the 2018 mechanism nearly a decade later")
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 1.0, f"{value:.1f}%", ha="center", fontsize=8)
    ax.text(2, 10.5, "largest valid Task 42 loss;\nits forward affinity is opposed", ha="center", fontsize=8)
    ax.text(3, 23.0, "target selected", ha="center", fontsize=8)
    save(fig, "fig06_ten_year_synthesis")


def main() -> None:
    data = rows()
    figure_1()
    figure_2(data)
    figure_3(data)
    figure_4(data)
    figure_5(data)
    figure_6(data)


if __name__ == "__main__":
    main()
