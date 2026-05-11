# -*- coding: utf-8 -*-
"""
Generate a dedicated CMC-vs-DRRC comparison figure.

The figure intentionally separates:
  - CMC local project measurements
  - DRRC literature-reported baseline values

This avoids overwriting the existing multi-baseline chart workflow.
"""

import csv
import io
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

BASE_DIR = Path(__file__).resolve().parent.parent
METRICS_CSV = BASE_DIR / "metrics" / "cmc_vs_drrc_summary.csv"
OUT_DIR = Path(__file__).resolve().parent
OUT_DIR.mkdir(exist_ok=True)


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Microsoft YaHei", "sans-serif"],
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})


CMC_COLOR = "#2563EB"
CMC_LIGHT = "#60A5FA"
DRRC_COLOR = "#F97316"
DRRC_LIGHT = "#FDBA74"
GRID_COLOR = "#D4D4D8"
TEXT_COLOR = "#18181B"


def load_rows():
    with METRICS_CSV.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def get_value(rows, system, item):
    for row in rows:
        if row["system"] == system and row["item"] == item:
            return row["value"]
    raise KeyError(f"Missing row for {system=} {item=}")


def get_density_rows(rows):
    wanted = [
        ("CMC", "Chu Shi Biao (Traditional, A direct)"),
        ("CMC", "Chu Shi Biao (Traditional, B t2s)"),
        ("CMC", "Bei Ying (Simplified)"),
        ("DRRC", "Coding potential"),
        ("DRRC", "Net information density"),
    ]
    labels = []
    values = []
    colors = []
    for system, item in wanted:
        labels.append(item)
        values.append(float(get_value(rows, system, item)))
        colors.append(CMC_COLOR if system == "CMC" else DRRC_COLOR)
    return labels, values, colors


def draw_density_panel(ax, rows):
    labels, values, colors = get_density_rows(rows)
    x = np.arange(len(labels))
    bars = ax.bar(x, values, color=colors, edgecolor="#111827", linewidth=0.8)

    ax.set_title("(a) Storage Density / Coding Potential")
    ax.set_ylabel("bits/nt")
    ax.set_xticks(x)
    ax.set_xticklabels([
        "CMC\nTrad A",
        "CMC\nTrad B",
        "CMC\nSimplified",
        "DRRC\nCoding",
        "DRRC\nNet",
    ])
    ax.grid(axis="y", linestyle=":", linewidth=0.6, color=GRID_COLOR, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(0, max(values) * 1.24)

    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.04,
            f"{value:.3f}" if value < 2 else f"{value:.2f}",
            ha="center",
            va="bottom",
            fontsize=8.5,
            color=TEXT_COLOR,
            fontweight="bold",
        )


def draw_constraint_panel(ax, rows):
    metrics = ["Global GC", "Max homopolymer"]
    cmc_vals = [
        float(get_value(rows, "CMC", "Global GC")),
        float(get_value(rows, "CMC", "Max homopolymer")),
    ]
    drrc_vals = [
        float(get_value(rows, "DRRC", "Global GC")),
        float(get_value(rows, "DRRC", "Max homopolymer")),
    ]

    x = np.arange(len(metrics))
    width = 0.34
    bars1 = ax.bar(x - width / 2, cmc_vals, width, label="CMC", color=CMC_LIGHT, edgecolor="#1D4ED8", linewidth=0.8)
    bars2 = ax.bar(x + width / 2, drrc_vals, width, label="DRRC", color=DRRC_LIGHT, edgecolor="#EA580C", linewidth=0.8)

    ax.set_title("(b) Core Sequence Constraints")
    ax.set_xticks(x)
    ax.set_xticklabels(["Global GC (%)", "Max homopolymer (nt)"])
    ax.grid(axis="y", linestyle=":", linewidth=0.6, color=GRID_COLOR, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="upper right", framealpha=0.9)
    ax.set_ylim(0, 60)

    for bars in (bars1, bars2):
        for bar in bars:
            value = bar.get_height()
            label = f"{value:.1f}" if value >= 10 else f"{int(value)}"
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                value + 1.0,
                label,
                ha="center",
                va="bottom",
                fontsize=8.5,
                color=TEXT_COLOR,
            )


def add_notes(fig):
    note = (
        "CMC bars use local project measurements. DRRC bars use literature-reported values from Bioinformatics (2024), "
        "so this figure is for method positioning rather than strict same-input rerun comparison."
    )
    fig.text(0.5, 0.02, note, ha="center", va="bottom", fontsize=8.5, color="#3F3F46")


def main():
    rows = load_rows()
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6))
    fig.suptitle("CMC vs DRRC: Project Results and Literature Baseline", fontsize=13, fontweight="bold", y=0.98)

    draw_density_panel(axes[0], rows)
    draw_constraint_panel(axes[1], rows)
    add_notes(fig)

    fig.tight_layout(rect=[0, 0.07, 1, 0.93], w_pad=2.4)

    svg_path = OUT_DIR / "fig_cmc_vs_drrc_comparison.svg"
    png_path = OUT_DIR / "fig_cmc_vs_drrc_comparison.png"
    fig.savefig(svg_path, format="svg")
    fig.savefig(png_path, format="png")
    plt.close(fig)

    print(f"[OK] {svg_path}")
    print(f"[OK] {png_path}")


if __name__ == "__main__":
    main()
