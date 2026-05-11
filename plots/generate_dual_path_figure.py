#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate fig_dual_path_preprocessing.png
========================================
Horizontal dual-path diagram (clean text style, no cell borders).
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

# -- Fonts -----------------------------------------------------------------
CJK = "Microsoft YaHei"

plt.rcParams.update({
    "font.family":        "sans-serif",
    "font.sans-serif":    ["Microsoft YaHei", "SimHei", "Arial",
                           "Helvetica", "DejaVu Sans"],
    "font.size":          10,
    "axes.unicode_minus": False,
    "savefig.dpi":        600,
    "savefig.bbox":       "tight",
    "savefig.pad_inches": 0.02,
})

# -- Colours ---------------------------------------------------------------
TXT     = "#202124"
SUB     = "#5F6368"
A_E     = "#1E8E3E"
A_BG    = "#F4FAF6"
A_HI    = "#C8E6C9"
B_E     = "#E37400"
B_BG    = "#FFFBF0"
B_HI    = "#FFE0B2"
BLUE    = "#1967D2"
RED     = "#D93025"
VIOLET  = "#7C4DFF"


# -- Helpers ---------------------------------------------------------------
def rbox(ax, x, y, w, h, fc, ec, lw=1.0, pad=0.006, zorder=2):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle=f"round,pad={pad}",
        fc=fc, ec=ec, lw=lw, zorder=zorder,
        transform=ax.transAxes, clip_on=False))


def arr(ax, x1, y1, x2, y2, color="#3C4043", lw=1.2, ms=13):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle="-|>",
        connectionstyle="arc3,rad=0",
        color=color, lw=lw, mutation_scale=ms, zorder=4,
        transform=ax.transAxes, clip_on=False))


def draw_chars(ax, cx, y, items, gap, fs, hi_idx, hi_bg, text_color=TXT):
    """Draw Chinese characters as plain text with optional soft highlight."""
    n = len(items)
    x0 = cx - (n - 1) * gap / 2
    xs = []
    for i, ch in enumerate(items):
        xi = x0 + i * gap
        xs.append(xi)
        if i in hi_idx:
            rbox(ax, xi - gap * 0.35, y - 0.035, gap * 0.7, 0.07,
                 fc=hi_bg, ec=hi_bg, lw=0, pad=0.002, zorder=3)
        ax.text(xi, y, ch, ha="center", va="center",
                fontsize=fs, color=text_color, fontname=CJK,
                transform=ax.transAxes, zorder=5)
    return xs


def draw_codes(ax, cx, y, items, gap, fs, hi_idx, hi_bg, text_color=TXT):
    """Draw Wubi codes as monospace text with optional soft highlight."""
    n = len(items)
    x0 = cx - (n - 1) * gap / 2
    xs = []
    for i, code in enumerate(items):
        xi = x0 + i * gap
        xs.append(xi)
        if i in hi_idx:
            rbox(ax, xi - gap * 0.38, y - 0.030, gap * 0.76, 0.06,
                 fc=hi_bg, ec=hi_bg, lw=0, pad=0.002, zorder=3)
        ax.text(xi, y, code, ha="center", va="center",
                fontsize=fs, color=text_color, fontfamily="monospace",
                transform=ax.transAxes, zorder=5)
    return xs


# -- Main ------------------------------------------------------------------
def main():
    fig, ax = plt.subplots(figsize=(14, 5.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Vertical positions
    y_title  = 0.955
    y_header = 0.890
    y_a      = 0.695
    y_input  = 0.490
    y_b      = 0.285
    y_leg    = 0.075

    # Column centres
    col_lab    = 0.065      # left labels
    col_char   = 0.21       # character group
    col_wubi   = 0.50       # wubi code group
    col_stream = 0.84       # codestream

    # Item spacing
    ch_gap   = 0.040        # between characters
    cd_gap   = 0.050        # between codes

    diff = {2, 3}

    # ── Data ──────────────────────────────────────────────────────
    trad   = ["先", "帝", "創", "業", "未", "半"]
    simp   = ["先", "帝", "创", "业", "未", "半"]
    cod_a  = ["tfhp", "upmh", "wgbj", "ogud", "fii", "ufj"]
    cod_b  = ["tfhp", "upmh", "wbjh", "ogd",  "fii", "ufj"]
    str_a  = "tfhp6upmh6wgbj6ogud6fii6ufj"
    str_b  = "tfhp6upmh6wbjh6ogd6fii6ufj"

    # =================================================================
    # TITLE
    # =================================================================
    ax.text(0.5, y_title,
            "Traditional Chinese Input: Wubi Preprocessing \u2014 Dual-Path Example",
            ha="center", va="center", fontsize=16, fontweight="bold",
            color=TXT, transform=ax.transAxes)

    # ── Stage headers ────────────────────────────────────────────
    for cx, lab in [(col_char, "Character Processing"),
                    (col_wubi, "Wubi-86 Encoding"),
                    (col_stream, "Codestream")]:
        ax.text(cx, y_header, lab, ha="center", va="center",
                fontsize=11.5, fontweight="bold", color=BLUE,
                transform=ax.transAxes)
        ax.plot([cx - 0.068, cx + 0.068],
                [y_header - 0.018, y_header - 0.018],
                color=BLUE, lw=0.8, alpha=0.35, transform=ax.transAxes)

    # =================================================================
    # INPUT TEXT  (centre row)
    # =================================================================
    ax.text(col_lab + 0.035, y_input + 0.065, "Input Text",
            ha="center", va="center", fontsize=11, fontweight="bold",
            color=BLUE, transform=ax.transAxes)

    xs_in = draw_chars(ax, col_lab + 0.035, y_input,
                       trad, ch_gap, 17, diff, "#DCEAFE", BLUE)

    ax.text(col_lab + 0.035, y_input - 0.060,
            "\u7e41\u4f53\u4e2d\u6587 (Traditional)",
            ha="center", va="center", fontsize=8.5, color=BLUE,
            transform=ax.transAxes, fontname=CJK)

    # =================================================================
    # PATH A  ─ green (Direct)
    # =================================================================
    rbox(ax, 0.005, y_a - 0.08, 0.96, 0.16,
         fc=A_BG, ec=A_E, lw=1.1, pad=0.006, zorder=1)

    # Label cluster (far left)
    lx_a = 0.035
    ax.text(lx_a, y_a + 0.025, "Path A",
            ha="center", va="center", fontsize=12, fontweight="bold",
            color=A_E, transform=ax.transAxes)
    ax.text(lx_a, y_a - 0.008, "Direct",
            ha="center", va="center", fontsize=8.5,
            color=A_E, transform=ax.transAxes)
    ax.text(lx_a, y_a - 0.042,
            "\u2605 Primary", ha="center", va="center",
            fontsize=7.5, fontweight="bold", color="white",
            transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.20", fc=A_E, ec=A_E, lw=0.6))

    # Characters
    xs_ac = draw_chars(ax, col_char, y_a, trad, ch_gap, 18, diff, A_HI)
    ax.text(col_char, y_a + 0.065, "trad_mapping.csv lookup",
            ha="center", va="center", fontsize=8.5, color=A_E,
            transform=ax.transAxes)
    ax.text(col_char, y_a - 0.060, "(characters unchanged)",
            ha="center", va="center", fontsize=7.5, color=SUB,
            transform=ax.transAxes)

    # Wubi codes
    xs_aw = draw_codes(ax, col_wubi, y_a, cod_a, cd_gap, 11.5, diff, A_HI)

    # Single group arrow:  characters --> codes
    arr(ax, xs_ac[-1] + 0.022, y_a, xs_aw[0] - 0.025, y_a,
        color=A_E, lw=1.0, ms=12)

    # Single group arrow:  codes --> codestream
    arr(ax, xs_aw[-1] + 0.028, y_a, col_stream - 0.105, y_a,
        color=A_E, lw=1.0, ms=12)

    # Codestream box
    ax.text(col_stream, y_a, str_a,
            ha="center", va="center", fontsize=8, fontweight="bold",
            color=A_E, fontfamily="monospace", zorder=8,
            transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.30", fc=A_BG, ec=A_E, lw=1.0))
    ax.text(col_stream, y_a - 0.058, 'delimiter: "6"',
            ha="center", va="center", fontsize=7.5, color=VIOLET,
            fontweight="bold", transform=ax.transAxes)

    # =================================================================
    # PATH B  ─ orange (T->S)
    # =================================================================
    rbox(ax, 0.005, y_b - 0.08, 0.96, 0.16,
         fc=B_BG, ec=B_E, lw=1.1, pad=0.006, zorder=1)

    lx_b = 0.035
    ax.text(lx_b, y_b + 0.025, "Path B",
            ha="center", va="center", fontsize=12, fontweight="bold",
            color=B_E, transform=ax.transAxes)
    ax.text(lx_b, y_b - 0.008, "T\u2192S",
            ha="center", va="center", fontsize=8.5,
            color=B_E, transform=ax.transAxes)
    ax.text(lx_b, y_b - 0.042, "Alternative",
            ha="center", va="center", fontsize=7.2, fontweight="bold",
            color=B_E, transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.20", fc="#FFF3E0",
                      ec=B_E, lw=0.6))

    # Characters
    xs_bc = draw_chars(ax, col_char, y_b, simp, ch_gap, 18, diff, B_HI)
    ax.text(col_char, y_b + 0.065, "OpenCC t2s conversion",
            ha="center", va="center", fontsize=8.5, color=B_E,
            transform=ax.transAxes)
    ax.text(col_char, y_b - 0.060,
            "\u5275\u2192\u521b  \u696d\u2192\u4e1a",
            ha="center", va="center", fontsize=7.5, color=B_E,
            transform=ax.transAxes, fontname=CJK)

    # Wubi codes
    xs_bw = draw_codes(ax, col_wubi, y_b, cod_b, cd_gap, 11.5, diff, B_HI)

    # Group arrows
    arr(ax, xs_bc[-1] + 0.022, y_b, xs_bw[0] - 0.025, y_b,
        color=B_E, lw=1.0, ms=12)
    arr(ax, xs_bw[-1] + 0.028, y_b, col_stream - 0.105, y_b,
        color=B_E, lw=1.0, ms=12)

    # Codestream box
    ax.text(col_stream, y_b, str_b,
            ha="center", va="center", fontsize=8, fontweight="bold",
            color=B_E, fontfamily="monospace", zorder=8,
            transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.30", fc=B_BG, ec=B_E, lw=1.0))
    ax.text(col_stream, y_b - 0.058, 'delimiter: "6"',
            ha="center", va="center", fontsize=7.5, color=VIOLET,
            fontweight="bold", transform=ax.transAxes)

    # =================================================================
    # Branch arrows:  Input --> Path A / Path B
    # =================================================================
    bx = 0.015
    arr(ax, bx, y_input + 0.035, bx, y_a - 0.065,
        color=A_E, lw=1.5, ms=13)
    arr(ax, bx, y_input - 0.035, bx, y_b + 0.065,
        color=B_E, lw=1.5, ms=13)

    # =================================================================
    # Vertical diff-guides + not-equal signs
    # =================================================================
    x_d1 = xs_aw[2]
    x_d2 = xs_aw[3]
    for xp in (x_d1, x_d2):
        ax.plot([xp, xp], [y_b - 0.065, y_a + 0.065],
                color=RED, lw=0.8, ls=(0, (2, 2)), alpha=0.35,
                transform=ax.transAxes)
        ax.text(xp, y_input, "\u2260",
                color=RED, fontsize=11, fontweight="bold",
                ha="center", va="center", transform=ax.transAxes,
                bbox=dict(boxstyle="circle,pad=0.12",
                          fc="#FFF7F7", ec="#EA4335", lw=0.9))

    # =================================================================
    # Legend
    # =================================================================
    rbox(ax, 0.17, y_leg - 0.013, 0.016, 0.026,
         fc=A_HI, ec=A_E, lw=0.9, pad=0.0015)
    ax.text(0.190, y_leg,
            "Path A: Traditional character (via trad_mapping.csv)",
            ha="left", va="center", fontsize=8.3, color=TXT,
            transform=ax.transAxes)

    rbox(ax, 0.56, y_leg - 0.013, 0.016, 0.026,
         fc=B_HI, ec=B_E, lw=0.9, pad=0.0015)
    ax.text(0.580, y_leg,
            "Path B: Simplified character (via OpenCC t2s)",
            ha="left", va="center", fontsize=8.3, color=TXT,
            transform=ax.transAxes)

    # =================================================================
    # Save
    # =================================================================
    out = Path(__file__).resolve().parent
    fig.savefig(out / "fig_dual_path_preprocessing.png")
    fig.savefig(out / "fig_dual_path_preprocessing.svg")
    fig.savefig(out / "fig_dual_path_preprocessing.pdf")
    plt.close(fig)
    print("Saved fig_dual_path_preprocessing (.png / .svg / .pdf)")


if __name__ == "__main__":
    main()
