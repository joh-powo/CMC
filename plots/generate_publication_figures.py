#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Publication-Quality Figures for CMC DNA Storage Manuscript
==========================================================
Generates 6 comparison figures following Nature/Science journal standards:
  - Fig 1: Storage Density
  - Fig 2: GC Content & Local Range
  - Fig 3: Maximum Homopolymer Length
  - Fig 4: Undesired Motif Content
  - Fig 5: Minimum Hamming Distance
  - Fig 6: Minimum Free Energy (MFE)

Style: Clean, minimal, high-contrast, Nature-compatible
"""

import sys, io, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import FancyBboxPatch
import numpy as np
from pathlib import Path

# ═══════════════════════════════════════════════════════════════
# Publication Style Configuration (Nature / Science standard)
# ═══════════════════════════════════════════════════════════════

plt.rcParams.update({
    # Typography — use Arial (Nature standard) with fallback
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'Microsoft YaHei', 'SimHei', 'DejaVu Sans'],
    'font.size': 8,
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'axes.unicode_minus': False,

    # Axes style — clean, thin
    'axes.linewidth': 0.6,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'xtick.minor.size': 1.5,
    'ytick.minor.size': 1.5,
    'xtick.direction': 'out',
    'ytick.direction': 'out',

    # Grid — very subtle
    'axes.grid': False,
    'grid.alpha': 0.2,
    'grid.linewidth': 0.4,
    'grid.linestyle': '-',

    # Figure quality
    'figure.dpi': 150,
    'savefig.dpi': 600,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,

    # Legend
    'legend.frameon': True,
    'legend.framealpha': 0.92,
    'legend.edgecolor': '#CCCCCC',
    'legend.fancybox': False,
    'legend.borderpad': 0.4,
    'legend.handlelength': 1.2,
    'legend.handletextpad': 0.4,
    'legend.columnspacing': 0.8,
})

# ═══════════════════════════════════════════════════════════════
# Color Palette — Nature-inspired, colorblind-safe
# ═══════════════════════════════════════════════════════════════
# Based on Wong (2011) Nature Methods colorblind-safe palette

COLORS = {
    'CMC':          '#0072B2',   # Blue (strong, primary)
    'YYC':          '#E69F00',   # Orange
    'DNA Fountain': '#009E73',   # Bluish-green
    'HEDGES':       '#CC79A7',   # Reddish-purple
    'DRRC':         '#D55E00',   # Vermillion
}

HATCHES = {
    'CMC':          None,
    'YYC':          None,
    'DNA Fountain': None,
    'HEDGES':       None,
    'DRRC':         None,
}

METHOD_ORDER = ['CMC', 'YYC', 'DNA Fountain', 'HEDGES', 'DRRC']
FILES = ['繁体版出師表', '背影']
FILE_LABELS = {'繁体版出師表': 'Chu Shi Biao\n(Traditional)', '背影': 'Bei Ying\n(Simplified)'}
FILE_LABELS_SHORT = {'繁体版出師表': 'Chu Shi Biao', '背影': 'Bei Ying'}

OUT_DIR = Path(__file__).resolve().parent
OUT_DIR.mkdir(exist_ok=True)

# ═══════════════════════════════════════════════════════════════
# Authoritative Data (from run_full_experiment.py RESULTS)
# ═══════════════════════════════════════════════════════════════

RESULTS = [
    {'method': 'CMC', 'file': '繁体版出師表', 'scheme': 'B方案(繁转简)',
     'density': 1.9422, 'gc_global': 50.0, 'gc_min': 45.0, 'gc_max': 55.0,
     'homopolymer': 2, 'motif': 0.0, 'hamming': 76,
     'mfe_mean': -30.6, 'mfe_min': -43.9, 'mfe_max': -21.2},
    {'method': 'CMC', 'file': '繁体版出師表', 'scheme': 'A方案(繁体直编)',
     'density': 1.9363, 'gc_global': 50.0, 'gc_min': 45.0, 'gc_max': 55.0,
     'homopolymer': 2, 'motif': 0.0, 'hamming': 81,
     'mfe_mean': -31.1, 'mfe_min': -43.8, 'mfe_max': -22.0},
    {'method': 'CMC', 'file': '背影', 'scheme': '简体存储',
     'density': 2.1734, 'gc_global': 50.0, 'gc_min': 45.0, 'gc_max': 55.0,
     'homopolymer': 2, 'motif': 0.0, 'hamming': 80,
     'mfe_mean': -27.8, 'mfe_min': -31.2, 'mfe_max': -24.3},
    {'method': 'YYC', 'file': '繁体版出師表', 'scheme': '',
     'density': 1.2368, 'gc_global': 54.9548, 'gc_min': 12.50, 'gc_max': 100.0,
     'homopolymer': 4, 'motif': 12.7673, 'hamming': 56,
     'mfe_mean': -13.0, 'mfe_min': -23.5, 'mfe_max': 1.7},
    {'method': 'YYC', 'file': '背影', 'scheme': '',
     'density': 1.2170, 'gc_global': 57.6966, 'gc_min': 22.22, 'gc_max': 100.0,
     'homopolymer': 4, 'motif': 13.8078, 'hamming': 41,
     'mfe_mean': -15.79, 'mfe_min': -26.9, 'mfe_max': 1.5},
    {'method': 'DNA Fountain', 'file': '繁体版出師表', 'scheme': '',
     'density': 1.0416, 'gc_global': 50.9695, 'gc_min': 34.0, 'gc_max': 68.0,
     'homopolymer': 3, 'motif': 36.8248, 'hamming': 17,
     'mfe_mean': -11.4447, 'mfe_min': -21.5, 'mfe_max': -1.6},
    {'method': 'DNA Fountain', 'file': '背影', 'scheme': '',
     'density': 1.0407, 'gc_global': 51.2859, 'gc_min': 32.0, 'gc_max': 68.0,
     'homopolymer': 3, 'motif': 36.5730, 'hamming': 17,
     'mfe_mean': -12.4023, 'mfe_min': -23.9, 'mfe_max': -1.4},
    {'method': 'HEDGES', 'file': '繁体版出師表', 'scheme': '',
     'density': 0.1573, 'gc_global': 50.1451, 'gc_min': 16.67, 'gc_max': 75.0,
     'homopolymer': 4, 'motif': 13.0405, 'hamming': 4,
     'mfe_mean': -88.5888, 'mfe_min': -106.0, 'mfe_max': -73.3},
    {'method': 'HEDGES', 'file': '背影', 'scheme': '',
     'density': 0.2729, 'gc_global': 50.5595, 'gc_min': 16.67, 'gc_max': 83.33,
     'homopolymer': 4, 'motif': 13.0876, 'hamming': 4,
     'mfe_mean': -89.6482, 'mfe_min': -111.5, 'mfe_max': -73.7},
    {'method': 'DRRC', 'file': '繁体版出師表', 'scheme': '',
     'density': 0.1594, 'gc_global': 50.3584, 'gc_min': 42.67, 'gc_max': 56.0,
     'homopolymer': 3, 'motif': 38.9244, 'hamming': 14,
     'mfe_mean': -3.437, 'mfe_min': -13.3, 'mfe_max': 0.5},
    {'method': 'DRRC', 'file': '背影', 'scheme': '',
     'density': 0.2373, 'gc_global': 50.3358, 'gc_min': 43.33, 'gc_max': 56.67,
     'homopolymer': 3, 'motif': 38.7277, 'hamming': 13,
     'mfe_mean': -2.9685, 'mfe_min': -12.0, 'mfe_max': 1.1},
]


def get_val(method, file_name, key, scheme_filter=None):
    """Get metric value for a method+file. For CMC 出師表, default to B方案."""
    for r in RESULTS:
        if r['method'] == method and r['file'] == file_name:
            if method == 'CMC':
                if scheme_filter and r['scheme'] != scheme_filter:
                    continue
                if not scheme_filter:
                    if file_name == '繁体版出師表' and r['scheme'] != 'B方案(繁转简)':
                        continue
                    if file_name == '背影' and r['scheme'] != '简体存储':
                        continue
            return r[key]
    return None


def add_significance_bracket(ax, x1, x2, y, text, lw=0.5, fontsize=6):
    """Add a bracket between two bars (for noting best performance)."""
    h = (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.015
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=lw, color='#333333')
    ax.text((x1 + x2) / 2, y + h * 1.2, text, ha='center', va='bottom',
            fontsize=fontsize, color='#333333', fontstyle='italic')


def save_figure(fig, name):
    """Save in both PNG (600 dpi) and SVG."""
    fig.savefig(OUT_DIR / f'{name}.png')
    fig.savefig(OUT_DIR / f'{name}.svg')
    plt.close(fig)
    print(f"  saved {name}.png/svg")


# ═══════════════════════════════════════════════════════════════
# Figure 1: Storage Density Comparison
# ═══════════════════════════════════════════════════════════════

def fig1_storage_density():
    print("\n[Fig 1] Storage Density ...")
    fig, axes = plt.subplots(1, 2, figsize=(6.0, 2.8), sharey=True,
                             gridspec_kw={'wspace': 0.08})

    for idx, fname in enumerate(FILES):
        ax = axes[idx]
        labels, vals, colors = [], [], []

        for m in METHOD_ORDER:
            if m == 'CMC' and fname == '繁体版出師表':
                for sl, sf in [('CMC (Direct)', 'A方案(繁体直编)'),
                               ('CMC (T→S)', 'B方案(繁转简)')]:
                    v = get_val('CMC', fname, 'density', sf)
                    labels.append(sl)
                    vals.append(v)
                    colors.append(COLORS['CMC'])
            else:
                v = get_val(m, fname, 'density')
                if v is not None:
                    labels.append(m)
                    vals.append(v)
                    colors.append(COLORS[m])

        x = np.arange(len(labels))
        bars = ax.bar(x, vals, width=0.62, color=colors,
                      edgecolor='white', linewidth=0.8, zorder=3)

        # Value annotations
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.03,
                    f'{val:.2f}', ha='center', va='bottom', fontsize=6,
                    fontweight='bold', color='#333333')

        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=6.5)
        ax.set_title(FILE_LABELS_SHORT[fname], fontsize=9, fontweight='bold', pad=6)

        if idx == 0:
            ax.set_ylabel('Storage Density (bits/nt)', fontsize=8)

        ax.set_ylim(0, max(vals) * 1.22)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(6, integer=False))
        ax.grid(axis='y', alpha=0.15, linewidth=0.4, zorder=0)
        ax.set_axisbelow(True)

        # Highlight CMC as best
        for i, (l, v) in enumerate(zip(labels, vals)):
            if 'CMC' in l:
                bars[i].set_edgecolor(COLORS['CMC'])
                bars[i].set_linewidth(1.2)

    fig.suptitle('(a) Storage Density Comparison', fontsize=10, fontweight='bold',
                 y=1.02, x=0.5)
    save_figure(fig, 'fig1_storage_density')


# ═══════════════════════════════════════════════════════════════
# Figure 2: GC Content with Local Range (Dumbbell / Lollipop)
# ═══════════════════════════════════════════════════════════════

def fig2_gc_content():
    print("\n[Fig 2] GC Content ...")
    from matplotlib.lines import Line2D

    fig, axes = plt.subplots(1, 2, figsize=(6.0, 3.2), sharey=True,
                             gridspec_kw={'wspace': 0.08})

    for idx, fname in enumerate(FILES):
        ax = axes[idx]
        labels, gc_globals, gc_mins, gc_maxs, colors = [], [], [], [], []

        for m in METHOD_ORDER:
            gc_g = get_val(m, fname, 'gc_global')
            gc_lo = get_val(m, fname, 'gc_min')
            gc_hi = get_val(m, fname, 'gc_max')
            if gc_g is None:
                continue
            labels.append(m)
            gc_globals.append(gc_g)
            gc_mins.append(gc_lo)
            gc_maxs.append(gc_hi)
            colors.append(COLORS[m])

        x = np.arange(len(labels))

        for i in range(len(labels)):
            # Range bar (vertical line from min to max)
            ax.plot([x[i], x[i]], [gc_mins[i], gc_maxs[i]], color=colors[i],
                    linewidth=3.0, solid_capstyle='round', alpha=0.30, zorder=2)
            # Global GC dot
            ax.scatter(x[i], gc_globals[i], color=colors[i], s=55, zorder=4,
                       edgecolors='white', linewidths=0.8)
            # Annotate global value
            offset_y = 3.5 if gc_globals[i] < 60 else -5
            ax.annotate(f'{gc_globals[i]:.1f}%', (x[i], gc_globals[i]),
                        textcoords='offset points', xytext=(0, offset_y),
                        ha='center', fontsize=5.5, color=colors[i], fontweight='bold')

        # Ideal line
        ax.axhline(y=50.0, color='#999999', linestyle='--', linewidth=0.8,
                   alpha=0.5, zorder=1)
        if idx == 1:
            ax.text(len(labels) - 0.3, 51.5, 'Ideal 50%',
                    fontsize=5.5, color='#999999', fontstyle='italic', ha='right')

        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=25, ha='right', fontsize=6.5)
        ax.set_title(FILE_LABELS_SHORT[fname], fontsize=9, fontweight='bold', pad=6)

        if idx == 0:
            ax.set_ylabel('GC Content (%)', fontsize=8)

        ax.set_ylim(0, 108)
        ax.yaxis.set_major_locator(mticker.MultipleLocator(20))
        ax.yaxis.set_minor_locator(mticker.MultipleLocator(10))
        ax.grid(axis='y', alpha=0.12, linewidth=0.4, zorder=0)
        ax.set_axisbelow(True)

    # Shared legend (upper-right of right panel)
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#666666',
               markersize=5, label='Global GC'),
        Line2D([0], [0], color='#666666', linewidth=3, alpha=0.30,
               label='Local range (min–max)'),
    ]
    axes[1].legend(handles=legend_elements, loc='upper right', fontsize=6,
                   borderaxespad=0.5)

    fig.suptitle('(b) GC Content — Global Mean with Local Range',
                 fontsize=10, fontweight='bold', y=1.02)
    save_figure(fig, 'fig2_gc_distribution')


# ═══════════════════════════════════════════════════════════════
# Figure 3: Maximum Homopolymer Length
# ═══════════════════════════════════════════════════════════════

def fig3_homopolymer():
    print("\n[Fig 3] Homopolymer ...")
    fig, axes = plt.subplots(1, 2, figsize=(6.0, 2.6), sharey=True,
                             gridspec_kw={'wspace': 0.08})

    for idx, fname in enumerate(FILES):
        ax = axes[idx]
        labels, vals, colors = [], [], []

        for m in METHOD_ORDER:
            v = get_val(m, fname, 'homopolymer')
            if v is not None:
                labels.append(m)
                vals.append(v)
                colors.append(COLORS[m])

        x = np.arange(len(labels))
        bars = ax.bar(x, vals, width=0.58, color=colors,
                      edgecolor='white', linewidth=0.8, zorder=3)

        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.08,
                    str(int(val)), ha='center', va='bottom', fontsize=7,
                    fontweight='bold', color='#333333')

        # Constraint line
        ax.axhline(y=2, color=COLORS['CMC'], linestyle=':', linewidth=0.8,
                   alpha=0.5, zorder=1)
        if idx == 1:
            ax.text(len(labels) - 0.5, 2.1, 'CMC limit = 2',
                    fontsize=5.5, color=COLORS['CMC'], fontstyle='italic',
                    ha='right')

        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=25, ha='right', fontsize=6.5)
        ax.set_title(FILE_LABELS_SHORT[fname], fontsize=9, fontweight='bold', pad=6)

        if idx == 0:
            ax.set_ylabel('Max Homopolymer Length (nt)', fontsize=8)

        ax.set_ylim(0, max(vals) * 1.3)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(6, integer=True))
        ax.grid(axis='y', alpha=0.15, linewidth=0.4, zorder=0)
        ax.set_axisbelow(True)

    fig.suptitle('(c) Maximum Homopolymer Length', fontsize=10, fontweight='bold',
                 y=1.02)
    save_figure(fig, 'fig4_homopolymer_comparison')


# ═══════════════════════════════════════════════════════════════
# Figure 4: Undesired Motif Content
# ═══════════════════════════════════════════════════════════════

def fig4_motif():
    print("\n[Fig 4] Motif Content ...")
    fig, axes = plt.subplots(1, 2, figsize=(6.0, 2.6), sharey=True,
                             gridspec_kw={'wspace': 0.08})

    for idx, fname in enumerate(FILES):
        ax = axes[idx]
        labels, vals, colors = [], [], []

        for m in METHOD_ORDER:
            v = get_val(m, fname, 'motif')
            if v is not None:
                labels.append(m)
                vals.append(v)
                colors.append(COLORS[m])

        x = np.arange(len(labels))
        bars = ax.bar(x, vals, width=0.58, color=colors,
                      edgecolor='white', linewidth=0.8, zorder=3)

        for bar, val in zip(bars, vals):
            label_txt = f'{val:.1f}%' if val > 0 else '0%'
            y_pos = bar.get_height() + 0.4 if val > 0 else 0.4
            ax.text(bar.get_x() + bar.get_width() / 2, y_pos,
                    label_txt, ha='center', va='bottom', fontsize=6,
                    fontweight='bold', color='#333333')

        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=25, ha='right', fontsize=6.5)
        ax.set_title(FILE_LABELS_SHORT[fname], fontsize=9, fontweight='bold', pad=6)

        if idx == 0:
            ax.set_ylabel('Undesired Motif Content (%)', fontsize=8)

        ax.set_ylim(0, max(vals) * 1.2)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(6))
        ax.grid(axis='y', alpha=0.15, linewidth=0.4, zorder=0)
        ax.set_axisbelow(True)

    fig.suptitle('(d) Undesired Motif Content', fontsize=10, fontweight='bold',
                 y=1.02)
    save_figure(fig, 'fig3_motif_comparison')


# ═══════════════════════════════════════════════════════════════
# Figure 5: Minimum Hamming Distance
# ═══════════════════════════════════════════════════════════════

def fig5_hamming():
    print("\n[Fig 5] Hamming Distance ...")
    fig, axes = plt.subplots(1, 2, figsize=(6.0, 2.8), sharey=True,
                             gridspec_kw={'wspace': 0.08})

    for idx, fname in enumerate(FILES):
        ax = axes[idx]
        labels, vals, colors = [], [], []

        for m in METHOD_ORDER:
            if m == 'CMC' and fname == '繁体版出師表':
                for sl, sf in [('CMC (Direct)', 'A方案(繁体直编)'),
                               ('CMC (T→S)', 'B方案(繁转简)')]:
                    v = get_val('CMC', fname, 'hamming', sf)
                    labels.append(sl)
                    vals.append(v)
                    colors.append(COLORS['CMC'])
            else:
                v = get_val(m, fname, 'hamming')
                if v is not None:
                    labels.append(m)
                    vals.append(v)
                    colors.append(COLORS[m])

        x = np.arange(len(labels))
        bars = ax.bar(x, vals, width=0.58, color=colors,
                      edgecolor='white', linewidth=0.8, zorder=3)

        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.8,
                    str(int(val)), ha='center', va='bottom', fontsize=6.5,
                    fontweight='bold', color='#333333')

        # Highlight CMC bars
        for i, l in enumerate(labels):
            if 'CMC' in l:
                bars[i].set_edgecolor(COLORS['CMC'])
                bars[i].set_linewidth(1.2)

        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=6.5)
        ax.set_title(FILE_LABELS_SHORT[fname], fontsize=9, fontweight='bold', pad=6)

        if idx == 0:
            ax.set_ylabel('Min. Hamming Distance', fontsize=8)

        ax.set_ylim(0, max(vals) * 1.2)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(6, integer=True))
        ax.grid(axis='y', alpha=0.15, linewidth=0.4, zorder=0)
        ax.set_axisbelow(True)

    fig.suptitle('(e) Minimum Hamming Distance (120 nt fragments)',
                 fontsize=10, fontweight='bold', y=1.02)
    save_figure(fig, 'fig5_hamming_comparison')


# ═══════════════════════════════════════════════════════════════
# Figure 6: MFE Comparison (with error bars)
# ═══════════════════════════════════════════════════════════════

def fig6_mfe():
    print("\n[Fig 6] MFE ...")
    fig, axes = plt.subplots(1, 2, figsize=(6.0, 2.8),
                             gridspec_kw={'wspace': 0.22})

    for idx, fname in enumerate(FILES):
        ax = axes[idx]
        labels, means, lo_err, hi_err, colors = [], [], [], [], []

        for m in METHOD_ORDER:
            v_mean = get_val(m, fname, 'mfe_mean')
            v_min = get_val(m, fname, 'mfe_min')
            v_max = get_val(m, fname, 'mfe_max')
            if v_mean is not None:
                labels.append(m)
                means.append(v_mean)
                # Error bar: distance from mean to min/max
                lo_err.append(abs(v_mean - v_min))
                hi_err.append(abs(v_max - v_mean))
                colors.append(COLORS[m])

        x = np.arange(len(labels))
        bars = ax.bar(x, means, width=0.58, color=colors,
                      edgecolor='white', linewidth=0.8, zorder=3,
                      yerr=[lo_err, hi_err], capsize=3,
                      error_kw={'linewidth': 0.8, 'capthick': 0.8,
                                'ecolor': '#555555', 'zorder': 5})

        # Value annotations (positioned below bars for negative values)
        for i, (bar, val) in enumerate(zip(bars, means)):
            y_pos = val - lo_err[i] - 2.5
            ax.text(bar.get_x() + bar.get_width() / 2, y_pos,
                    f'{val:.1f}', ha='center', va='top', fontsize=5.5,
                    fontweight='bold', color='#333333')

        ax.axhline(y=0, color='#AAAAAA', linewidth=0.6, zorder=1)

        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=25, ha='right', fontsize=6.5)
        ax.set_title(FILE_LABELS_SHORT[fname], fontsize=9, fontweight='bold', pad=6)

        if idx == 0:
            ax.set_ylabel('MFE (kcal/mol)', fontsize=8)

        ax.yaxis.set_major_locator(mticker.MaxNLocator(7))
        ax.grid(axis='y', alpha=0.15, linewidth=0.4, zorder=0)
        ax.set_axisbelow(True)

    fig.suptitle('(f) Minimum Free Energy (MFE) — Mean ± Range',
                 fontsize=10, fontweight='bold', y=1.02)
    save_figure(fig, 'fig6_mfe_comparison')


# ═══════════════════════════════════════════════════════════════
# Figure 7 (Bonus): Radar / Spider Chart — Normalized Comparison
# ═══════════════════════════════════════════════════════════════

def fig7_radar():
    """Radar chart showing normalized metrics for each method (出師表 only)."""
    print("\n[Fig 7] Radar chart ...")

    categories = ['Storage\nDensity', 'GC\nStability', 'Homopolymer\nControl',
                  'Motif\nAvoidance', 'Hamming\nDistance', 'MFE\nStability']
    N = len(categories)

    # Normalize: higher = better, 0-1 scale
    fname = '繁体版出師表'
    raw_data = {}
    for m in METHOD_ORDER:
        density = get_val(m, fname, 'density') or 0
        gc_dev = abs((get_val(m, fname, 'gc_global') or 50) - 50)  # lower is better
        homo = get_val(m, fname, 'homopolymer') or 0  # lower is better
        motif = get_val(m, fname, 'motif') or 0  # lower is better
        hamming = get_val(m, fname, 'hamming') or 0  # higher is better
        gc_range = (get_val(m, fname, 'gc_max') or 50) - (get_val(m, fname, 'gc_min') or 50)

        raw_data[m] = [density, gc_range, homo, motif, hamming, abs(get_val(m, fname, 'mfe_mean') or 0)]

    # Normalize each dimension
    norm_data = {}
    for m in METHOD_ORDER:
        r = raw_data[m]
        # density: higher is better → linear 0-1
        all_d = [raw_data[mm][0] for mm in METHOD_ORDER]
        d_norm = r[0] / max(all_d) if max(all_d) > 0 else 0

        # gc_range: lower is better → invert
        all_gcr = [raw_data[mm][1] for mm in METHOD_ORDER]
        gcr_norm = 1 - (r[1] / max(all_gcr)) if max(all_gcr) > 0 else 1

        # homo: lower is better → invert
        all_h = [raw_data[mm][2] for mm in METHOD_ORDER]
        h_norm = 1 - ((r[2] - min(all_h)) / (max(all_h) - min(all_h))) if max(all_h) > min(all_h) else 1

        # motif: lower is better → invert
        all_mo = [raw_data[mm][3] for mm in METHOD_ORDER]
        mo_norm = 1 - (r[3] / max(all_mo)) if max(all_mo) > 0 else 1

        # hamming: higher is better
        all_ham = [raw_data[mm][4] for mm in METHOD_ORDER]
        ham_norm = r[4] / max(all_ham) if max(all_ham) > 0 else 0

        # MFE: closer to 0 is better (less secondary structure) → invert abs
        all_mfe = [raw_data[mm][5] for mm in METHOD_ORDER]
        mfe_norm = 1 - (r[5] / max(all_mfe)) if max(all_mfe) > 0 else 1

        norm_data[m] = [d_norm, gcr_norm, h_norm, mo_norm, ham_norm, mfe_norm]

    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles += angles[:1]  # close the polygon

    fig, ax = plt.subplots(figsize=(4.0, 4.0), subplot_kw=dict(polar=True))

    for m in METHOD_ORDER:
        values = norm_data[m] + [norm_data[m][0]]
        ax.plot(angles, values, 'o-', linewidth=1.2, markersize=3.5,
                color=COLORS[m], label=m, alpha=0.85)
        ax.fill(angles, values, alpha=0.08, color=COLORS[m])

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_thetagrids(np.degrees(angles[:-1]), categories, fontsize=6.5)

    ax.set_ylim(0, 1.05)
    ax.set_yticks([0.25, 0.50, 0.75, 1.0])
    ax.set_yticklabels(['0.25', '0.50', '0.75', '1.00'], fontsize=5.5,
                        color='#888888')
    ax.yaxis.grid(True, linewidth=0.3, alpha=0.3)
    ax.xaxis.grid(True, linewidth=0.3, alpha=0.3)
    ax.spines['polar'].set_linewidth(0.4)

    ax.legend(loc='upper right', bbox_to_anchor=(1.28, 1.12), fontsize=6.5,
              framealpha=0.9, edgecolor='#CCCCCC')

    fig.suptitle('Normalized Performance Comparison\n(Chu Shi Biao)',
                 fontsize=9, fontweight='bold', y=1.0)
    plt.tight_layout()
    save_figure(fig, 'fig_radar_comparison')


# ═══════════════════════════════════════════════════════════════
# Figure 8 (Bonus): GC Fluctuation Range — Grouped Horizontal Bar
# ═══════════════════════════════════════════════════════════════

def fig8_gc_fluctuation():
    """Horizontal bar chart split into two panels — one per file."""
    print("\n[Fig 8] GC fluctuation range ...")
    fig, axes = plt.subplots(1, 2, figsize=(6.0, 2.6), sharex=True,
                             gridspec_kw={'wspace': 0.12})

    for idx, fname in enumerate(FILES):
        ax = axes[idx]
        labels, vals, colors = [], [], []

        for m in METHOD_ORDER:
            gc_min = get_val(m, fname, 'gc_min')
            gc_max = get_val(m, fname, 'gc_max')
            if gc_min is not None and gc_max is not None:
                fluct = gc_max - gc_min
                labels.append(m)
                vals.append(fluct)
                colors.append(COLORS[m])

        x = np.arange(len(labels))
        bars = ax.bar(x, vals, width=0.58, color=colors,
                      edgecolor='white', linewidth=0.8, zorder=3)

        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1.0,
                    f'{val:.1f}%', ha='center', va='bottom', fontsize=6,
                    fontweight='bold', color='#444444')

        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=25, ha='right', fontsize=6.5)
        ax.set_title(FILE_LABELS_SHORT[fname], fontsize=9, fontweight='bold', pad=6)

        if idx == 0:
            ax.set_ylabel('GC Fluctuation Range (%)', fontsize=8)

        ax.set_ylim(0, max(vals) * 1.2)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(6))
        ax.grid(axis='y', alpha=0.15, linewidth=0.4, zorder=0)
        ax.set_axisbelow(True)

    fig.suptitle('GC Fluctuation Range (Lower is Better)',
                 fontsize=10, fontweight='bold', y=1.02)
    save_figure(fig, 'fig_gc_fluctuation')


# ═══════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("=" * 60)
    print("  Generating Publication-Quality Figures")
    print("  Style: Nature / Science journal standard")
    print("=" * 60)

    fig1_storage_density()
    fig2_gc_content()
    fig3_homopolymer()
    fig4_motif()
    fig5_hamming()
    fig6_mfe()
    fig7_radar()
    fig8_gc_fluctuation()

    print("\n" + "=" * 60)
    print("  All figures generated successfully!")
    print(f"  Output: {OUT_DIR}")
    print("=" * 60)
