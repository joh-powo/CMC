# -*- coding: utf-8 -*-
"""
Publication-quality comparison charts for CMC paper.
Generates:
  1. Simplified Chinese text (Beiying) - 4-panel bar chart
  2. Traditional Chinese text (Chu Shi Biao) - 4-panel bar chart
  3. Radar chart - Multi-dimensional overview (supplementary)
  4. GC content comparison with ideal range band (supplementary)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from pathlib import Path
import sys
import io

# fix GBK stdout on Windows
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# -- Global style --
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'Microsoft YaHei', 'sans-serif'],
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'axes.titleweight': 'bold',
    'xtick.labelsize': 9.5,
    'ytick.labelsize': 9.5,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'axes.edgecolor': '#333333',
})

OUT = Path(__file__).parent
OUT.mkdir(exist_ok=True)

# -- Data --
methods = ['CMC', 'YYC', 'DNA Fountain', 'HEDGES']

# Traditional Chinese (Chu Shi Biao)
trad = {
    'density':     [1.94, 0.16, 0.97, 0.15],
    'gc':          [50.0, 50.76, 50.97, 45.64],
    'homopolymer': [2, 49, 3, 8],
    'hamming':     [76, 41, 17, 4],
}

# Simplified Chinese (Beiying)
simp = {
    'density':     [2.17, 0.24, 0.98, 0.26],
    'gc':          [50.0, 50.89, 51.29, 46.88],
    'homopolymer': [2, 64, 3, 7],
    'hamming':     [80, 47, 17, 4],
}

# Palette: CMC highlighted, others muted
COLORS = ['#2563EB', '#94A3B8', '#F59E0B', '#6366F1']
EDGE_COLORS = ['#1D4ED8', '#64748B', '#D97706', '#4F46E5']
HATCHES = ['///', '', '...', 'xxx']


def draw_four_panel(data, sample_label, filename):
    """Draw a single 4-panel bar chart figure."""
    fig, axes = plt.subplots(2, 2, figsize=(8.5, 7))
    fig.suptitle('Encoding Performance Comparison - ' + sample_label,
                 fontsize=14, fontweight='bold', y=0.98)

    x = np.arange(len(methods))
    bar_w = 0.55

    panels = [
        ('density',     '(a) Information Density',       'bits/nt',   None,  None),
        ('gc',          '(b) GC Content',                '%',         40,    55),
        ('homopolymer', '(c) Max Homopolymer Length',     'nt',        None,  None),
        ('hamming',     '(d) Min Hamming Distance',       '',          None,  None),
    ]

    for ax, (key, title, ylabel, ymin, ymax) in zip(axes.flat, panels):
        vals = data[key]
        bars = ax.bar(x, vals, bar_w,
                      color=COLORS, edgecolor=EDGE_COLORS,
                      linewidth=0.8, zorder=3)

        # Apply hatches
        for bar, h in zip(bars, HATCHES):
            bar.set_hatch(h)

        # Value labels
        for i, (bar, v) in enumerate(zip(bars, vals)):
            offset = max(vals) * 0.03
            fontw = 'bold' if i == 0 else 'normal'
            color = '#1D4ED8' if i == 0 else '#333333'
            fmt = '{:.2f}'.format(v) if isinstance(v, float) and v < 10 else str(v)
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + offset,
                    fmt, ha='center', va='bottom', fontsize=8.5,
                    fontweight=fontw, color=color)

        # GC panel: add ideal range band
        if key == 'gc':
            ax.axhspan(45, 55, color='#10B981', alpha=0.08, zorder=1)
            ax.axhline(50, color='#10B981', linestyle='--', linewidth=0.8,
                       alpha=0.6, zorder=2, label='Ideal 50%')
            ax.legend(loc='lower left', framealpha=0.7, edgecolor='#cccccc')

        ax.set_xticks(x)
        ax.set_xticklabels(methods, fontweight='medium')
        ax.set_ylabel(ylabel, fontweight='medium')
        ax.set_title(title, pad=8, fontsize=11)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', linestyle=':', linewidth=0.4, alpha=0.6, zorder=0)

        if ymin is not None:
            ax.set_ylim(bottom=ymin)
        if ymax is not None:
            ax.set_ylim(top=ymax + (ymax - (ymin or 0)) * 0.15)

        # Add padding at top for labels
        cur_ymin, cur_ymax = ax.get_ylim()
        ax.set_ylim(cur_ymin, cur_ymax + (cur_ymax - cur_ymin) * 0.12)

    fig.tight_layout(rect=[0, 0, 1, 0.94], h_pad=3.0, w_pad=2.5)
    fig.savefig(OUT / filename, format='svg')
    fig.savefig(OUT / filename.replace('.svg', '.png'), format='png')
    plt.close(fig)
    print('  [OK] ' + filename)


def draw_radar(filename):
    """Radar chart normalising each metric to [0, 1] where bigger = better."""
    categories = ['Information\nDensity', 'GC Balance', 'Homopolymer\nControl',
                  'Hamming\nDistance']
    N = len(categories)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles += angles[:1]

    # Average of two samples
    raw = {}
    for m, td, sd in zip(methods,
                         zip(trad['density'], trad['gc'], trad['homopolymer'], trad['hamming']),
                         zip(simp['density'], simp['gc'], simp['homopolymer'], simp['hamming'])):
        raw[m] = [(a + b) / 2 for a, b in zip(td, sd)]

    # Normalise
    def norm(idx, invert=False):
        col = [raw[m][idx] for m in methods]
        mn, mx = min(col), max(col)
        rng = mx - mn if mx != mn else 1
        out = [(v - mn) / rng for v in col]
        if invert:
            out = [1 - v for v in out]
        return out

    density_n  = norm(0)
    gc_dev     = [abs(raw[m][1] - 50.0) for m in methods]
    gc_mn, gc_mx = min(gc_dev), max(gc_dev)
    gc_rng = gc_mx - gc_mn if gc_mx != gc_mn else 1
    gc_n       = [1 - (v - gc_mn) / gc_rng for v in gc_dev]
    homo_n     = norm(2, invert=True)
    hamming_n  = norm(3)

    fig, ax = plt.subplots(figsize=(5.5, 5.5), subplot_kw=dict(polar=True))
    fig.suptitle('Multi-Metric Performance Radar (Average of Two Samples)',
                 fontsize=12, fontweight='bold', y=0.97)

    radar_colors = ['#2563EB', '#94A3B8', '#F59E0B', '#6366F1']
    markers = ['o', 's', 'D', '^']

    for i, m in enumerate(methods):
        values = [density_n[i], gc_n[i], homo_n[i], hamming_n[i]]
        values += values[:1]
        lw = 2.2 if i == 0 else 1.4
        ax.plot(angles, values, 'o-', linewidth=lw, label=m,
                color=radar_colors[i], marker=markers[i], markersize=6)
        ax.fill(angles, values, alpha=0.08 if i != 0 else 0.15,
                color=radar_colors[i])

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=9.5, fontweight='medium')
    ax.set_ylim(0, 1.15)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(['0.25', '0.50', '0.75', '1.00'], fontsize=7.5, color='#666')
    ax.grid(linewidth=0.4, linestyle=':', alpha=0.6)
    ax.legend(loc='lower right', bbox_to_anchor=(1.25, -0.05),
              framealpha=0.85, edgecolor='#cccccc')

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(OUT / filename, format='svg')
    fig.savefig(OUT / filename.replace('.svg', '.png'), format='png')
    plt.close(fig)
    print('  [OK] ' + filename)


def draw_gc_band(filename):
    """GC content grouped bar chart with ideal range band for both samples."""
    fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=True)
    fig.suptitle('Global GC Content Comparison', fontsize=13, fontweight='bold', y=0.99)

    x = np.arange(len(methods))
    bar_w = 0.5

    for ax, data, label in [(axes[0], trad, 'Traditional Chinese\n(Chu Shi Biao)'),
                            (axes[1], simp, 'Simplified Chinese\n(Beiying)')]:
        ax.axhspan(45, 55, color='#10B981', alpha=0.10, zorder=1,
                   label='Ideal range (45-55%)')
        ax.axhline(50, color='#10B981', linestyle='--', linewidth=0.9,
                   alpha=0.5, zorder=2)

        bars = ax.bar(x, data['gc'], bar_w,
                      color=COLORS, edgecolor=EDGE_COLORS,
                      linewidth=0.8, zorder=3)
        for bar, h in zip(bars, HATCHES):
            bar.set_hatch(h)

        for i, (bar, v) in enumerate(zip(bars, data['gc'])):
            fw = 'bold' if i == 0 else 'normal'
            clr = '#1D4ED8' if i == 0 else '#333333'
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.25,
                    '{:.1f}%'.format(v), ha='center', va='bottom', fontsize=8.5,
                    fontweight=fw, color=clr)

        ax.set_xticks(x)
        ax.set_xticklabels(methods, fontweight='medium')
        ax.set_xlabel(label, fontsize=10, fontweight='medium')
        ax.set_ylim(40, 56)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', linestyle=':', linewidth=0.4, alpha=0.6, zorder=0)

    axes[0].set_ylabel('GC Content (%)', fontweight='medium')
    axes[0].legend(loc='lower left', fontsize=8, framealpha=0.7, edgecolor='#cccccc')

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(OUT / filename, format='svg')
    fig.savefig(OUT / filename.replace('.svg', '.png'), format='png')
    plt.close(fig)
    print('  [OK] ' + filename)


if __name__ == '__main__':
    print('Generating charts ...')
    draw_four_panel(simp, 'Simplified Chinese Text (Beiying)',
                    'fig_simplified_comparison.svg')
    draw_four_panel(trad, 'Traditional Chinese Text (Chu Shi Biao)',
                    'fig_traditional_comparison.svg')
    draw_radar('fig_radar_comparison.svg')
    draw_gc_band('fig_gc_comparison.svg')
    print('Done.')
