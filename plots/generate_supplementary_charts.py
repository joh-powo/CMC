# -*- coding: utf-8 -*-
"""
Supplementary charts for CMC paper:
  1. Local GC distribution line chart (24nt window, similar to paper Fig.5)
  2. Encoding efficiency waterfall chart (text -> wubi -> compressed -> DNA)
  3. Hamming distance heatmap (pairwise distance matrix)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path
import sys, io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'Microsoft YaHei', 'sans-serif'],
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'axes.titleweight': 'bold',
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.8,
    'axes.edgecolor': '#333333',
})

OUT = Path(__file__).parent

# ── Load DNA sequences ────────────────────────────────────────────
dna_trad = (OUT / 'dna_trad.txt').read_text().strip()
dna_simp = (OUT / 'dna_simp.txt').read_text().strip()

# ── Hamming distance data (from MATLAB runs) ─────────────────────
# Traditional (Chu Shi Biao) - 6 fragments x 120nt
hamming_trad = np.zeros((6, 6))
pairs_trad = {
    (0,1):89, (0,2):97, (0,3):84, (0,4):89, (0,5):85,
    (1,2):92, (1,3):90, (1,4):100,(1,5):81,
    (2,3):86, (2,4):84, (2,5):99,
    (3,4):76, (3,5):80,
    (4,5):83,
}
for (i,j), v in pairs_trad.items():
    hamming_trad[i,j] = v
    hamming_trad[j,i] = v

# Simplified (Beiying) - 6 fragments x 120nt
hamming_simp = np.zeros((6, 6))
pairs_simp = {
    (0,1):83, (0,2):88, (0,3):86, (0,4):84, (0,5):90,
    (1,2):85, (1,3):91, (1,4):89, (1,5):88,
    (2,3):87, (2,4):88, (2,5):80,
    (3,4):80, (3,5):91,
    (4,5):89,
}
for (i,j), v in pairs_simp.items():
    hamming_simp[i,j] = v
    hamming_simp[j,i] = v

# Stage data (from MATLAB output)
# Format: (stage_label, bits_value)
# Original text: chars * 16 bits per char (all Chinese)
# Wubi: wubi_len * 8 bits (ASCII)
# Compressed: compressed_bytes * 8 bits
# DNA: total_nt * 2 bits (each nt stores ~2 bits theoretically)
stages_trad = {
    'chars': 752, 'wubi_len': 2215, 'compressed_bytes': 1060, 'dna_nt': 6195
}
stages_simp = {
    'chars': 1305, 'wubi_len': 3784, 'compressed_bytes': 1639, 'dna_nt': 9607
}


# ════════════════════════════════════════════════════════════════════
# Chart 1: Local GC Distribution (120nt sliding window, step=24nt)
# ════════════════════════════════════════════════════════════════════
def compute_local_gc_sliding(dna, window=120, step=24):
    """Compute GC content with a sliding window for smoother curves."""
    gc_values = []
    positions = []
    idx = 0
    for start in range(0, len(dna) - window + 1, step):
        segment = dna[start:start + window]
        gc = sum(1 for b in segment if b in 'GCgc') / len(segment) * 100
        gc_values.append(gc)
        idx += 1
        positions.append(idx)
    return positions, gc_values


def draw_local_gc(filename):
    """Local GC distribution line chart (120nt sliding window, step=24nt)."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharey=True)
    fig.suptitle('Local GC Content Distribution (120 nt sliding window, step = 24 nt)',
                 fontsize=13, fontweight='bold', y=1.0)

    samples = [
        (axes[0], dna_trad, 'Traditional Chinese (Chu Shi Biao)', '#2563EB'),
        (axes[1], dna_simp, 'Simplified Chinese (Beiying)', '#2563EB'),
    ]

    for ax, dna, title, color in samples:
        pos, gc = compute_local_gc_sliding(dna, window=120, step=24)

        # Show up to 60 points
        n_show = min(60, len(pos))
        pos = pos[:n_show]
        gc = gc[:n_show]

        # Ideal band
        ax.axhspan(45, 55, color='#10B981', alpha=0.10, zorder=0,
                   label='Ideal range (45-55%)')
        ax.axhline(50, color='#10B981', linestyle='--', linewidth=0.8,
                   alpha=0.5, zorder=1)

        # CMC line
        ax.plot(pos, gc, '-o', color=color, linewidth=1.5,
                markersize=3.5, markeredgewidth=0.4, markeredgecolor='white',
                zorder=3, label='CMC')

        # Fill between deviation
        ax.fill_between(pos, 50, gc, alpha=0.10, color=color, zorder=2)

        ax.set_xlabel('Sliding Window Position', fontweight='medium')
        ax.set_title(title, fontsize=11, pad=8)
        ax.set_ylim(36, 65)
        ax.set_xlim(0.5, n_show + 0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', linestyle=':', linewidth=0.4, alpha=0.5)
        ax.legend(loc='upper right', fontsize=8, framealpha=0.7)

        # Stats annotation
        mean_gc = np.mean(gc)
        std_gc = np.std(gc)
        ax.text(0.03, 0.06,
                'Mean: {:.1f}%\nStd: {:.1f}%\nRange: {:.1f}-{:.1f}%'.format(
                    mean_gc, std_gc, min(gc), max(gc)),
                transform=ax.transAxes, fontsize=8,
                verticalalignment='bottom',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                         edgecolor='#cccccc', alpha=0.85))

    axes[0].set_ylabel('GC Content (%)', fontweight='medium')

    fig.tight_layout()
    fig.savefig(OUT / filename, format='svg')
    fig.savefig(OUT / filename.replace('.svg', '.png'), format='png')
    plt.close(fig)
    print('  [OK] ' + filename)


# ════════════════════════════════════════════════════════════════════
# Chart 2: Encoding Efficiency Waterfall
# ════════════════════════════════════════════════════════════════════
def draw_waterfall(filename):
    """Waterfall chart showing data volume at each encoding stage."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    fig.suptitle('CMC Encoding Pipeline - Data Volume at Each Stage',
                 fontsize=14, fontweight='bold', y=1.01)

    datasets = [
        (axes[0], stages_trad, 'Traditional Chinese (Chu Shi Biao)'),
        (axes[1], stages_simp, 'Simplified Chinese (Beiying)'),
    ]

    stage_labels = [
        'Original\nText',
        'Wubi\nEncoding',
        'Brotli\nCompression',
        'DNA\nSequence'
    ]

    bar_colors = ['#64748B', '#3B82F6', '#10B981', '#2563EB']
    arrow_color = '#94A3B8'

    for ax, st, title in datasets:
        # Convert to bits for comparison
        text_bits = st['chars'] * 16           # Chinese chars @ 16 bits
        wubi_bits = st['wubi_len'] * 8         # ASCII chars @ 8 bits
        comp_bits = st['compressed_bytes'] * 8 # bytes @ 8 bits
        dna_bits  = st['dna_nt'] * 2           # nucleotides (encoding domain)

        values = [text_bits, wubi_bits, comp_bits, dna_bits]
        x = np.arange(len(stage_labels))

        bars = ax.bar(x, values, 0.55, color=bar_colors,
                      edgecolor=['#475569', '#1D4ED8', '#059669', '#1D4ED8'],
                      linewidth=0.8, zorder=3)

        # Value labels
        for i, (bar, v) in enumerate(zip(bars, values)):
            # Show both bits and human-readable
            if i == 0:
                label = '{:,} bits\n({} chars)'.format(v, st['chars'])
            elif i == 1:
                ratio = v / text_bits * 100
                label = '{:,} bits\n({} chars, {:.0f}%)'.format(v, st['wubi_len'], ratio)
            elif i == 2:
                ratio = v / text_bits * 100
                label = '{:,} bits\n({} B, {:.0f}%)'.format(v, st['compressed_bytes'], ratio)
            else:
                label = '{:,} bits\n({} nt)'.format(v, st['dna_nt'])

            fontw = 'bold' if i == 3 else 'normal'
            clr = '#1D4ED8' if i == 3 else '#333333'
            ax.text(bar.get_x() + bar.get_width()/2,
                    bar.get_height() + max(values) * 0.02,
                    label, ha='center', va='bottom', fontsize=7.5,
                    fontweight=fontw, color=clr, linespacing=1.3)

        # Arrows between bars
        for i in range(len(values) - 1):
            change = values[i+1] - values[i]
            pct = change / values[i] * 100
            sign = '+' if change > 0 else ''
            mid_y = max(values) * 0.02
            ax.annotate('{}{:.0f}%'.format(sign, pct),
                       xy=(i + 0.5, mid_y), fontsize=7.5,
                       ha='center', va='bottom', color='#EF4444' if change > 0 else '#10B981',
                       fontweight='bold')

        ax.set_xticks(x)
        ax.set_xticklabels(stage_labels, fontsize=9, fontweight='medium')
        ax.set_ylabel('Data Volume (bits)', fontweight='medium')
        ax.set_title(title, fontsize=11, pad=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', linestyle=':', linewidth=0.4, alpha=0.5, zorder=0)

        # Padding
        cur_ymin, cur_ymax = ax.get_ylim()
        ax.set_ylim(0, cur_ymax * 1.25)

        # Storage density annotation
        density = text_bits / st['dna_nt']
        ax.text(0.97, 0.97,
                'Density: {:.2f} bits/nt'.format(density),
                transform=ax.transAxes, fontsize=9, fontweight='bold',
                ha='right', va='top', color='#1D4ED8',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='#EFF6FF',
                         edgecolor='#93C5FD', alpha=0.9))

    fig.tight_layout()
    fig.savefig(OUT / filename, format='svg')
    fig.savefig(OUT / filename.replace('.svg', '.png'), format='png')
    plt.close(fig)
    print('  [OK] ' + filename)


# ════════════════════════════════════════════════════════════════════
# Chart 3: Hamming Distance Heatmap
# ════════════════════════════════════════════════════════════════════
def draw_heatmap(filename):
    """Pairwise Hamming distance heatmap for 6 fragments."""
    import matplotlib.gridspec as gridspec

    fig = plt.figure(figsize=(13, 6))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 0.05], wspace=0.35)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    cax = fig.add_subplot(gs[0, 2])

    fig.suptitle('Pairwise Hamming Distance Between 120 nt Fragments',
                 fontsize=14, fontweight='bold', y=1.02)

    datasets = [
        (ax1, hamming_trad, 'Traditional Chinese (Chu Shi Biao)', 76),
        (ax2, hamming_simp, 'Simplified Chinese (Beiying)', 80),
    ]

    frag_labels = ['F1', 'F2', 'F3', 'F4', 'F5', 'F6']
    im = None

    for ax, matrix, title, min_hd in datasets:
        # Mask diagonal
        mask = np.eye(6, dtype=bool)
        display = np.copy(matrix)
        display[mask] = np.nan

        im = ax.imshow(display, cmap='YlGnBu', vmin=70, vmax=105,
                       aspect='equal', interpolation='nearest')

        # Annotate cells
        for i in range(6):
            for j in range(6):
                if i == j:
                    ax.text(j, i, '-', ha='center', va='center',
                            fontsize=16, color='#94A3B8')
                else:
                    val = int(matrix[i, j])
                    is_min = (val == min_hd)
                    ax.text(j, i, str(val), ha='center', va='center',
                            fontsize=15,
                            fontweight='bold' if is_min else 'normal',
                            color='#DC2626' if is_min else 'black')

        ax.set_xticks(range(6))
        ax.set_yticks(range(6))
        ax.set_xticklabels(frag_labels, fontweight='medium')
        ax.set_yticklabels(frag_labels, fontweight='medium')
        ax.set_title(title + '\n(Min HD = {})'.format(min_hd),
                     fontsize=11, pad=10)

        # Grid
        for edge in np.arange(-0.5, 6, 1):
            ax.axhline(edge, color='white', linewidth=1.5)
            ax.axvline(edge, color='white', linewidth=1.5)

    # Colorbar in dedicated axis
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label('Hamming Distance', fontweight='medium')
    cbar.outline.set_edgecolor('#cccccc')

    fig.subplots_adjust(top=0.88, bottom=0.08, left=0.06, right=0.93)
    fig.savefig(OUT / filename, format='svg')
    fig.savefig(OUT / filename.replace('.svg', '.png'), format='png')
    plt.close(fig)
    print('  [OK] ' + filename)


# ── Main ──────────────────────────────────────────────────────────
if __name__ == '__main__':
    print('Generating supplementary charts ...')
    draw_local_gc('fig_local_gc_distribution.svg')
    draw_waterfall('fig_encoding_waterfall.svg')
    draw_heatmap('fig_hamming_heatmap.svg')
    print('Done.')
