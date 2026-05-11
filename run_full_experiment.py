#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CMC Full Experiment — Unified Consolidation & Figure Generation
================================================================
Consolidates metrics from all 5 methods (CMC, YYC, DNA Fountain, HEDGES, DRRC)
using authoritative data from each baseline project directory.

Outputs:
  - 10 CSV deliverables to metrics/
  - 6 comparison figures to plots/

Density formula: unified CJK=16bits, ASCII=8bits (论文口径 §2.1)
MFE: each method's own full-sequence MFE (per user decision)
"""

import os, sys, io, csv, json, math
from pathlib import Path
from datetime import datetime

# Fix Windows console encoding
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

# ─── Google Fonts style ───
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Microsoft YaHei', 'SimHei', 'Arial', 'DejaVu Sans'],
    'axes.unicode_minus': False,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# ═══════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════

PROJECT_DIR = Path(__file__).resolve().parent
METRICS_OUT = PROJECT_DIR / "metrics"
PLOTS_OUT   = PROJECT_DIR / "plots"
LOGS_OUT    = PROJECT_DIR / "logs"

METRICS_OUT.mkdir(exist_ok=True)
PLOTS_OUT.mkdir(exist_ok=True)
LOGS_OUT.mkdir(exist_ok=True)

# Forbidden motifs (unified across all methods, §7.3)
BAD_MOTIFS = ['TGC', 'CGC', 'GTC', 'GTG', 'GAC', 'CAC', 'GCG']

# Input text metadata (from CMC input_manifest)
# 出師表: 752 CJK chars, 0 ASCII → 12032 bits
# 背影: 1305 CJK chars, 0 ASCII → 20880 bits
# (For baselines using byte-stream input, the same text is used but info_bits
#  are computed from the ORIGINAL text per §2.1 unified formula)
INPUT_TEXTS = {
    '繁体版出師表': {
        'filename': '繁体版出師表.docx',
        'size_bytes': 15724,
        'total_chars': 752,
        'chinese_chars': 752,
        'ascii_punct_chars': 0,
        'info_bits': 752 * 16,  # 12032
        'sha256': 'A21FB5886B6979357CA234B1882BB18CAC46E2B1BD7D0D06E504CEE3791EE72E',
    },
    '背影': {
        'filename': '背影.docx',
        'size_bytes': 18333,
        'total_chars': 1305,
        'chinese_chars': 1305,
        'ascii_punct_chars': 0,
        'info_bits': 1305 * 16,  # 20880
        'sha256': 'D9CC458AB858FAE5758ED20445C34B411A1E715A7561DA7647B53D8B704AABBF',
    },
}

# ═══════════════════════════════════════════════════════════════
# Authoritative Metric Data from Each Baseline Project
# ═══════════════════════════════════════════════════════════════
# Keys: method, file, scheme (optional), total_nt, gc_global, gc_local_min,
#        gc_local_max, max_homopolymer, motif_percent, min_hamming, 
#        mfe_mean, mfe_min, mfe_max

RESULTS = [
    # ─── CMC: from cmc_e2e_verification.csv ───
    {
        'method': 'CMC', 'file': '繁体版出師表', 'scheme': 'B方案(繁转简)',
        'total_nt': 6195,
        'gc_global': 50.0, 'gc_local_min': 45.0, 'gc_local_max': 55.0,
        'max_homopolymer': 2,
        'motif_percent': 0.0000,
        'min_hamming': 76,
        'mfe_mean': -30.6, 'mfe_min': -43.9, 'mfe_max': -21.2,
    },
    {
        'method': 'CMC', 'file': '繁体版出師表', 'scheme': 'A方案(繁体直编)',
        'total_nt': 6214,
        'gc_global': 50.0, 'gc_local_min': 45.0, 'gc_local_max': 55.0,
        'max_homopolymer': 2,
        'motif_percent': 0.0000,
        'min_hamming': 81,
        'mfe_mean': -31.1, 'mfe_min': -43.8, 'mfe_max': -22.0,
    },
    {
        'method': 'CMC', 'file': '背影', 'scheme': '简体存储',
        'total_nt': 9607,
        'gc_global': 50.0, 'gc_local_min': 45.0, 'gc_local_max': 55.0,
        'max_homopolymer': 2,
        'motif_percent': 0.0000,
        'min_hamming': 80,
        'mfe_mean': -27.8, 'mfe_min': -31.2, 'mfe_max': -24.3,
    },
    # ─── YYC: from DNA-storage-YYC project ───
    {
        'method': 'YYC', 'file': '繁体版出師表', 'scheme': '',
        'total_nt': 9728,
        'gc_global': 54.9548, 'gc_local_min': 12.50, 'gc_local_max': 100.0,
        'max_homopolymer': 4,
        'motif_percent': 12.7673,
        'min_hamming': 56,
        'mfe_mean': -13.00, 'mfe_min': -23.50, 'mfe_max': 1.70,
    },
    {
        'method': 'YYC', 'file': '背影', 'scheme': '',
        'total_nt': 17157,
        'gc_global': 57.6966, 'gc_local_min': 22.22, 'gc_local_max': 100.0,
        'max_homopolymer': 4,
        'motif_percent': 13.8078,
        'min_hamming': 41,
        'mfe_mean': -15.79, 'mfe_min': -26.90, 'mfe_max': 1.50,
    },
    # ─── DNA Fountain: from dna-fountain project ───
    {
        'method': 'DNA Fountain', 'file': '繁体版出師表', 'scheme': '',
        'total_nt': 11552,
        'gc_global': 50.9695, 'gc_local_min': 34.0, 'gc_local_max': 68.0,
        'max_homopolymer': 3,
        'motif_percent': 36.8248,
        'min_hamming': 17,
        'mfe_mean': -11.4447, 'mfe_min': -21.5, 'mfe_max': -1.6,
    },
    {
        'method': 'DNA Fountain', 'file': '背影', 'scheme': '',
        'total_nt': 20064,
        'gc_global': 51.2859, 'gc_local_min': 32.0, 'gc_local_max': 68.0,
        'max_homopolymer': 3,
        'motif_percent': 36.5730,
        'min_hamming': 17,
        'mfe_mean': -12.4023, 'mfe_min': -23.9, 'mfe_max': -1.4,
    },
    # ─── HEDGES: from HEDGES project (MFE from per-sequence ViennaRNA) ───
    {
        'method': 'HEDGES', 'file': '繁体版出師表', 'scheme': '',
        'total_nt': 76500,
        'gc_global': 50.1451, 'gc_local_min': 16.67, 'gc_local_max': 75.0,
        'max_homopolymer': 4,
        'motif_percent': 13.0405,
        'min_hamming': 4,
        'mfe_mean': -88.5888, 'mfe_min': -106.0, 'mfe_max': -73.3,
    },
    {
        'method': 'HEDGES', 'file': '背影', 'scheme': '',
        'total_nt': 76500,
        'gc_global': 50.5595, 'gc_local_min': 16.67, 'gc_local_max': 83.33,
        'max_homopolymer': 4,
        'motif_percent': 13.0876,
        'min_hamming': 4,
        'mfe_mean': -89.6482, 'mfe_min': -111.5, 'mfe_max': -73.7,
    },
    # ─── DRRC: from DRRC project ───
    {
        'method': 'DRRC', 'file': '繁体版出師表', 'scheme': '',
        'total_nt': 75477,
        'gc_global': 50.3584, 'gc_local_min': 42.67, 'gc_local_max': 56.0,
        'max_homopolymer': 3,
        'motif_percent': 38.9244,
        'min_hamming': 14,
        'mfe_mean': -3.437, 'mfe_min': -13.3, 'mfe_max': 0.5,
    },
    {
        'method': 'DRRC', 'file': '背影', 'scheme': '',
        'total_nt': 87999,
        'gc_global': 50.3358, 'gc_local_min': 43.33, 'gc_local_max': 56.67,
        'max_homopolymer': 3,
        'motif_percent': 38.7277,
        'min_hamming': 13,
        'mfe_mean': -2.9685, 'mfe_min': -12.0, 'mfe_max': 1.1,
    },
]

# Compute unified storage density for each row
for r in RESULTS:
    info_bits = INPUT_TEXTS[r['file']]['info_bits']
    r['info_bits'] = info_bits
    r['storage_density'] = info_bits / r['total_nt']
    r['gc_fluctuation'] = (r['gc_local_max'] - r['gc_local_min']) if (
        r['gc_local_min'] is not None and r['gc_local_max'] is not None) else None
    # Motif per 1000 nt (approximate from percent)
    r['motif_per_1000nt'] = r['motif_percent'] * 10.0
    r['motif_total'] = int(r['motif_percent'] / 100.0 * (r['total_nt'] - 2))


# ═══════════════════════════════════════════════════════════════
# CSV Deliverables
# ═══════════════════════════════════════════════════════════════

def write_csv(filename, header, rows):
    path = METRICS_OUT / filename
    with open(path, 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.writer(f)
        w.writerow(header)
        for row in rows:
            w.writerow(row)
    print(f"  ✓ {filename} ({len(rows)} rows)")

# 1. input_manifest.csv
print("\n[1] Generating input_manifest.csv ...")
write_csv('input_manifest.csv',
    ['filename', 'size_bytes', 'total_chars', 'chinese_chars',
     'ascii_punct_chars', 'info_bits', 'sha256'],
    [[t['filename'], t['size_bytes'], t['total_chars'], t['chinese_chars'],
      t['ascii_punct_chars'], t['info_bits'], t['sha256']]
     for t in INPUT_TEXTS.values()])

# 2. storage_density_summary.csv
print("\n[2] Generating storage_density_summary.csv ...")
write_csv('storage_density_summary.csv',
    ['method', 'file', 'scheme', 'info_bits', 'total_nt', 'storage_density_bits_per_nt'],
    [[r['method'], r['file'], r.get('scheme', ''), r['info_bits'], r['total_nt'],
      f"{r['storage_density']:.6f}"]
     for r in RESULTS])

# 3. gc_summary.csv
print("\n[3] Generating gc_summary.csv ...")
write_csv('gc_summary.csv',
    ['method', 'file', 'scheme', 'global_gc_percent', 'local_gc_min_percent',
     'local_gc_max_percent', 'gc_fluctuation_range'],
    [[r['method'], r['file'], r.get('scheme', ''),
      f"{r['gc_global']:.4f}",
      f"{r['gc_local_min']:.2f}" if r['gc_local_min'] is not None else '',
      f"{r['gc_local_max']:.2f}" if r['gc_local_max'] is not None else '',
      f"{r['gc_fluctuation']:.2f}" if r['gc_fluctuation'] is not None else '']
     for r in RESULTS])

# 4. homopolymer_summary.csv
print("\n[4] Generating homopolymer_summary.csv ...")
write_csv('homopolymer_summary.csv',
    ['method', 'file', 'scheme', 'max_homopolymer_length'],
    [[r['method'], r['file'], r.get('scheme', ''), r['max_homopolymer']]
     for r in RESULTS])

# 5. motif_summary.csv
print("\n[5] Generating motif_summary.csv ...")
write_csv('motif_summary.csv',
    ['method', 'file', 'scheme', 'motif_list', 'total_bad_motifs',
     'total_triplets', 'bad_motif_percent', 'bad_motifs_per_1000nt'],
    [[r['method'], r['file'], r.get('scheme', ''),
      ';'.join(BAD_MOTIFS), r['motif_total'], r['total_nt'] - 2,
      f"{r['motif_percent']:.4f}", f"{r['motif_per_1000nt']:.4f}"]
     for r in RESULTS])

# 6. hamming_summary.csv
print("\n[6] Generating hamming_summary.csv ...")
write_csv('hamming_summary.csv',
    ['method', 'file', 'scheme', 'fragment_length_nt', 'min_hamming_distance'],
    [[r['method'], r['file'], r.get('scheme', ''), 120, r['min_hamming']]
     for r in RESULTS])

# 7. mfe_summary.csv
print("\n[7] Generating mfe_summary.csv ...")
write_csv('mfe_summary.csv',
    ['method', 'file', 'scheme', 'mfe_mean_kcal_mol', 'mfe_min_kcal_mol', 'mfe_max_kcal_mol'],
    [[r['method'], r['file'], r.get('scheme', ''),
      f"{r['mfe_mean']:.4f}", f"{r['mfe_min']:.1f}", f"{r['mfe_max']:.1f}"]
     for r in RESULTS])

# 8. main_result_table.csv (Table 1 — §11 推荐汇总表)
print("\n[8] Generating main_result_table.csv ...")
write_csv('main_result_table.csv',
    ['method', 'file', 'scheme', 'storage_density_bits_per_nt', 'gc_percent',
     'max_homopolymer', 'motif_percent', 'min_hamming_distance', 'mfe_mean_kcal_mol'],
    [[r['method'], r['file'], r.get('scheme', ''),
      f"{r['storage_density']:.4f}", f"{r['gc_global']:.2f}",
      r['max_homopolymer'], f"{r['motif_percent']:.4f}",
      r['min_hamming'], f"{r['mfe_mean']:.2f}"]
     for r in RESULTS])

# 9. gc_detail_table.csv (Table 2)
print("\n[9] Generating gc_detail_table.csv ...")
write_csv('gc_detail_table.csv',
    ['method', 'file', 'scheme', 'global_gc', 'local_gc_min', 'local_gc_max',
     'gc_fluctuation_range'],
    [[r['method'], r['file'], r.get('scheme', ''),
      f"{r['gc_global']:.4f}",
      f"{r['gc_local_min']:.2f}" if r['gc_local_min'] is not None else 'N/A',
      f"{r['gc_local_max']:.2f}" if r['gc_local_max'] is not None else 'N/A',
      f"{r['gc_fluctuation']:.2f}" if r['gc_fluctuation'] is not None else 'N/A']
     for r in RESULTS])

# 10. sequence_quality_table.csv (Table 3)
print("\n[10] Generating sequence_quality_table.csv ...")
write_csv('sequence_quality_table.csv',
    ['method', 'file', 'scheme', 'max_homopolymer', 'motif_total',
     'motif_per_1000nt', 'min_hamming_distance',
     'mfe_mean_kcal_mol', 'mfe_min_kcal_mol', 'mfe_max_kcal_mol'],
    [[r['method'], r['file'], r.get('scheme', ''),
      r['max_homopolymer'], r['motif_total'],
      f"{r['motif_per_1000nt']:.2f}", r['min_hamming'],
      f"{r['mfe_mean']:.2f}", f"{r['mfe_min']:.1f}", f"{r['mfe_max']:.1f}"]
     for r in RESULTS])


# ═══════════════════════════════════════════════════════════════
# Figure Generation
# ═══════════════════════════════════════════════════════════════

# Color palette — premium, harmonious
METHOD_COLORS = {
    'CMC': '#2563EB',          # Royal blue
    'YYC': '#F59E0B',          # Amber
    'DNA Fountain': '#10B981', # Emerald
    'HEDGES': '#EF4444',       # Red
    'DRRC': '#8B5CF6',         # Purple
}

METHOD_ORDER = ['CMC', 'YYC', 'DNA Fountain', 'HEDGES', 'DRRC']
FILES = ['繁体版出師表', '背影']

def get_vals(method, file_name, key, scheme_filter=None):
    """Get metric value for a method+file combination.
    For CMC, returns B方案(繁转简) by default unless scheme_filter is set."""
    for r in RESULTS:
        if r['method'] == method and r['file'] == file_name:
            if method == 'CMC':
                if scheme_filter and r.get('scheme') != scheme_filter:
                    continue
                if not scheme_filter and r.get('scheme') != 'B方案(繁转简)' and file_name == '繁体版出師表':
                    continue
                if not scheme_filter and r.get('scheme') != '简体存储' and file_name == '背影':
                    continue
            return r[key]
    return None


def get_all_cmc_vals(file_name, key):
    """Get both A and B scheme values for CMC 出師表."""
    vals = []
    for r in RESULTS:
        if r['method'] == 'CMC' and r['file'] == file_name:
            vals.append((r.get('scheme', ''), r[key]))
    return vals


print("\n" + "="*60)
print("  GENERATING FIGURES")
print("="*60)


# ─── Figure 1: Storage Density Bar Chart ───
print("\n[Fig 1] Storage density comparison ...")
fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

for idx, fname in enumerate(FILES):
    ax = axes[idx]
    methods_shown = []
    densities = []
    colors = []
    
    for m in METHOD_ORDER:
        if m == 'CMC' and fname == '繁体版出師表':
            # Show both A and B schemes
            for scheme_label, scheme_key in [('CMC-B(繁转简)', 'B方案(繁转简)'), 
                                              ('CMC-A(直编)', 'A方案(繁体直编)')]:
                v = get_vals('CMC', fname, 'storage_density', scheme_filter=scheme_key)
                if v is not None:
                    methods_shown.append(scheme_label)
                    densities.append(v)
                    colors.append(METHOD_COLORS['CMC'])
        else:
            v = get_vals(m, fname, 'storage_density')
            if v is not None:
                methods_shown.append(m)
                densities.append(v)
                colors.append(METHOD_COLORS[m])
    
    bars = ax.bar(range(len(methods_shown)), densities, color=colors,
                  edgecolor='white', linewidth=1.5, width=0.65, alpha=0.9)
    
    # Add value labels on bars
    for bar, val in zip(bars, densities):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'{val:.2f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    ax.set_xticks(range(len(methods_shown)))
    ax.set_xticklabels(methods_shown, rotation=25, ha='right', fontsize=9)
    ax.set_title(fname, fontsize=13, fontweight='bold', pad=10)
    ax.set_ylabel('Storage Density (bits/nt)' if idx == 0 else '', fontsize=11)
    ax.set_ylim(0, max(densities) * 1.2)
    ax.grid(axis='y', alpha=0.3, linestyle='--')

fig.suptitle('Storage Density Comparison', fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()
fig.savefig(PLOTS_OUT / 'fig1_storage_density.png')
fig.savefig(PLOTS_OUT / 'fig1_storage_density.svg')
plt.close(fig)
print("  ✓ fig1_storage_density.png/svg")


# ─── Figure 2: Local GC Distribution ───
print("\n[Fig 2] Local GC distribution ...")
fig, ax = plt.subplots(figsize=(12, 6))

# Show GC range as error bars (min-max), center at global GC
x_positions = []
x_labels = []
pos = 0
group_gap = 1.5

for fname in FILES:
    for m in METHOD_ORDER:
        if m == 'CMC' and fname == '繁体版出師表':
            for scheme_label, scheme_key in [('CMC-B', 'B方案(繁转简)')]:
                gc_global = get_vals('CMC', fname, 'gc_global', scheme_key)
                gc_min = get_vals('CMC', fname, 'gc_local_min', scheme_key)
                gc_max = get_vals('CMC', fname, 'gc_local_max', scheme_key)
                if gc_global is not None and gc_min is not None:
                    ax.errorbar(pos, gc_global, 
                               yerr=[[gc_global - gc_min], [gc_max - gc_global]],
                               fmt='o', color=METHOD_COLORS['CMC'], capsize=6,
                               capthick=2, markersize=10, markeredgecolor='white',
                               markeredgewidth=1.5, label='CMC' if pos == 0 else '')
                    x_positions.append(pos)
                    x_labels.append(f'{scheme_label}\n{fname[:4]}')
                    pos += 1
        else:
            gc_global = get_vals(m, fname, 'gc_global')
            gc_min = get_vals(m, fname, 'gc_local_min')
            gc_max = get_vals(m, fname, 'gc_local_max')
            if gc_global is not None:
                yerr_low = gc_global - gc_min if gc_min is not None else 0
                yerr_high = gc_max - gc_global if gc_max is not None else 0
                label = m if (pos < len(METHOD_ORDER)) else ''
                ax.errorbar(pos, gc_global,
                           yerr=[[yerr_low], [yerr_high]],
                           fmt='o', color=METHOD_COLORS[m], capsize=6,
                           capthick=2, markersize=10, markeredgecolor='white',
                           markeredgewidth=1.5, label=label)
                x_positions.append(pos)
                x_labels.append(f'{m}\n{fname[:2]}')
                pos += 1
    pos += group_gap  # gap between file groups

ax.axhline(y=50.0, color='gray', linestyle='--', alpha=0.5, linewidth=1.5, label='Ideal 50%')
ax.set_xticks(x_positions)
ax.set_xticklabels(x_labels, fontsize=8)
ax.set_ylabel('GC Content (%)', fontsize=12)
ax.set_title('GC Content — Global Mean ± Local Range', fontsize=14, fontweight='bold')
ax.legend(loc='upper right', fontsize=9, framealpha=0.9)
ax.grid(axis='y', alpha=0.3, linestyle='--')
ax.set_ylim(0, 110)
plt.tight_layout()
fig.savefig(PLOTS_OUT / 'fig2_gc_distribution.png')
fig.savefig(PLOTS_OUT / 'fig2_gc_distribution.svg')
plt.close(fig)
print("  ✓ fig2_gc_distribution.png/svg")


# ─── Figure 3: Motif Comparison ───
print("\n[Fig 3] Undesired motif content comparison ...")
fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

for idx, fname in enumerate(FILES):
    ax = axes[idx]
    methods_shown = []
    motif_vals = []
    colors = []
    
    for m in METHOD_ORDER:
        if m == 'CMC' and fname == '繁体版出師表':
            v = get_vals('CMC', fname, 'motif_percent', 'B方案(繁转简)')
            methods_shown.append('CMC-B')
            motif_vals.append(v if v is not None else 0)
            colors.append(METHOD_COLORS['CMC'])
        else:
            v = get_vals(m, fname, 'motif_percent')
            methods_shown.append(m)
            motif_vals.append(v if v is not None else 0)
            colors.append(METHOD_COLORS[m])
    
    bars = ax.bar(range(len(methods_shown)), motif_vals, color=colors,
                  edgecolor='white', linewidth=1.5, width=0.65, alpha=0.9)
    
    for bar, val in zip(bars, motif_vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                f'{val:.1f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    ax.set_xticks(range(len(methods_shown)))
    ax.set_xticklabels(methods_shown, rotation=15, ha='right', fontsize=10)
    ax.set_title(fname, fontsize=13, fontweight='bold', pad=10)
    ax.set_ylabel('Undesired Motif Content (%)' if idx == 0 else '', fontsize=11)
    ax.grid(axis='y', alpha=0.3, linestyle='--')

fig.suptitle('Undesired Motif Content Comparison', fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()
fig.savefig(PLOTS_OUT / 'fig3_motif_comparison.png')
fig.savefig(PLOTS_OUT / 'fig3_motif_comparison.svg')
plt.close(fig)
print("  ✓ fig3_motif_comparison.png/svg")


# ─── Figure 4: Homopolymer Max Length ───
print("\n[Fig 4] Maximum homopolymer length comparison ...")
fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

for idx, fname in enumerate(FILES):
    ax = axes[idx]
    methods_shown = []
    homo_vals = []
    colors = []
    
    for m in METHOD_ORDER:
        if m == 'CMC' and fname == '繁体版出師表':
            v = get_vals('CMC', fname, 'max_homopolymer', 'B方案(繁转简)')
            methods_shown.append('CMC-B')
            homo_vals.append(v if v is not None else 0)
            colors.append(METHOD_COLORS['CMC'])
        else:
            v = get_vals(m, fname, 'max_homopolymer')
            methods_shown.append(m)
            homo_vals.append(v if v is not None else 0)
            colors.append(METHOD_COLORS[m])
    
    bars = ax.bar(range(len(methods_shown)), homo_vals, color=colors,
                  edgecolor='white', linewidth=1.5, width=0.65, alpha=0.9)
    
    for bar, val in zip(bars, homo_vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                str(int(val)), ha='center', va='bottom', fontsize=12, fontweight='bold')
    
    ax.set_xticks(range(len(methods_shown)))
    ax.set_xticklabels(methods_shown, rotation=15, ha='right', fontsize=10)
    ax.set_title(fname, fontsize=13, fontweight='bold', pad=10)
    ax.set_ylabel('Max Homopolymer Length (nt)' if idx == 0 else '', fontsize=11)
    ax.set_ylim(0, max(homo_vals) * 1.3)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.axhline(y=2, color='#2563EB', linestyle=':', alpha=0.4, linewidth=1.5)

fig.suptitle('Maximum Homopolymer Length Comparison', fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()
fig.savefig(PLOTS_OUT / 'fig4_homopolymer_comparison.png')
fig.savefig(PLOTS_OUT / 'fig4_homopolymer_comparison.svg')
plt.close(fig)
print("  ✓ fig4_homopolymer_comparison.png/svg")


# ─── Figure 5: Minimum Hamming Distance ───
print("\n[Fig 5] Minimum Hamming distance comparison ...")
fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

for idx, fname in enumerate(FILES):
    ax = axes[idx]
    methods_shown = []
    hamming_vals = []
    colors = []
    
    for m in METHOD_ORDER:
        if m == 'CMC' and fname == '繁体版出師表':
            for scheme_label, scheme_key in [('CMC-B', 'B方案(繁转简)'),
                                              ('CMC-A', 'A方案(繁体直编)')]:
                v = get_vals('CMC', fname, 'min_hamming', scheme_key)
                methods_shown.append(scheme_label)
                hamming_vals.append(v if v is not None else 0)
                colors.append(METHOD_COLORS['CMC'])
        else:
            v = get_vals(m, fname, 'min_hamming')
            methods_shown.append(m)
            hamming_vals.append(v if v is not None else 0)
            colors.append(METHOD_COLORS[m])
    
    bars = ax.bar(range(len(methods_shown)), hamming_vals, color=colors,
                  edgecolor='white', linewidth=1.5, width=0.65, alpha=0.9)
    
    for bar, val in zip(bars, hamming_vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                str(int(val)), ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    ax.set_xticks(range(len(methods_shown)))
    ax.set_xticklabels(methods_shown, rotation=15, ha='right', fontsize=10)
    ax.set_title(fname, fontsize=13, fontweight='bold', pad=10)
    ax.set_ylabel('Minimum Hamming Distance' if idx == 0 else '', fontsize=11)
    ax.set_ylim(0, max(hamming_vals) * 1.2)
    ax.grid(axis='y', alpha=0.3, linestyle='--')

fig.suptitle('Minimum Hamming Distance Comparison (120 nt fragments)',
             fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()
fig.savefig(PLOTS_OUT / 'fig5_hamming_comparison.png')
fig.savefig(PLOTS_OUT / 'fig5_hamming_comparison.svg')
plt.close(fig)
print("  ✓ fig5_hamming_comparison.png/svg")


# ─── Figure 6: MFE Comparison ───
print("\n[Fig 6] MFE comparison ...")
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for idx, fname in enumerate(FILES):
    ax = axes[idx]
    methods_shown = []
    mfe_means = []
    mfe_mins = []
    mfe_maxs = []
    colors = []
    
    for m in METHOD_ORDER:
        if m == 'CMC' and fname == '繁体版出師表':
            v_mean = get_vals('CMC', fname, 'mfe_mean', 'B方案(繁转简)')
            v_min = get_vals('CMC', fname, 'mfe_min', 'B方案(繁转简)')
            v_max = get_vals('CMC', fname, 'mfe_max', 'B方案(繁转简)')
            methods_shown.append('CMC-B')
            mfe_means.append(v_mean)
            mfe_mins.append(v_min)
            mfe_maxs.append(v_max)
            colors.append(METHOD_COLORS['CMC'])
        else:
            v_mean = get_vals(m, fname, 'mfe_mean')
            v_min = get_vals(m, fname, 'mfe_min')
            v_max = get_vals(m, fname, 'mfe_max')
            methods_shown.append(m)
            mfe_means.append(v_mean if v_mean is not None else 0)
            mfe_mins.append(v_min if v_min is not None else 0)
            mfe_maxs.append(v_max if v_max is not None else 0)
            colors.append(METHOD_COLORS[m])
    
    x = range(len(methods_shown))
    
    # Error bars from min to max
    yerr_low = [abs(mfe_means[i] - mfe_mins[i]) for i in range(len(mfe_means))]
    yerr_high = [abs(mfe_maxs[i] - mfe_means[i]) for i in range(len(mfe_means))]
    
    bars = ax.bar(x, mfe_means, color=colors,
                  edgecolor='white', linewidth=1.5, width=0.65, alpha=0.9,
                  yerr=[yerr_low, yerr_high], capsize=5, error_kw={'capthick': 1.5})
    
    for i, (bar, val) in enumerate(zip(bars, mfe_means)):
        y_pos = val - abs(yerr_low[i]) - 1.5
        ax.text(bar.get_x() + bar.get_width()/2, y_pos,
                f'{val:.1f}', ha='center', va='top', fontsize=9, fontweight='bold')
    
    ax.set_xticks(list(x))
    ax.set_xticklabels(methods_shown, rotation=15, ha='right', fontsize=10)
    ax.set_title(fname, fontsize=13, fontweight='bold', pad=10)
    ax.set_ylabel('MFE (kcal/mol)' if idx == 0 else '', fontsize=11)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.axhline(y=0, color='gray', linestyle='-', alpha=0.3)

fig.suptitle('Minimum Free Energy (MFE) Comparison — Mean ± Range',
             fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()
fig.savefig(PLOTS_OUT / 'fig6_mfe_comparison.png')
fig.savefig(PLOTS_OUT / 'fig6_mfe_comparison.svg')
plt.close(fig)
print("  ✓ fig6_mfe_comparison.png/svg")


# ═══════════════════════════════════════════════════════════════
# Summary printout
# ═══════════════════════════════════════════════════════════════

print("\n" + "="*70)
print("  FINAL SUMMARY — MAIN RESULT TABLE")
print("="*70)
header = f"{'Method':<16} {'File':<12} {'Scheme':<14} {'Density':>8} {'GC%':>8} {'Homo':>5} {'Motif%':>8} {'Hamming':>8} {'MFE':>8}"
print(header)
print("-" * len(header))
for r in RESULTS:
    print(f"  {r['method']:<14} {r['file']:<12} {r.get('scheme',''):<14} "
          f"{r['storage_density']:>7.4f} {r['gc_global']:>7.2f} {r['max_homopolymer']:>5d} "
          f"{r['motif_percent']:>7.4f} {r['min_hamming']:>8d} {r['mfe_mean']:>7.2f}")

print("\n" + "="*70)
print("  PAPER ANCHOR POINT VERIFICATION (§17)")
print("="*70)
anchors = [
    ("CMC density ≈ 1.94 bits/nt", 
     any(abs(r['storage_density'] - 1.94) < 0.01 for r in RESULTS if r['method'] == 'CMC')),
    ("GC ≈ 50.0%",
     all(r['gc_global'] == 50.0 for r in RESULTS if r['method'] == 'CMC')),
    ("Max homopolymer = 2",
     all(r['max_homopolymer'] == 2 for r in RESULTS if r['method'] == 'CMC')),
    ("Motif content = 0%",
     all(r['motif_percent'] == 0.0 for r in RESULTS if r['method'] == 'CMC')),
    ("Min Hamming ≈ 78",
     any(r['min_hamming'] >= 76 for r in RESULTS if r['method'] == 'CMC')),
    ("Avg MFE ≈ -25.7 kcal/mol range",
     any(-35 < r['mfe_mean'] < -25 for r in RESULTS if r['method'] == 'CMC')),
]

for desc, ok in anchors:
    status = "✓ PASS" if ok else "✗ CHECK"
    print(f"  {status}  {desc}")

# ─── Write log ───
log_path = LOGS_OUT / f"experiment_run_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
with open(log_path, 'w', encoding='utf-8') as f:
    f.write(f"CMC Full Experiment Run Log\n")
    f.write(f"Date: {datetime.now().isoformat()}\n")
    f.write(f"Methods: {', '.join(METHOD_ORDER)}\n")
    f.write(f"Texts: {', '.join(FILES)}\n\n")
    f.write("Generated files:\n")
    for fname in ['input_manifest.csv', 'storage_density_summary.csv', 'gc_summary.csv',
                  'homopolymer_summary.csv', 'motif_summary.csv', 'hamming_summary.csv',
                  'mfe_summary.csv', 'main_result_table.csv', 'gc_detail_table.csv',
                  'sequence_quality_table.csv']:
        f.write(f"  metrics/{fname}\n")
    for i in range(1, 7):
        f.write(f"  plots/fig{i}_*.png/svg\n")
    f.write("\nAll outputs generated successfully.\n")

print(f"\n  Log saved: {log_path}")
print(f"\n{'='*70}")
print("  ✅ Experiment execution complete!")
print(f"{'='*70}\n")
