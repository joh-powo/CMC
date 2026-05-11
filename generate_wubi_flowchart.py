#!/usr/bin/env python3
"""
Publication-quality Wubi Preprocessing Flowchart (English version)
for the CMC DNA storage paper.
Wider layout, generous whitespace, Nature/Science aesthetic.
"""

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = ['Arial', 'Microsoft YaHei', 'SimHei']
matplotlib.rcParams['axes.unicode_minus'] = False

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patheffects as pe

# ── Palette ──
P = {
    'bg':       '#FFFFFF',
    'inp_f':    '#DCEEFB',  'inp_e':    '#4A90D9',
    'proc_f':   '#FFF4E0',  'proc_e':   '#E8963E',
    'mod_f':    '#E2F4E6',  'mod_e':    '#4CAF6E',
    'out_f':    '#F0E4F7',  'out_e':    '#9B59B6',
    'arrow':    '#4A5568',
    'txt':      '#1A1A2E',
    'sub':      '#6B7B8D',
    'hi':       '#D93025',
    'ann_f':    '#FFF9E6',  'ann_e':    '#E0A800',
    'dat_f':    '#F0F4F8',  'dat_e':    '#8FA4B8',
    'strip_f':  '#F7F9FC',  'strip_e':  '#C4CDD6',
}


def _box(ax, x, y, w, h, title, subtitle=None,
         fc='#DCEEFB', ec='#4A90D9', fs=13, subfs=9.5,
         lw=1.5, bold=True, radius='round,pad=0.18'):
    patch = FancyBboxPatch((x - w/2, y - h/2), w, h,
                           boxstyle=radius, fc=fc, ec=ec, lw=lw, zorder=2)
    patch.set_path_effects([
        pe.withSimplePatchShadow(offset=(1, -1),
                                 shadow_rgbFace='#C8C8C8', alpha=0.20)])
    ax.add_patch(patch)
    wt = 'bold' if bold else 'normal'
    if subtitle:
        ax.text(x, y + 0.22, title, ha='center', va='center',
                fontsize=fs, fontweight=wt, color=P['txt'], zorder=3)
        ax.text(x, y - 0.22, subtitle, ha='center', va='center',
                fontsize=subfs, color=P['sub'], zorder=3, style='italic')
    else:
        ax.text(x, y, title, ha='center', va='center',
                fontsize=fs, fontweight=wt, color=P['txt'], zorder=3)


def _arr(ax, x1, y1, x2, y2, col=None, lw=1.6, cs='arc3,rad=0'):
    col = col or P['arrow']
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle='-|>',
        connectionstyle=cs, color=col, lw=lw,
        mutation_scale=15, zorder=1))


def _note(ax, x, y, w, h, t1, t2=None, fc=None, ec=None, ls='--'):
    fc, ec = fc or P['ann_f'], ec or P['ann_e']
    p = FancyBboxPatch((x - w/2, y - h/2), w, h,
                       boxstyle='round,pad=0.10',
                       fc=fc, ec=ec, lw=1.0, ls=ls, zorder=2)
    ax.add_patch(p)
    if t2:
        ax.text(x, y + 0.13, t1, ha='center', va='center',
                fontsize=9, fontweight='bold', color=P['txt'], zorder=3)
        ax.text(x, y - 0.13, t2, ha='center', va='center',
                fontsize=7.5, color=P['sub'], zorder=3, style='italic',
                fontfamily='monospace')
    else:
        ax.text(x, y, t1, ha='center', va='center',
                fontsize=9, color=P['txt'], zorder=3)


def _char_strip(ax, x, y, w, h, chars, highlights=None, fs=13):
    """Horizontal strip of example Chinese characters."""
    n = len(chars)
    iw = w / n
    highlights = highlights or {}
    for i, ch in enumerate(chars):
        ix = x - w/2 + iw * i + iw/2
        bg = highlights.get(i, P['strip_f'])
        p = FancyBboxPatch((ix - iw/2 + 0.03, y - h/2), iw - 0.06, h,
                           boxstyle='round,pad=0.04',
                           fc=bg, ec=P['strip_e'], lw=0.6, zorder=2)
        ax.add_patch(p)
        ax.text(ix, y, ch, ha='center', va='center',
                fontsize=fs, color=P['txt'], zorder=3)


def _code_strip(ax, x, y, w, h, codes, fs=12):
    """Horizontal strip of Wubi codes with delimiter highlighting."""
    n = len(codes)
    iw = w / n
    for i, c in enumerate(codes):
        ix = x - w/2 + iw * i + iw/2
        bg = '#E8DEF8' if c == '6' else '#EDE7F6'
        ec = P['out_e'] if c == '6' else '#C4B5D6'
        lw = 1.2 if c == '6' else 0.6
        p = FancyBboxPatch((ix - iw/2 + 0.02, y - h/2), iw - 0.04, h,
                           boxstyle='round,pad=0.04',
                           fc=bg, ec=ec, lw=lw, zorder=2)
        ax.add_patch(p)
        wt = 'bold' if c == '6' else 'normal'
        ax.text(ix, y, c, ha='center', va='center',
                fontsize=fs, color=P['txt'], zorder=3,
                fontweight=wt, fontfamily='monospace')


def main():
    fig, ax = plt.subplots(figsize=(12, 18))
    fig.patch.set_facecolor(P['bg'])
    ax.set_facecolor(P['bg'])
    ax.set_xlim(-6.5, 6.5)
    ax.set_ylim(-2, 20.5)
    ax.set_aspect('equal')
    ax.axis('off')

    # ─── Title ───
    ax.text(0, 20.0, 'Wubi Preprocessing Pipeline',
            ha='center', fontsize=22, fontweight='bold', color=P['txt'])
    ax.text(0, 19.45, 'Chinese Text → Wubi Codestream for CMC DNA Storage Encoding',
            ha='center', fontsize=12, color=P['sub'], style='italic')

    # ── Layout constants ──
    CX = 0
    BW = 5.0    # box width
    BH = 1.05   # box height
    GAP = 2.3   # vertical gap between box centers

    y = 18.0    # start y

    # =============== STEP 1: Input ===============
    _box(ax, CX, y, BW, BH,
         'Step 1: Chinese Text Input',
         'Accept Traditional or Simplified Chinese text',
         fc=P['inp_f'], ec=P['inp_e'])
    _char_strip(ax, CX, y - 0.80, 4.2, 0.45,
                ['先', '帝', '創', '業', '未', '半'])
    ax.text(CX + 2.6, y - 0.80, '(Traditional)',
            fontsize=8, color=P['sub'], ha='left', va='center', style='italic')

    _arr(ax, CX, y - 1.12, CX, y - GAP + BH/2 + 0.08)

    # =============== STEP 2: OpenCC ===============
    y -= GAP
    _box(ax, CX, y, BW, BH,
         'Step 2: Traditional → Simplified',
         'OpenCC conversion (t2s), triggered only for traditional input',
         fc=P['proc_f'], ec=P['proc_e'])

    # Side note: OpenCC
    nx = CX + 4.8
    _note(ax, nx, y + 0.2, 2.6, 0.58,
          'OpenCC Engine', 'conversion="t2s"')
    _arr(ax, nx - 1.3, y + 0.2, CX + BW/2, y + 0.05,
         col=P['ann_e'], lw=1.0, cs='arc3,rad=-0.12')

    _arr(ax, CX, y - BH/2 - 0.08, CX, y - GAP + BH/2 + 0.08)

    # =============== STEP 3: Simplified Text ===============
    y -= GAP
    _box(ax, CX, y, BW, BH,
         'Step 3: Simplified Chinese Text',
         'Unified simplified form for Wubi lookup',
         fc=P['inp_f'], ec=P['inp_e'])
    _char_strip(ax, CX, y - 0.80, 4.2, 0.45,
                ['先', '帝', '创', '业', '未', '半'],
                highlights={2: '#FFECB3', 3: '#FFECB3'})
    # Change highlight
    ax.annotate('', xy=(CX - 0.55, y - 1.10),
                xytext=(CX - 0.55, y - 0.83),
                arrowprops=dict(arrowstyle='->', color=P['hi'], lw=1.0))
    ax.text(CX - 0.30, y - 1.12, '創→创  業→业',
            fontsize=8.5, color=P['hi'], ha='left', va='center')

    _arr(ax, CX, y - 1.22, CX, y - GAP + BH/2 + 0.08)

    # =============== STEP 4: Wubi-86 Lookup ===============
    y -= GAP + 0.2
    _box(ax, CX, y, BW, BH,
         'Step 4: Wubi-86 Dictionary Lookup',
         'Map each character to its Wubi-86 code (a–z)',
         fc=P['mod_f'], ec=P['mod_e'])

    # Data source
    _note(ax, CX + 4.8, y, 2.6, 0.58,
          'character.csv', '~190K entries',
          fc=P['dat_f'], ec=P['dat_e'])
    _arr(ax, CX + 4.8 - 1.3, y, CX + BW/2, y,
         col=P['dat_e'], lw=1.0)

    _arr(ax, CX, y - BH/2 - 0.08, CX, y - GAP + BH/2 + 0.08)

    # =============== STEP 5: Punctuation Stripping ===============
    y -= GAP
    _box(ax, CX, y, BW, BH,
         'Step 5: Punctuation Stripping',
         'Remove punctuation; record positions for lossless recovery',
         fc=P['proc_f'], ec=P['proc_e'])

    _note(ax, CX - 4.8, y + 0.15, 2.6, 0.58,
          'punctMap', 'positions[] + chars{}')
    _arr(ax, CX - 4.8 + 1.3, y + 0.15, CX - BW/2, y + 0.05,
         col=P['ann_e'], lw=1.0, cs='arc3,rad=0.12')

    _arr(ax, CX, y - BH/2 - 0.08, CX, y - GAP + BH/2 + 0.08)

    # =============== STEP 6: Collision Disambiguation ===============
    y -= GAP
    _box(ax, CX, y, BW, BH,
         'Step 6: Collision Disambiguation',
         'Resolve one-code-multi-character ambiguities',
         fc=P['proc_f'], ec=P['proc_e'])

    _note(ax, CX - 4.8, y + 0.15, 2.6, 0.58,
          'collisionMap', 'positions[] + indices[]')
    _arr(ax, CX - 4.8 + 1.3, y + 0.15, CX - BW/2, y + 0.05,
         col=P['ann_e'], lw=1.0, cs='arc3,rad=0.12')

    # Example
    ax.text(CX + 3.2, y - 0.25,
            'e.g.  jghh → 量 / 里',
            fontsize=9, color=P['sub'], ha='left', va='center', style='italic')

    _arr(ax, CX, y - BH/2 - 0.08, CX, y - GAP + BH/2 + 0.08)

    # =============== STEP 7: Output ===============
    y -= GAP
    _box(ax, CX, y, BW, BH,
         'Step 7: Wubi Codestream Output',
         'Alphabet codes (a–z) separated by delimiter "6"',
         fc=P['out_f'], ec=P['out_e'])

    # Codestream example
    codes = list('tfhp') + ['6'] + list('jghh') + ['6'] + list('wtaj')
    _code_strip(ax, CX, y - 0.82, 5.2, 0.42, codes)
    ax.text(CX + 3.1, y - 0.82, '← delimiter "6"',
            fontsize=9, color=P['out_e'], ha='left', va='center',
            fontweight='bold')

    # Scope bracket: wubi_encoder.m (steps 4–6)
    bx = CX - BW/2 - 0.65
    y4 = 18.0 - GAP * 3 - 0.2   # step 4 y
    y6 = y4 - GAP * 2            # step 6 y
    bt = y4 + BH/2 + 0.15
    bb = y6 - BH/2 - 0.15
    ax.plot([bx, bx - 0.2, bx - 0.2, bx],
            [bt, bt, bb, bb],
            color=P['mod_e'], lw=2.0, zorder=1, solid_capstyle='round')
    ax.text(bx - 0.30, (bt + bb) / 2, 'wubi_encoder.m',
            fontsize=9.5, color=P['mod_e'], ha='right', va='center',
            rotation=90, fontfamily='monospace', fontweight='bold')

    # Downstream pointer
    _arr(ax, CX, y - 1.12, CX, y - 1.75)
    ax.text(CX, y - 2.05,
            '→  Brotli Compression  →  Constrained Mixed-Radix DNA Encoding',
            ha='center', fontsize=11, color=P['sub'], style='italic')

    # ── Save ──
    base = r'e:\研\李孜淇\Projects\CMRC\plots\fig_wubi_preprocessing'
    for ext in ('png', 'svg', 'pdf'):
        fig.savefig(f'{base}.{ext}', dpi=300, bbox_inches='tight',
                    facecolor=fig.get_facecolor())
        print(f'Saved: {base}.{ext}')
    plt.close(fig)


if __name__ == '__main__':
    main()
