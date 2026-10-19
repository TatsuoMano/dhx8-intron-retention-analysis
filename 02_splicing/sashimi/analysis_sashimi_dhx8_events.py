#!/usr/bin/env python3
"""
Sashimi-style plots for DHX8-knockdown splicing events (Supplementary Fig. S8A-E).

Reconstruction of the LOST dedicated DHX8 sashimi script. Custom coverage
(samtools depth) + junction arcs (STAR SJ.out.tab), matplotlib Arc/FancyArrow.
Redrawn from v2 BAMs (siCont_2 corrected 2026-07-06) for provenance consistency.

Panels (mm10 / GRCm38):
  S8A  Smarcc2  A5SS  (chr10:+)
  S8B  Sf1      A3SS  (chr19:+)
  S8C  Dpp7     SE    (chr2:-)   inclusion vs skipping
  S8D  Rpl4     IR    (chr9:+)   retained intron
  S8E  Ccdc71   CE?   (chr9:+)   MAJIQ-called 175bp putative CE (actually IR)

Coverage is mean RPKM-scaled read depth per group (siControl n=3, siDHX8 n=3).
Junction numbers = group-mean STAR unique-junction reads at the relevant donor/acceptor.
"""

import subprocess
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Rectangle
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')

# ---- Arial / editable-PDF settings (project standard) --------------------
matplotlib.rcParams['font.family'] = ['Arial', 'Helvetica', 'Liberation Sans',
                                      'Nimbus Sans', 'DejaVu Sans']
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['axes.unicode_minus'] = False

SAMTOOLS = '/usr/bin/samtools'
BAM_DIR = '/mnt/HDD2/Sakai_RNA-seq/output/mapping'
OUT_DIR = '/mnt/HDD2/unified_results/sashimi_plots'

CTRL = ['siCont_1', 'siCont_2', 'siCont_3']   # corrected siCont_2 (v2)
KD   = ['siDHX8_1', 'siDHX8_2', 'siDHX8_3']

CTRL_COLOR = '#4C72B0'
KD_COLOR   = '#C44E52'

def bam(sample):
    return f'{BAM_DIR}/{sample}_Aligned.sortedByCoord.out.bam'

def sj(sample):
    return f'{BAM_DIR}/{sample}_SJ.out.tab'

# --------------------------------------------------------------------------
def get_coverage(sample, chrom, start, end):
    """Per-base depth over [start,end] (1-based, inclusive) via samtools depth."""
    region = f'{chrom}:{start}-{end}'
    cmd = [SAMTOOLS, 'depth', '-a', '-r', region, bam(sample)]
    cov = np.zeros(end - start + 1)
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    for line in res.stdout.strip().split('\n'):
        if not line:
            continue
        f = line.split('\t')
        pos = int(f[1]) - start
        if 0 <= pos < len(cov):
            cov[pos] = int(f[2])
    return cov

def get_junctions(sample, chrom, start, end):
    """Return dict {(intron_start, intron_end): unique_reads} within window.
    STAR SJ.out.tab: col1 chrom, col2 intron start (1-based), col3 intron end,
    col7 uniquely-mapping reads across the junction."""
    out = {}
    try:
        with open(sj(sample)) as fh:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) < 7 or f[0] != chrom:
                    continue
                s, e, u = int(f[1]), int(f[2]), int(f[6])
                if u <= 0:
                    continue
                if s >= start and e <= end:
                    out[(s, e)] = u
    except FileNotFoundError:
        pass
    return out

def group_coverage(samples, chrom, start, end):
    return np.mean([get_coverage(s, chrom, start, end) for s in samples], axis=0)

def junc_mean(samples, chrom, jstart, jend, tol=3):
    """Group-mean reads for a junction with donor/acceptor near (jstart,jend)."""
    vals = []
    for s in samples:
        jj = get_junctions(s, chrom, jstart - 10000, jend + 10000)
        best = 0
        for (a, b), u in jj.items():
            if abs(a - jstart) <= tol and abs(b - jend) <= tol:
                best = max(best, u)
        vals.append(best)
    return float(np.mean(vals)) if vals else 0.0

def draw_arc(ax, x0, x1, ymax, count, color, lw=2.5, up=True, label=True,
             frac=0.55, fs=8):
    mid = (x0 + x1) / 2.0
    width = abs(x1 - x0)
    h = ymax * frac
    base = 0.0
    arc = Arc((mid, base), width, 2 * h, angle=0, theta1=0, theta2=180,
              color=color, lw=lw, zorder=6)
    ax.add_patch(arc)
    if label:
        ax.text(mid, base + h * 0.96, f'{count:.0f}', ha='center', va='bottom',
                fontsize=fs, color=color, fontweight='bold', zorder=7)

def style_cov_ax(ax, cov, start, end, color, label, ymax=None):
    x = np.arange(start, end + 1)
    ax.fill_between(x, cov, color=color, alpha=0.35, lw=0)
    ax.plot(x, cov, color=color, lw=0.6)
    if ymax is None:
        ymax = max(cov.max() * 1.35, 1)
    ax.set_xlim(start, end)
    ax.set_ylim(0, ymax)
    ax.set_ylabel(label, fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelbottom=False)
    return ymax

def save(fig, base):
    fig.savefig(f'{OUT_DIR}/{base}.pdf', bbox_inches='tight')
    fig.savefig(f'{OUT_DIR}/{base}.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved {base}.pdf / .png')

# ==========================================================================
# S8A  Smarcc2  A5SS  chr10:+
#   common upstream exon ends 128487933 (rMATS shortES/longES start)
#   short 5'SS donor end = 128488017 ; long 5'SS donor end = 128488359 (alt 342bp)
#   downstream acceptor exon starts ~128488895 (flankingES)
# ==========================================================================
def panel_smarcc2():
    print('S8A Smarcc2 (A5SS chr10:+)')
    chrom = 'chr10'
    PS, PE = 128487500, 128490400
    common_s, common_e = 128487933 - 130, 128488017   # upstream exon w/ common region
    alt_s, alt_e = 128488017, 128488359               # alt 342bp (short->long donor)
    down_s, down_e = 128488895, 128490174             # downstream exon
    # donor of short = 128488017 -> acceptor 128488896 (SJ 128488018-128488895)
    # donor of long  = 128488359 -> acceptor 128488896 (SJ 128488360-128488895)
    ctrl_short = junc_mean(CTRL, chrom, 128488018, 128488895)
    ctrl_long  = junc_mean(CTRL, chrom, 128488360, 128488895)
    kd_short   = junc_mean(KD,   chrom, 128488018, 128488895)
    kd_long    = junc_mean(KD,   chrom, 128488360, 128488895)
    print(f'   short 5SS  ctrl={ctrl_short:.0f} kd={kd_short:.0f}')
    print(f'   long  5SS  ctrl={ctrl_long:.0f}  kd={kd_long:.0f}')
    c_cov = group_coverage(CTRL, chrom, PS, PE)
    k_cov = group_coverage(KD,   chrom, PS, PE)

    fig, (a1, a2, ag) = plt.subplots(
        3, 1, figsize=(7.2, 4.6),
        gridspec_kw={'height_ratios': [1, 1, 0.28]}, sharex=True)
    ymax = max(c_cov.max(), k_cov.max()) * 1.35
    for ax, cov, col, lab, cs, cl in [
            (a1, c_cov, CTRL_COLOR, 'siControl (n=3)\nRPKM', ctrl_short, ctrl_long),
            (a2, k_cov, KD_COLOR,   'siDHX8 (n=3)\nRPKM',    kd_short, kd_long)]:
        style_cov_ax(ax, cov, PS, PE, col, lab, ymax=ymax)
        ax.axvspan(alt_s, alt_e, color='#F0A030', alpha=0.18, zorder=0)  # alt region
        draw_arc(ax, alt_s, down_s, ymax, cs, '#7B68B6', lw=3.0, frac=0.50)   # short 5SS
        draw_arc(ax, alt_e, down_s, ymax, cl, '#2E8B57', lw=2.2, frac=0.68)   # long 5SS
    # gene model
    ag.set_ylim(0, 1); ag.set_xlim(PS, PE); ag.axis('off')
    ag.add_patch(Rectangle((common_s, 0.3), common_e - common_s, 0.4,
                 facecolor='#3B6FB0', edgecolor='k', lw=0.6))
    ag.add_patch(Rectangle((alt_s, 0.3), alt_e - alt_s, 0.4,
                 facecolor='#F0A030', edgecolor='k', lw=0.6))
    ag.add_patch(Rectangle((down_s, 0.3), down_e - down_s, 0.4,
                 facecolor='#3B6FB0', edgecolor='k', lw=0.6))
    ag.plot([common_e, alt_s], [0.5, 0.5], color='k', lw=0.7)
    ag.text((alt_s+alt_e)/2, 0.02, 'alt. (342 bp)', ha='center', va='top',
            fontsize=7, color='#B06000')
    ag.text(PS+50, 0.85, 'Smarcc2 (+ strand)', fontsize=7.5, style='italic')
    a1.legend(handles=[mpatches.Patch(color='#2E8B57', label='Long 5′SS'),
                       mpatches.Patch(color='#7B68B6', label='Short 5′SS')],
              fontsize=7, loc='upper right', framealpha=0.9)
    a1.set_title('Smarcc2 — A5SS  (dPSI = −0.309, FDR = 5.0×10⁻⁵)',
                 fontsize=10, fontweight='bold')
    fig.tight_layout()
    save(fig, 'sashimi_a5ss_smarcc2')
    return {'ctrl_short': ctrl_short, 'ctrl_long': ctrl_long,
            'kd_short': kd_short, 'kd_long': kd_long}

# ==========================================================================
# S8B  Sf1  A3SS  chr19:+
#   upstream exon ends ~128488... no: SJ shows donor 6374608 -> acceptors
#   long 3'SS acceptor  = 6374718 (SJ 6374609-6374718 ... exon) ; then exon->6375657
#   Represent proximal (short) vs distal (long) 3'SS at the alt-3' exon.
#   Use rMATS-consistent alt exon 6374718..6374898 window; downstream common exon 6375657-6377203
# ==========================================================================
def panel_sf1():
    print('S8B Sf1 (A3SS chr19:+)')
    chrom = 'chr19'
    PS, PE = 6373900, 6377400
    up_s, up_e = 6374548, 6374608          # upstream exon (donor 6374608)
    altA_s = 6374718                       # long 3'SS acceptor
    altB_s = 6374898                       # short 3'SS acceptor
    alt_e  = 6375657                       # then exon extends toward common
    down_s, down_e = 6375657, 6377203      # downstream common exon
    # junctions: donor 6374608 -> long acceptor 6374718 (SJ 6374609-6374718)
    #            acceptor-exon donor 6374898 -> 6375657 (SJ 6374899-6375657)
    ctrl_long  = junc_mean(CTRL, chrom, 6374609, 6374718)
    ctrl_short = junc_mean(CTRL, chrom, 6374899, 6375657)
    kd_long    = junc_mean(KD,   chrom, 6374609, 6374718)
    kd_short   = junc_mean(KD,   chrom, 6374899, 6375657)
    print(f'   long 3SS (6374718)  ctrl={ctrl_long:.0f} kd={kd_long:.0f}')
    print(f'   short/exon donor    ctrl={ctrl_short:.0f} kd={kd_short:.0f}')
    c_cov = group_coverage(CTRL, chrom, PS, PE)
    k_cov = group_coverage(KD,   chrom, PS, PE)

    fig, (a1, a2, ag) = plt.subplots(
        3, 1, figsize=(7.2, 4.6),
        gridspec_kw={'height_ratios': [1, 1, 0.28]}, sharex=True)
    ymax = max(c_cov.max(), k_cov.max()) * 1.35
    for ax, cov, col, lab, cl, cs in [
            (a1, c_cov, CTRL_COLOR, 'siControl (n=3)\nRPKM', ctrl_long, ctrl_short),
            (a2, k_cov, KD_COLOR,   'siDHX8 (n=3)\nRPKM',    kd_long, kd_short)]:
        style_cov_ax(ax, cov, PS, PE, col, lab, ymax=ymax)
        ax.axvspan(altA_s, altB_s, color='#F0A030', alpha=0.18, zorder=0)
        draw_arc(ax, up_e, altA_s, ymax, cl, '#2E8B57', lw=2.4, frac=0.60)   # long 3SS
        draw_arc(ax, altB_s, down_s, ymax, cs, '#7B68B6', lw=3.0, frac=0.50) # exon donor
    ag.set_ylim(0, 1); ag.set_xlim(PS, PE); ag.axis('off')
    ag.add_patch(Rectangle((up_s, 0.3), up_e - up_s, 0.4,
                 facecolor='#3B6FB0', edgecolor='k', lw=0.6))
    ag.add_patch(Rectangle((altA_s, 0.3), altB_s - altA_s, 0.4,
                 facecolor='#F0A030', edgecolor='k', lw=0.6))
    ag.add_patch(Rectangle((down_s, 0.3), down_e - down_s, 0.4,
                 facecolor='#3B6FB0', edgecolor='k', lw=0.6))
    ag.plot([up_e, altA_s], [0.5, 0.5], color='k', lw=0.7)
    ag.text((altA_s+altB_s)/2, 0.02, 'alt. 3′SS', ha='center', va='top',
            fontsize=7, color='#B06000')
    ag.text(PS+50, 0.85, 'Sf1 (+ strand)', fontsize=7.5, style='italic')
    a1.legend(handles=[mpatches.Patch(color='#2E8B57', label='Long 3′SS'),
                       mpatches.Patch(color='#7B68B6', label='Short 3′SS')],
              fontsize=7, loc='upper right', framealpha=0.9)
    a1.set_title('Sf1 — A3SS  (dPSI = +0.195, FDR = 1.8×10⁻¹²)',
                 fontsize=10, fontweight='bold')
    fig.tight_layout()
    save(fig, 'sashimi_a3ss_sf1')
    return {'ctrl_long': ctrl_long, 'ctrl_short': ctrl_short,
            'kd_long': kd_long, 'kd_short': kd_short}

# ==========================================================================
# S8C  Dpp7  SE  chr2:-
#   skipped exon 25353607-25353726 ; upstream 25353359-25353516 ; downstream 25354232-25354293
#   inclusion (up)  SJ 25353517-25353607 ; inclusion (dn) SJ 25353727-25354232
#   skipping        SJ 25353517-25354232
# ==========================================================================
def panel_dpp7():
    print('S8C Dpp7 (SE chr2:-)')
    chrom = 'chr2'
    PS, PE = 25353200, 25354450
    up_s, up_e = 25353359, 25353516
    ex_s, ex_e = 25353607, 25353726       # skipped exon
    dn_s, dn_e = 25354232, 25354293
    ctrl_inc_up = junc_mean(CTRL, chrom, 25353517, 25353607)
    ctrl_inc_dn = junc_mean(CTRL, chrom, 25353727, 25354232)
    ctrl_skip   = junc_mean(CTRL, chrom, 25353517, 25354232)
    kd_inc_up   = junc_mean(KD,   chrom, 25353517, 25353607)
    kd_inc_dn   = junc_mean(KD,   chrom, 25353727, 25354232)
    kd_skip     = junc_mean(KD,   chrom, 25353517, 25354232)
    print(f'   ctrl inc_up={ctrl_inc_up:.0f} inc_dn={ctrl_inc_dn:.0f} skip={ctrl_skip:.0f}')
    print(f'   kd   inc_up={kd_inc_up:.0f} inc_dn={kd_inc_dn:.0f} skip={kd_skip:.0f}')
    c_cov = group_coverage(CTRL, chrom, PS, PE)
    k_cov = group_coverage(KD,   chrom, PS, PE)

    fig, (a1, a2, ag) = plt.subplots(
        3, 1, figsize=(7.2, 4.6),
        gridspec_kw={'height_ratios': [1, 1, 0.28]}, sharex=True)
    ymax = max(c_cov.max(), k_cov.max()) * 1.35
    for ax, cov, col, lab, iu, idn, sk in [
            (a1, c_cov, CTRL_COLOR, 'siControl (n=3)\nRPKM', ctrl_inc_up, ctrl_inc_dn, ctrl_skip),
            (a2, k_cov, KD_COLOR,   'siDHX8 (n=3)\nRPKM',    kd_inc_up, kd_inc_dn, kd_skip)]:
        style_cov_ax(ax, cov, PS, PE, col, lab, ymax=ymax)
        ax.axvspan(ex_s, ex_e, color='#F0A030', alpha=0.16, zorder=0)
        draw_arc(ax, up_e, ex_s, ymax, iu, col, lw=2.0, frac=0.42)          # inclusion up
        draw_arc(ax, ex_e, dn_s, ymax, idn, col, lw=2.0, frac=0.42)         # inclusion dn
        draw_arc(ax, up_e, dn_s, ymax, sk, '#E8890C', lw=2.6, frac=0.70)    # skipping
    ag.set_ylim(0, 1); ag.set_xlim(PS, PE); ag.axis('off')
    for (s, e, cfill) in [(up_s, up_e, '#3B6FB0'), (ex_s, ex_e, '#E8890C'), (dn_s, dn_e, '#3B6FB0')]:
        ag.add_patch(Rectangle((s, 0.3), e - s, 0.4, facecolor=cfill, edgecolor='k', lw=0.6))
    ag.plot([up_e, ex_s], [0.5, 0.5], color='k', lw=0.7)
    ag.plot([ex_e, dn_s], [0.5, 0.5], color='k', lw=0.7)
    ag.text((ex_s+ex_e)/2, 0.02, 'skipped exon', ha='center', va='top',
            fontsize=7, color='#B06000')
    ag.text(PS+30, 0.85, 'Dpp7 (− strand)', fontsize=7.5, style='italic')
    a1.set_title('Dpp7 — Exon skipping in DHX8 KD (dPSI = +0.450, FDR = 0.0)',
                 fontsize=10, fontweight='bold')
    fig.tight_layout()
    save(fig, 'sashimi_se_dpp7')
    return {'ctrl_inc_up': ctrl_inc_up, 'ctrl_inc_dn': ctrl_inc_dn, 'ctrl_skip': ctrl_skip,
            'kd_inc_up': kd_inc_up, 'kd_inc_dn': kd_inc_dn, 'kd_skip': kd_skip}

# ==========================================================================
# S8D  Rpl4  IR  chr9:+
#   exon3 = 64175194-64175300 ; exon4 = 64175584-64175722
#   retained intron = 64175301-64175583 (283 bp)
#   flanking splice junction SJ 64175301-64175583
# ==========================================================================
def panel_rpl4():
    print('S8D Rpl4 (IR chr9:+)')
    chrom = 'chr9'
    PS, PE = 64175000, 64175900
    ex3_s, ex3_e = 64175194, 64175300
    intr_s, intr_e = 64175301, 64175583
    ex4_s, ex4_e = 64175584, 64175722
    ctrl_j = junc_mean(CTRL, chrom, 64175301, 64175583)
    kd_j   = junc_mean(KD,   chrom, 64175301, 64175583)
    print(f'   flanking junction  ctrl={ctrl_j:.0f} kd={kd_j:.0f}')
    c_cov = group_coverage(CTRL, chrom, PS, PE)
    k_cov = group_coverage(KD,   chrom, PS, PE)

    fig, (a1, a2, ag) = plt.subplots(
        3, 1, figsize=(7.2, 4.6),
        gridspec_kw={'height_ratios': [1, 1, 0.28]}, sharex=True)
    ymax = max(c_cov.max(), k_cov.max()) * 1.35
    for ax, cov, col, lab, jc in [
            (a1, c_cov, CTRL_COLOR, 'siControl (n=3)\nRPKM', ctrl_j),
            (a2, k_cov, KD_COLOR,   'siDHX8 (n=3)\nRPKM',    kd_j)]:
        style_cov_ax(ax, cov, PS, PE, col, lab, ymax=ymax)
        ax.axvspan(intr_s, intr_e, color='#F5B0B0', alpha=0.30, zorder=0)
        draw_arc(ax, intr_s, intr_e, ymax, jc, col, lw=2.6, frac=0.75)
    ag.set_ylim(0, 1); ag.set_xlim(PS, PE); ag.axis('off')
    ag.add_patch(Rectangle((ex3_s, 0.3), ex3_e - ex3_s, 0.4, facecolor='#333333', edgecolor='k', lw=0.6))
    ag.add_patch(Rectangle((ex4_s, 0.3), ex4_e - ex4_s, 0.4, facecolor='#333333', edgecolor='k', lw=0.6))
    ag.add_patch(Rectangle((intr_s, 0.32), intr_e - intr_s, 0.36, facecolor='none',
                 edgecolor='#C44E52', lw=1.2, ls='--'))
    ag.text((intr_s+intr_e)/2, 0.02, 'retained intron (283 bp)', ha='center', va='top',
            fontsize=7, color='#C44E52')
    ag.text(PS+20, 0.85, 'Rpl4 (+ strand)', fontsize=7.5, style='italic')
    a1.set_title('Rpl4 — IR in DHX8 KD (logFC = 4.47, baseMean = 16156)',
                 fontsize=10, fontweight='bold')
    fig.tight_layout()
    save(fig, 'sashimi_ir_rpl4')
    return {'ctrl_flank': ctrl_j, 'kd_flank': kd_j}

# ==========================================================================
# S8E  Ccdc71  pseudo-CE / IR  chr9:+
#   E1 = 108460527-108460617 ; putative CE (MAJIQ) = 108461420-108461594 (175bp)
#   E2 = 108462938-108465938
#   canonical junction E1->E2  SJ 108460618-108462937
# ==========================================================================
def panel_ccdc71():
    print('S8E Ccdc71 (pseudo-CE/IR chr9:+)')
    chrom = 'chr9'
    PS, PE = 108460300, 108463400
    e1_s, e1_e = 108460527, 108460617
    ce_s, ce_e = 108461420, 108461594     # 175 bp MAJIQ putative CE
    e2_s, e2_e = 108462938, 108463400     # show start of E2
    ctrl_can = junc_mean(CTRL, chrom, 108460618, 108462937)
    kd_can   = junc_mean(KD,   chrom, 108460618, 108462937)
    print(f'   canonical E1->E2 junction  ctrl={ctrl_can:.0f} kd={kd_can:.0f}')
    c_cov = group_coverage(CTRL, chrom, PS, PE)
    k_cov = group_coverage(KD,   chrom, PS, PE)

    fig, (a1, a2, ag) = plt.subplots(
        3, 1, figsize=(7.2, 4.6),
        gridspec_kw={'height_ratios': [1, 1, 0.28]}, sharex=True)
    ymax = max(c_cov.max(), k_cov.max()) * 1.35
    for ax, cov, col, lab, jc in [
            (a1, c_cov, CTRL_COLOR, 'siControl (n=3)\nRPKM', ctrl_can),
            (a2, k_cov, KD_COLOR,   'siDHX8 (n=3)\nRPKM',    kd_can)]:
        style_cov_ax(ax, cov, PS, PE, col, lab, ymax=ymax)
        ax.axvspan(e1_e, e2_s, color='#F5B0B0', alpha=0.16, zorder=0)          # intron/IR
        ax.axvspan(ce_s, ce_e, color='#F0A030', alpha=0.40, zorder=0)          # 175bp CE window
        draw_arc(ax, e1_e, e2_s, ymax, jc, col, lw=2.6, frac=0.78)
    ag.set_ylim(0, 1); ag.set_xlim(PS, PE); ag.axis('off')
    ag.add_patch(Rectangle((e1_s, 0.3), e1_e - e1_s, 0.4, facecolor='#3B6FB0', edgecolor='k', lw=0.6))
    ag.add_patch(Rectangle((ce_s, 0.3), ce_e - ce_s, 0.4, facecolor='#E8890C', edgecolor='k', lw=0.6))
    ag.add_patch(Rectangle((e2_s, 0.3), e2_e - e2_s, 0.4, facecolor='#3B6FB0', edgecolor='k', lw=0.6))
    ag.plot([e1_e, e2_s], [0.5, 0.5], color='#888888', lw=0.7)
    ag.text((ce_s+ce_e)/2, 0.9, 'CE (175 bp)', ha='center', va='bottom',
            fontsize=7, color='#B06000')
    ag.annotate('', xy=(e2_s, 0.05), xytext=(e1_e, 0.05),
                arrowprops=dict(arrowstyle='<->', color='#C44E52', lw=0.9))
    ag.text((e1_e+e2_s)/2, 0.02, 'IR (intron retention)', ha='center', va='top',
            fontsize=7, color='#C44E52')
    ag.text(PS+20, 0.85, 'Ccdc71 (chr9, + strand)', fontsize=7.5, style='italic')
    a1.set_title('Ccdc71 — DHX8 KD: IR + cryptic exon? (MAJIQ dPSI = +0.52)',
                 fontsize=10, fontweight='bold')
    fig.tight_layout()
    save(fig, 'sashimi_ce_ccdc71_dhx8')
    return {'ctrl_canonical': ctrl_can, 'kd_canonical': kd_can}

# ==========================================================================
if __name__ == '__main__':
    print('=== DHX8 sashimi reconstruction (v2 BAMs, corrected siCont_2) ===')
    summary = {}
    summary['S8A_Smarcc2'] = panel_smarcc2()
    summary['S8B_Sf1']     = panel_sf1()
    summary['S8C_Dpp7']    = panel_dpp7()
    summary['S8D_Rpl4']    = panel_rpl4()
    summary['S8E_Ccdc71']  = panel_ccdc71()
    print('\n=== JUNCTION-COUNT SUMMARY (group means, v2) ===')
    for k, v in summary.items():
        print(k, {kk: round(vv, 1) for kk, vv in v.items()})
