#!/usr/bin/env python3
"""RECOVERED SCRIPT (2026-09-17). Source: the ssh heredoc (python3 - <<'PYEOF' ... PYEOF on ncnp-linux) recorded in
~/.claude-mem/logs/claude-mem-2026-03-31.log at 2026-03-31 18:36:38 JST, i.e. the command that wrote
/mnt/HDD2/unified_results/sashimi_plots/sashimi_ce_unc13a_gse126543.png (the file placed in SuppFigures_TM2.ai as
figures_for_final/SupFig07/sfig_sashimi_ce_unc13a.png; pixel correlation 1.000). The script was never saved to disk
(this is exactly the /tmp-script loss described in ~/.claude/CLAUDE.md). Body below is verbatim from the log.
Runs on ncnp-linux (pysam, BAMs of GSE126543, GTF). Re-rendered at the v2 cell size through regen_runner.py.
"""
import pysam
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict
import subprocess
import re as regex_module

BAM_DIR = "/mnt/HDD2/sorted_nuclei_GSE126543/mapping"
# CORRECT: sample_sheet TDP43_positive = TDP-43(+) = normal
NORMAL_BAMS = [f"{BAM_DIR}/SRR857{i}/SRR857{i}.sorted.bam" for i in [1937,1939,1941,1944,1947,1949,1951]]
# CORRECT: sample_sheet TDP43_negative = TDP-43(-) = disease
DISEASE_BAMS = [f"{BAM_DIR}/SRR857{i}/SRR857{i}.sorted.bam" for i in [1938,1940,1942,1945,1948,1950,1952]]
GTF = "/mnt/HDD2/hg38.ncbiRefSeq.gtf"

def get_exons(gene, chrom, rs, re_val):
    cmd = f"grep -w '{gene}' {GTF} | grep -w exon | grep -w '{chrom}'"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    tx_exons = defaultdict(list)
    for line in result.stdout.strip().split('\n'):
        if not line: continue
        fields = line.split('\t')
        if len(fields) < 9: continue
        s, e = int(fields[3]), int(fields[4])
        m = regex_module.search(r'transcript_id "([^"]+)"', fields[8])
        if m: tx_exons[m.group(1)].append((s, e))
    best = []; best_c = 0
    for tid, tex in tx_exons.items():
        c = sum(1 for s, e in tex if s < re_val and e > rs)
        if c > best_c:
            best_c = c
            best = [(s, e) for s, e in tex if s < re_val + 500 and e > rs - 500]
    return sorted(set(best))

def get_data(bams, chrom, start, end):
    length = end - start
    coverage = np.zeros(length)
    junctions = defaultdict(int)
    total_reads = 0
    for bam_path in bams:
        try:
            with pysam.AlignmentFile(bam_path, "rb") as bam:
                for s in bam.get_index_statistics():
                    total_reads += s.mapped
                for read in bam.fetch(chrom, start, end):
                    if read.is_secondary or read.is_supplementary: continue
                    for bs, be in read.get_blocks():
                        b0, b1 = max(bs-start,0), min(be-start,length)
                        if b0 < b1: coverage[b0:b1] += 1
                    if read.cigartuples:
                        ref_pos = read.reference_start
                        for op, l in read.cigartuples:
                            if op == 3:
                                js, je = ref_pos, ref_pos + l
                                if start <= js <= end and start <= je <= end:
                                    junctions[(js, je)] += 1
                            if op in [0,2,3,7,8]: ref_pos += l
        except: pass
    rpm = total_reads / 1e6 if total_reads > 0 else 1
    coverage /= rpm
    return coverage, dict(junctions)

def draw_arc(ax, x1, x2, count, y_base, color, max_h, lw=None):
    mid = (x1+x2)/2
    t = np.linspace(0, np.pi, 80)
    x_arr = x1 + (x2-x1)*t/np.pi
    y = y_base + max_h*np.sin(t)
    if lw is None: lw = np.clip(count/30, 1, 5)
    ax.plot(x_arr, y, color=color, linewidth=lw, alpha=0.8, solid_capstyle='round')
    ax.text(mid, y_base+max_h*1.08, str(count), ha='center', va='bottom',
            fontsize=9, color=color, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.1', facecolor='white', edgecolor='none', alpha=0.8))

chrom = 'chr19'; rs = 17639000; re_val = 17646000
ce_s = 17642414; ce_e = 17642541

normal_cov, normal_junc = get_data(NORMAL_BAMS, chrom, rs, re_val)
disease_cov, disease_junc = get_data(DISEASE_BAMS, chrom, rs, re_val)
exons = get_exons('UNC13A', chrom, rs, re_val)

# Verify: CE junction = 17641556-17642413, normal junction = 17641556-17642844
print(f"TDP-43(+) normal: CE junc(17641556-17642413)={normal_junc.get((17641556,17642413),0)}, normal junc(17641556-17642844)={normal_junc.get((17641556,17642844),0)}")
print(f"TDP-43(-) disease: CE junc(17641556-17642413)={disease_junc.get((17641556,17642413),0)}, normal junc(17641556-17642844)={disease_junc.get((17641556,17642844),0)}")

fig, axes = plt.subplots(3, 1, figsize=(14, 8), height_ratios=[3, 3, 1],
                         gridspec_kw={'hspace': 0.08})

max_cov = max(normal_cov.max(), disease_cov.max())
if max_cov == 0: max_cov = 10
arc_base = max_cov * 0.55
arc_h = max_cov * 0.35
x = np.arange(len(normal_cov))
ce_x_s = ce_s - rs
ce_x_e = ce_e - rs

for idx, (cov, juncs, label, color) in enumerate([
    (normal_cov, normal_junc, 'TDP-43(+): normal', '#2166AC'),
    (disease_cov, disease_junc, 'TDP-43(-): disease (nuclear TDP-43 loss)', '#D6604D'),
]):
    ax = axes[idx]
    ax.fill_between(x, 0, cov, color=color, alpha=0.3)
    ax.plot(x, cov, color=color, linewidth=0.5)
    ax.axvspan(ce_x_s, ce_x_e, alpha=0.2, color='#FFAA00', zorder=0)
    
    for (js,je), cnt in sorted(juncs.items(), key=lambda xv:-xv[1])[:8]:
        if cnt < 3: continue
        touches_ce = (ce_s-20 < js < ce_e+20) or (ce_s-20 < je < ce_e+20)
        jcolor = '#FF6600' if touches_ce else color
        draw_arc(ax, js-rs, je-rs, cnt, arc_base, jcolor, arc_h)
    
    ax.set_ylim(0, max_cov * 1.5)
    ax.set_xlim(0, len(cov))
    ax.set_ylabel(label, fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelbottom=False)

axes[0].set_title('UNC13A — Cryptic Exon (GSE126543, FANS)', fontsize=13, fontweight='bold')

ax = axes[2]
y_mid = 0.5; eh = 0.35
for es, ee in exons:
    ds = max(es, rs) - rs
    de = min(ee, re_val) - rs
    if ds < de:
        ax.add_patch(plt.Rectangle((ds, y_mid-eh/2), de-ds, eh,
            facecolor='#4393C3', edgecolor='black', linewidth=0.8))
sorted_exons = sorted(exons)
for i in range(len(sorted_exons)-1):
    e1 = min(sorted_exons[i][1], re_val) - rs
    e2 = max(sorted_exons[i+1][0], rs) - rs
    if e1 < e2:
        ax.plot([e1, e2], [y_mid, y_mid], 'k-', linewidth=0.8)
ax.add_patch(plt.Rectangle((ce_x_s, y_mid-eh/2), ce_x_e-ce_x_s, eh,
    facecolor='#FFAA00', edgecolor='#CC7700', linewidth=1.5, alpha=0.8))
ax.text((ce_x_s+ce_x_e)/2, y_mid+eh/2+0.05, 'CE (127bp)', ha='center', va='bottom',
        fontsize=8, color='#CC7700', fontweight='bold')
ax.set_xlim(0, len(normal_cov)); ax.set_ylim(-0.1, 1.1); ax.axis('off')
ax.text(len(normal_cov)*0.5, -0.05, 'UNC13A (chr19)', ha='center', fontsize=9, fontstyle='italic')

plt.savefig('/mnt/HDD2/unified_results/sashimi_plots/sashimi_ce_unc13a_gse126543.png',
            dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved!")