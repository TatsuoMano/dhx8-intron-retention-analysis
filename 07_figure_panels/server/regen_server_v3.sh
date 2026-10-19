#!/bin/bash
# v3 panels, server scripts (ncnp-linux). Working dir: /mnt/HDD2/unified_results/figregen_2026-09-29 (project dir, not /tmp).
# Same rules as regen_local_v3.sh. Outputs: panels/ (copied back to Analysis/2026-09-29_illustrator_v3/panels).
set -u
cd /mnt/HDD2/unified_results/figregen_2026-09-29
export PYTHONPATH=$HOME/bin
UR=/mnt/HDD2/unified_results
R="python3 regen_runner.py"
K="--keep-look --keep-ticks"
run() { $R "$@" 2>&1 | grep -E "regen_runner ->|Error|Traceback|no figure" ; }
~/bin/notify_telegram.sh "[DHX8-Sakai] figure regen v3: started (30 panels)" >/dev/null 2>&1
run $UR/_fig3a_volcano_sq.py --figsize 2.165 2.2 $K --axes-rect 0.2 0.17 0.97 0.97 --out panels/fig3a_volcano
run $UR/_fig3b_c2cp_v5.py --figsize 4.45 3.4 $K --out panels/fig3b_gsea_c2cp
run $UR/_fig5a_ir_volcano.py --figsize 2.69 2.8 $K --out panels/fig5a_ir_volcano   # about square axes (TM2: landscape); the full height of the B block (3.4 in) stretches the plot
# Fig 5B: four graphs of the same size with the same axes rectangle
for i in 0 1 2; do run $UR/_fig5_rnaquant_v7_v2.py --figsize 1.9 1.64 --panel $i $K --ylim 0 1.75 --axes-rect 0.2 0.22 0.97 0.87 --out panels/fig5b${i}_rna_quantity; done   # y range: the 7-pt legend sits above the peaks
run $UR/_fig5_rnaquant_v7_v2.py --figsize 1.9 1.64 --panel 3 $K --xlim -0.5 2.95 --ylim -3.4 5.2 --axes-rect 0.2 0.22 0.97 0.87 --out panels/fig5b3_rna_quantity   # room for the "-0.008" label inside the frame
run $UR/_fig6_venns_v2.py --figsize 2.165 1.86 --select fig6a_venn_khong $K --out panels/fig6a_venn_khong
run $UR/_fig6_venns_v2.py --figsize 2.165 1.86 --select fig6d_venn_namkoong $K --out panels/fig6b_venn_namkoong_as
run $UR/_fig6_venns_v2.py --figsize 2.165 1.86 --select namkoong_hs $K --out panels/fig6c_venn_namkoong_hs
run $UR/_fig6_venns_v2.py --figsize 2.165 1.86 --select fig6c_venn_moon $K --out panels/fig6d_venn_ren
# S12A/B: same scale and axes rectangle; counts on the centre line as in TM2 (the server script lifts them off it): the overlap count
# centred on the lens, the SG-only count at the right of its circle, so that "3" and "35" do not read as "335"; label of the pre-v2fix script corrected
V="--layout none --axes-rect 0 0 1 0.87 --xlim -0.59 1.396 --ylim -1.066 0.575"     # 72.5 pt per unit: IR circle 81 pt, room for "(Markmiller 2018)"
run $UR/_s11_venns_v2.py --figsize 2.0 1.9 --select youn --text-replace "(1,150)" "(1,149)" $K $V --text-move "^3$" 0.06 -0.14 --text-move "^35$" -0.105 0.06 --out panels/sfigS12a_venn_youn
run $UR/_s11_venns_v2.py --figsize 2.0 1.9 --select markmiller --text-replace "(1,150)" "(1,149)" $K $V --text-move "^9$" 0.062 -0.14 --text-move "^26$" -0.115 0.06 --out panels/sfigS12b_venn_markmiller
run $UR/_sfig7_gsea_v2.py --figsize 3.95 3.3 --select gobp $K --out panels/sfigS8a_gsea_gobp
run $UR/_sfig7_gsea_v2.py --figsize 2.66 3.3 --select hallmark $K --out panels/sfigS8b_gsea_hallmark
# S9: "Depth" (the script averages samtools depth), FDR exponents as mathtext superscripts in Arial (10⁻⁵ has no glyphs in Arial), junction counts lifted off their arcs; C: rMATS reports FDR 0 -> "FDR < 10⁻¹⁶" (author decision 2026-10-01)
run $UR/analysis_sashimi_dhx8_events.py --figsize 3.3 2.2 --cap 7 $K --text-replace "RPKM" "Depth" --text-replace "FDR = 0.0" "FDR < 10⁻¹⁶" --sup-mathtext bold --sup-min-pt 5.5 --sup-raise 0.33 --text-offset-pt "^(42|91|248|245|120|124|8|77|85|4028|2839|23|395)$" 0 2.5 --text-offset-pt "^(93|69)$" -10 1 --text-offset-pt "^573$" 0 6 --text-offset-pt "^CE \(175" 18 0 --text-offset-pt "^Ccdc71 \(chr9" -14 0 --text-offset-pt "^IR \(intron" 0 -2.5 --title-regex "\s+\(" "\n(" --out-dir panels/sashimi   # titles in two lines: at 7 pt the longest is wider than the panel
run sashimi_ce_gse126543_recovered.py --figsize 3.3 2.2 --cap 7 $K --select unc13a --cwd $UR --text-replace "TDP-43(-): disease (nuclear TDP-43 loss)" "TDP-43(−)\ndisease" --text-replace "TDP-43(+): normal" "TDP-43(+)\nnormal" --text-offset-pt "^2714$" 0 8 --text-offset-pt "^3183$" 0 8 --text-offset-pt "^42$" 0 16 --out panels/sfigS9f_sashimi_ce_unc13a
# S10: A and B rendered separately; leader lines end at the top of the gene labels
run $UR/majiq/_bespoke_ce_fdr_arial.py --figsize 4.7 3.3 --panel 0 --title-off $K --ylim -2.45 2.9 --leader-to-edge --out panels/sfigS10a_ce_candidates
run $UR/majiq/_bespoke_ce_fdr_arial.py --figsize 1.9 3.3 --panel 1 --title-off $K --text-offset-pt "^FDR=0\.1$" 0 4 --text-replace "GSE126543 (TDP-43+/-)" "GSE126543\n(TDP-43+/-)" --out panels/sfigS10b_ce_fdr
run $UR/sg_analysis/gsea_dhx8/sg_permutation_test_lenmatched.py --figsize 6.73 3.1 --select sg_permutation_test $K --legend-out --out panels/sfigS11a_permutation
run $UR/sg_analysis/gsea_dhx8/sg_permutation_test_lenmatched.py --figsize 6.73 2.75 --select diagnostics $K --legend-out --out panels/sfigS11b_permutation_diagnostics
ls panels/*.pdf panels/sashimi/*.pdf | wc -l
~/bin/notify_telegram.sh "[DHX8-Sakai] figure regen v3: ALL COMPLETE ($(ls panels/*.pdf panels/sashimi/*.pdf 2>/dev/null | wc -l) PDFs)" >/dev/null 2>&1
echo V3_DONE
