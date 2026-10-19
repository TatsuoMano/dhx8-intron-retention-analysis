#!/bin/bash
# v3 panels, local scripts: every computational panel re-rendered at the size it has in the v3 layout (TM2 arrangement at
# 174 mm with 1.5 mm page margins, content width 6.73 in; layouts changed only where 7-pt text does not fit), with --keep-look
# (line widths, markers and text scaled with the figure as in TM2, then text floors 7 pt / numeric ticks 6 pt) and --keep-ticks.
# Original scripts are run unchanged through regen_runner.py.
set -u
cd "$(dirname "$0")"
R="python3.13 ../2026-09-17_figure_regen/regen_runner.py"; A=..; FF=../../Manuscript/Figures/figures_for_final; O=panels
K="--keep-look --keep-ticks"
run() { $R "$@" 2>&1 | grep -E "regen_runner ->|Error|Traceback|no figure" ; }
# Fig 3C: legend below in one column (inside, at 7 pt, it covers GUGUGU); same left/right axes edges as 3A
run $A/2026-07-07_transite_gene_lists/plot_fig3cd_ugrepeat.py --figsize 2.165 2.7 --select fig3c_kmer $K --label-wrap 30 --legend-out --legend-ncol 1 --text-offset-pt "^adj\. P" -16 2 --axes-rect 0.2 0.36 0.97 0.91 --out $O/fig3c_kmer_volcano
run $A/2026-07-07_transite_gene_lists/plot_fig3cd_ugrepeat.py --figsize 3.9 1.5 --select fig3d_motif $K --sci-notation --sup-min-pt 5.5 --sup-raise 0.33 --out $O/fig3d_motif_table   # adj. P as 3.3×10⁻¹⁵ (superscript 5.5 pt), not 3.3e-15; 1.5 in tall (bottom level with 3C) so that the exponents clear the header rule; --sup-raise 0.33 = Illustrator's superscript shift (mathtext's default 0.7 puts them 4.9 pt up on Arial)
run $A/2026-06-26_figure_fontfix/_fig4e_summary_v2.py --figsize 5.6 1.85 $K --cap 7 --title-off --out $O/fig4b_integrated_summary
run $A/2026-07-23_fig6e_fontfix/replot_fig6e_biancon_khongmatch.py --figsize 2.165 1.86 --select high_khongmatch $K --out $O/fig6e1_venn_biancon_high
run $A/2026-07-23_fig6e_fontfix/replot_fig6e_biancon_khongmatch.py --figsize 2.165 1.86 --select medhigh_khongmatch $K --out $O/fig6e2_venn_biancon_medhigh
run $A/2026-05-14_relabel_moon_to_ren/regenerate_fig6e_table.py --figsize 6.73 2.0 $K --cap 7 --text-remove "^Fig 6F" --out $O/fig6f_sg_overlap_table   # TM2 clips this title away
# Fig 7B/7C: same axes top and bottom; 7C legend below (at 7 pt it covers bars at every corner), x range widened so that "-0.64" clears the axis
run $A/2026-06-26_figure_fontfix/_fig7a_phase_response_v4.py --figsize 3.9 3.6 --panel 0 $K --axes-rect 0.13 0.19 0.98 0.87 --out $O/fig7b_phase_nes
run $A/2026-06-26_figure_fontfix/_fig7a_phase_response_v4.py --figsize 2.7 3.6 --panel 1 $K --legend-out --xlim -1.3 2.0 --axes-rect 0.27 0.19 0.96 0.87 --out $O/fig7c_helicases
run $A/2026-06-26_figure_fontfix/_fig7b_encode_v2.py --figsize 4.6 1.5 $K --cap 7 --title-off --out $O/fig7d_encode_table
# Fig 8A and S14B/C graphs (heights chosen so that the axes keep about the aspect they have in TM2); p values as 7.3×10⁻⁵ in italic (the labels are italic); interaction graphs: y to 1.06 so that the two-line label clears the error bar of the right grey point (at 1.03 it crossed it)
run $A/2026-06-26_figure_fontfix/_fig8a_eye_v2.py --figsize 2.3 2.0 --cap 7 --select R2_main_interaction $K --sci-notation --sup-mathtext italic --sup-min-pt 5.5 --sup-raise 0.33 --ylim 0.75 1.06 --out $O/fig8a2_eye_interaction_R2
run $A/2026-06-26_figure_fontfix/_fig8a_eye_v2.py --figsize 2.3 2.0 --cap 7 --select R2_main $K --sci-notation --sup-mathtext italic --sup-min-pt 5.5 --sup-raise 0.33 --clear-xaxis --ylim 0 1.5 --text-move "^interaction" 0 0.1 --text-offset-pt "^\*{3,4}$" 0 -0.8 --out $O/fig8a1_eye_bars_R2
run $A/2026-06-26_figure_fontfix/_fig8a_eye_v2.py --figsize 2.43 2.3 --cap 7 --select R3_interaction $K --sci-notation --sup-mathtext italic --sup-min-pt 5.5 --sup-raise 0.33 --ylim 0.75 1.06 --out $O/sfigS14c_eye_interaction_R3
run $A/2026-06-26_figure_fontfix/_fig8a_eye_v2.py --figsize 2.43 2.3 --cap 7 --select R3 $K --sci-notation --sup-mathtext italic --sup-min-pt 5.5 --sup-raise 0.33 --clear-xaxis --ylim 0 1.5 --text-move "^interaction" 0 0.1 --text-offset-pt "^\*{3,4}$" 0 -0.8 --out $O/sfigS14b_eye_bars_R3
# S1: taller, so that the two-line labels of HB and CTD fit when rotated and the motif labels clear the axis
run $FF/SupFig13/plot_dhx8_domain.py --figsize 6.73 2.4 $K --text-replace "helicase-deficient" "helicase-impaired" --text-rotate "^(HB|CTD)" 90 --text-offset-pt "^I \(GKT\)" -7 0 --text-offset-pt "^III " 6 0 --text-offset-pt "^II \(DEAH\)" -3 0 --out $O/sfigS1_dhx8_domain
run $A/2026-06-26_figure_fontfix/_s11c_gocc.py --figsize 2.49 1.9 $K --cap 7 --axes-rect 0.6 0.3 0.9 0.85 --text-offset-pt "^GO:CC" -24 0 --tick-replace "cytoplasmic ribonucleoprotein granule" "cytoplasmic ribonucleoprotein\ngranule" --out $O/sfigS12c_gocc
run $A/2026-07-11_s12_finalize/replot_s12_v1style.py --figsize 6.73 2.3 $K --cap 7 --legend-loc "upper right" --axes-ylim 2 0 1.08 --out $O/sfigS13_motif_features   # C: room for the 7-pt legend above the histogram
