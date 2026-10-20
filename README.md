# Analysis code: DHX8 helicase dysfunction promotes intron retention in stress granule-enriched transcripts and TDP-43-associated foci formation

This repository holds the analysis and plotting code for the article

> Sakai R.†, Mano T.†, Fujikake N., Kabuta C., Hashimoto T., Kabuta T. DHX8 helicase dysfunction promotes intron retention in stress granule-enriched transcripts and TDP-43-associated foci formation. *Nucleic Acids Research* (submitted, 2026). †These authors contributed equally.

Correspondence: Ryohei Sakai (ryohei-sakai@nms.ac.jp) and Tomohiro Kabuta (kabuta@ncnp.go.jp). The computational analyses were done by Tatsuo Mano.

The scripts are the versions that produced the reported results and the figure panels of the submitted version. They are kept as they were run, with their original file names. Superseded versions are not included. Data, count matrices, alignments and figures are not included (see [Data](#data)).

## Layout

| Folder | Contents | Ran on |
|---|---|---|
| `01_preprocessing/` | Trimming, STAR alignment, featureCounts gene counts, strandedness check (Supplementary Table S1), download and alignment of the comparison datasets | server |
| `02_splicing/` | Intron and exon feature definition, Welch pre-step and limma eBayes intron retention (IR), rMATS (DHX8 KD, TDP-43 datasets, ENCODE), MAJIQ and cryptic-exon (CE) calling with the mirror-method FDR, sashimi plots | server |
| `03_expression_and_pathways/` | DESeq2, fgsea (C2:CP, GO:BP, Hallmark), spliceosome-phase gene sets and phase GSEA | server |
| `04_stress_granule_enrichment/` | Fisher overlap of IR genes and DEG-UP genes with SG transcriptomes and proteomes, Fig. 6F table, length x expression matched permutation test, g:Profiler GO:CC | server |
| `05_sequence_features_and_motifs/` | RBP motif density, k-mer, GC and G-run analysis of retained introns; Transite 3′ UTR k-mer analysis | server |
| `06_statistics_wetlab/` | Recomputation of the statistics of the GraphPad Prism graphs, per-cell analysis of Supplementary Fig. S2B, Drosophila eye 2 x 2 Welch–James analysis, Brown–Forsythe test, fly-level sensitivity analysis | local (macOS) |
| `07_figure_panels/` | Plotting scripts for the computational panels, the rendering wrapper `regen_runner.py` and the exact render commands (`server/regen_server_v3.sh`, `local/regen_local_v3.sh`, `server/rerun_s9_2026-10-01.sh`, `server/rerun_s9_2026-10-20.sh`) | server and local |
| `lib/figprov.py` | Figure provenance helper imported by the plotting scripts and by `regen_runner.py` (writes PDF + PNG and a JSON sidecar with input and script hashes) | both |
| `environment/` | Software versions of the environments the scripts ran in | |

"Server" is the authors' Linux workstation (Ubuntu, data under `/mnt/HDD2` and reference files under `/mnt/NAS/database`). "Local" is a macOS workstation; local scripts read the GraphPad Prism files and intermediate tables from the authors' Dropbox project folder (`/Users/tatsuomano/Dropbox/Research/DHX8-Sakai/...`).

**Absolute paths.** The scripts contain the absolute paths of the authors' machines (`/mnt/HDD2/...`, `/mnt/NAS/database/...`, `/home/tatsuomano-linux2/...`, `/Users/tatsuomano/...`). They have to be changed to local paths before the scripts are run elsewhere. Several shell wrappers also call `~/bin/notify_telegram.sh`, a private job-completion notifier that is not part of this repository; the calls can be deleted or replaced by `true`. MAJIQ needs its own academic licence file (`--license`), which is not included. The shell scripts are bash scripts and need bash 4 or later (`recount_tdp43_datasets.sh` needs 4.3 for `declare -n`); the `/bin/bash` 3.2 of macOS is too old.

## Results, figure panels and scripts

Paths are relative to the repository root. "Inputs" names the files as they are on the authors' server or workstation.

### Pre-processing and expression

| Reported result / panel | Script | Inputs |
|---|---|---|
| Read trimming of the DHX8 KD libraries (Trimmomatic PE) | `01_preprocessing/DHX8_KD/011_trim_qc.sh` (with `000_sample_list.txt`, `001_config.txt`); re-trim of siCont_2: `01_preprocessing/DHX8_KD/_v2_stage1_rerun.sh` step 1 | raw FASTQ (GEO GSE337869) |
| STAR alignment of the DHX8 KD libraries to mm10 | `01_preprocessing/DHX8_KD/050_mapping.sh`, `070_samtools.sh`; siCont_2: `_v2_stage1b_rerun.sh` | trimmed FASTQ, STAR mm10 index |
| featureCounts gene counts (`-s 2`) | `01_preprocessing/DHX8_KD/_v2_stage1b_rerun.sh` step 3 | `*_sorted.bam`, `mm10.ncbiRefSeq.gtf` |
| Strandedness check, Supplementary Table S1 | `01_preprocessing/strandedness/build_strand_table.sh`, `parse_strand_table.py` | featureCounts summaries for `-s 0/1/2` |
| GSE65973 download and STAR alignment | `01_preprocessing/GSE65973/01_download.sh`, `02_star_mapping.sh`, `auto_pipeline.sh` | GEO GSE65973, `sample_sheet.tsv` |
| GSE126543 download, STAR 2-pass alignment, counts | `01_preprocessing/GSE126543/sorted_nuclei_GSE126543_pipeline.sh` (steps 1–3) | GEO GSE126543 (SRA), `sample_sheet.tsv` |
| featureCounts recounts of GSE65973 (`-s 1`; gene, intron, flanking exon) and GSE126543 (`-s 0`; gene, intron, flanking exon) | `01_preprocessing/recount_tdp43_datasets.sh` (reconstructed; by default (`--check`) it compares each command with the header of its count file and runs nothing; `--run` runs featureCounts) | dataset BAMs, `mm10.ncbiRefSeq.gtf`, `hg38.ncbiRefSeq.gtf`, intron and exon SAF files |
| ENCODE knockdown BAMs | `01_preprocessing/encode/download_all.sh` | ENCODE accessions in `bam_metadata.tsv` |
| DESeq2, DHX8 KD and GSE65973; Fig. 3A input | `03_expression_and_pathways/deseq2_mouse.R` | `counts_s2/gene_counts.txt`; GSE65973 `counts_s1/gene_counts.txt` |
| DESeq2, GSE126543 | `03_expression_and_pathways/deseq2_unified.R` (GSE126543 section) | `sorted_nuclei_GSE126543/counts/gene_counts.txt` |
| DESeq2 table in CSV layout for downstream scripts | `03_expression_and_pathways/deseq2_tsv_to_csv.sh` (reconstructed, see below) | `deseq2_dhx8_kd_vs_ctrl.tsv` |
| fgsea on the Wald statistic (C2:CP, GO:BP, Hallmark); Fig. 3B, Supplementary Fig. S8A, S8B | `03_expression_and_pathways/fgsea_unified.R` | `deseq2_*.tsv`, MSigDB GMT files |
| Supplementary Table S2 (phase membership with DHX8 KD log2FC, padj, direction, and `Fig7B_subset`, membership of the curated subset whose mean log2FC is printed in Fig. 7B) | `03_expression_and_pathways/build_table_s2.py` (reconstructed; `Fig7B_subset` added 2026-10-20) | `spliceosome_phase_genesets.gmt`, `deseq2_dhx8_kd_vs_ctrl.tsv` (current table, after the siCont_2 correction) |
| Spliceosome-phase GSEA (Supplementary Table S2 phases); Fig. 7B NES | `03_expression_and_pathways/fgsea_spliceosome_phase.R` | `deseq2_results.csv`, `spliceosome_phase_genesets.gmt` (included) |
| Fig. 3A DEG volcano | `07_figure_panels/server/_fig3a_volcano_sq.py` | `deseq2_dhx8_kd_vs_ctrl.tsv` |
| Fig. 3B C2:CP bar plot | `07_figure_panels/server/_fig3b_c2cp_v5.py` | `gsea/fgsea_c2cp_dhx8.tsv` |
| Supplementary Fig. S8A, S8B GO:BP and Hallmark bar plots | `07_figure_panels/server/_sfig7_gsea_v2.py` | `gsea/fgsea_gobp_dhx8.tsv`, `fgsea_hallmark_dhx8.tsv` |
| Fig. 7B per-phase mean log2FC of the curated 83-factor subset (27 E/A, 16 B, 6 Bact/NTC, 30 regulatory, 2 post-catalytic, 2 catalytic): +0.408, +0.369, +0.342, +0.260, −0.011, +0.011 | `03_expression_and_pathways/fig7b_phase_mean_log2fc.py` | `deseq2_dhx8_kd_vs_ctrl.tsv`; `spliceosome_phase_genesets.gmt` (all-member means, printed for comparison only) |
| Fig. 7B phase NES bars with the mean log2FC labels; Fig. 7C helicase log2FC | `07_figure_panels/local/_fig7a_phase_response_v4.py`. The script draws both panels in one figure and takes no arguments; `--panel` is an option of `regen_runner.py`, which keeps axes 0 (NES, Fig. 7B) or axes 1 (helicases, Fig. 7C) and removes the other (`regen_local_v3.sh`) | values written in the script: NES from `fgsea_spliceosome_phase.tsv`, mean log2FC from `fig7b_phase_mean_log2fc.py`, helicase log2FC from `deseq2_dhx8_kd_vs_ctrl.tsv`. The script header says the mean log2FC comes from the GMT gene sets; the printed values are the curated-subset means. The header is left as it was run, so that the script hash matches the provenance sidecar of the panels |

### Splicing

| Reported result / panel | Script | Inputs |
|---|---|---|
| Intron and flanking-exon features (SAF, 280,767 introns) | `02_splicing/intron_retention/generate_intron_gtf.py`, `step1_per_intron_count.py` (`run_pipeline.sh` step 1) | `mm10.ncbiRefSeq.gtf` |
| Welch pre-step and tested-intron set (36,190 introns) | `02_splicing/intron_retention/run_strandfix_dhx8.sh` → `step2_identify_retained_introns.py`, `validation_step1_limma.R` | DHX8 KD BAMs, SAF files |
| limma eBayes IR: 1,258 retained introns in 1,150 genes; Fig. 4B IR row, Fig. 5A | `02_splicing/intron_retention/_v2_limma_correct.sh` → `run_limma_strandfix_v2.R` | `counts_s2/{intron,exon}_counts.txt`, `intron_metadata.tsv`, `ir_results_all_dhx8.tsv` |
| Intron coordinates and length joined to the IR tables | `02_splicing/intron_retention/join_ir_meta.py` (reconstructed), called by `04_stress_granule_enrichment/_v2_sg_fix.sh`, `_v2_sg_fix2.sh` | `limma_*_dhx8_s2.tsv`, `intron_metadata.tsv` |
| limma IR for GSE65973 and GSE126543 (0 retained introns; Fig. 4B) | `02_splicing/intron_retention/limma_ir_mouse.R` (GSE65973), `limma_ir_unified.R` (GSE126543; hg38 features from `step1_hg38_intron_saf.py`) | intron and exon counts of each dataset |
| rMATS, DHX8 KD (SE 4,204, A3SS 599, A5SS 458; Fig. 4B, Fig. 5B) | `02_splicing/rmats/_v2_rmats_dhx8.sh` | DHX8 KD BAMs, `mm10.ncbiRefSeq.gtf` |
| SE direction: 3,596 of 4,204 significant SE events (86%) with control PSI > KD PSI, exact binomial P = 10^−512.4 | `02_splicing/rmats/se_direction_binomial.py` | `rmats/dhx8_v2/SE.MATS.JCEC.txt`, `summary.txt` |
| rMATS, GSE65973 and GSE126543 (Fig. 4B) | `02_splicing/rmats/rmats_tdp43_datasets_2026-03-28.sh` (reconstructed) | dataset BAMs, sample sheets |
| rMATS, ENCODE spliceosome-factor knockdowns (Fig. 7D) | `02_splicing/rmats/run_encode_splicing.sh`, `run_remaining.sh` (RBFOX2, RBM39); the summary block at the end of `run_encode_splicing.sh` filters on fixed awk columns that are not FDR and IncLevelDifference in these files and was not used (see `encode_ri_percent_fig7d.py`) | ENCODE BAMs, Ensembl GRCh38.109 GTF |
| Fig. 7D SE and RI counts and RI% = RI/(SE + A3SS + A5SS + RI), FDR < 0.05 and \|IncLevelDifference\| > 0.1, MXE excluded | `02_splicing/rmats/encode_ri_percent_fig7d.py` | `encode_rbp_kd/splicing/rmats_<factor>/*.MATS.JCEC.txt` |
| MAJIQ build, psi-coverage, deltapsi for DHX8 KD (22 de novo junctions, 12 CE; Fig. 4B CE row) | `02_splicing/majiq/_v2_majiq_dhx8.sh` (build), `_v2_motif_majiq.sh` (MAJIQ block), `analysis_deltapsi_integration.py` | `majiq/dhx8/experiments.tsv`, GENCODE vM25 GFF3 |
| MAJIQ for GSE126543 and GSE65973 (build, psi-coverage) | `02_splicing/majiq/run_majiq_robust.sh` | `majiq/<dataset>/experiments.tsv`, GENCODE v38 / vM25 GFF3 |
| MAJIQ deltapsi with splice-graph annotation for GSE126543, GSE65973 (and GSE27218): `deltapsi_annotated.tsv` | `02_splicing/majiq/deltapsi_annotated_2026-03-29.sh` (reconstructed; by default (`--check`) it compares the command below with the command recorded in each output and runs nothing; `--run` runs deltapsi) | `<dataset>_psicov.zarr`, `build/splicegraph.zarr` |
| Supplementary Fig. S10A, S10B CE candidates and mirror-method FDR | `02_splicing/majiq/_bespoke_ce_fdr_arial.py` | `majiq/all_datasets_significant_events.tsv` |
| Supplementary Fig. S9A–S9E sashimi plots (Smarcc2, Sf1, Dpp7, Rpl4, Ccdc71) | `02_splicing/sashimi/analysis_sashimi_dhx8_events.py` | DHX8 KD BAMs (`samtools depth`), STAR `SJ.out.tab` |
| Supplementary Fig. S9F UNC13A cryptic exon (GSE126543) | `02_splicing/sashimi/sashimi_ce_gse126543_recovered.py` | GSE126543 BAMs, `hg38.ncbiRefSeq.gtf` |
| Fig. 4B event counts: rMATS SE/A3SS/A5SS (and MXE, RI) at FDR < 0.05, limma retained and released introns and genes, MAJIQ de novo junctions, CE candidates and strong CEs, for the three datasets | `02_splicing/count_events_fig4b.py` | rMATS `*.MATS.JCEC.txt`, `limma_ir_results_*.tsv`, `majiq/<dataset>/deltapsi_annotated.tsv` |
| Fig. 4B integrated splicing summary table | `07_figure_panels/local/_fig4e_summary_v2.py` | values written in the script, from `count_events_fig4b.py` |
| Fig. 5A IR volcano | `07_figure_panels/server/_fig5a_ir_volcano.py` | `limma_ir_results_dhx8.tsv` |
| Fig. 5B RNA quantity by event type (Welch t-test, Bonferroni x3, Cohen's d) | `07_figure_panels/server/_fig5_rnaquant_v7_v2.py` (the statistics are computed in this script) | `deseq2_dhx8_kd_vs_ctrl.tsv`, `limma_retained_dhx8_s2.tsv`, `rmats/dhx8_v2/{SE,A3SS,A5SS}.MATS.JCEC.txt` |
| 95% CI of Cohen's d for Fig. 5B (IR 0.179–0.300; large-sample analytic interval, a percentile bootstrap is printed for comparison) | `02_splicing/ir_rna_level_effect_size.py` | as Fig. 5B |
| Fig. 7D ENCODE RI proportion table | `07_figure_panels/local/_fig7b_encode_v2.py` | values written in the script, from `encode_ri_percent_fig7d.py` |

### Stress granule enrichment

| Reported result / panel | Script | Inputs |
|---|---|---|
| Fisher overlap of IR genes and DEG-UP genes with SG transcriptomes and proteomes (universe 12,829; Namkoong dataset-specific universe) | `04_stress_granule_enrichment/analysis_sg_transcriptome_gsea.py` (run by `_v2_sg_fix2.sh`) | `deseq2_results.csv`, `limma_{retained,ir_results}_dhx8.tsv`, published SG tables (Supplementary Table S5) |
| Fig. 6F table values (IR OR vs DEG-UP OR) | `04_stress_granule_enrichment/_fig6f_table_v2.py`; drawn by `07_figure_panels/local/regenerate_fig6e_table.py` (values written in the script) | as above |
| Length x expression matched permutation test (10,000 permutations); Supplementary Fig. S11A, S11B | `04_stress_granule_enrichment/sg_permutation_test_lenmatched.py` (`run_sg_perm.sh`, `_v2_sg.sh`) | `limma_retained_dhx8.tsv`, `deseq2_dhx8_kd_vs_ctrl.tsv`, `mm10_gene_length.tsv`, Khong and Biancon sets |
| g:Profiler GO:CC of the IR genes (g:SCS) | `04_stress_granule_enrichment/gocc_v2.py` (reconstructed) | `limma_retained_dhx8_s2.tsv` |
| Fig. 6A–6D Venns (Khong, Namkoong AS, Namkoong HS, Ren); Supplementary Fig. S12A, S12B (Youn, Markmiller) | `07_figure_panels/server/_fig6_venns_v2.py`, `_s11_venns_v2.py` | `sg_rna_ir_overlap.tsv` |
| Fig. 6E Biancon Venns | `07_figure_panels/local/replot_fig6e_biancon_khongmatch.py` | values written in the script, from `sg_rna_ir_overlap.tsv` |
| Supplementary Fig. S12C GO:CC bar chart | `07_figure_panels/local/_s11c_gocc.py` | values written in the script, from `gocc_ir_v2_gprofiler.tsv` |

### Sequence features and motifs

| Reported result / panel | Script | Inputs |
|---|---|---|
| RBP motif density (14 motifs, Mann–Whitney U, Bonferroni), 4/5/6-mer chi-squared, GC and G-run of 1,258 retained vs 3x length-matched introns (KS D = 0.011) | `05_sequence_features_and_motifs/analysis_motif_enrichment.py` | `limma_retained_dhx8.tsv`, `limma_ir_results_dhx8.tsv`, mm10 FASTA |
| Supplementary Fig. S13 motif, GC and length panels | `07_figure_panels/local/replot_s12_v1style.py` | `motif_enrichment_results.tsv`, `s12_per_intron_features_v2.tsv`, `s12_length_arrays_v2.tsv`, `s12_ks_v2.txt` |
| Transite foreground and background tables: one coding RefSeq transcript per gene (NM_ preferred, then XM_; longest by total exon length); DEG-UP log2FC > 0.5 set (1,545), strict DEG-UP log2FC > 1 (459 rows, 451 with a 3′ UTR), MEF-expressed background (11,972; 11,900 with a 3′ UTR) | `05_sequence_features_and_motifs/transite/build_transite_tables.py` (reconstructed) | `mm10.ncbiRefSeq.gtf`, `deseq2_dhx8_kd_vs_ctrl.tsv`, `mm10_3utr_seqs.rds` |
| mm10 3′ UTR sequences for Transite | `05_sequence_features_and_motifs/transite/extract_3utr.R` | `mm10.ncbiRefSeq.gtf`, mm10 FASTA |
| Transite hexamer TSMA, strict DEG-UP (n = 451) vs MEF-expressed background and vs all mm10 3′ UTRs | `05_sequence_features_and_motifs/transite/run_tsma_strict.R` | `deseq2_dhx8_kd_vs_ctrl.tsv`, `foreground_up_only_v2.tsv`, `background_all_expressed_v2.tsv`, `mm10_3utr_seqs.rds` |
| Transite length-matched background (KS D = 0.016, P > 0.99) | `05_sequence_features_and_motifs/transite/run_tsma_lenmatch.R` | `foreground_up_strict_v2.tsv`, `background_all_expressed_v2.tsv`, `mm10_3utr_seqs.rds` |
| Fig. 3C k-mer volcano, Fig. 3D hexamer/RBP table | `07_figure_panels/local/plot_fig3cd_ugrepeat.py` | `fig3c_kmer_strict_mef.tsv` |

### Wet-lab statistics and the Drosophila experiments

| Reported result / panel | Script | Inputs |
|---|---|---|
| Fig. 1C printed P values (two-way ANOVA + Šídák) | `06_statistics_wetlab/prism_stats.py` (`fig1c`) | `Fig1c.pzfx` |
| Fig. 1E adjusted P values (one-way ANOVA + Tukey; check of the Prism values) | `06_statistics_wetlab/prism_stats.py` (`fig1e`) | `Fig1f. P-EIF2a定量_240112 .pzfx` |
| Fig. 8B TDP-43 and TDP-35 P values (one-way ANOVA + Dunnett; check of the Prism values) | `06_statistics_wetlab/prism_stats.py` (`fig8b`) | `Fig8b4. Fujikake-_定量データ_250716.pzfx` |
| Supplementary Fig. S7A, S7B unpaired t-tests (check of the Prism values) | `06_statistics_wetlab/prism_stats.py` (`s7`) | `SupFig12a._DHX8KD_tdp43_actb_gapdh まとめ_250306.pzfx`, `SupFig12b&c_231222 定量.pzfx` |
| Supplementary Fig. S2B per-cell statistics (Wilcoxon signed-rank, Welch) | `06_statistics_wetlab/analysis_s2b_percell.py` | `SupFig02a. 2SGs定量_230609IC.pzfx` |
| GraphPad Prism file reader used by the two scripts above | `06_statistics_wetlab/prism_panels.py` | |
| Fig. 8A pigmented-area bars and interaction plot, Supplementary Fig. S14B, S14C; Welch–James 2 x 2 interaction and simple-effect Welch t-tests for both RNAi lines | `06_statistics_wetlab/_fig8a_eye_v2.py` | `Fig8a.2 260514_analysis.pzfx` (`pigmantation` table) |
| Brown–Forsythe test of the eye groups (line #2 P ≈ 1 × 10⁻⁵) | `06_statistics_wetlab/brown_forsythe_eye.py` (reconstructed) | `Fig8a.2 260514_analysis.pzfx` |
| Conventional two-way ANOVA for the principal line, and fly-level sensitivity of the eye tests (eyes averaged per fly over all pairings; not reported in the article) | `06_statistics_wetlab/fly_level_sensitivity/analysis_fly_level_sensitivity.py` | `Fig8a.2 260514_analysis.pzfx`, collaborator spreadsheet `260508-2.xlsx` (cross-check) |

`_fig8a_eye_v2.py` is the figure version of the original eye analysis (`analysis_drosophila_eye_twoway.py`, not included). The statistics code is the same and both write a byte-identical `eye_twoway_summary_stats.tsv`.

### Other panels

| Panel | Script | Inputs |
|---|---|---|
| Supplementary Fig. S1 DHX8 domain architecture | `07_figure_panels/local/plot_dhx8_domain.py` | coordinates written in the script (UniProt Q14562) |

### Panels without code

These panels are micrographs, blots, drawings or graphs made in GraphPad Prism 9, and have no analysis code: Fig. 1A, 1B, 1D, 1F; Fig. 1C and 1E graphs (Prism); Fig. 2A–2E (Fig. 2E line scan measured in ImageJ and plotted in Prism); Fig. 4A (schematic); Fig. 7A (BioRender); Fig. 8A photographs; Fig. 8B blots and graphs (Prism); Fig. 8C (BioRender); Supplementary Fig. S2A (time-lapse frames), S2B graph (Prism), S3, S4, S5A, S5B (ImageJ line scan in Prism), S6, S7 (Prism graphs and blots) and S14A (photographs). Image measurements (foci counts, band densitometry, eye pigmentation, line scans) were made by hand in ImageJ without macros. The Prism files (`.pzfx`) hold the measurements and the Prism analyses.

Three values printed on the assembled figures differ from the original Prism output: the Fig. 1C P values were replaced by the values from `prism_stats.py` (0.0014 and 0.0005), the Prism P values of Supplementary Fig. S2B were removed (the per-cell statistics in the legend come from `analysis_s2b_percell.py`), and the first-lane label of Fig. 8B reads "GFP".

## How the figure panels were rendered

Every computational panel was rendered by running its plotting script unchanged through `07_figure_panels/regen_runner.py`. The wrapper intercepts `savefig`, resizes the selected figure to its final size in the layout, sets Arial, clamps text to at least 7 pt (numeric ticks 6 pt), applies the listed text edits and saves the panel with `lib/figprov.py` (PDF, PNG and a provenance sidecar). It does no analysis. The exact commands, with all arguments, are in

- `07_figure_panels/server/regen_server_v3.sh` (server panels: Fig. 3A, 3B, 5A, 5B, 6A–6D, S8, S9, S10, S11, S12A, S12B);
- `07_figure_panels/server/rerun_s9_2026-10-01.sh` (re-render of Supplementary Fig. S9A–S9E for revision 2) and `07_figure_panels/server/rerun_s9_2026-10-20.sh` (the final re-render: the same command with the S9E title corrected to MAJIQ dPSI = +0.53 by `--text-replace`);
- `07_figure_panels/local/regen_local_v3.sh` (local panels: Fig. 3C, 3D, 4B, 6E, 6F, 7B, 7C, 7D, 8A, S1, S12C, S13, S14B, S14C).

Example:

```bash
export PYTHONPATH=$PWD/lib      # figprov.py
python3 07_figure_panels/regen_runner.py 07_figure_panels/server/_fig3a_volcano_sq.py \
    --figsize 2.165 2.2 --keep-look --keep-ticks --axes-rect 0.2 0.17 0.97 0.97 --out panels/fig3a_volcano
```

The scripts insert `~/bin` into `sys.path` to import `figprov`; set `PYTHONPATH` to `lib/` (as above) or copy `lib/figprov.py` to `~/bin`. The SHA-256 of every plotting script in this repository equals the script hash recorded in the provenance sidecar of the panel it produced. The panels were then placed in Adobe Illustrator; the layout tooling is not included because it copies micrographs and collaborators' artwork from private Illustrator files and does no computation.

## Reconstructed files

Nine steps were run as interactive commands or from `/tmp` and were not kept as files. Their code was recovered verbatim from the session logs or from the wrapper that wrote them, and is marked `RECONSTRUCTED` in its header. The two reconstructed shell scripts that would overwrite server outputs (`recount_tdp43_datasets.sh`, `deltapsi_annotated_2026-03-29.sh`) only check by default and run only with `--run`:

| File | Origin |
|---|---|
| `02_splicing/intron_retention/join_ir_meta.py` | heredoc in `04_stress_granule_enrichment/_v2_sg_fix.sh` (ran as `/tmp/join_ir_meta.py`) |
| `02_splicing/rmats/rmats_tdp43_datasets_2026-03-28.sh` | interactive rMATS commands of 2026-03-28 (Telegram calls removed) |
| `03_expression_and_pathways/deseq2_tsv_to_csv.sh` | awk one-liner of 2026-07-06 |
| `04_stress_granule_enrichment/gocc_v2.py` | `/tmp/gocc_v2.py` of 2026-07-07 |
| `01_preprocessing/recount_tdp43_datasets.sh` | interactive featureCounts commands of 2026-03-27 (GSE126543) and 2026-03-28 (GSE65973); equal to the command lines in the count-file headers |
| `02_splicing/majiq/deltapsi_annotated_2026-03-29.sh` | interactive MAJIQ deltapsi command of 2026-03-29; equal to the command recorded in each `deltapsi_annotated.tsv` |
| `03_expression_and_pathways/build_table_s2.py` | `/tmp/build_s3.py` of 2026-06-03 (output paths turned into arguments; the `Fig7B_subset` column was added on 2026-10-20) |
| `05_sequence_features_and_motifs/transite/build_transite_tables.py` | `/tmp/gen_transite_tables2.py` of 2026-07-08 and the strict and 3′ UTR-mapped subsets of 2026-07-13 |
| `06_statistics_wetlab/brown_forsythe_eye.py` | written for this release; recomputes the Brown–Forsythe P values reported in the Methods (line #2 P = 9.9 × 10⁻⁶) |

## Reproduction checks (2026-10-20)

The scripts added on 2026-10-20 were run against the server data without writing to it (`ssh ncnp-linux 'python3 -' < script.py`, `ssh ncnp-linux 'bash -s' < script.sh`, which runs the shell scripts in their default `--check` mode). The printed numbers are compared with the revision-3 text, tables and figure panels. The differences found in the first run were corrected in revision 3 (listed after the table); after the corrections every value below agrees.

| Script | Printed | In the article | Result |
|---|---|---|---|
| `02_splicing/count_events_fig4b.py` | rMATS SE/A3SS/A5SS/MXE: DHX8 KD 4,204/599/458/450, GSE65973 204/42/36/28, GSE126543 667/170/86/177; limma: DHX8 KD 36,190 tested, 1,258 retained introns in 1,150 genes, 75 released; GSE65973 and GSE126543 0 retained; MAJIQ de novo/CE/strong CE: DHX8 KD 22/12/1 (Ccdc71, dpsi_mean +0.53), GSE126543 160/24/9 (UNC13A +0.60, STMN2 +0.64), GSE65973 2/0/0 | Fig. 4B, Results, Methods (36,190 tested introns), Fig. 5A and Supplementary Fig. S9E legends | reproduced |
| `02_splicing/rmats/se_direction_binomial.py` | 3,596 of 4,204 (85.5%) in the skipping direction; two-sided exact binomial P = 10^−512.4 (one-sided 10^−512.7) | 86%, P < 10⁻⁵¹² | reproduced |
| `02_splicing/ir_rna_level_effect_size.py` | IR n = 1,149 vs 11,680, mean log2FC +0.156, d = +0.240, Bonferroni P = 1.9 × 10⁻¹⁹, analytic 95% CI 0.179–0.300 (bootstrap 0.189–0.289); SE d = +0.023, P = 0.60; A3SS/A5SS d = −0.020, P = 1.00 | d = +0.240, 95% CI 0.179–0.300 (Methods: large-sample normal approximation), P = 1.9 × 10⁻¹⁹; SE +0.023, 0.60; A3SS/A5SS −0.020, 1.00 | reproduced |
| `02_splicing/rmats/encode_ri_percent_fig7d.py` | SE/RI/RI%: PRPF8 1,756/1,306/35.05, EFTUD2 1,173/546/24.45, SNRNP200 1,900/729/20.87, PRPF6 737/308/22.38, RBFOX2 648/250/22.10, RBM39 5,016/1,601/19.33 | Fig. 7D and Results: 35.1, 24.5, 20.9, 22.4, 22.1, 19.3 | reproduced (all six OK) |
| `03_expression_and_pathways/fig7b_phase_mean_log2fc.py` | curated-subset means E/A +0.4079, B +0.3691, Bact/NTC +0.3416, regulatory +0.2600, post-catalytic −0.0110, catalytic +0.0110 (current DESeq2) | Fig. 7B labels +0.408, +0.369, +0.342, +0.260, −0.011, +0.011 | reproduced |
| `03_expression_and_pathways/build_table_s2.py` | 164 gene–phase rows, 163 genes, phases 41/28/24/10/20/41; UP 122, ns 33, DOWN 8, NA 1; `Fig7B_subset` = yes in 83 rows (UP 57, ns 19, DOWN 7) | Supplementary Table S2 (revision 3) | output byte-identical to `Manuscript/Tables/SupTableS2/Supplementary_Table_S2.tsv` |
| `05_sequence_features_and_motifs/transite/build_transite_tables.py` | 652 strict DEG-UP genes → 459 with a coding transcript → 451 with a 3′ UTR; background 12,829 genes → 11,972 coding → 11,900 with a 3′ UTR | n = 451; background "the 11,900 expressed genes with a coding transcript and a mapped 3′ UTR" | reproduced; all six tables byte-identical to the server and Dropbox copies |
| `01_preprocessing/recount_tdp43_datasets.sh` | six commands equal to the featureCounts header of their output | | reproduced (commands only; featureCounts was not re-run, since the server is read-only here) |
| `02_splicing/majiq/deltapsi_annotated_2026-03-29.sh` | commands equal to the command recorded in `deltapsi_annotated.tsv` of GSE126543, GSE65973 and GSE27218 | | reproduced (commands only) |

Corrected in revision 3 after the first run of these checks:

- Supplementary Table S2. The log2FC, padj and direction columns of the earlier table were byte-identical to a build from the DESeq2 table before the siCont_2 correction (`/mnt/HDD2/_v1_corrupt_snapshot_20260706/downstream_results/deseq2_dhx8_kd_vs_ctrl.tsv`; UP 118, ns 37, DOWN 8, NA 1). The table was rebuilt from the current DESeq2; Ppie, Syf2, Prpf18 and Rbm8a moved from ns to UP. The `Fig7B_subset` column was added.
- MAJIQ dPSI. Ccdc71 +0.52 became +0.53 (Results, Supplementary Fig. S9E legend, and the S9E panel title, which `analysis_sashimi_dhx8_events.py` writes as +0.52 and `rerun_s9_2026-10-20.sh` replaces at render time). UNC13A 0.62 and STMN2 0.65, which were the TDP-43(−) PSI (0.616 and 0.652), became the dPSI 0.60 and 0.64.
- GSE65973. "No detectable de novo junctions" became two significant de novo junctions (Ppfia2 and Cnksr2), neither of which meets the CE criteria.
- Fig. 7D and Results. PRPF8 35.0% and EFTUD2 24.4% became 35.1% and 24.5% (the panel values in `_fig7b_encode_v2.py` were corrected and the panel re-rendered).
- Methods. The MAJIQ significance rule now includes the mean |dPSI| ≥ 0.1 floor that the code applies (the code takes probability_changing ≥ 0.95, with MAJIQ's default changing threshold of 0.2; the de novo, CE and strong-CE counts are the same with P > 0.95 alone); the Transite background, the Cohen's d interval, the tested-intron rule and the feature counting (featureCounts, not htseq-count) are described as in the code; the software versions are those listed below.

## Software

Versions installed in the environments where the scripts ran (queried on 2026-10-19, and recorded in the figure provenance sidecars). Since revision 3 the Methods and Supplementary Table S3 give these versions for the tools they cite.

| Software | Version | Used for |
|---|---|---|
| Trimmomatic | 0.39 | read trimming |
| STAR | 2.7.10a | alignment |
| Subread featureCounts | 2.0.3 | gene, intron and exon counts |
| samtools | 1.13 | BAM sorting, indexing, depth |
| R | 4.3.2 | |
| DESeq2 | 1.42.1 | differential expression |
| limma / edgeR | 3.58.1 / 4.0.16 | IR (limma eBayes) |
| fgsea | 1.26.0 | GSEA |
| MSigDB gene sets | v2024.1 (Mm, Hs), as named in `fgsea_unified.R` | GSEA |
| transite | 1.20.1 | 3′ UTR k-mer analysis |
| GenomicFeatures / Biostrings / Rsamtools | 1.54.4 / 2.70.3 / 2.18.0 | 3′ UTR extraction |
| rMATS (rmats-turbo) | 4.3.0 (conda env `rmats`; the 2026-03-28 runs used a source build in `/mnt/HDD2/rmats-turbo`) | splicing events |
| MAJIQ | 3.0.19.dev1, a development build (`3.0.19.dev1+gde48918a`) | LSVs, de novo junctions |
| Python | 3.8.15, 3.10.12, 3.13.11 (server); 3.13.12 (local) | see the next table |
| scipy | 1.10.1, 1.15.3, 1.17.0 (server); 1.17.1 (local) | statistical tests |
| statsmodels | 0.14.6 (local; the same version is installed on the server) | ANOVA, Tukey and Šídák tests in `06_statistics_wetlab/prism_stats.py` |
| pysam | 0.23.3 | intron sequences (`analysis_motif_enrichment.py`), BAM reading (`sashimi_ce_gse126543_recovered.py`) |
| g:Profiler client (gprofiler-official) | 1.0.0 (web service queried 2026-07-07) | GO:CC |
| GraphPad Prism | 9 | wet-lab graphs and statistics |
| ImageJ (NIH) | | image measurements |

Python environments:

| Environment | Python | Packages | Scripts |
|---|---|---|---|
| server, conda `clipper3` (miniconda3) | 3.8.15 | pysam 0.23.3, numpy 1.24.4, pandas 2.0.3, scipy 1.10.1 | `analysis_motif_enrichment.py` |
| server, miniconda3 base | 3.13.11 | numpy 2.4.2, pandas 3.0.1, scipy 1.17.0 | `sg_permutation_test_lenmatched.py` via `run_sg_perm.sh` |
| server, system `python3` | 3.10.12 | matplotlib 3.10.8, matplotlib-venn 1.1.2, numpy 2.2.6, pandas 2.3.3, scipy 1.15.3, statsmodels 0.14.6, seaborn 0.13.2, openpyxl 3.1.5, pysam 0.23.3, gprofiler-official 1.0.0 | other server Python scripts and all server panels |
| local (macOS) | 3.13.12 | matplotlib 3.10.8, matplotlib-venn 1.1.2, numpy 2.4.3, pandas 3.0.1, scipy 1.17.1, statsmodels 0.14.6, openpyxl | `06_statistics_wetlab/`, local panels |

`environment/` lists the same versions as pip and R requirement files.

## Data

- DHX8 knockdown RNA-seq (this study): GEO **GSE337869** (private until publication).
- dNLS-TDP-43 mouse cortex: GEO **GSE65973**.
- TDP-43-positive and -negative sorted neuronal nuclei: GEO **GSE126543**.
- ENCODE shRNA knockdowns (Supplementary Table S3), knockdown experiment and its control experiment: PRPF8 ENCSR998MZP (control ENCSR856ZRV), SNRNP200 ENCSR003LSA (ENCSR104ABF), EFTUD2 ENCSR620OKS (ENCSR689PHN), PRPF6 ENCSR529QEZ (ENCSR491FOC) in HepG2; RBFOX2 ENCSR336DFS (ENCSR667PLJ), RBM39 ENCSR678WOA (ENCSR164MUK) in K562. The file accessions of the BAMs used (two knockdown and two control replicates per factor) are in `01_preprocessing/encode/bam_metadata.tsv`.
- Published SG transcriptome and proteome tables: Supplementary Table S5 (supplementary files of the cited studies; GSE99304 for Khong et al.).
- Reference annotation: UCSC mm10 and hg38 genomes and `ncbiRefSeq` GTFs; Ensembl GRCh38.109 (ENCODE rMATS); GENCODE vM25 and v38 GFF3 (MAJIQ).
- The GraphPad Prism files with the wet-lab measurements are not in this repository and are available from the corresponding authors.

## Licence

MIT (see `LICENSE`).
