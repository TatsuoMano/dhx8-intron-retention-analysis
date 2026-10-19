# RECONSTRUCTED FILE (code release, 2026-10-19).
# This script ran on the server as /tmp/join_ir_meta.py (2026-07-06/07), written by the heredoc in
# 04_stress_granule_enrichment/_v2_sg_fix.sh and called again by _v2_sg_fix2.sh. The /tmp copy no
# longer exists; the body below is copied verbatim from that heredoc.
# Usage: python3 join_ir_meta.py <limma_*_dhx8_s2.tsv> <intron_metadata.tsv> <out_13col.tsv>
import pandas as pd, sys
ir = pd.read_csv(sys.argv[1], sep='\t')
meta = pd.read_csv(sys.argv[2], sep='\t')
addc = ['chr','start','end','strand','intron_length']
m = ir.merge(meta[['intron_id']+addc], on='intron_id', how='left')
cols = ['intron_id','logFC','avg_expr','t_stat','pvalue','adj_p','chr','start','end','strand','gene_name','intron_length','retained']
m = m[[c for c in cols if c in m.columns]]
nmiss = int(m['intron_length'].isna().sum())
m.to_csv(sys.argv[3], sep='\t', index=False)
print(f"  joined rows={len(m)} intron_length_missing={nmiss}")
