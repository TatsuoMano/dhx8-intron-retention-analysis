# RECONSTRUCTED FILE (code release, 2026-10-19).
# This script ran on the server as /tmp/gocc_v2.py on 2026-07-07 14:37 JST and wrote
# gocc_ir_v2_gprofiler.tsv (the values drawn in Supplementary Fig. S12C). The /tmp copy no longer
# exists; the body below is copied verbatim from the session log that recorded the command.
# Requires the g:Profiler web service (gprofiler-official 1.0.0); results depend on the g:Profiler
# database release at the time of the query.
import pandas as pd
from gprofiler import GProfiler
df=pd.read_csv("/mnt/HDD2/TDP43_comparison/validation_20260226/limma_retained_dhx8_s2.tsv",sep="\t")
genes=sorted(df["gene_name"].dropna().unique().tolist())
print("n IR genes:",len(genes))
gp=GProfiler(return_dataframe=True)
res=gp.profile(organism="mmusculus",query=genes,sources=["GO:CC"],
               user_threshold=0.05,significance_threshold_method="g_SCS",
               no_evidences=True)
print("n sig GO:CC terms:",len(res))
res=res.sort_values("p_value")
cols=["native","name","p_value","term_size","query_size","intersection_size"]
# show the granule terms explicitly + top 12
key=res[res["name"].str.contains("granule|stress",case=False,na=False)]
print("=== GRANULE / STRESS terms ===")
print(key[cols].to_string(index=False))
print("=== TOP 12 by p ===")
print(res[cols].head(12).to_string(index=False))
res[cols].to_csv("/mnt/HDD2/unified_results/sg_analysis/gsea_dhx8/gocc_ir_v2_gprofiler.tsv",sep="\t",index=False)
print("SAVED gocc_ir_v2_gprofiler.tsv")
