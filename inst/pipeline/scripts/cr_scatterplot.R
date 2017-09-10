suppressPackageStartupMessages(library(GetoptLong))

GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       "chr=s", "chromosome name",
	       "fdr=f", "cutoff for fdr",
	       "methdiff=f", "cutoff for methylation difference",
	       head = "Visualize correlation between epigenomic signals by Gviz",
	       foot = "",
	       script_name = "-e \"epik.cmd::epik()\" cr_scatterplot")

library(epik)
load_epik_config(config)
initialize_project_directory(PROJECT_DIR)

sig_cr = readRDS(qq("@{PROJECT_DIR}/rds_cr/cr_reduce_@{chr}_fdr_@{fdr}_methdiff_@{methdiff}.rds"))

gn = extract_field_from_gencode(GTF_FILE, level = "gene", primary_key = "gene_id", field = "gene_name")

sig_cr$gene_symbol = gn[sig_cr$gene_id]

pdf(qq("@{PROJECT_DIR}/image/cr_scatterplot_@{chr}.pdf"))
cr_scatterplot(sig_cr, EXPR, text_column = c("gene_symbol", "gene_tss_dist", "merged_windows", "corr_p"))
dev.off()
