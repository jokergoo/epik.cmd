
suppressPackageStartupMessages(library(GetoptLong))

GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       "fdr=f", "cutoff for fdr",
	       "methdiff=f", "cutoff for methylation difference",
	       head = "Visualize correlation between epigenomic signals by EnrichedHeatmap",
	       foot = "",
	       script_name = "-e \"epik.cmd::epik()\" cr_enriched_heatmap")

library(epik)
load_epik_config(config)
initialize_project_directory(PROJECT_DIR)


cr = readRDS(qq("@{PROJECT_DIR}/rds/cr_all.rds"))

ha = HeatmapAnnotation(df = SAMPLE, col = COLOR)

pdf(qq("@{PROJECT_DIR}/image/cr_enriched_heatmap_at_cgi_fdr_@{fdr}_methdiff_@{methdiff}.pdf"), width = 10 + length(MARKS)*3, height = 10)
cr_enriched_heatmap_at_cgi(cr, TXDB, EXPR, GENOMIC_FEATURES$cgi, marks = MARKS,
	fdr_cutoff = fdr, meth_diff_cutoff = methdiff, type = "neg", expr_annotation = ha)
cr_enriched_heatmap_at_cgi(cr, TXDB, EXPR, GENOMIC_FEATURES$cgi, marks = MARKS,
	fdr_cutoff = fdr, meth_diff_cutoff = methdiff, type = "pos", expr_annotation = ha)
dev.off()

pdf(qq("@{PROJECT_DIR}/image/cr_enriched_heatmap_at_tss_fdr_@{fdr}_methdiff_@{methdiff}.pdf"), width = 10 + length(MARKS)*3, height = 10)
cr_enriched_heatmap_at_tss(cr, TXDB, EXPR, GENOMIC_FEATURES$cgi, marks = MARKS,
	fdr_cutoff = fdr, meth_diff_cutoff = methdiff, type = "neg", expr_annotation = ha)
cr_enriched_heatmap_at_tss(cr, TXDB, EXPR, GENOMIC_FEATURES$cgi, marks = MARKS,
	fdr_cutoff = fdr, meth_diff_cutoff = methdiff, type = "pos", expr_annotation = ha)
dev.off()

pdf(qq("@{PROJECT_DIR}/image/cr_enriched_heatmap_at_gene_fdr_@{fdr}_methdiff_@{methdiff}.pdf"), width = 10 + length(MARKS)*3, height = 10)
cr_enriched_heatmap_at_gene(cr, TXDB, EXPR, GENOMIC_FEATURES$cgi, marks = MARKS,
	k = 1, expr_annotation = ha)
cr_enriched_heatmap_at_gene(cr, TXDB, EXPR, GENOMIC_FEATURES$cgi, marks = MARKS,
	k = 4, expr_annotation = ha)
dev.off()
