
suppressPackageStartupMessages(library(GetoptLong))

GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       "fdr=f", "cutoff for fdr",
	       "methdiff=f", "cutoff for methylation difference",
	       head = "Visualize correlation between epigenomic signals by Gviz",
	       foot = "",
	       script_name = "-e \"epik.cmd::epik()\" cr_gviz")

library(epik)
load_epik_config(config)
initialize_project_directory(PROJECT_DIR)

cr = readRDS(qq("@{PROJECT_DIR}/rds_cr/cr_all.rds"))
	
pdf("@{PROJECT_DIR}/image/cr_hilbert_curve.pdf", width = 8, height = 8)
cr_hilbert_curve(cr_all)
dev.off()


library(ComplexHeatmap)
ha = HeatmapAnnotation(df = SAMPLE, col = COLOR, show_annotation_name = TRUE)

pdf("@{PROJECT_DIR}/image/cr_enriched_heatmap.pdf", width = 12, height = 18)
cr = cr_enriched_heatmap(cr, TXDB, EXPR, ha)
dev.off()

pdf("@{PROJECT_DIR}/image/sig_cr_compare_cutoff.pdf", width = 12, height = 8)
sig_cr_compare_cutoff(cr, TXDB)
dev.off()
