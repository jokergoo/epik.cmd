
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

cr_all = NULL
for(chr in CHROMOSOME) {
	cr = readRDS(qq("@{PROJECT_DIR}/rds_cr/cr_reduce_@{chr}_fdr_@{fdr}_methdiff_@{methdiff}.rds"))
	cr_all = cr_concatenate(cr_all, cr)
}

library(ComplexHeatmap)
ha = HeatmapAnnotation(df = SAMPLE, col = COLOR, show_annotation_name = TRUE)

pdf(qq("@{PROJECT_DIR}/image/sig_cr_heatmap_fdr_@{fdr}_methdiff_@{methdiff}.pdf"), width = 12, height = 15)
sig_cr_heatmap(cr_all, TXDB, EXPR, ha = ha,
	gf_list = list(CGI = GENOMIC_FEATURES$cgi, shore = GENOMIC_FEATURES$cgi_shore))
dev.off()
