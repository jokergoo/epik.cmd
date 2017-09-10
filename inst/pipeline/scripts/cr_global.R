
suppressPackageStartupMessages(library(GetoptLong))

subgroup = NULL
GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       "subgroup=s{1,}", "subgroups for test",
	       head = "Visualize correlation between epigenomic signals by Gviz",
	       foot = "",
	       script_name = "-e \"epik.cmd::epik()\" cr_gviz")

library(epik)
load_epik_config(config)


if(is.null(subgroup)) {
	subgroup = unique(SAMPLE$subgroup)
}

SAMPLE = SAMPLE[SAMPLE$subgroup %in% subgroup, , drop = FALSE]
unique_subgroup = unique(SAMPLE$subgroup)

subgroup_str = paste(sort(unique_subgroup), collapse = "_")

cr = readRDS(qq("@{PROJECT_DIR}/rds_cr/@{subgroup_str}/cr_@{chr}_@{subgroup_str}.rds"))

dir.create(qq("@{PROJECT_DIR}/image/@{subgroup_str}"), showWarnings = FALSE)

pdf("@{PROJECT_DIR}/image/@{subgroup_str}/cr_hilbert_curve_@{subgroup_str}.pdf", width = 8, height = 8)
cr_hilbert_curve(cr_all)
dev.off()

library(ComplexHeatmap)
ha = HeatmapAnnotation(df = SAMPLE, col = COLOR, show_annotation_name = TRUE)

pdf("@{PROJECT_DIR}/image/@{subgroup_str}/cr_enriched_heatmap_@{subgroup_str}.pdf", width = 12, height = 18)
cr = cr_enriched_heatmap(cr, TXDB, EXPR, ha)
dev.off()

pdf("@{PROJECT_DIR}/image/@{subgroup_str}/sig_cr_compare_cutoff_@{subgroup_str}.pdf", width = 12, height = 8)
sig_cr_compare_cutoff(cr, TXDB)
dev.off()
