
suppressPackageStartupMessages(library(GetoptLong))

subgroup = NULL
# email = NULL
GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       "subgroup=s{1,}", "subgroups for test",
	       # "email=s", "email registed on David web service",
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

cr_all = readRDS(qq("@{PROJECT_DIR}/rds_cr/@{subgroup_str}/cr_@{CHROMOSOME[1]}_@{subgroup_str}.rds"))
for(chr in CHROMOSOME[-1]) {
	cr = readRDS(qq("@{PROJECT_DIR}/rds_cr/@{subgroup_str}/cr_@{chr}_@{subgroup_str}.rds"))
	cr_all = cr_concatenate(cr_all, cr)
}

dir.create(qq("@{PROJECT_DIR}/image/@{subgroup_str}"), showWarnings = FALSE)

pdf(qq("@{PROJECT_DIR}/image/@{subgroup_str}/cr_hilbert_curve_@{subgroup_str}.pdf"), width = 8, height = 8)
cr_hilbert_curve(cr_all)
dev.off()

library(ComplexHeatmap)
ha = HeatmapAnnotation(df = SAMPLE[, "subgroup", drop = FALSE], col = COLOR, show_annotation_name = TRUE)

pdf(qq("@{PROJECT_DIR}/image/@{subgroup_str}/cr_enriched_heatmap_@{subgroup_str}.pdf"), width = 12, height = 14)
cr_all = cr_enriched_heatmap(cr_all, TXDB, EXPR, ha)
dev.off()

saveRDS(cr_all, file = qq("@{PROJECT_DIR}/rds_cr/@{subgroup_str}/cr_all_@{subgroup_str}.rds"))

pdf(qq("@{PROJECT_DIR}/image/@{subgroup_str}/sig_cr_compare_cutoff_@{subgroup_str}.pdf"), width = 12, height = 8)
sig_cr_compare_cutoff(cr_all, TXDB)
dev.off()

pdf(qq("@{PROJECT_DIR}/image/@{subgroup_str}/cr_enriched_heatmap_only_sig_@{subgroup_str}.pdf"), width = 10, height = 10)
for(fdr_cutoff in c(0.1, 0.05, 0.01)) {
	for(meth_diff_cutoff in c(0, 0.1, 0.2, 0.3)) {
		qqcat("making enriched heatmap only for sig CRs (fdr = @{fdr_cutoff}, meth_diff_cutoff = @{meth_diff_cutoff})\n")
		sig_cr_enriched_heatmap(cr_all, TXDB, fdr_cutoff = fdr_cutoff, meth_diff_cutoff = meth_diff_cutoff)		
	}
}
dev.off()

pdf(qq("@{PROJECT_DIR}/image/@{subgroup_str}/cr_genes_gtrellis_@{subgroup_str}.pdf"), width = 14, height = 12)
cr_genes_gtrellis(cr_all, TXDB, EXPR)
dev.off()


# ## functional enrichment
# if(!is.null(email)) {
# 	pdf(qq("@{PROJECT_DIR}/image/@{subgroup_str}/cr_genes_enriched_functions_@{subgroup_str}.pdf"), width = 14, height = 12)
# 	try(cr_genes_function_enrichment(cr_all, email))
# 	dev.off()
# }
