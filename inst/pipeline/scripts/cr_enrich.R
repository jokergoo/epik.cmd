
suppressPackageStartupMessages(library(GetoptLong))

subgroup = NULL
fdr = NULL
methdiff = NULL
GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       "subgroup=s{1,}", "subgroups for test",
	       "fdr=f", "cutoff for fdr",
	       "methdiff=f", "cutoff for methylation difference",
	       head = "Visualize correlation between epigenomic signals by EnrichedHeatmap",
	       foot = "",
	       script_name = "-e \"epik.cmd::epik()\" cr_enriched_heatmap")

library(epik)
load_epik_config(config)


if(is.null(subgroup)) {
	subgroup = unique(SAMPLE$subgroup)
}

SAMPLE = SAMPLE[SAMPLE$subgroup %in% subgroup, , drop = FALSE]
unique_subgroup = unique(SAMPLE$subgroup)
subgroup_str = paste(sort(unique_subgroup), collapse = "_")

if(is.null(fdr) || is.null(methdiff)) {
	files = file.list(qq("@{PROJECT_DIR}/rds_cr/@{subgroup_str}/"))
	files_reduce = grep("cr_reduce_.*?_fdr_(.*?)_methdiff_(.*?)\\.rds", files, value = TRUE)
	if(length(files_reduce)) {
		fdr_detected = as.numeric(unique(gsub("cr_reduce_.*?_fdr_(.*?)_methdiff_.*?\\.rds", files_reduce)))
		methdiff_detected = as.numeric(unique(gsub("cr_reduce_.*?_fdr_.*?_methdiff_(.*?)\\.rds", files_reduce)))

		if(length(fdr_detected) == 1 && length(methdiff_detected) == 1) {
			fdr = fdr_detected
			methdiff = methdiff_detected
		} else {
			stop("`fdr` and `methdiff` should be specified.")
		}
	} else {
		stop("`fdr` and `methdiff` should be specified.")
	}
}



cr_all = NULL
for(chr in CHROMOSOME) {
	cr = readRDS(qq("@{PROJECT_DIR}/rds_cr/@{subgroup_str}/cr_@{chr}.rds"))
	cr_all = cr_concatenate(cr_all, cr)
}

dir.create(qq("@{PROJECT_DIR}/image/@{subgroup_str}"), showWarnings = FALSE)

ha = HeatmapAnnotation(df = SAMPLE, col = COLOR)

pdf(qq("@{PROJECT_DIR}/image/@{subgroup_str}/cr_enriched_heatmap_at_cgi_fdr_@{fdr}_methdiff_@{methdiff}.pdf"), width = 10 + length(MARKS)*3, height = 10)
cr_enriched_heatmap_at_cgi(cr_all, TXDB, EXPR, GENOMIC_FEATURES$cgi, marks = MARKS,
	fdr_cutoff = fdr, meth_diff_cutoff = methdiff, type = "neg", expr_annotation = ha)
cr_enriched_heatmap_at_cgi(cr_all, TXDB, EXPR, GENOMIC_FEATURES$cgi, marks = MARKS,
	fdr_cutoff = fdr, meth_diff_cutoff = methdiff, type = "pos", expr_annotation = ha)
dev.off()

pdf(qq("@{PROJECT_DIR}/image/@{subgroup_str}/cr_enriched_heatmap_at_tss_fdr_@{fdr}_methdiff_@{methdiff}.pdf"), width = 10 + length(MARKS)*3, height = 10)
cr_enriched_heatmap_at_tss(cr_all, TXDB, EXPR, GENOMIC_FEATURES$cgi, marks = MARKS,
	fdr_cutoff = fdr, meth_diff_cutoff = methdiff, type = "neg", expr_annotation = ha)
cr_enriched_heatmap_at_tss(cr_all, TXDB, EXPR, GENOMIC_FEATURES$cgi, marks = MARKS,
	fdr_cutoff = fdr, meth_diff_cutoff = methdiff, type = "pos", expr_annotation = ha)
dev.off()

pdf(qq("@{PROJECT_DIR}/image/@{subgroup_str}/cr_enriched_heatmap_at_gene_fdr_@{fdr}_methdiff_@{methdiff}.pdf"), width = 10 + length(MARKS)*3, height = 10)
cr_enriched_heatmap_at_gene(cr_all, TXDB, EXPR, GENOMIC_FEATURES$cgi, marks = MARKS,
	k = 1, expr_annotation = ha)
cr_enriched_heatmap_at_gene(cr_all, TXDB, EXPR, GENOMIC_FEATURES$cgi, marks = MARKS,
	k = 4, expr_annotation = ha)
dev.off()
