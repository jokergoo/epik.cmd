
suppressPackageStartupMessages(library(GetoptLong))

subgroup = NULL
fdr = NULL
methdiff = NULL
where = "tss"
GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       "subgroup=s{1,}", "subgroups for test",
	       "fdr=f", "cutoff for fdr",
	       "methdiff=f", "cutoff for methylation difference",
	       "where=s", "target regions to enrich. Possible values are 'tss', 'gene' and 'tss_cgi'",
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
	files = list.files(qq("@{PROJECT_DIR}/rds_cr/@{subgroup_str}/"))
	files_reduce = grep("cr_reduce_.*?_fdr_(.*?)_methdiff_(.*?)\\.rds", files, value = TRUE)
	if(length(files_reduce)) {
		fdr_detected = as.numeric(unique(gsub("cr_reduce_.*?_fdr_(.*?)_methdiff_.*?\\.rds", "\\1", files_reduce)))
		methdiff_detected = as.numeric(unique(gsub("cr_reduce_.*?_fdr_.*?_methdiff_(.*?)_.*\\.rds", "\\1", files_reduce)))


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

# cr_all = NULL
# for(chr in CHROMOSOME) {
# 	cr = readRDS(qq("@{PROJECT_DIR}/rds_cr/@{subgroup_str}/cr_@{chr}_@{subgroup_str}.rds"))
# 	cr_all = cr_concatenate(cr_all, cr)
# }

if(!file.exists(qq("@{PROJECT_DIR}/rds_cr/@{subgroup_str}/cr_all_@{subgroup_str}.rds"))) {
	stop("You need to run 'Rscript -e \"epik.cmd:epik()\" cr_global' first.")
}
cr_all = readRDS(qq("@{PROJECT_DIR}/rds_cr/@{subgroup_str}/cr_all_@{subgroup_str}.rds"))
# metadata(cr_all)$cr_param$km = NULL

dir.create(qq("@{PROJECT_DIR}/image/@{subgroup_str}"), showWarnings = FALSE)

library(ComplexHeatmap)
ha = HeatmapAnnotation(df = SAMPLE[, "subgroup", drop = FALSE], col = COLOR)

if(where == "tss_cgi" || where == "cgi") {
	pdf(qq("@{PROJECT_DIR}/image/@{subgroup_str}/cr_enriched_heatmap_at_tss_cgi_fdr_@{fdr}_methdiff_@{methdiff}_@{subgroup_str}.pdf"), width = 10 + length(MARKS)*3, height = 10)
	cat("enrich neg_CR to tss-CGIs\n")
	cr_enriched_heatmap_at_tss_cgi(cr_all, TXDB, EXPR, GENOMIC_FEATURE_LIST$cgi, marks = MARKS,
		fdr_cutoff = fdr, meth_diff_cutoff = methdiff, type = "neg", expr_ha = ha)
	cat("enrich pos_CR to tss-CGIs\n")
	cr_enriched_heatmap_at_tss_cgi(cr_all, TXDB, EXPR, GENOMIC_FEATURE_LIST$cgi, marks = MARKS,
		fdr_cutoff = fdr, meth_diff_cutoff = methdiff, type = "pos", expr_ha = ha)
	dev.off()

} else if(where == "tss") {
	pdf(qq("@{PROJECT_DIR}/image/@{subgroup_str}/cr_enriched_heatmap_at_tss_fdr_@{fdr}_methdiff_@{methdiff}_@{subgroup_str}.pdf"), width = 10 + length(MARKS)*3, height = 10)
	cat("enrich neg_CR to tss\n")
	cr_enriched_heatmap_at_tss(cr_all, TXDB, EXPR, GENOMIC_FEATURE_LIST$cgi, marks = MARKS,
		fdr_cutoff = fdr, meth_diff_cutoff = methdiff, type = "neg", expr_ha = ha)
	cat("enrich pos_CR to tss\n")
	cr_enriched_heatmap_at_tss(cr_all, TXDB, EXPR, GENOMIC_FEATURE_LIST$cgi, marks = MARKS,
		fdr_cutoff = fdr, meth_diff_cutoff = methdiff, type = "pos", expr_ha = ha)
	dev.off()
} else if(where == "gene") {
	pdf(qq("@{PROJECT_DIR}/image/@{subgroup_str}/cr_enriched_heatmap_at_gene_fdr_@{fdr}_methdiff_@{methdiff}_@{subgroup_str}.pdf"), width = 10 + length(MARKS)*3, height = 10)
	cat("enrich CR to genes\n")
	cr_enriched_heatmap_at_gene(cr_all, TXDB, EXPR, GENOMIC_FEATURE_LIST$cgi, K = 1, marks = MARKS, expr_ha = ha)
	cr_enriched_heatmap_at_gene(cr_all, TXDB, EXPR, GENOMIC_FEATURE_LIST$cgi, K = 2, marks = MARKS, expr_ha = ha)
	cr_enriched_heatmap_at_gene(cr_all, TXDB, EXPR, GENOMIC_FEATURE_LIST$cgi, K = 3, marks = MARKS, expr_ha = ha)
	cr_enriched_heatmap_at_gene(cr_all, TXDB, EXPR, GENOMIC_FEATURE_LIST$cgi, K = 4, marks = MARKS, expr_ha = ha)
	dev.off()
}
