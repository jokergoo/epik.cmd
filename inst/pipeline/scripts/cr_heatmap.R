
suppressPackageStartupMessages(library(GetoptLong))

subgroup = NULL
fdr = NULL
methdiff = NULL
GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       "subgroup=s{1,}", "subgroups for test",
	       "fdr=f", "cutoff for fdr",
	       "methdiff=f", "cutoff for methylation difference",
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
	cr = readRDS(qq("@{PROJECT_DIR}/rds_cr/@{subgroup_str}/cr_reduce_@{chr}_fdr_@{fdr}_methdiff_@{methdiff}.rds"))
	cr_all = cr_concatenate(cr_all, cr)
}

library(ComplexHeatmap)
ha = HeatmapAnnotation(df = SAMPLE, col = COLOR, show_annotation_name = TRUE)

dir.create(qq("@{PROJECT_DIR}/image/@{subgroup_str}"), showWarnings = FALSE)

pdf(qq("@{PROJECT_DIR}/image/@{subgroup_str}/cr_heatmap_fdr_@{fdr}_methdiff_@{methdiff}_@{subgroup_str}.pdf"), width = 12, height = 15)
sig_cr_heatmap(cr_all, TXDB, EXPR, ha = ha,
	gf_list = list(CGI = GENOMIC_FEATURES$cgi, shore = GENOMIC_FEATURES$cgi_shore))
dev.off()
