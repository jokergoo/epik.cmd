suppressPackageStartupMessages(library(GetoptLong))

subgroup = NULL
fdr = NULL
methdiff = NULL
GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       "chr=s", "chromosome name",
	       "subgroup=s{1,}", "subgroups for test",
	       "fdr=f", "cutoff for fdr",
	       "methdiff=f", "cutoff for methylation difference",
	       head = "Visualize correlation between epigenomic signals by Gviz",
	       foot = "",
	       script_name = "-e \"epik.cmd::epik()\" cr_scatterplot")

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

sig_cr = readRDS(qq("@{PROJECT_DIR}/rds_cr/@{subgroup_str}/cr_reduce_@{chr}_fdr_@{fdr}_methdiff_@{methdiff}_@{subgroup_str}.rds"))

gn = extract_field_from_gencode(GTF_FILE, level = "gene", primary_key = "gene_id", field = "gene_name")

sig_cr$gene_symbol = gn[sig_cr$gene_id]

dir.create(qq("@{PROJECT_DIR}/image/@{subgroup_str}/scatterplot"), recursive = TRUE,  showWarnings = FALSE)

pdf(qq("@{PROJECT_DIR}/image/@{subgroup_str}/scatterplot/cr_scatterplot_@{chr}_fdr_@{fdr}_methdiff_@{methdiff}_@{subgroup_str}.pdf"))
cr_scatterplot(sig_cr, EXPR, text_column = c("gene_symbol", "gene_tss_dist", "merged_windows", "corr_p"))
dev.off()
