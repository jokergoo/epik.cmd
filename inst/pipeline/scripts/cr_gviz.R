
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

dir.create(qq("@{PROJECT_DIR}/image/@{subgroup_str}/gviz/@{chr}"), recursive = TRUE, showWarnings = FALSE)

if(length(MARKS)) {
	peak_list = lapply(MARKS, get_peak_list, chr = chr)
	names(peak_list) = MARKS
} else {
	peak_list = NULL
}

tx_list = transcriptsBy(TXDB, "gene")
gn = extract_field_from_gencode(GTF_FILE, level = "gene", primary_key = "gene_id", field = "gene_name")

cr_param = metadata(sig_cr)
subgroup = cr_param$subgroup
subgroup_level = unique(subgroup)

for(gi in unique(sig_cr$gene_id)) {
	# (sm*st+0.2*ns+4)*0.5+8 + 0.12*nt
	sm = length(MARKS)
	st = length(subgroup_level)
	ns = nrow(SAMPLE)
	nt = length(tx_list[[gi]])
	pdf(qq("@{PROJECT_DIR}/image/@{subgroup_str}/gviz/@{chr}/cr_gviz_@{gi}_@{gn[gi]}_fdr_@{fdr}_methdiff_@{methdiff}_@{subgroup_str}.pdf"), width = 8, height = ((sm*st+0.2*ns+4)*0.5+8 + 0.12*nt)*0.5)
	cr_gviz(sig_cr, gi, EXPR, TXDB, gf_list = list(CGI = GENOMIC_FEATURE_LIST$cgi), hm_list = peak_list, title = qq("@{gi} (@{gn[gi]}), @{as.vector(strand(tx_list[[gi]]))[1]} strand"))
	dev.off()
}
