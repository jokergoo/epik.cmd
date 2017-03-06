
suppressPackageStartupMessages(library(GetoptLong))

GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       "chr=s", "chromosome name",
	       "fdr=f", "cutoff for fdr",
	       "methdiff=f", "cutoff for methylation difference",
	       head = "Visualize correlation between epigenomic signals by Gviz",
	       foot = "",
	       script_name = "-e \"epik.cmd::epik()\" cr_gviz")

library(epik)
load_config(config, use_std_dir = TRUE)

sig_cr = readRDS(qq("@{PROJECT_DIR}/rds_cr/cr_reduce_@{chr}_fdr_@{fdr}_methdiff_@{methdiff}.rds"))

dir.create(qq("@{PROJECT_DIR}/gviz/@{chr}"), recursive = TRUE, showWarnings = FALSE)

if(length(MARK)) {
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
	pdf(qq("@{PROJECT_DIR}/gviz/@{chr}/cr_gviz_@{gi}_@{gn[gi]}.pdf"), width = 8, height = ((sm*st+0.2*ns+4)*0.5+8 + 0.12*nt)*0.5)
	cr_gviz(sig_cr, gi, EXPR, TXDB, gf_list = list(CGI = GENOMIC_FEATURES$cgi), hm_list = peak_list, title = qq("@{gi} (@{gn[gi]})"))
	dev.off()
}
