
suppressPackageStartupMessages(library(GetoptLong))

chr = "chr1"
fdr = 0.05
methdiff = 0.1
subgroup = NULL
GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	        "chr=s", "A single chromosome name",
	        "subgroup=s{1,}", "subgroups for test",
	        "fdr=f", "cutoff for FDRs",
	        "methdiff=f", "cutoff for methylation difference",
	       head = "Filter significant CRs and merge overlapping CRs",
	       foot = "",
	       script_name = "-e \"epik.cmd::epik()\" cr_reduce")

library(epik)
load_epik_config(config)

if(is.null(subgroup)) {
	subgroup = unique(SAMPLE$subgroup)
}

SAMPLE = SAMPLE[SAMPLE$subgroup %in% subgroup, , drop = FALSE]
unique_subgroup = unique(SAMPLE$subgroup)

subgroup_str = paste(sort(unique_subgroup), collapse = "_")

cr = readRDS(qq("@{PROJECT_DIR}/rds_cr/@{subgroup_str}/cr_@{chr}_@{subgroup_str}.rds"))

cr = get_sig_cr(cr, fdr = fdr, meth_diff = methdiff)
cr2 = cr_reduce(cr, TXDB, EXPR)
saveRDS(cr2, file = qq("@{PROJECT_DIR}/rds_cr/@{subgroup_str}/cr_reduce_@{chr}_fdr_@{fdr}_methdiff_@{methdiff}.rds"))
