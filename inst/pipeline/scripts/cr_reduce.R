
suppressPackageStartupMessages(library(GetoptLong))

chr = "chr1"
fdr = 0.05
methdiff = 0.1
GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	        "chr=s", "A single chromosome name",
	        "fdr=f", "cutoff for FDRs",
	        "methdiff=f", "cutoff for methylation difference",
	       head = "Filter significant CRs and merge overlapping CRs",
	       foot = "",
	       script_name = "-e \"epik.cmd::epik()\" cr_reduce")

library(epik)
load_config(config, use_std_dir = TRUE)

cr = readRDS(qq("@{PROJECT_DIR}/rds_cr/cr_@{chr}.rds"))

cr = get_sig_cr(cr, fdr = fdr, meth_diff = methdiff)
cr2 = cr_reduce(cr, TXDB, EXPR)
saveRDS(cr2, file = qq("@{PROJECT_DIR}/rds_cr/cr_reduce_@{chr}_fdr_@{fdr}_methdiff_@{methdiff}.rds"))
