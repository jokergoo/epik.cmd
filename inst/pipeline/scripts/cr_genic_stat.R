
suppressPackageStartupMessages(library(GetoptLong))

GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       "fdr=f", "cutoff for fdr",
	       "methdiff=f", "cutoff for methylation difference",
	       head = "Genic statistics",
	       foot = "",
	       script_name = "-e \"epik.cmd::epik()\" cr_genic_stat")

library(epik)
load_config(config, use_std_dir = TRUE)

cr_all = NULL
for(chr in CHROMOSOME) {
	cr = readRDS(qq("@{PROJECT_DIR}/rds_cr/cr_reduce_@{chr}_fdr_@{fdr}_methdiff_@{methdiff}.rds"))
	cr_all = cr_concatenate(cr_all, cr)
}

pdf(qq("@{PROJECT_DIR}/image/cr_genic_stat.pdf"), width = 15, height = 4)
cr_genic_stat(cr_all, TXDB)
dev.off()
