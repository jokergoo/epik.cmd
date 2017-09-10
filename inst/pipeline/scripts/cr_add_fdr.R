
suppressPackageStartupMessages(library(GetoptLong))

GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       head = "Calculate FDRs",
	       foot = "CRs from all chromosomes are concatenated first and FDRs for the correlation test
	              as well for the ANOVA test (if there are more than one subgroups) are calculated. 
	              In the end, the RDS files are updated for each chromosome.",
	       script_name = "-e \"epik.cmd::epik()\" cr_add_fdr")

library(epik)
load_epik_config(config)
initialize_project_directory(PROJECT_DIR)

cr_all = NULL
for(chr in CHROMOSOME) {
	cr = readRDS(qq("@{PROJECT_DIR}/rds_cr/cr_@{chr}.rds"))
	cr_all = cr_concatenate(cr_all, cr)
}

cr_all = cr_add_fdr_column(cr_all)
saveRDS(cr_all, file = qq("@{PROJECT_DIR}/rds_cr/cr_all.rds"))
for(chr in CHROMOSOME) {
	saveRDS(cr_all[seqnames(cr_all) == chr], file = qq("@{PROJECT_DIR}/rds_cr/cr_@{chr}.rds"))
}
