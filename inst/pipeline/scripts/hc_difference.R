
suppressPackageStartupMessages(library(GetoptLong))

GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       "subgroup1=s", "Name of subgroup 1",
	       "subgroup2=s", "Name of subgroup 2",
	       head = "Visualize global difference by Hilbert curve",
	       foot = "",
	       script_name = "-e \"epik.cmd::epik()\" hc_difference")

library(epik)
load_config(config, use_std_dir = TRUE)

if(!(subgroup1 %in% SAMPLE$subgroup)) {
	stop("Cannot find ", subgroup1, " in `SAMPLE$subgroup")
}
if(!(subgroup2 %in% SAMPLE$subgroup)) {
	stop("Cannot find ", subgroup2, " in `SAMPLE$subgroup")
}

pdf(qq("@{PROJECT_DIR}/image/hc_difference_methylation_@{subgroup1}_vs_@{subgroup2}.pdf"), width = 9, height = 8)
gr_meth = hilbert_curve_methylation_difference(subgroup = SAMPLE$subgroup, comparison = c(subgroup1, subgroup2))
dev.off()

if(!is.null(MARKS)) {
	pdf(qq("@{PROJECT_DIR}/image/hc_difference_histone_mark_@{subgroup1}_vs_@{subgroup2}.pdf"), width = 9, height = 8)
	gr_list = lapply(MARKS, function(mk) hilbert_curve_chipseq_difference(mk, subgroup = SAMPLE$subgroup, comparison = c(subgroup1, subgroup2))
	names(gr_list) = MARKS
	dev.off()

	gr_meth$diff = gr_meth$mean_group1 - gr_meth$mean_group2
	gr_list = lapply(gr_list, function(gr) {
		gr$diff = gr$mean_group1 - gr$mean_group2
		gr
	})

	pdf(qq("@{PROJECT_DIR}/image/general_chipseq_association_to_methylation_@{subgroup1}_vs_@{subgroup2}.pdf"), width = 12, height = 10)
	general_chipseq_association_to_methylation(gr_list, gr_meth)
	dev.off()

	pdf(qq("@{PROJECT_DIR}/image/general_chipseq_association_@{subgroup1}_vs_@{subgroup2}.pdf"), width = 12, height = 10)
	general_chipseq_association(gr_list, q = seq(0.1, 0.9, by = 0.1))
	dev.off()
}

