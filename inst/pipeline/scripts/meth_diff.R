 
suppressPackageStartupMessages(library(GetoptLong))
windowsize = 1000
gf = NULL
subgroup = NULL
GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       "windowsize=i", "window size to segment each region",
	       "subgroup=s{1,}", "subgroups for test",
	       "gf=s{1,}", "names of genomic features",
	       head = "Differential methylation in genomic features",
	       foot = "",
	       script_name = "-e \"epik.cmd::epik()\" diff_meth_in_genomic_features")

library(epik)
load_epik_config(config)

if(!is.null(gf)) {
	GENOMIC_FEATURE_LIST = GENOMIC_FEATURE_LIST[intersect(names(GENOMIC_FEATURE_LIST), gf)]
}


if(is.null(subgroup)) {
	subgroup = unique(SAMPLE$subgroup)
}

SAMPLE = SAMPLE[SAMPLE$subgroup %in% subgroup, , drop = FALSE]
if(length(unique(SAMPLE$subgroup)) <= 1) {
	stop("Number of subgroups is less than 2.")
}

unique_subgroup = unique(SAMPLE$subgroup)

sample_id = rownames(SAMPLE)

library(ComplexHeatmap)
ha = HeatmapAnnotation(df = SAMPLE, col = COLOR, show_annotation_name = TRUE)

library(EnrichedHeatmap)
GENOMIC_FEATURE_LIST = lapply(GENOMIC_FEATURE_LIST, function(gr) {
	gr = makeWindows(gr, w = windowsize, short.keep = TRUE)
	gr[width(gr) >= windowsize/2]
})

gr_list = get_mean_methylation_in_genomic_features(sample_id, chromosome = CHROMOSOME, genomic_features = GENOMIC_FEATURE_LIST)
gr_name = names(gr_list)
# saveRDS(gr_list, file = qq("@{OUTPUT_DIR}/rds/mean_meth_gr_list.rds"))

subgroup_str = paste(sort(unique_subgroup), collapse = "_")

pdf(qq("@{PROJECT_DIR}/image/heatmap_diff_methylation_in_genomic_features_@{subgroup_str}.pdf"), width = 16, height = 16)
for(i in seq_along(gr_list)) {
	qqcat("making heatmap for @{gr_name[i]}\n")
	heatmap_diff_methylation_in_genomic_features(gr_list[[i]], subgroup = SAMPLE$subgroup, column_title = gr_name[i], ha = ha)
}
dev.off()
