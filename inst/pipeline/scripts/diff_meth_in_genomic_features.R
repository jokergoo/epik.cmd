
suppressPackageStartupMessages(library(GetoptLong))
GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       head = "Differential methylation in genomic features",
	       foot = "",
	       script_name = "-e \"epik.cmd::epik()\" diff_meth_in_genomic_features")

library(epik)
load_config(config, use_std_dir = TRUE)

sample_id = rownames(SAMPLE)

library(ComplexHeatmap)
ha = HeatmapAnnotation(df = SAMPLE, col = COLOR, show_annotation_name = TRUE)

gr_list = get_mean_methylation_in_genomic_features(sample_id, chromosome = CHROMOSOME, genomic_features = GENOMIC_FEATURE_LIST)
gr_name = names(gr_list)
saveRDS(gr_list, file = qq("@{OUTPUT_DIR}/rds/mean_meth_gr_list.rds"))

pdf(qq("@{PROJECT_DIR}/image/heatmap_diff_methylation_in_genomic_features.pdf"), width = 16, height = 16)
for(i in seq_along(gr_list)) {
	qqcat("making heatmap for @{gr_name[i]}\n")
	heatmap_diff_methylation_in_genomic_features(gr_list[[i]], subgroup = SAMPLE$subgroup, column_title = gr_name[i], ha = ha)
}
dev.off()
