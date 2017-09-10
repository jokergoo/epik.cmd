
suppressPackageStartupMessages(library(GetoptLong))
GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       head = "Make plots for global methylation patterns.",
	       foot = "",
	       script_name = "-e \"epik.cmd::epik()\" general_methylation")

library(epik)
load_config(config)

sample_id = rownames(SAMPLE)


#### density plots ####

library(ComplexHeatmap)
ha = HeatmapAnnotation(df = SAMPLE, col = COLOR, show_annotation_name = TRUE)

pdf(qq("@{PROJECT_DIR}/image/global_methylation_distribution.pdf"), width = 8, height = 8)
methylation_global_distribution(sample_id, subgroup = SAMPLE$subgroup, ha = ha)
dev.off()
pdf(qq("@{PROJECT_DIR}/image/global_methylation_distribution_cgi.pdf"), width = 8, height = 8)
methylation_global_distribution(sample_id, subgroup = SAMPLE$subgroup, ha = ha, background = GENOMIC_FEATURE_LIST$cgi)
dev.off()
pdf(qq("@{PROJECT_DIR}/image/global_methylation_distribution_cgi_shore.pdf"), width = 8, height = 8)
methylation_global_distribution(sample_id, subgroup = SAMPLE$subgroup, ha = ha, background = GENOMIC_FEATURE_LIST$cgi_shore)
dev.off()


#### QC plots #####

pdf(qq("@{PROJECT_DIR}/image/basic_WGBS_qc.pdf"), width = 9, height = 6)
methylation_qcplot(sample_id, chromosome = CHROMOSOME)
dev.off()
pdf(qq("@{PROJECT_DIR}/image/basic_WGBS_qc_in_cgi.pdf"), width = 9, height = 6)
methylation_qcplot(sample_id, chromosome = CHROMOSOME, background = GENOMIC_FEATURE_LIST$cgi)
dev.off()
pdf(qq("@{PROJECT_DIR}/image/basic_WGBS_qc_in_cgi_shore.pdf"), width = 9, height = 6)
methylation_qcplot(sample_id, chromosome = CHROMOSOME, background = GENOMIC_FEATURE_LIST$cgi_shore)
dev.off()


#### gtrellis plot ####

pdf(qq("@{PROJECT_DIR}/image/methylation_gtrellis.pdf"), width = 10, height = 12)
methylation_gtrellis_multiple_samples(sample_id, subgroup = SAMPLE$subgroup, nrow = 3, compact = TRUE)
dev.off()
