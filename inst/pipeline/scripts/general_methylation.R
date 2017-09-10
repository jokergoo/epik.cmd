
suppressPackageStartupMessages(library(GetoptLong))
GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       head = "Calculate correlation between methylation and gene expression",
	       foot = "",
	       script_name = "-e \"epik.cmd::epik()\" general_methylation")

library(epik)
load_epik_config(config)
initialize_project_directory(PROJECT_DIR)

library(ComplexHeatmap)

sample_id = rownames(SAMPLE)

ha = HeatmapAnnotation(df = SAMPLE, col = COLOR, show_annotation_name = TRUE)

cat("general methylation distribution whole genome-wide\n")
pdf(qq("@{PROJECT_DIR}/image/general_methylation_distribution.pdf"), width = 8, height = 8)
global_methylation_distribution(sample_id = sample_id, subgroup = SAMPLE$subgroup, ha = ha, plot_cov = TRUE)
dev.off()

cat("general methylation distribution in cpg islands\n")
pdf(qq("@{PROJECT_DIR}/image/general_methylation_distribution_cgi.pdf"), width = 8, height = 8)
global_methylation_distribution(sample_id = sample_id, subgroup = SAMPLE$subgroup, 
	background = GENOMIC_FEATURE_LIST$cgi, ha = ha, plot_cov = TRUE)
dev.off()


cat("general methylation distribution in cgi shores\n")
pdf(qq("@{PROJECT_DIR}/image/general_methylation_distribution_cgi_shores.pdf"), width = 8, height = 8)
shore = GENOMIC_FEATURE_LIST$cgi_shore
global_methylation_distribution(sample_id = sample_id, subgroup = SAMPLE$subgroup, 
	background = shore, ha = ha, plot_cov = TRUE)
dev.off()

# cat("general methylation distribution in neither cgi nor shores\n")
# pdf(qq("@{PROJECT_DIR}/image/general_methylation_distribution_neither_cgi_nor_shores.pdf"), width = 10, height = 10)
# chromInfo = getChromInfoFromUCSC(GENOME)
# chromInfo = chromInfo[chromInfo$chrom %in% CHROMOSOME, ]
# chromGr = GRanges(chromInfo$chrom, ranges = IRanges(rep(1, nrow(chromInfo)), chromInfo$length))
# complement = setdiff(chromGr, union(GENOMIC_FEATURE_LIST$cgi, GENOMIC_FEATURE_LIST$cgi_shore))
# global_methylation_distribution(sample_id = sample_id, subgroup = SAMPLE$subgroup, chromosome = CHROMOSOME, 
# 	background = complement, ha = ha)
# dev.off()

pdf(qq("@{PROJECT_DIR}/image/basic_WGBS_qc.pdf"), width = 9, height = 6)
wgbs_qcplot(sample_id, chromosome = CHROMOSOME)
dev.off()

pdf(qq("@{PROJECT_DIR}/image/basic_WGBS_qc_in_cgi.pdf"), width = 9, height = 6)
wgbs_qcplot(sample_id, chromosome = CHROMOSOME, background = GENOMIC_FEATURE_LIST$cgi)
dev.off()

pdf(qq("@{PROJECT_DIR}/image/basic_WGBS_qc_in_cgi_shore.pdf"), width = 9, height = 6)
wgbs_qcplot(sample_id, chromosome = CHROMOSOME, background = GENOMIC_FEATURE_LIST$cgi_shore)
dev.off()

# pdf(qq("@{PROJECT_DIR}/image/basic_WGBS_qc_neither_cgi_nor_shores.pdf"), width = 9, height = 6)
# wgbs_qcplot(sample_id, chromosome = CHROMOSOME, background = complement)
# dev.off()


pdf(qq("@{PROJECT_DIR}/image/methylation_gtrellis.pdf"), width = 10, height = 12)
gtrellis_methylation_for_multiple_samples(sample_id, subgroup = SAMPLE$subgroup, nrow = 3, compact = TRUE)
dev.off()
