suppressPackageStartupMessages(library(GetoptLong))

chr = "chr1"
extend = 5000
windowsize = 6
windowstep = 3
maxwindow = 10000
subgroup = NULL
GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       "chr=s", "A single chromosome name",
	       "subgroup=s{1,}", "subgroups for test",
	       "extend=i", "Extension of the gene model, both upstream and downstream. (5kb by default)",
	       "windowsize=i", "Number CpG sites in a window. (6bp by default)",
	       "windowstep=i", "Step of the sliding window, measured in number of CpG sites. (3bp by default)",
	       "maxwindow=i", "Maximum width of a window. (10kb by default)",
	       head = "Calculate correlation between methylation and gene expression",
	       foot = "",
	       script_name = "-e \"epik.cmd::epik()\" correlated_regions")

library(epik)
load_config(config)

if(is.null(subgroup)) {
	subgroup = unique(SAMPLE$subgroup)
}

SAMPLE = SAMPLE[SAMPLE$subgroup %in% subgroup, , drop = FALSE]
if(length(unique(SAMPLE$subgroup)) <= 1) {
	stop("Number of subgroups is less than 2.")
}

unique_subgroup = unique(SAMPLE$subgroup)

sample_id = rownames(SAMPLE)

if(is.null(EXPR)) stop("'EXPR' should be defined.")
if(is.null(TXDB)) stop("'TXDB' should be defined.")

if(!chr %in% CHROMOSOME) stop("'chr' should be included in 'CHROMOSOME'.")

cr = correlated_regions(sample_id, EXPR, TXDB, chr = chr, subgroup = SAMPLE$subgroup, col = COLOR$subgroup,
	window_size = windowsize, window_step = windowstep, max_width = maxwindow, genome = GENOME)

subgroup_str = paste(sort(unique_subgroup), collapse = "_")
dir.create(qq("@{PROJECT_DIR}/rds_cr/@{subgroup_str}"), showWarnings = FALSE)

saveRDS(cr, file = qq("@{PROJECT_DIR}/rds_cr/@{subgroup_str}/cr_@{chr}_@{subgroup_str}.rds"))
