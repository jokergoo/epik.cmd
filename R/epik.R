
# == title
# Run pre-defiend scripts
#
# == details
# There are several R scripts which can be run directly. The path of all scripts can be obtained by 
# 
#    dir(system.file("pipeline", "scripts", package = "epik.cmd"), pattern = "\\.R$")
#
# You can either directly run these R scripts by:
#
#    Rscript full-path-of-cmd.R [options]
#
# or use the short command:
#
#    Rscript -e "epik.cmd::epik()" cmd [options]
#
# to get the list of commands:
#
#    Rscript -e "epikc.cmd::epik()"
#
# For each cmd, use ``Rscript -e "epik.cmd::epik()" cmd --help`` to get help.
#
# Basically all scripts need ``--config`` option which corresponds to a configuration R file
# that defines how to import data. A template for the configuration file is at:
#
#   system.file("pipeline", "pipeline_configure_template.R", package = "epik.cmd")
#
# Available commands (``cmd``) are:
#
# -chromatin_states_transitions visualize chromatin states transitions by chord diagram
# -correlated_enriched visualize enrichment of correlated regions on tss/cgi/tfbs/enhancers
# -correlated_regions find regions in which methylation is correlated to expression of the associated gene
# -correlated_regions_downstream visualize statistics of correlated regions, visualize genome-wide distribution of correlated regions by Hilbert curve
# -correlated_regions_filter only keep correlation regions with significant correlations. Subgroup specificity for each class is calculated if needed. 
# -correlated_regions_gviz visualize correlated regions and other associated information by Gviz package
# -correlated_regions_reduce reduce correlated regions
# -differential_methylation_in_cgi_and_shore visualize differentially methylated regions in cgi and shores
# -differential_methylation_in_genomic_features visualize differentially methylated regions in a set of genomic features.
# -general_methylation_distribution use heatmap to visualize methylation distribution
# -methylation_subtype_classification_in_cgi_and_shore classify subgroups by methylation in cgi and shores
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
epik = function() {

	all_cmds = dir(system.file("pipeline", "scripts", package = "epik.cmd"), pattern = "\\.R$")
	all_cmds = gsub("\\.R$", "", all_cmds)
	
	msg = 'Usage: Rscript -e "epik.cmd::epik()" cmd [options]\n\nAvailable cmd:\n\n'
	msg = paste0(msg, qq("  @{all_cmds}\n"), "\n", 
		"For each cmd, use 'Rscript -e \"epik.cmd::epik()\" cmd --help' to get help.\n\n")

	x = commandArgs(trailingOnly = TRUE)
	R_binary = file.path(R.home("bin"), "R")

	if(length(x) == 0) {
		cat(msg)
		if(interactive()) {
			warning("It should only be used under command line.")
			return(invisible(NULL))
		} else {
			q(save = "no")
		}
	}

	if(!x[1] %in% all_cmds) {
		cat("'", x, "' is not supported.\n\n", sep = "")
		cat(msg)
		if(interactive()) {
			warning("It should only be used under command line.")
			return(invisible(NULL))
		} else {
			q(save = "no")
		}
	}
	
	GetoptLong:::source(qq("@{system.file('pipeline', 'scripts', package = 'epik.cmd')}/@{x[1]}.R"), argv = paste(x[-1], collapse=' '))
}
