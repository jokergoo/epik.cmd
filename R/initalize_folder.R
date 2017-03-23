
# == title
# Initialize project directories
#
# == param
# -project_dir path of the project directory
#
# == value
# No value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
initialize_project_directory = function(project_dir) {
	dir.create(project_dir, showWarnings = FALSE)
	dir.create(paste0(project_dir, "/image"), showWarnings = FALSE)
	dir.create(paste0(project_dir, "/gviz"), showWarnings = FALSE)
	dir.create(paste0(project_dir, "/rds"), showWarnings = FALSE)
	dir.create(paste0(project_dir, "/rds_cr"), showWarnings = FALSE)
	dir.create(paste0(project_dir, "/temp"), showWarnings = FALSE)
	dir.create(paste0(project_dir, "/enriched_heatmap"), showWarnings = FALSE)

	message("Following folders created:")
	message(qq("  @{project_dir}/image/ for general plots"))
	message(qq("  @{project_dir}/gviz/ for Gviz plots"))
	message(qq("  @{project_dir}/rds/ for rds files"))
	message(qq("  @{project_dir}/rds_cr/ for CR-related rds files"))
	message(qq("  @{project_dir}/temp/ for temporary files"))
	message(qq("  @{project_dir}/enriched_heatmap/ for enriched heatmaps"))
}

