
suppressPackageStartupMessages(library(GetoptLong))
min1 = 0.5
min2 = 0.5
methdiff = 0
equalscale = TRUE
GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       "min1=f", "If there are multiple samples in the group 1, it is possible that 
	                  a segment has more than one states asigned to it. If the recurrency 
	                  of each state is relatively low, it means there is no one dominant 
	                  state for this segment and it should be removed. This argument controls 
	                  the minimal value for the recurrency of states in a given segment. 
	                  The value is a percent. (0.5 by default)",
	       "min2=f", "Same as `min1`, but for samples in group 2. (0.5 by default)",
	       "methdiff=f", "The segments for which the methylation difference between two groups 
	                     is less than this value are removed. (0 by default)",
	       "equalscale!", "When there are more than one comparisons, should all chord diagrams 
	                      have same scale? (TRUE by default)",
	       head = "Visualize chromatin states transitions in two subgroups",
	       foot = "To successfully run it, `chipseq_hooks$chromHMM()` must be set beforehand in 
	               the configuration file.

	               The orders of chromatin states in the plot can be 
	               set by `CHROMATIN_STATES_ORDER` (a vector of state names) and the color can 
	               be set by `CHROMATIN_STATES_COLOR` (a vector of colors for which names are 
	               state names) in the configuration file.",
	       script_name = "-e \"epik.cmd::epik()\" chromatin_states_transitions")

library(epik)
load_config(config, use_std_dir = TRUE)

sample_id = rownames(SAMPLE)
chromHMM_list = get_chromHMM_list(sample_id)

# not all samples in `sample_id` has chromHMM data
l = !sapply(chromHMM_list, is.null)
sample_id = sample_id[l]
chromHMM_list = chromHMM_list[l]
subgroup = SAMPLE$subgroup[l]

if(length(sample_id)) stop("no sample detected.")
if(length(unique(subgroup)) <= 1) stop("less than 2 subgroups.")

subgroup_level = unique(subgroup)
nc = length(subgroup_level)

# make pair-wise comparison
mat_list = list()
for(i in 1:(nc-1)) {
	for(j in (i+1):nc) {

		message(qq("making chromatin transition matrix between @{subgroup_level[i]} and @{subgroup_level[j]}"))

		gr_list_1 = chromHMM_list[subgroup == subgroup_level[i]]
		gr_list_2 = chromHMM_list[subgroup == subgroup_level[j]]
		min1 = min1 * length(gr_list_1)
		min2 = min2 * length(gr_list_2)

		mat = make_transition_matrix_from_chromHMM(gr_list_1, gr_list_2, window = window, min_1 = min1, min_2 = min_2, meth_diff = methdiff)

		if(exists("CHROMATIN_STATES_ORDER")) {
			mat = mat[CHROMATIN_STATES_ORDER, CHROMATIN_STATES_ORDER]
		}

		mat_list[[paste0(i, "_", j)]] = mat
	}
}

if(equalscale) {
	im = which.max(sapply(mat_list, sum))[1]
	max_mat1 = mat_list[[im]]

	im = which.max(sapply(mat_list, function(m) {
		cn = intersect(rownames(m), colnames(m))
		for(i in cn) {
			m[i, i] = 0
		}
		sum(m)
	}))[1]
	max_mat2 = mat_list[[im]]
} else {
	max_mat1 = NULL
	max_mat2 = NULL
}


if(!exists("CHROMATIN_STATES_COLOR")) {
	CHROMATIN_STATES_COLOR = NULL
}

pdf(qq("@{PROJECT_DIR}/image/chromatin_states_transitions.pdf", width = 8, height = 8))
for(i in 1:(nc-1)) {
	for(j in (i+1):nc) {

		if(nrow(mat) > 6) {
			legend_position = c("bottomleft", "bottomright")
		} else {
			legend_position = "bottomleft"
		}

		mat = mat_list[[paste0(i, "_", j)]]

		chromatin_states_transition_chord_diagram(mat, max_mat = max_mat1, group_names = subgroup_level[c(i, j)], 
			remove_unchanged_transition = FALSE, legend_position = legend_position, state_col = CHROMATIN_STATES_COLOR)
		chromatin_states_transition_chord_diagram(mat, max_mat = max_mat2, group_names = subgroup_level[c(i, j)], 
			remove_unchanged_transition = TRUE, legend_position = legend_position, state_col = CHROMATIN_STATES_COLOR)
	}
}
dev.off()
