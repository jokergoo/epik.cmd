
suppressPackageStartupMessages(library(GetoptLong))
min1 = 0.5
min2 = 0.5
methdiff = 0
GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       "subgroup1=s", "name of subgroup1",
	       "subgroup2=s", "name of subgroup2",
	       "min1=f", "If there are multiple samples in the group 1, it is possible that 
	                  a segment has more than one states asigned to it. If the recurrency 
	                  of each state is relatively low, it means there is no one dominant 
	                  state for this segment and it should be removed. This argument controls 
	                  the minimal value for the recurrency of states in a given segment. 
	                  The value is a percent. (0.5 by default)",
	       "min2=f", "Same as `min1`, but for samples in group 2. (0.5 by default)",
	       "methdiff=f", "The segments for which the methylation difference between two groups 
	                     is less than this value are removed. (0 by default)",
	       head = "Visualize chromatin states transitions in two subgroups",
	       foot = "To successfully run it, `chipseq_hooks$chromHMM()` must be set beforehand in 
	               the configuration file.

	               The orders of chromatin states in the plot can be 
	               set by `CHROMATIN_STATES_ORDER` (a vector of state names) and the color can 
	               be set by `CHROMATIN_STATES_COLOR` (a vector of colors for which names are 
	               state names) in the configuration file.",
	       script_name = "-e \"epik.cmd::epik()\" chromatin_states_transitions")

library(epik)
load_epik_config(config)

if(!(subgroup1 %in% SAMPLE$subgroup)) {
	stop("Cannot find ", subgroup1, " in `SAMPLE$subgroup. ")
}
if(!(subgroup2 %in% SAMPLE$subgroup)) {
	stop("Cannot find ", subgroup2, " in `SAMPLE$subgroup")
}

SAMPLE = SAMPLE[SAMPLE$subgroup %in% c(subgroup1, subgroup2), , drop = FALSE]

sample_id = rownames(SAMPLE)
chromHMM_list = get_chromHMM_list(sample_id)

# not all samples in `sample_id` has chromHMM data
l = !sapply(chromHMM_list, is.null)
chromHMM_list = chromHMM_list[l]
SAMPLE = SAMPLE[names(chromHMM_list), , drop = FALSE]
subgroup = SAMPLE$subgroup

if(length(chromHMM_list) <= 1) stop("no sample detected.")
if(length(unique(subgroup)) <= 1) stop("less than 2 subgroups.")

message(qq("making chromatin transition matrix between @{subgroup1} and @{subgroup2}"))

gr_list_1 = chromHMM_list[subgroup == subgroup1]
gr_list_2 = chromHMM_list[subgroup == subgroup2]

mat = make_transition_matrix_from_chromHMM(gr_list_1, gr_list_2, min_1 = min1, min_2 = min2, meth_diff = methdiff)

if(exists("CHROMATIN_STATES_ORDER")) {
	mat = mat[CHROMATIN_STATES_ORDER, CHROMATIN_STATES_ORDER]
}

if(!exists("CHROMATIN_STATES_COLOR")) {
	CHROMATIN_STATES_COLOR = NULL
}

pdf(qq("@{PROJECT_DIR}/image/chromatin_states_transitions_@{subgroup1}_vs_@{subgroup2}.pdf"), width = 8, height = 8)

if(nrow(mat) > 6) {
	legend_position = c("bottomleft", "bottomright")
} else {
	legend_position = "bottomleft"
}

chromatin_states_transition_chord_diagram(mat, group_names = c(subgroup1, subgroup2), 
	remove_unchanged_transition = FALSE, legend_position = legend_position, state_col = CHROMATIN_STATES_COLOR)
chromatin_states_transition_chord_diagram(mat, group_names = c(subgroup1, subgroup2), 
	remove_unchanged_transition = TRUE, legend_position = legend_position, state_col = CHROMATIN_STATES_COLOR)

dev.off()
