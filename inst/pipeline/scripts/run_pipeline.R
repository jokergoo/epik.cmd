
suppressPackageStartupMessages(library(GetoptLong))

enforce = FALSE
execute = TRUE
prefix = ""
email = NULL
submitby = "qsub"
GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       "enforce!", "Enforce to run all jobs (by default FALSE)",
	       "execute!", "Whether submit jobs to clusters (by default FALSE)",
	       "email=s", "Email address the messages should be sent to",
	       "submitby=s", "which job system. (by default qsub)",
	       "prefix=s", "job prefix",
	       head = "",
	       foot = "",
	       script_name = "-e \"epik.cmd::epik()\" run_pipeline")


library(epik.cmd)
epik_pipeline(config, prefix = prefix, email = email, enforce = enforce, submit_by = submitby)
