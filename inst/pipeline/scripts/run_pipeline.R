
suppressPackageStartupMessages(library(GetoptLong))

GetoptLong("config=s", "A configuration R script. Check the help page of `load_config()` 
	                   function to find out how to properly set one.",
	       "enforce", "Enforce to run all jobs (by default FALSE)",
	       "execute", "Whether submit jobs to clusters (by default FALSE)",
	       "email=s", "Email address the messages should be sent to",
	       "submitby=s", "which job system. (by default qsub)"
	       head = "",
	       foot = "",
	       script_name = "-e \"epik.cmd::epik()\" run_pipeline")


library(epik.cmd)
epik_pipeline()
