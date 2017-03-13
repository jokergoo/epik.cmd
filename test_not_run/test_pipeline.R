
x = single_job("ls", name = "test", walltime="1:00:00", mem="1G", 
	nodes = 1, email = "z.gu@dkfz.de", tmpdir = ".", execute = FALSE)
x = single_job("ls -l", name = "test2", walltime="1:00:00", mem="1G", 
	nodes = 1, email = "z.gu@dkfz.de", tmpdir = ".", dependency = x, execute = FALSE)
