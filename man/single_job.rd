\name{single_job}
\alias{single_job}
\title{
Submit a single job to the pipeline
}
\description{
Submit a single job to the pipeline
}
\usage{
single_job(cmd, name, walltime = "1:00:00", mem = "1G", nodes = 1,
    email = NULL, dependency = NULL, enforce = FALSE, tmpdir = tempdir(), submit_by = "qsub",
    execute = TRUE)
}
\arguments{

  \item{cmd}{bash commands}
  \item{name}{name of the job}
  \item{walltime}{running time for the job}
  \item{mem}{memory for the job}
  \item{nodes}{nodes for the job}
  \item{email}{email the messages should be sent to}
  \item{dependency}{dependency of upstream jobs. The value is a vector of job IDs.}
  \item{enforce}{enforce to rerun the job}
  \item{tmpdir}{temp directory where the bash scripts are put}
  \item{submit_by}{which job system}
  \item{execute}{whether submit the job}

}
\value{
Job id
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
