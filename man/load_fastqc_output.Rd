% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_functions.R
\name{load_fastqc_output}
\alias{load_fastqc_output}
\title{Format FastQC Output}
\usage{
load_fastqc_output(files, sample_names, ...)
}
\arguments{
\item{files}{a character vector of file names (including absolute or relative
file paths) to parse. Files should all end with \verb{_fastqc.zip}, as output by
FastQC.}

\item{sample_names}{a character vector the same length as \code{files} with the
values to put in the \code{sample} column of each combined data frame. Sample
names should be unique to account for the fact that there are potentially
two read files per specimen.}

\item{...}{additional arguments to pass to \code{\link[fastqcr:qc_read_collection]{fastqcr::qc_read_collection()}},
for example different values for \code{modules} or \code{verbose}.}
}
\value{
a named list of data.frames containing stats for all samples listed
in \code{files}. There will be one data frame per FastQC "module". See the
\code{modules} argument of \code{\link[fastqcr:qc_read_collection]{fastqcr::qc_read_collection()}} for the list of data
available.
}
\description{
Loads a folder of fastqc.zip files output by FastQC, combines data for all
samples, and cleans up the final data. This is a wrapper around
\code{\link[fastqcr:qc_read_collection]{fastqcr::qc_read_collection()}}.
}
\details{
This function uses \code{\link[fastqcr:qc_read_collection]{fastqcr::qc_read_collection()}} to parse all of the zip
files and combine them into data.frames. Some of the data is then cleaned up
to fix column typing, naming, and/or missing Percentage columns. The final
output will be a list of data.frames corresponding to the modules requested
from \code{\link[fastqcr:qc_read_collection]{fastqcr::qc_read_collection()}} (by default, all modules).

\strong{Note:} This function can take anywhere from several seconds to over an hour, depending
on how many files there are.
}
\examples{
\dontrun{
# Fastqc output for 2 samples (example1 and example2) that had both R1 and R2
# fastq files.
files <- c(
  "data/example1_1_fastqc.zip", "data/example1_2_fastqc.zip",
  "data/example2_1_fastqc.zip", "data/example2_2_fastqc.zip"
)
sample_names <- c("example1_R1", "example1_R2", "example2_R1", "example2_R2")

fastqc_data <- load_fastqc_output(files, sample_names)

# Only get certain modules
fastqc_data_short <- load_fastqc_output(
  files, sample_names,
  modules = c("Summary", "Basic Statistics", "Overrepresented sequences")
)
}
}
\seealso{
\code{\link[fastqcr:qc_read_collection]{fastqcr::qc_read_collection()}}, \code{\link[fastqcr:qc_read]{fastqcr::qc_read()}}
}
