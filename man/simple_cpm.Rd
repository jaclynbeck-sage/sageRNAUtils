% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization_functions.R
\name{simple_cpm}
\alias{simple_cpm}
\title{Simple CPM Operation on a Counts Matrix}
\usage{
simple_cpm(data, library_size = NULL, size_factors = NULL)
}
\arguments{
\item{data}{a matrix or matrix-like object where rows are genes and columns
are samples. All values in \code{data} should be >= 0.}

\item{library_size}{(optional) a numeric vector of library sizes for each
column in \code{data}. If \code{NULL}, \code{library_size} will be \code{colSums(data)}.
Defaults to \code{NULL}.}

\item{size_factors}{(optional) a numeric vector of additional size factors to
include in the library size, for example TMM factors. Defaults to \code{NULL}.}
}
\value{
an object the same type and shape as data with all counts transformed
to CPM
}
\description{
Transforms count data into counts per million (CPM) by dividing by library
size and multiplying by 1 million.
}
\details{
This function is named \code{simple_cpm} to avoid name collision with
\code{\link[edgeR:cpm]{edgeR::cpm()}}, which has more advanced configuration options.
}
\examples{
counts <- round(matrix(runif(1000, min = 0, max = 100), ncol = 10))
cpm_data <- simple_cpm(counts)
all(colSums(cpm_data) == 1e6)

# With size factors
tmm <- edgeR::normLibSizes(counts)
cpm_data <- simple_cpm(counts, size_factors = tmm)
}
\seealso{
\code{\link[=simple_log2norm]{simple_log2norm()}}
}
