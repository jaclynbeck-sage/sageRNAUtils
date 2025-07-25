% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization_functions.R
\name{log2_cpm_to_counts}
\alias{log2_cpm_to_counts}
\title{Convert log2-CPM Normalized Values Back to Counts}
\usage{
log2_cpm_to_counts(data, library_size, size_factors = NULL, pseudocount = 0.5)
}
\arguments{
\item{data}{a matrix or matrix-like object where rows are genes and columns
are samples. All values in \code{data} should on the log2 scale and normalized
as in \code{\link[=simple_log2norm]{simple_log2norm()}}.}

\item{library_size}{a numeric vector describing the library size of each
sample in \code{data}. The length of \code{library_size} must match the number of
columns in \code{data}.}

\item{size_factors}{(optional) a numeric vector of additional size factors to
include in the library size, for example TMM factors. Defaults to \code{NULL}.}

\item{pseudocount}{The pseudocount that was used in the original call to
\code{\link[=simple_log2norm]{simple_log2norm()}}. Defaults to 0.5.}
}
\value{
an object the same type and shape as \code{data} where log2-CPM values
have been converted back to integer counts. Any non-integer values will be
rounded off.
}
\description{
Reverses the operation performed by \code{\link[=simple_log2norm]{simple_log2norm()}} to convert log2-CPM
normalized values back to integer counts.
}
\examples{
# Make a matrix of normalized data
counts <- round(matrix(runif(1000, min = 0, max = 100), ncol = 10))
log_data <- simple_log2norm(counts, pseudocount = 1)

# Reverse the operation
counts2 <- log2_cpm_to_counts(log_data,
                              library_size = colSums(counts),
                              pseudocount = 1)
all(counts == counts2)
}
\seealso{
\code{\link[=simple_log2norm]{simple_log2norm()}}, \code{\link[=cpm_to_counts]{cpm_to_counts()}}
}
