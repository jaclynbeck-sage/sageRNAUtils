#' Simple CPM Operation on a Counts Matrix
#'
#' Transforms count data into counts per million (CPM) by dividing by library
#' size and multiplying by 1 million.
#'
#' This function is named `simple_cpm` to avoid confusion with [edgeR::cpm()],
#' which handles pseudocounts differently and has the ability to incorporate
#' TMM values into the library size.
#'
#' @param data a matrix or matrix-like object where rows are genes and columns
#'   are samples. All values in `data` should be >= 0.
#'
#' @return an object the same type and shape as data with all counts transformed
#'   to CPM
#' @export
#'
#' @seealso [simple_lognorm()]
#'
#' @examples
#' counts <- round(matrix(runif(1000, min = 0, max = 1e8), nrow = 100))
#' cpm_data <- simple_cpm(counts)
#' all(colSums(cpm_data) == 1e6)
simple_cpm <- function(data) {
  if (inherits(data, "Matrix")) {
    sweep(data, 2, Matrix::colSums(data), "/") * 1e6
  } else {
    sweep(data, 2, colSums(data), "/") * 1e6
  }
}


#' Simple log2-Normalization of a Counts Matrix
#'
#' Normalizes a counts matrix by transforming to CPM values and taking
#' log2(cpm + pseudocount)
#'
#' @param data a matrix or matrix-like object where rows are genes and columns
#'   are samples. Data should be un-normalized counts. All values in `data`
#'   should be >= 0.
#' @param pseudocount a pseudocount to add to every CPM value when taking the
#'   log2, to avoid taking log2(0). Defaults to 0.5
#'
#' @return an object the same type and shape as data with all counts transformed
#'   to log2(CPM)
#' @export
#'
#' @seealso [simple_cpm()]
#'
#' @examples
#' # Example using a pseudocount of 1, which makes all log2 values positive and
#' # count values of 0 remain 0.
#' counts <- round(matrix(runif(1000, min = 0, max = 1e8), nrow = 100))
#' log_data <- simple_lognorm(counts, pseudocount = 1)
simple_lognorm <- function(data, pseudocount = 0.5) {
  if (pseudocount == 1 & inherits(data, "sparseMatrix")) {
    # Preserves sparsity
    sparse_data <- simple_cpm(data)
    sparse_data@x <- log2(sparse_data@x + pseudocount)
    sparse_data
  } else {
    log2(simple_cpm(data) + pseudocount)
  }
}


#' Convert CPM Values Back to Counts
#'
#' @param library_size a numeric vector describing the library size of each
#'   sample in `data`. The length of `library_size` must match the number of
#'   columns in `data`. Typically library size is the sum of all counts in a
#'   sample, but this could also be modified with TMM or other size factors as
#'   well.
#' @inheritParams simple_cpm
#'
#' @return an object the same type and shape as `data` where CPM values have
#'   been converted back to integer counts. Any non-integer values will be
#'   rounded off.
#' @export
#'
#' @seealso [simple_cpm()]
#'
#' @examples
#' # Make a matrix of CPM data
#' counts <- round(matrix(runif(1000, min = 0, max = 1e8), nrow = 100))
#' cpm_data <- simple_cpm(counts)
#'
#' # Reverse the operation
#' counts2 <- cpm_to_counts(cpm_data, library_size = colSums(counts))
#' all(counts == counts2)
cpm_to_counts <- function(data, library_size) {
  round(sweep(data, 2, library_size, "*") / 1e6)
}
