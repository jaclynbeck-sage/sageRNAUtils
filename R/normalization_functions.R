#' Simple CPM Operation on a Counts Matrix
#'
#' Transforms count data into counts per million (CPM) by dividing by library
#' size and multiplying by 1 million.
#'
#' This function is named `simple_cpm` to avoid name collision with
#' [edgeR::cpm()], which has more advanced configuration options.
#'
#' @param data a matrix or matrix-like object where rows are genes and columns
#'   are samples. All values in `data` should be >= 0.
#' @param library_size (optional) a numeric vector of library sizes for each
#'   column in `data`. If `NULL`, `library_size` will be `colSums(data)`.
#'   Defaults to `NULL`.
#' @param size_factors (optional) a numeric vector of additional size factors to
#'   include in the library size, for example TMM factors. Defaults to `NULL`.
#'
#' @return an object the same type and shape as data with all counts transformed
#'   to CPM
#' @export
#'
#' @seealso [simple_lognorm()]
#'
#' @examples
#' counts <- round(matrix(runif(1000, min = 0, max = 100), ncol = 10))
#' cpm_data <- simple_cpm(counts)
#' all(colSums(cpm_data) == 1e6)
#'
#' # With size factors
#' tmm <- edgeR::normLibSizes(counts)
#' cpm_data <- simple_cpm(counts, size_factors = tmm)
simple_cpm <- function(data, library_size = NULL, size_factors = NULL) {
  if (is.null(library_size)) {
    # Use diferent functions for sparse vs dense matrices
    if (inherits(data, "Matrix")) {
      library_size <- Matrix::colSums(data)
    } else {
      library_size <- colSums(data)
    }
  }

  if (length(library_size) != ncol(data)) {
    stop("The length of 'library_size' doesn't match the number of columns in 'data'.")
  }

  if (!is.null(size_factors)) {
    if (length(size_factors) != length(library_size)) {
      stop("'library_size' and 'size_factors' are not the same length.")
    }
    library_size <- library_size * size_factors
  }

  sweep(data, 2, library_size, "/") * 1e6
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
#' @inheritParams simple_cpm
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
#' counts <- round(matrix(runif(1000, min = 0, max = 100), ncol = 10))
#' log_data <- simple_lognorm(counts, pseudocount = 1)
simple_lognorm <- function(data, library_size = NULL, size_factors = NULL, pseudocount = 0.5) {
  if (pseudocount == 1 & inherits(data, "sparseMatrix")) {
    # Preserves sparsity
    sparse_data <- simple_cpm(data, library_size = library_size,
                              size_factors = size_factors)
    sparse_data@x <- log2(sparse_data@x + pseudocount)
    return(sparse_data)
  } else {
    return(log2(simple_cpm(data, library_size = library_size, size_factors = size_factors)
                + pseudocount))
  }
}


#' Convert CPM Values Back to Counts
#'
#' Multiplies CPM values by library size / 1e6 to obtain the original counts
#' values.
#'
#' @param library_size a numeric vector describing the library size of each
#'   sample in `data`. The length of `library_size` must match the number of
#'   columns in `data`.
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
#' counts <- round(matrix(runif(1000, min = 0, max = 100), ncol = 10))
#' cpm_data <- simple_cpm(counts)
#'
#' # Reverse the operation
#' counts2 <- cpm_to_counts(cpm_data, library_size = colSums(counts))
#' all(counts == counts2)
cpm_to_counts <- function(data, library_size, size_factors = NULL) {
  if (length(library_size) != ncol(data)) {
    stop("The length of 'library_size' doesn't match the number of columns in 'data'.")
  }

  if (!is.null(size_factors)) {
    if (length(size_factors) != length(library_size)) {
      stop("'library_size' and 'size_factors' are not the same length.")
    }
    library_size <- library_size * size_factors
  }

  round(sweep(data, 2, library_size, "*") / 1e6)
}


#' Convert log2-CPM Normalized Values Back to Counts
#'
#' Reverses the operation performed by [simple_lognorm()] to convert log2-CPM
#' normalized values back to integer counts.
#'
#' @inheritParams cpm_to_counts
#' @param data a matrix or matrix-like object where rows are genes and columns
#'   are samples. All values in `data` should on the log2 scale and normalized
#'   as in [simple_lognorm()].
#' @param pseudocount The pseudocount that was used in the original call to
#'   [simple_lognorm()]. Defaults to 0.5.
#'
#' @return an object the same type and shape as `data` where log2-CPM values
#'   have been converted back to integer counts. Any non-integer values will be
#'   rounded off.
#' @export
#' @seealso [simple_lognorm()], [cpm_to_counts()]
#'
#' @examples
#' # Make a matrix of normalized data
#' counts <- round(matrix(runif(1000, min = 0, max = 100), ncol = 10))
#' log_data <- simple_lognorm(counts, pseudocount = 1)
#'
#' # Reverse the operation
#' counts2 <- log_cpm_to_counts(log_data,
#'                              library_size = colSums(counts),
#'                              pseudocount = 1)
#' all(counts == counts2)
log_cpm_to_counts <- function(data, library_size, size_factors = NULL, pseudocount = 0.5) {
  if (pseudocount == 1 & inherits(data, "sparseMatrix")) {
    sparse_data <- data
    sparse_data@x <- 2^sparse_data@x - pseudocount
    cpm_to_counts(sparse_data, library_size = library_size, size_factors = size_factors)
  } else {
    cpm_to_counts(2^data - pseudocount, library_size = library_size, size_factors = size_factors)
  }
}


#' Convert edgeR-style Log-CPM Values to Counts
#'
#' Reverses the operation performed by [edgeR::cpm()] with argument `log =
#' TRUE`, which is slightly different than the simple `log2(CPM + pseudocount)`
#' operation.
#'
#' @details
#' When edgeR calculates log-CPM, it does the following:
#'
#' `log_cpm = [log(counts + prior) - log(adj_lib_size) + log(1e6)] / log(2)`
#'
#' with:
#'
#' `prior = pseudocount * lib_size / avg(lib_size)`
#'
#' `adj_lib_size = lib_size + 2 * prior`
#'
#' where `lib_size` is either `colSums(counts)`, `colSums(counts) *
#' tmm_factors`, or some other pre-calculated offset value.
#'
#' To reverse this operation, we calculate `prior` and `adj_lib_size` for all
#' samples as above and undo the normalization:
#'
#' `counts = e^[log_cpm * log(2) + log(adj_lib_size) - log(1e6)] - prior`
#'
#' @inheritParams cpm_to_counts
#' @param data a matrix or matrix-like object where rows are genes and columns
#'   are samples. Data should be on the log2-scale, exactly as generated by a
#'   call to [edgeR::cpm()] with argument `log = TRUE`.
#' @param prior_count (optional) the `prior.count` argument that was used in the
#' original call to [edgeR::cpm()]. Defaults to 2, which is the same default in
#' the `edgeR` function.
#'
#' @return an object the same type and shape as `data` where edgeR-style log-CPM
#'   values have been converted back to integer counts. Any non-integer values
#'   will be rounded off.
#' @export
#'
#' @seealso [edgeR::cpm()]
#'
#' @examples
#' counts <- round(matrix(runif(1000, min = 0, max = 100), ncol = 10))
#' log_cpm <- edgeR::cpm(counts, log = TRUE)
#' new_counts <- edger_log_cpm_to_counts(log_cpm, library_size = colSums(counts))
#' all(new_counts == counts)
#'
#' # With size factors
#' dge <- edgeR::DGEList(counts)
#' dge <- edgeR::normLibSizes(dge)
#' log_cpm <- edgeR::cpm(dge, log = TRUE)
#' new_counts <- edger_log_cpm_to_counts(log_cpm,
#'                                       library_size = dge$samples$lib.size,
#'                                       size_factors = dge$samples$norm.factors)
#' all(new_counts == counts)
edger_log_cpm_to_counts <- function(data, library_size, size_factors = NULL, prior_count = 2) {
  if (length(library_size) != ncol(data)) {
    stop("The length of 'library_size' doesn't match the number of columns in 'data'.")
  }

  if (!is.null(size_factors)) {
    if (length(size_factors) != length(library_size)) {
      stop("'library_size' and 'size_factors' are not the same length.")
    }
    library_size <- library_size * size_factors
  }

  prior <- prior_count * library_size / mean(library_size)
  adj_lib_size <- library_size + 2 * prior

  counts <- exp(sweep(data * log(2), 2, log(adj_lib_size) - log(1e6), "+"))
  counts <- round(sweep(counts, 2, prior, "-"))

  return(counts)
}


#' Convert CQN-Normalized Values to Counts
#'
#' Takes normalized values generated by [cqn::cqn()], which are on the
#' log2-scale, and converts them back to "corrected" counts.
#'
#' @details
#' [cqn::cqn()] normalizes the data prior to GC content adjustment as follows:
#'
#' `log_data = log2(counts + 1) - log2(lib_size / 10^6)`
#'
#' where `lib_size` is the library size of each sample, by default `colSums(counts)`
#' in the `cqn` function.
#'
#' To reverse this operation on the output of `cqn` to get counts we do:
#'
#' `counts = 2^(cqn_data + log2(lib_size / 10^6)) - 1`
#'
#' @param data a matrix or matrix-like object where rows are genes and columns
#'   are samples. Data should be on the log2-scale, exactly as generated by a
#'   call to [cqn::cqn()].
#' @param library_size a numeric vector describing the library size of each
#'   sample in `data`, which should match what was input to the original call to
#'   [cqn::cqn()]. Typically, `library_size` will be equal to the `sizeFactors`
#'   argument if it was supplied to `cqn`, or to the sum of counts for each
#'   sample if `sizeFactors` was left as the default (`NULL`). The length of
#'   `library_size` must match the number of columns in `data`.
#'
#' @return an object the same type and shape as `data` where cqn-normalized
#'   values have been converted back to integer counts. Any non-integer values
#'   will be rounded off, and values less than 0 will be set to 0.
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate fake counts, GC content, and lengths for each gene
#' counts <- round(matrix(runif(1000, min = 0, max = 100), ncol = 10))
#' gc_content <- abs(rnorm(100, mean = 0.5, sd = 0.25))
#' gene_length <- abs(round(rnorm(100, mean = 200, sd = 100)))
#'
#' # Run CQN
#' cqn_data <- cqn::cqn(counts, x = gc_content, lengths = gene_length)
#' cqn_log <- cqn_data$y + cqn_data$offset
#'
#' # Get corrected counts
#' corrected_counts <- cqn_to_counts(cqn_log, library_size = colSums(counts))
#' }
cqn_to_counts <- function(data, library_size) {
  # TODO accept "cqn" objects too
  counts <- sweep(data, 2, log2(library_size / 1e6), "+")
  counts <- round(2^counts - 1)
  counts[counts < 0] <- 0

  return(counts)
}
