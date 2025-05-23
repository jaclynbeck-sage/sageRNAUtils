# simple_cpm -------------------------------------------------------------------

test_that("simple_cpm gives CPM", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  counts_cpm <- simple_cpm(counts)

  expect_true(inherits(counts_cpm, "matrix"))
  expect_identical(colSums(counts_cpm), c(1e6, 1e6))
  expect_identical(counts_cpm, counts / 15 * 1e6)
})

test_that("simple_cpm retains matrix row and column names", {
  genes <- paste0("gene", 1:6)
  counts <- matrix(data = c(0:5, 5:0), nrow = 6,
                   dimnames = list(genes, c("s1", "s2")))
  counts_cpm <- simple_cpm(counts)

  expect_identical(colnames(counts_cpm), c("s1", "s2"))
  expect_identical(rownames(counts_cpm), genes)
})

test_that("simple_cpm works with data.frames", {
  genes <- paste0("gene", 1:6)
  counts <- data.frame(s1 = 0:5, s2 = 5:0, row.names = genes)
  counts_cpm <- simple_cpm(counts)

  expect_true(inherits(counts_cpm, "data.frame"))
  expect_identical(colSums(counts_cpm), c(s1 = 1e6, s2 = 1e6))
  expect_identical(counts_cpm, counts / 15 * 1e6)
})

test_that("simple_cpm works with sparse matrices", {
  genes <- paste0("gene", 1:16)
  counts <- matrix(c(0:5, rep(0, 10), 5:0, rep(0, 10)), nrow = 16,
                   dimnames = list(genes, c("s1", "s2")))
  counts <- Matrix::Matrix(counts, sparse = TRUE)

  counts_cpm <- simple_cpm(counts)

  expect_true(inherits(counts_cpm, "CsparseMatrix"))
  expect_identical(Matrix::colSums(counts_cpm), c(s1 = 1e6, s2 = 1e6))
  expect_identical(counts_cpm, counts / 15 * 1e6)
})

test_that("simple_cpm works with non-integers", {
  counts <- matrix(data = c(0:5, 5:0) + 0.3, nrow = 6)
  counts_cpm <- simple_cpm(counts)

  # Uses expect_equal instead of expect_identical because of floating point precision
  expect_equal(colSums(counts_cpm), c(1e6, 1e6))
  expect_identical(counts_cpm, counts / 16.8 * 1e6)
})

test_that("simple_cpm works with different library_size", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  library_size <- c(10, 20)

  counts_cpm <- simple_cpm(counts, library_size = library_size)
  expected <- sweep(counts, 2, library_size, "/") * 1e6

  expect_identical(counts_cpm, expected)
})

test_that("simple_cpm works with different size_factors", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  size_factors <- c(0.75, 1.5)

  counts_cpm <- simple_cpm(counts, size_factors = size_factors)
  expected <- sweep(counts, 2, colSums(counts) * size_factors, "/") * 1e6

  expect_identical(counts_cpm, expected)
})

test_that("simple_cpm works with different library_size and size_factors", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  library_size <- c(10, 20)
  size_factors <- c(0.75, 1.5)

  counts_cpm <- simple_cpm(counts, library_size = library_size, size_factors = size_factors)
  expected <- sweep(counts, 2, library_size * size_factors, "/") * 1e6

  expect_identical(counts_cpm, expected)
})

test_that("simple_cpm checks for correct lengths", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  library_size <- c(10, 20, 30)
  size_factors <- c(0.75, 1.5, 1, 1)

  expect_error(simple_cpm(counts, library_size = library_size),
               regexp = "The length of 'library_size'")

  expect_error(simple_cpm(counts, size_factors = size_factors),
               regexp = "'library_size' and 'size_factors'")

  expect_error(simple_cpm(counts, library_size[1:2], size_factors = size_factors),
               regexp = "'library_size' and 'size_factors'")
})


# simple_lognorm ---------------------------------------------------------------

test_that("simple_lognorm returns log2-normalized values", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  counts_log <- simple_lognorm(counts)

  counts_cpm <- counts / 15 * 1e6

  expect_true(inherits(counts_log, "matrix"))
  expect_identical(counts_log, log2(counts_cpm + 0.5))
})

test_that("simple_lognorm works with CPM instead of counts input", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  counts_cpm <- counts / 15 * 1e6
  counts_log <- simple_lognorm(counts_cpm)

  expect_identical(counts_log, log2(counts_cpm + 0.5))
})

test_that("simple_lognorm can use a different pseudocount", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  counts_log <- simple_lognorm(counts, pseudocount = 1)

  counts_cpm <- counts / 15 * 1e6

  expect_identical(counts_log, log2(counts_cpm + 1))
})

test_that("simple_lognorm works with data.frames", {
  genes <- paste0("gene", 1:6)
  counts <- data.frame(s1 = 0:5, s2 = 5:0, row.names = genes)
  counts_log <- simple_lognorm(counts)

  counts_cpm <- counts / 15 * 1e6

  expect_true(inherits(counts_log, "data.frame"))
  expect_identical(counts_log, log2(counts_cpm + 0.5))
})

test_that("simple_lognorm works with sparse matrices, pseudocount = 0.5", {
  genes <- paste0("gene", 1:16)
  counts <- matrix(c(0:5, rep(0, 10), 5:0, rep(0, 10)), nrow = 16,
                   dimnames = list(genes, c("s1", "s2")))
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  counts_cpm <- counts / 15 * 1e6

  counts_log <- simple_lognorm(counts)

  expect_true(inherits(counts_log, "dgeMatrix"))
  expect_identical(counts_log, log2(counts_cpm + 0.5))
})

test_that("simple_lognorm preserves sparsity with pseudocount = 1", {
  genes <- paste0("gene", 1:16)
  counts <- matrix(c(0:5, rep(0, 10), 5:0, rep(0, 10)), nrow = 16,
                   dimnames = list(genes, c("s1", "s2")))
  counts <- Matrix::Matrix(counts, sparse = TRUE)

  counts_log <- simple_lognorm(counts, pseudocount = 1)

  expected <- Matrix::Matrix(log2(counts / 15 * 1e6 + 1), sparse = TRUE)

  expect_true(inherits(counts_log, "CsparseMatrix"))
  expect_identical(counts_log, expected)
})

test_that("simple_lognorm retains matrix row and column names", {
  genes <- paste0("gene", 1:6)
  counts <- matrix(c(0:5, 5:0), nrow = 6, dimnames = list(genes, c("s1", "s2")))
  counts_log <- simple_lognorm(counts)

  expect_identical(colnames(counts_log), c("s1", "s2"))
  expect_identical(rownames(counts_log), genes)
})

test_that("simple_lognorm works with different library_size", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  library_size <- c(10, 20)

  counts_log <- simple_lognorm(counts, library_size = library_size)
  expected <- log2(sweep(counts, 2, library_size, "/") * 1e6 + 0.5)

  expect_identical(counts_log, expected)
})

test_that("simple_lognorm works with different size_factors", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  size_factors <- c(0.75, 1.5)

  counts_log <- simple_lognorm(counts, size_factors = size_factors)
  expected <- log2(sweep(counts, 2, colSums(counts) * size_factors, "/") * 1e6 + 0.5)

  expect_identical(counts_log, expected)
})

test_that("simple_lognorm works with different library_size and size_factors", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  library_size <- c(10, 20)
  size_factors <- c(0.75, 1.5)

  counts_log <- simple_lognorm(counts, library_size = library_size, size_factors = size_factors)
  expected <- log2(sweep(counts, 2, library_size * size_factors, "/") * 1e6 + 0.5)

  expect_identical(counts_log, expected)
})

test_that("simple_lognorm checks for correct lengths", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  library_size <- c(10, 20, 30)
  size_factors <- c(0.75, 1.5, 1, 1)

  expect_error(simple_lognorm(counts, library_size = library_size),
               regexp = "The length of 'library_size'")

  expect_error(simple_lognorm(counts, size_factors = size_factors),
               regexp = "'library_size' and 'size_factors'")

  expect_error(simple_lognorm(counts, library_size[1:2], size_factors = size_factors),
               regexp = "'library_size' and 'size_factors'")
})



# cpm_to_counts ----------------------------------------------------------------

test_that("cpm_to_counts returns counts", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  counts_cpm <- counts / 15 * 1e6
  new_counts <- cpm_to_counts(counts_cpm, library_size = colSums(counts))

  expect_true(inherits(new_counts, "matrix"))

  # Using expect_equal instead of expect_identical because new_counts will be
  # a matrix of doubles and counts is a matrix of integers
  expect_equal(new_counts, counts)
})

test_that("cpm_to_counts rounds non-integers", {
  counts <- matrix(data = c(0:5, 5:0) + 0.3, nrow = 6)
  counts_cpm <- counts / 16.8 * 1e6
  new_counts <- cpm_to_counts(counts_cpm, library_size = colSums(counts))

  # Using expect_equal instead of expect_identical because new_counts will be
  # a matrix of doubles and counts is a matrix of integers
  expect_equal(new_counts, round(counts))
})

test_that("cpm_to_counts works on data.frames", {
  genes <- paste0("gene", 1:6)
  counts <- data.frame(s1 = 0:5, s2 = 5:0, row.names = genes)
  counts_cpm <- counts / 15 * 1e6
  new_counts <- cpm_to_counts(counts_cpm, library_size = colSums(counts))

  expect_true(inherits(new_counts, "data.frame"))

  # Using expect_equal instead of expect_identical because new_counts will be
  # a matrix of doubles and counts is a matrix of integers
  expect_equal(new_counts, round(counts))
})

test_that("cpm_to_counts works on sparse matrices", {
  genes <- paste0("gene", 1:16)
  counts <- matrix(c(0:5, rep(0, 10), 5:0, rep(0, 10)), nrow = 16,
                   dimnames = list(genes, c("s1", "s2")))
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  counts_cpm <- counts / 15 * 1e6

  new_counts <- cpm_to_counts(counts_cpm, library_size = Matrix::colSums(counts))

  expect_true(inherits(new_counts, "CsparseMatrix"))
  expect_identical(new_counts, counts)
})

test_that("cpm_to_counts retains matrix row and column names", {
  genes <- paste0("gene", 1:6)
  counts <- matrix(data = c(0:5, 5:0), nrow = 6,
                   dimnames = list(genes, c("s1", "s2")))
  counts_cpm <- counts / 15 * 1e6
  new_counts <- cpm_to_counts(counts_cpm, library_size = colSums(counts))

  expect_identical(colnames(new_counts), c("s1", "s2"))
  expect_identical(rownames(new_counts), genes)
})

test_that("cpm_to_counts uses different library_size", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  library_size <- c(10, 20)

  counts_cpm <- sweep(counts, 2, library_size, "/") * 1e6
  new_counts <- cpm_to_counts(counts_cpm, library_size = library_size)

  # Using expect_equal instead of expect_identical because new_counts will be
  # a matrix of doubles and counts is a matrix of integers
  expect_equal(new_counts, counts)
})

test_that("cpm_to_counts uses different size_factors", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  size_factors <- c(0.75, 1.5)

  counts_cpm <- sweep(counts, 2, colSums(counts) * size_factors, "/") * 1e6
  new_counts <- cpm_to_counts(counts_cpm, library_size = colSums(counts),
                              size_factors = size_factors)

  # Using expect_equal instead of expect_identical because new_counts will be
  # a matrix of doubles and counts is a matrix of integers
  expect_equal(new_counts, counts)
})

test_that("cpm_to_counts uses different library_size and size_factors", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  library_size <- c(10, 20)
  size_factors <- c(0.75, 1.5)

  counts_cpm <- sweep(counts, 2, library_size * size_factors, "/") * 1e6
  new_counts <- cpm_to_counts(counts_cpm, library_size = library_size,
                              size_factors = size_factors)

  # Using expect_equal instead of expect_identical because new_counts will be
  # a matrix of doubles and counts is a matrix of integers
  expect_equal(new_counts, counts)
})

test_that("cpm_to_counts checks for correct lengths", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  library_size <- c(10, 20, 30)
  size_factors <- c(0.75, 1.5, 1, 1)

  counts_cpm <- sweep(counts, 2, colSums(counts), "/") * 1e6

  expect_error(cpm_to_counts(counts_cpm, library_size = library_size),
               regexp = "The length of 'library_size'")

  expect_error(cpm_to_counts(counts_cpm, library_size = colSums(counts),
                             size_factors = size_factors),
               regexp = "'library_size' and 'size_factors'")
})


# log_cpm_to_counts ----------------------------------------------------------------

test_that("log_cpm_to_counts returns counts", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  counts_log <- log2(counts / 15 * 1e6 + 0.5)
  new_counts <- log_cpm_to_counts(counts_log, library_size = colSums(counts))

  expect_true(inherits(new_counts, "matrix"))

  # Using expect_equal instead of expect_identical because new_counts will be
  # a matrix of doubles and counts is a matrix of integers
  expect_equal(new_counts, counts)
})

test_that("log_cpm_to_counts rounds non-integers", {
  counts <- matrix(data = c(0:5, 5:0) + 0.3, nrow = 6)
  counts_log <- log2(counts / 16.8 * 1e6 + 0.5)
  new_counts <- log_cpm_to_counts(counts_log, library_size = colSums(counts))

  # Using expect_equal instead of expect_identical because new_counts will be
  # a matrix of doubles and counts is a matrix of integers
  expect_equal(new_counts, round(counts))
})

test_that("log_cpm_to_counts works on data.frames", {
  genes <- paste0("gene", 1:6)
  counts <- data.frame(s1 = 0:5, s2 = 5:0, row.names = genes)
  counts_log <- log2(counts / 15 * 1e6 + 0.5)
  new_counts <- log_cpm_to_counts(counts_log, library_size = colSums(counts))

  expect_true(inherits(new_counts, "data.frame"))

  # Using expect_equal instead of expect_identical because new_counts will be
  # a matrix of doubles and counts is a matrix of integers
  expect_equal(new_counts, round(counts))
})

test_that("log_cpm_to_counts works on sparse matrices", {
  genes <- paste0("gene", 1:16)
  counts <- matrix(c(0:5, rep(0, 10), 5:0, rep(0, 10)), nrow = 16,
                   dimnames = list(genes, c("s1", "s2")))
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  counts_log <- counts
  counts_log@x <- log2(counts_log@x / 15 * 1e6 + 1)

  new_counts <- log_cpm_to_counts(counts_log,
                                  library_size = Matrix::colSums(counts),
                                  pseudocount = 1)

  expect_true(inherits(new_counts, "CsparseMatrix"))
  expect_identical(new_counts, counts)
})

test_that("log_cpm_to_counts retains matrix row and column names", {
  genes <- paste0("gene", 1:6)
  counts <- matrix(data = c(0:5, 5:0), nrow = 6,
                   dimnames = list(genes, c("s1", "s2")))
  counts_log <- log2(counts / 15 * 1e6 + 0.5)
  new_counts <- log_cpm_to_counts(counts_log, library_size = colSums(counts))

  expect_identical(colnames(new_counts), c("s1", "s2"))
  expect_identical(rownames(new_counts), genes)
})

test_that("log_cpm_to_counts uses different pseudocount", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  counts_log <- log2(counts / 15 * 1e6 + 2)
  new_counts <- log_cpm_to_counts(counts_log,
                                  library_size = colSums(counts),
                                  pseudocount = 2)
  expect_equal(new_counts, counts)
})

test_that("log_cpm_to_counts uses different library_size", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  library_size <- c(10, 20)

  counts_log <- log2(sweep(counts, 2, library_size, "/") * 1e6 + 0.5)
  new_counts <- log_cpm_to_counts(counts_log, library_size = library_size)

  # Using expect_equal instead of expect_identical because new_counts will be
  # a matrix of doubles and counts is a matrix of integers
  expect_equal(new_counts, counts)
})

test_that("log_cpm_to_counts uses different size_factors", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  size_factors <- c(0.75, 1.5)

  counts_log <- log2(sweep(counts, 2, colSums(counts) * size_factors, "/") * 1e6 + 0.5)
  new_counts <- log_cpm_to_counts(counts_log, library_size = colSums(counts),
                                  size_factors = size_factors)

  # Using expect_equal instead of expect_identical because new_counts will be
  # a matrix of doubles and counts is a matrix of integers
  expect_equal(new_counts, counts)
})

test_that("log_cpm_to_counts uses different library_size and size_factors", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  library_size <- c(10, 20)
  size_factors <- c(0.75, 1.5)

  counts_log <- log2(sweep(counts, 2, library_size * size_factors, "/") * 1e6 + 0.5)
  new_counts <- log_cpm_to_counts(counts_log, library_size = library_size,
                                  size_factors = size_factors)

  # Using expect_equal instead of expect_identical because new_counts will be
  # a matrix of doubles and counts is a matrix of integers
  expect_equal(new_counts, counts)
})

test_that("log_cpm_to_counts checks for correct lengths", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  library_size <- c(10, 20, 30)
  size_factors <- c(0.75, 1.5, 1, 1)

  counts_log <- log2(sweep(counts, 2, colSums(counts), "/") * 1e6 + 0.5)

  expect_error(log_cpm_to_counts(counts_log, library_size = library_size),
               regexp = "The length of 'library_size'")

  expect_error(log_cpm_to_counts(counts_log, library_size = colSums(counts),
                                 size_factors = size_factors),
               regexp = "'library_size' and 'size_factors'")
})


# edger_log_cpm_to_counts ------------------------------------------------------

test_that("edger_log_cpm_to_counts correctly returns counts", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  counts_log <- edgeR::cpm(counts, log = TRUE)

  new_counts <- edger_log_cpm_to_counts(counts_log, library_size = colSums(counts))

  expect_true(inherits(new_counts, "matrix"))
  expect_equal(new_counts, counts)
})

test_that("edger_log_cpm_to_counts rounds non-integers", {
  counts <- matrix(data = c(0:5, 5:0) + 0.3, nrow = 6)
  counts_log <- edgeR::cpm(counts, log = TRUE)

  new_counts <- edger_log_cpm_to_counts(counts_log, library_size = colSums(counts))

  # Using expect_equal instead of expect_identical because new_counts will be
  # a matrix of doubles and counts is a matrix of integers
  expect_equal(new_counts, round(counts))
})

test_that("edger_log_cpm_to_counts works on data.frames", {
  genes <- paste0("gene", 1:6)
  counts <- data.frame(s1 = 0:5, s2 = 5:0, row.names = genes)
  counts_log <- as.data.frame(edgeR::cpm(counts, log = TRUE))

  new_counts <- edger_log_cpm_to_counts(counts_log, library_size = colSums(counts))

  expect_true(inherits(new_counts, "data.frame"))

  # Using expect_equal instead of expect_identical because new_counts will be
  # a matrix of doubles and counts is a matrix of integers
  expect_equal(new_counts, round(counts))
})

test_that("edger_log_cpm_to_counts retains matrix row and column names", {
  genes <- paste0("gene", 1:6)
  counts <- matrix(data = c(0:5, 5:0), nrow = 6,
                   dimnames = list(genes, c("s1", "s2")))
  counts_log <- edgeR::cpm(counts, log = TRUE)

  new_counts <- edger_log_cpm_to_counts(counts_log, library_size = colSums(counts))

  expect_identical(colnames(new_counts), c("s1", "s2"))
  expect_identical(rownames(new_counts), genes)
})

test_that("edger_log_cpm_to_counts uses different prior_count", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  counts_log <- edgeR::cpm(counts, log = TRUE, prior.count = 0.1)
  new_counts <- edger_log_cpm_to_counts(counts_log,
                                        library_size = colSums(counts),
                                        prior_count = 0.1)
  expect_equal(new_counts, counts)
})

test_that("edger_log_cpm_to_counts uses different library_size", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  dge <- edgeR::DGEList(counts)
  dge$samples$lib.size <- c(10, 20)

  counts_log <- edgeR::cpm(counts, log = TRUE, lib.size = dge$samples$lib.size)
  new_counts <- edger_log_cpm_to_counts(counts_log,
                                        library_size = dge$samples$lib.size)
  expect_equal(new_counts, counts)
})

test_that("edger_log_cpm_to_counts uses different size_factors", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  dge <- edgeR::DGEList(counts)
  dge <- edgeR::normLibSizes(dge)

  counts_log <- edgeR::cpm(counts, log = TRUE)
  new_counts <- edger_log_cpm_to_counts(counts_log,
                                        library_size = dge$samples$lib.size,
                                        size_factors = dge$samples$norm.factors)
  expect_equal(new_counts, counts)
})

test_that("edger_log_cpm_to_counts uses different library_size and size_factors", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  dge <- edgeR::DGEList(counts)
  dge <- edgeR::normLibSizes(dge)
  dge$samples$lib.size <- c(10, 20)

  counts_log <- edgeR::cpm(counts, log = TRUE, lib.size = dge$samples$lib.size)
  new_counts <- edger_log_cpm_to_counts(counts_log,
                                        library_size = dge$samples$lib.size,
                                        size_factors = dge$samples$norm.factors)
  expect_equal(new_counts, counts)
})

test_that("edger_log_cpm_to_counts checks for correct lengths", {
  counts <- matrix(data = c(0:5, 5:0), nrow = 6)
  library_size <- c(10, 20, 30)
  size_factors <- c(0.75, 1.5, 1, 1)

  counts_log <- edgeR::cpm(counts, log = TRUE)

  expect_error(edger_log_cpm_to_counts(counts_log, library_size = library_size),
               regexp = "The length of 'library_size'")

  expect_error(log_cpm_to_counts(counts_log, library_size = colSums(counts),
                                 size_factors = size_factors),
               regexp = "'library_size' and 'size_factors'")
})


# TODO tests for cqn_to_counts
