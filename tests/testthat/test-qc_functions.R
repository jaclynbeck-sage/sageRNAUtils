# is_outlier_IQR ---------------------------------------------------------------

test_that("is_outlier_IQR finds no outliers", {
  data <- c(rep(10:20, 5))
  outliers <- is_outlier_IQR(data)

  expect_true(sum(outliers) == 0)
})

test_that("is_outlier_IQR finds upper outliers only", {
  data <- c(-5, -4, 10:20, 99, 100)
  outliers <- is_outlier_IQR(data, tail = "upper")

  expect_true(sum(outliers) == 2)
  expect_identical(data[outliers], c(99, 100))
})

test_that("is_outlier_IQR finds lower outliers only", {
  data <- c(-5, -4, 10:20, 99, 100)
  outliers <- is_outlier_IQR(data, tail = "lower")

  expect_true(sum(outliers) == 2)
  expect_identical(data[outliers], c(-5, -4))
})

test_that("is_outlier_IQR finds all outliers", {
  data <- c(-5, -4, 10:20, 99, 100)
  outliers <- is_outlier_IQR(data)

  expect_true(sum(outliers) == 4)
  expect_identical(data[outliers], c(-5, -4, 99, 100))
})

test_that("is_outlier_IQR uses IQR_mult argument", {
  data <- c(-5, -4, 10:20, 40, 41, 99, 100)

  outliers1.5 <- is_outlier_IQR(data, IQR_mult = 1.5)
  outliers3 <- is_outlier_IQR(data, IQR_mult = 3)

  expect_true(sum(outliers1.5) == 6)
  expect_true(sum(outliers3) == 2)
  expect_identical(data[outliers1.5], c(-5, -4, 40, 41, 99, 100))
  expect_identical(data[outliers3], c(99, 100))
})

test_that("is_outlier_IQR handles NA values", {
  data <- c(-5, -4, NA, 10:20, 99, NA, 100)
  outliers <- is_outlier_IQR(data)

  expect_true(sum(outliers) == 6)
  expect_identical(data[outliers], c(-5, -4, NA, 99, NA, 100))
})

test_that("is_outlier_IQR works on matrices", {
  data <- c(-5, -4, 10:20, 99, 100)
  data <- cbind(data, data + 1, data - 1)
  outliers <- is_outlier_IQR(data)

  expect_true(sum(outliers) == 12)
  expect_true(nrow(outliers) == nrow(data))
  expect_true(ncol(outliers) == ncol(data))
  expect_identical(data[outliers], c(-5, -4, 99, 100, -4, -3, 100, 101, -6, -5, 98, 99))
})


# is_outlier_SD ---------------------------------------------------------------

test_that("is_outlier_SD finds no outliers", {
  data <- c(rep(10:20, 5))
  outliers <- is_outlier_SD(data)

  expect_true(sum(outliers) == 0)
})

test_that("is_outlier_SD finds upper outliers only", {
  data <- c(-51, -50, rep(10:20, 10), 99, 100)
  outliers <- is_outlier_SD(data, tail = "upper")

  expect_true(sum(outliers) == 2)
  expect_identical(data[outliers], c(99, 100))
})

test_that("is_outlier_SD finds lower outliers only", {
  data <- c(-51, -50, rep(10:20, 10), 99, 100)
  outliers <- is_outlier_SD(data, tail = "lower")

  expect_true(sum(outliers) == 2)
  expect_identical(data[outliers], c(-51, -50))
})

test_that("is_outlier_SD finds all outliers", {
  data <- c(-51, -50, rep(10:20, 10), 99, 100)
  outliers <- is_outlier_SD(data)

  expect_true(sum(outliers) == 4)
  expect_identical(data[outliers], c(-51, -50, 99, 100))
})

test_that("is_outlier_SD uses n_sds argument", {
  data <- c(-51, -50, rep(10:20, 10), 99, 100)

  outliers4 <- is_outlier_SD(data, n_sds = 4)
  outliers5 <- is_outlier_SD(data, n_sds = 5)

  expect_true(sum(outliers4) == 4)
  expect_true(sum(outliers5) == 2)
  expect_identical(data[outliers4], c(-51, -50, 99, 100))
  expect_identical(data[outliers5], c(99, 100))
})

test_that("is_outlier_SD handles NA values", {
  data <- c(-51, -50, NA, rep(10:20, 10), 99, NA, 100)
  outliers <- is_outlier_SD(data)

  expect_true(sum(outliers) == 6)
  expect_identical(data[outliers], c(-51, -50, NA, 99, NA, 100))
})

test_that("is_outlier_SD works on matrices", {
  data <- c(-51, -50, rep(10:20, 10), 99, 100)
  data <- cbind(data, data + 1, data - 1)
  outliers <- is_outlier_SD(data)

  expect_true(sum(outliers) == 12)
  expect_true(nrow(outliers) == nrow(data))
  expect_true(ncol(outliers) == ncol(data))
  expect_identical(data[outliers], c(-51, -50, 99, 100, -50, -49, 100, 101, -52, -51, 98, 99))
})


# find_pca_outliers ------------------------------------------------------------

test_that("find_pca_outliers correctly calculates PCA and thresholds", {
  data <- matrix_for_pca_no_outliers()
  results <- find_pca_outliers(data)

  expected_pca <- prcomp(t(data), center = TRUE, scale. = TRUE)$x
  expected_thresh1 <- 4 * sd(expected_pca[, 1])
  expected_thresh2 <- 4 * sd(expected_pca[, 2])

  expect_identical(names(results), c("pca_df", "pc1_threshold", "pc2_threshold", "outliers"))
  expect_identical(results$pca_df[, colnames(expected_pca)], as.data.frame(expected_pca))
  expect_identical(results$pc1_threshold, expected_thresh1)
  expect_identical(results$pc2_threshold, expected_thresh2)
})

test_that("find_pca_outliers finds no outliers with no metadata", {
  data <- matrix_for_pca_no_outliers()

  results <- find_pca_outliers(data)

  expect_identical(names(results), c("pca_df", "pc1_threshold", "pc2_threshold", "outliers"))
  expect_length(results$outliers, 0)
  expect_identical(results$pca_df$sample, colnames(data))
})

test_that("find_pca_outliers finds no outliers with metadata with default sample column name", {
  data <- matrix_for_pca_no_outliers()

  metadata <- data.frame(specimenID = colnames(data),
                         extra_column = 1:ncol(data))

  results <- find_pca_outliers(data, metadata = metadata)

  expect_identical(names(results), c("pca_df", "pc1_threshold", "pc2_threshold", "outliers"))
  expect_length(results$outliers, 0)
  expect_true("specimenID" %in% colnames(results$pca_df))
  expect_identical(sort(results$pca_df$specimenID), sort(colnames(data)))
})

test_that("find_pca_outliers finds outliers with no metadata", {
  data <- matrix_for_pca_with_outliers()

  results <- find_pca_outliers(data)

  expect_identical(names(results), c("pca_df", "pc1_threshold", "pc2_threshold", "outliers"))
  expect_length(results$outliers, 2)
  expect_identical(results$pca_df$sample, colnames(data))
  expect_identical(results$outliers, c("sample_out1", "sample_out2"))
})

test_that("find_pca_outliers finds outliers with metadata with default sample column name", {
  data <- matrix_for_pca_with_outliers()

  metadata <- data.frame(specimenID = colnames(data),
                         extra_column = 1:ncol(data))

  results <- find_pca_outliers(data, metadata = metadata)

  expect_identical(names(results), c("pca_df", "pc1_threshold", "pc2_threshold", "outliers"))
  expect_length(results$outliers, 2)
  expect_true("specimenID" %in% colnames(results$pca_df))
  expect_true("extra_column" %in% colnames(results$pca_df))

  # pca_df might come out with the specimenIDs in a different sorted order when
  # metadata is used, so both the colnames and the extra column variable need to
  # be sorted that way too
  expected_extra_column <- dplyr::arrange(metadata, specimenID) |>
    dplyr::pull(extra_column)
  received_extra_column <- dplyr::arrange(results$pca_df, specimenID) |>
    dplyr::pull(extra_column)

  expect_identical(sort(results$pca_df$specimenID), sort(colnames(data)))
  expect_identical(received_extra_column, expected_extra_column)
  expect_identical(results$outliers, c("sample_out1", "sample_out2"))
})

test_that("find_pca_outliers uses metadata with different sample column name", {
  data <- matrix_for_pca_with_outliers()

  metadata <- data.frame(sample_ids = colnames(data),
                         extra_column = 1:ncol(data))

  results <- find_pca_outliers(data, metadata = metadata, sample_colname = "sample_ids")

  expect_true("sample_ids" %in% colnames(results$pca_df))
  expect_true("extra_column" %in% colnames(results$pca_df))
  expect_identical(sort(results$pca_df$sample_ids), sort(colnames(data)))

  expected_extra_column <- dplyr::arrange(metadata, sample_ids) |>
    dplyr::pull(extra_column)
  received_extra_column <- dplyr::arrange(results$pca_df, sample_ids) |>
    dplyr::pull(extra_column)

  expect_identical(received_extra_column, expected_extra_column)
})

test_that("find_pca_outliers correctly uses metadata that is out of order", {
  data <- matrix_for_pca_with_outliers()

  metadata <- data.frame(specimenID = colnames(data),
                         extra_column = 1:ncol(data)) |>
    dplyr::arrange(-extra_column) # Flip the order

  expected_extra_column <- dplyr::arrange(metadata, specimenID) |>
    dplyr::pull(extra_column)

  results <- find_pca_outliers(data, metadata = metadata)

  expect_identical(sort(results$pca_df$specimenID), sort(colnames(data)))

  expected_extra_column <- dplyr::arrange(metadata, specimenID) |>
    dplyr::pull(extra_column)
  received_extra_column <- dplyr::arrange(results$pca_df, specimenID) |>
    dplyr::pull(extra_column)

  expect_identical(received_extra_column, expected_extra_column)
})

test_that("find_pca_outliers finds no outliers with high n_sds argument", {
  data <- matrix_for_pca_with_outliers()
  results <- find_pca_outliers(data, n_sds = 10)

  expect_length(results$outliers, 0)
})

test_that("find_pca_outliers uses genes specified in gene_info", {
  data <- matrix_for_pca_with_outliers()

  # Ensures samples 21 and 22 won't be outliers unless we only use genes 1:10
  data[11:60, 21:22] <- data[11:60, 19:20]

  gene_info <- data.frame(ensembl_gene_id = paste0("gene", 1:nrow(data)),
                          gene_biotype = c(rep("protein_coding", 10), rep("non_coding", 50)),
                          chromosome_name = "2")

  results_all <- find_pca_outliers(data)
  results_sub <- find_pca_outliers(data, gene_info = gene_info)

  expect_length(results_all$outliers, 0)
  expect_length(results_sub$outliers, 2)
  expect_identical(results_sub$outliers, c("sample_out1", "sample_out2"))
})

test_that("find_pca_outliers can use gene_info with imperfect overlap with data", {
  data <- matrix_for_pca_with_outliers()

  # Half of data's genes exist in gene_info, plus gene_info has genes that don't
  # exist in data. Gene names are shuffled so they are not in the same order as
  # data's genes.
  gene_info <- data.frame(ensembl_gene_id = c(sample(rownames(data), 30, replace = FALSE),
                                              paste0("extra_gene", 1:10)),
                          gene_biotype = "protein_coding",
                          chromosome_name = "2")

  results <- find_pca_outliers(data)

  expect_false(all(gene_info$ensembl_gene_id %in% rownames(data)))
  expect_false(all(rownames(data) %in% gene_info$ensembl_gene_id))
  expect_length(results$outliers, 2)
  expect_identical(results$outliers, c("sample_out1", "sample_out2"))
})

test_that("find_pca_outliers checks for correct sample column name", {
  data <- matrix_for_pca_with_outliers()

  metadata <- data.frame(sample_ids = colnames(data),
                         extra_column = 1:ncol(data))

  expect_error(find_pca_outliers(data, metadata = metadata),
               regexp = "is not a valid column")
})

test_that("find_pca_outliers checks for missing samples in metadata", {
  data <- matrix_for_pca_with_outliers()

  metadata <- data.frame(specimenID = colnames(data)[1:10],
                         extra_column = 1:10)

  expect_error(find_pca_outliers(data, metadata = metadata),
               regexp = "missing samples")
})


# find_pca_outliers_by_group ---------------------------------------------------

test_that("find_pca_outliers_by_group splits samples into groups with no metadata", {
  data <- matrix_for_pca_with_groups()
  pca_group <- stringr::str_replace(colnames(data), "_.*", "")

  results_nosplit <- find_pca_outliers(data)
  results_split <- find_pca_outliers_by_group(data, pca_group)

  expect_identical(names(results_split), c("group_results", "outliers"))
  expect_length(results_nosplit$outliers, 0)
  expect_length(results_split$outliers, 4)
  expect_identical(results_split$outliers,
                   c("group1_sample_out1", "group1_sample_out2",
                     "group3_sample_out3", "group3_sample_out4"))
  expect_identical(names(results_split$group_results),
                   c("group1", "group2", "group3"))
})

test_that("find_pca_outliers_by_group splits samples by metadata column name", {
  data <- matrix_for_pca_with_groups()

  metadata <- data.frame(specimenID = colnames(data),
                         extra_column = 1:ncol(data),
                         group = stringr::str_replace(colnames(data), "_.*", ""))

  results_nosplit <- find_pca_outliers(data)
  results_split <- find_pca_outliers_by_group(data, pca_group = "group",
                                              metadata = metadata)

  expect_identical(names(results_split), c("group_results", "outliers"))
  expect_length(results_nosplit$outliers, 0)
  expect_length(results_split$outliers, 4)
  expect_identical(results_split$outliers,
                   c("group1_sample_out1", "group1_sample_out2",
                     "group3_sample_out3", "group3_sample_out4"))
  expect_identical(names(results_split$group_results),
                   c("group1", "group2", "group3"))

  expect_true(all(sapply(results_split$group_results, function(res) {
    "extra_column" %in% colnames(res$pca_df)
  })))

  expected_extra_column <- dplyr::arrange(metadata, specimenID) |>
    dplyr::pull(extra_column)
  received_extra_column <- lapply(results_split$group_results, function(res) {
    dplyr::arrange(res$pca_df, specimenID) |> dplyr::pull(extra_column)
  }) |>
    unlist() |>
    as.vector() # Remove names

  expect_identical(received_extra_column, expected_extra_column)
})

test_that("find_pca_outliers_by_group splits samples into groups with metadata column values", {
  data <- matrix_for_pca_with_groups()

  metadata <- data.frame(specimenID = colnames(data),
                         extra_column = 1:ncol(data),
                         group = stringr::str_replace(colnames(data), "_.*", ""))

  results_split <- find_pca_outliers_by_group(data, pca_group = metadata$group,
                                              metadata = metadata)

  expect_identical(names(results_split), c("group_results", "outliers"))
  expect_length(results_split$outliers, 4)
  expect_identical(results_split$outliers,
                   c("group1_sample_out1", "group1_sample_out2",
                     "group3_sample_out3", "group3_sample_out4"))
  expect_identical(names(results_split$group_results),
                   c("group1", "group2", "group3"))

  expect_true(all(sapply(results_split$group_results, function(res) {
    "extra_column" %in% colnames(res$pca_df)
  })))

  expected_extra_column <- dplyr::arrange(metadata, specimenID) |>
    dplyr::pull(extra_column)
  received_extra_column <- lapply(results_split$group_results, function(res) {
    dplyr::arrange(res$pca_df, specimenID) |> dplyr::pull(extra_column)
  }) |>
    unlist() |>
    as.vector() # Remove names

  expect_identical(received_extra_column, expected_extra_column)
})

test_that("find_pca_outliers_by_group checks for correct group column name", {
  data <- matrix_for_pca_with_groups()

  metadata <- data.frame(specimenID = colnames(data),
                         group = stringr::str_replace(colnames(data), "_.*", ""))

  expect_error(find_pca_outliers_by_group(data, pca_group = "bad_name",
                                          metadata = metadata),
               regexp = "is not a valid column")
})

test_that("find_pca_outliers_by_group checks for correct sample column name", {
  data <- matrix_for_pca_with_groups()

  metadata <- data.frame(sample_ids = colnames(data),
                         group = stringr::str_replace(colnames(data), "_.*", ""))

  expect_error(find_pca_outliers_by_group(data, pca_group = "group",
                                          metadata = metadata),
               regexp = "is not a valid column")
})

test_that("find_pca_outliers_by_group checks for missing samples in metadata", {
  data <- matrix_for_pca_with_groups()

  metadata <- data.frame(specimenID = colnames(data)[1:10],
                         group = stringr::str_replace(colnames(data)[1:10], "_.*", ""))

  expect_error(find_pca_outliers_by_group(data, pca_group = "group", metadata = metadata),
               regexp = "missing samples")
})

test_that("find_pca_outliers_by_group checks for missing samples in pca_group", {
  data <- matrix_for_pca_with_groups()
  pca_group <- stringr::str_replace(colnames(data)[1:10], "_.*", "")

  expect_error(find_pca_outliers_by_group(data, pca_group = pca_group),
               regexp = "number of samples")
})

test_that("find_pca_outliers_by_group ignores groups that are too small", {
  data <- matrix_for_pca_with_groups()
  pca_group <- stringr::str_replace(colnames(data), "_.*", "")
  pca_group[60:64] <- "small_group" # Removes outliers from group 3

  expect_message(
    {results <- find_pca_outliers_by_group(data, pca_group = pca_group)},
    regexp = "'small_group' is too small"
  )

  expect_length(results$group_results, 3)
  expect_length(results$outliers, 2)
  expect_identical(results$outliers, c("group1_sample_out1", "group1_sample_out2"))
  expect_identical(names(results$group_results), c("group1", "group2", "group3"))
})

test_that("find_pca_outliers_by_group returns empty lists when all groups are too small", {
  data <- matrix_for_pca_with_groups()[, 1:9]
  pca_group <- rep("group1", 9)

  expect_message(
    {results <- find_pca_outliers_by_group(data, pca_group = pca_group)},
    regexp = "'group1' is too small"
  )

  expect_length(results$group_results, 0)
  expect_null(results$outliers)
})

test_that("find_pca_outliers_by_group uses different min_group_size", {
  data <- matrix_for_pca_with_groups()[, 1:60] # Makes Group 3 have < 20 samples
  pca_group <- stringr::str_replace(colnames(data), "_.*", "")

  expect_message(
    {results <- find_pca_outliers_by_group(data, pca_group = pca_group,
                                           min_group_size = 20)},
    regexp = "'group3' is too small"
  )

  expect_length(results$group_results, 2)
  expect_length(results$outliers, 2)
  expect_identical(results$outliers, c("group1_sample_out1", "group1_sample_out2"))
  expect_identical(names(results$group_results), c("group1", "group2"))
})

# TODO test that n_sds and gene_info pass through to find_pca_outliers


# find_sex_mismatches ----------------------------------------------------------

test_that("find_sex_mismatches works with gene symbols", {
  data <- matrix_for_sex_mismatches_with_symbols()
  metadata <- data.frame(specimenID = colnames(data),
                         sex = c(rep("male", 10), rep("female", 10)))

  results <- find_sex_mismatches(metadata, data)

  expect_identical(names(results), c("sex_check_df", "y_expr_threshold", "mismatches"))
  expect_identical(results$y_expr_threshold, 2.0)
  expect_length(results$mismatches, 0)
})

test_that("find_sex_mismatches works with Ensembl IDs", {
  data <- matrix_for_sex_mismatches_with_ensembl_ids()
  metadata <- data.frame(specimenID = colnames(data),
                         sex = c(rep("male", 10), rep("female", 10)))

  results <- find_sex_mismatches(metadata, data)

  expect_identical(names(results), c("sex_check_df", "y_expr_threshold", "mismatches"))
  expect_identical(results$y_expr_threshold, 2.0)
  expect_length(results$mismatches, 0)
})

test_that("find_sex_mismatches works with versioned Ensembl IDs", {
  data <- matrix_for_sex_mismatches_with_ensembl_ids()
  rownames(data) <- paste0(rownames(data), ".", 1:nrow(data))

  metadata <- data.frame(specimenID = colnames(data),
                         sex = c(rep("male", 10), rep("female", 10)))

  results <- find_sex_mismatches(metadata, data)

  expect_identical(names(results), c("sex_check_df", "y_expr_threshold", "mismatches"))
  expect_identical(results$y_expr_threshold, 2.0)
  expect_length(results$mismatches, 0)
})

test_that("find_sex_mismatches correctly calculates mean Y expression", {
  data <- matrix_for_sex_mismatches_with_symbols()
  metadata <- data.frame(specimenID = colnames(data),
                         sex = c(rep("male", 10), rep("female", 10)))

  results <- find_sex_mismatches(metadata, data)

  # Expression values in sex_check_df should match the expression values of
  # the sex genes in data
  expect_equal(as.matrix(results$sex_check_df[, 3:7]),
               t(data[1:5, results$sex_check_df$specimenID]),
               ignore_attr = TRUE) # Ignore missing col/rownames

  mean_Y <- rowMeans(results$sex_check_df[4:7])
  expect_identical(results$sex_check_df$mean_Y, mean_Y)
})

test_that("find_sex_mismatches finds outliers", {
  data <- matrix_for_sex_mismatches_with_symbols()
  metadata <- data.frame(specimenID = colnames(data),
                         sex = c("female", rep("male", 9),
                                 rep("female", 9), "male"))

  results <- find_sex_mismatches(metadata, data)

  expect_length(results$mismatches, 2)
  expect_identical(results$mismatches, c("sample_F10", "sample_M1"))
})

test_that("find_sex_mismatches uses different sample column name", {
  data <- matrix_for_sex_mismatches_with_symbols()
  metadata <- data.frame(sample_ids = colnames(data),
                         sex = c("female", rep("male", 9),
                                 rep("female", 9), "male"))

  results <- find_sex_mismatches(metadata, data, sample_colname = "sample_ids")

  expect_length(results$mismatches, 2)
  expect_identical(results$mismatches, c("sample_F10", "sample_M1"))
})

test_that("find_sex_mismatches uses different sex column name", {
  data <- matrix_for_sex_mismatches_with_symbols()
  metadata <- data.frame(specimenID = colnames(data),
                         assigned_sex = c("female", rep("male", 9),
                                          rep("female", 9), "male"))

  results <- find_sex_mismatches(metadata, data, sex_colname = "assigned_sex")

  expect_length(results$mismatches, 2)
  expect_identical(results$mismatches, c("sample_F10", "sample_M1"))
})

test_that("find_sex_mismatches uses different Y expression threshold", {
  data <- matrix_for_sex_mismatches_with_symbols()
  metadata <- data.frame(specimenID = colnames(data),
                         sex = c("female", rep("male", 9),
                                 rep("female", 9), "male"))

  # Move the Y expression lower so that raising the threshold removes this
  # sample from the outliers list
  data[, 1] <- data[, 1] - 5

  results_2 <- find_sex_mismatches(metadata, data,
                                   y_expr_threshold = 2.0)
  results_8 <- find_sex_mismatches(metadata, data,
                                   y_expr_threshold = 8.0)

  expect_length(results_2$mismatches, 2)
  expect_length(results_8$mismatches, 1)
  expect_identical(results_8$mismatches, c("sample_F10"))
})

test_that("find_sex_mismatches checks for missing genes in data", {
  data <- matrix_for_sex_mismatches_with_symbols()
  data <- data[-2, ]

  metadata <- data.frame(specimenID = colnames(data),
                         sex = c(rep("male", 10), rep("female", 10)))

  expect_error(find_sex_mismatches(metadata, data),
               regexp = "missing the required genes")
})

test_that("find_sex_mismatches checks for correct sample column name", {
  data <- matrix_for_sex_mismatches_with_symbols()
  metadata <- data.frame(sample_ids = colnames(data),
                         sex = c(rep("male", 10), rep("female", 10)))

  expect_error(find_sex_mismatches(metadata, data),
               regexp = "not a valid column")
})

test_that("find_sex_mismatches checks for correct sex column name", {
  data <- matrix_for_sex_mismatches_with_symbols()
  metadata <- data.frame(specimenID = colnames(data),
                         assigned_sex = c(rep("male", 10), rep("female", 10)))

  expect_error(find_sex_mismatches(metadata, data),
               regexp = "not a valid column")
})

test_that("find_sex_mismatches checks for missing samples in metadata", {
  data <- matrix_for_sex_mismatches_with_symbols()
  metadata <- data.frame(specimenID = colnames(data)[1:5],
                         sex = rep("male", 5))

  expect_error(find_sex_mismatches(metadata, data),
               regexp = "missing samples")
})


# TODO fastqc read function tests
